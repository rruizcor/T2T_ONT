"""Join hg38 candidate sites with their Picard-lifted CHM13 positions.

For every site:
  * allele mapping hg38 -> CHM13 (swap / reverse-complement aware)
  * homopolymer proximity on BOTH references (ref and alt haplotype)
  * exclusion BEDs on both references
  * +/- window sequence comparison between the references (edlib global alignment):
    mismatches, insertions, deletions, and whether any difference lies in a homopolymer
Outputs a pairs table and per-reference CHROM\tPOS target lists for bcftools -R.
"""
import bisect
import re
import sys
from collections import Counter, defaultdict

import edlib
import pysam

sys.stderr = open(snakemake.log[0], "w")

W = int(snakemake.params.window)
HP = int(snakemake.params.hp_len)
HPF = int(snakemake.params.hp_flank)
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def read_vcf(path):
    rows = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = {}
            for kv in f[7].split(";"):
                if kv == ".":
                    continue
                k, _, v = kv.partition("=")
                info[k] = v if v else True
            rows[f[2]] = dict(chrom=f[0], pos=int(f[1]), ref=f[3].upper(), alt=f[4].upper(), info=info)
    return rows


def read_fai(path):
    order = {}
    lens = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            n, l = line.split("\t")[:2]
            order[n] = i
            lens[n] = int(l)
    return order, lens


def load_bed(path):
    if not path:
        return {}
    iv = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            iv[f[0]].append((int(f[1]), int(f[2])))
    out = {}
    for c, lst in iv.items():
        lst.sort()
        merged = []
        for s, e in lst:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = ([m[0] for m in merged], [m[1] for m in merged])
    return out


def in_bed(bed, chrom, pos0):
    if chrom not in bed:
        return False
    starts, ends = bed[chrom]
    i = bisect.bisect_right(starts, pos0) - 1
    return i >= 0 and pos0 < ends[i]


def hp_mask(seq):
    """Boolean list: position lies in a run of >= HP identical bases."""
    mask = [False] * len(seq)
    for m in re.finditer(r"(A+|C+|G+|T+)", seq):
        if m.end() - m.start() >= HP:
            for i in range(m.start(), m.end()):
                mask[i] = True
    return mask


def near_homopolymer(fa, chrom, pos1, alt, clen):
    """Run >= HP touching [pos-HPF, pos+HPF] on the ref OR alt haplotype."""
    s = max(0, pos1 - 1 - HPF - HP)
    e = min(clen, pos1 + HPF + HP)
    seq = fa.fetch(chrom, s, e).upper()
    i = pos1 - 1 - s
    lo, hi = max(0, i - HPF), min(len(seq), i + HPF + 1)
    for hap in (seq, seq[:i] + alt + seq[i + 1 :]):
        m = hp_mask(hap)
        if any(m[lo:hi]):
            return True
    return False


def window(fa, chrom, pos1, clen):
    s = max(0, pos1 - 1 - W)
    e = min(clen, pos1 + W)
    return fa.fetch(chrom, s, e).upper()


def compare(seq38, seq13):
    """Global alignment of the two reference windows. Returns mm, ins, del, diff_in_hp
    (ins/del relative to hg38: ins = bases present in CHM13 only)."""
    if seq38 == seq13:
        return 0, 0, 0, False
    res = edlib.align(seq13, seq38, mode="NW", task="path")
    mask = hp_mask(seq38)
    mm = ins = dele = 0
    hp_hit = False
    t = 0  # position in hg38 (target)
    for n, op in re.findall(r"(\d+)([=XID])", res["cigar"]):
        n = int(n)
        if op == "=":
            t += n
        elif op == "X":
            mm += n
            hp_hit |= any(mask[t : t + n])
            t += n
        elif op == "I":  # extra bases in query (CHM13)
            ins += n
            hp_hit |= any(mask[max(0, t - 1) : t + 1])
        elif op == "D":  # bases in hg38 missing from CHM13
            dele += n
            hp_hit |= any(mask[t : t + n])
            t += n
    return mm, ins, dele, hp_hit


h38 = read_vcf(snakemake.input.hg38)
l13 = read_vcf(snakemake.input.chm13)
fa38 = pysam.FastaFile(snakemake.input.fa_hg38)
fa13 = pysam.FastaFile(snakemake.input.fa_chm13)
ord38, len38 = read_fai(snakemake.input.fai_hg38)
ord13, len13 = read_fai(snakemake.input.fai_chm13)
bed38 = load_bed(snakemake.params.exclude_hg38)
bed13 = load_bed(snakemake.params.exclude_chm13)

counts = Counter()
counts["hg38_candidates"] = len(h38)
with open(snakemake.input.rej) as fh:
    for line in fh:
        if not line.startswith("#"):
            counts["liftover_rejected:" + line.split("\t")[6]] += 1

# several hg38 sites can lift to the same CHM13 position -> ambiguous, drop all
dest = Counter((r["chrom"], r["pos"]) for r in l13.values())

rows = []
for key, a in h38.items():
    b = l13.get(key)
    if b is None:
        continue
    counts["lifted"] += 1
    if dest[(b["chrom"], b["pos"])] > 1:
        counts["drop:many_to_one_liftover"] += 1
        continue
    swapped = "SwappedAlleles" in b["info"]
    revcomp = "ReverseComplementedAlleles" in b["info"]
    a_ref, a_alt = a["ref"], a["alt"]
    if revcomp:
        a_ref, a_alt = rc(a_ref), rc(a_alt)
    exp_ref, exp_alt = (a_alt, a_ref) if swapped else (a_ref, a_alt)
    if (b["ref"], b["alt"]) != (exp_ref, exp_alt):
        counts["drop:allele_inconsistent"] += 1
        continue
    if in_bed(bed13, b["chrom"], b["pos"] - 1):
        counts["drop:excluded_bed_chm13"] += 1
        continue
    if in_bed(bed38, a["chrom"], a["pos"] - 1):
        counts["drop:excluded_bed_hg38"] += 1
        continue
    if near_homopolymer(fa38, a["chrom"], a["pos"], a["alt"], len38[a["chrom"]]) or near_homopolymer(
        fa13, b["chrom"], b["pos"], b["alt"], len13[b["chrom"]]
    ):
        counts["drop:homopolymer"] += 1
        continue
    s38 = window(fa38, a["chrom"], a["pos"], len38[a["chrom"]])
    s13 = window(fa13, b["chrom"], b["pos"], len13[b["chrom"]])
    if revcomp:
        s13 = rc(s13)
    mm, ins, dele, hp_hit = compare(s38, s13)
    rows.append(
        dict(
            site_id=key,
            chrom_hg38=a["chrom"], pos_hg38=a["pos"], ref_hg38=a["ref"], alt_hg38=a["alt"],
            chrom_chm13=b["chrom"], pos_chm13=b["pos"], ref_chm13=b["ref"], alt_chm13=b["alt"],
            # the hg38 REF / ALT alleles, as written on the CHM13 strand
            hg38ref_on_chm13=a_ref, hg38alt_on_chm13=a_alt,
            swapped=int(swapped), revcomp=int(revcomp),
            prox_mismatch=mm, prox_ins=ins, prox_del=dele,
            prox_diff=int(mm + ins + dele > 0), prox_diff_in_hp=int(hp_hit),
            cand_vaf=a["info"].get("CAND_VAF", "."),
            clnsig=a["info"].get("CLNSIG", "."), gene=a["info"].get("GENE", "."),
            clnid=a["info"].get("CLNID", "."),
        )
    )
    counts["kept"] += 1
    if swapped:
        counts["kept:swapped"] += 1
    if mm + ins + dele:
        counts["kept:proximal_diff"] += 1

cols = list(rows[0].keys()) if rows else [
    "site_id", "chrom_hg38", "pos_hg38", "ref_hg38", "alt_hg38", "chrom_chm13", "pos_chm13",
    "ref_chm13", "alt_chm13", "hg38ref_on_chm13", "hg38alt_on_chm13", "swapped", "revcomp",
    "prox_mismatch", "prox_ins", "prox_del", "prox_diff", "prox_diff_in_hp", "cand_vaf",
    "clnsig", "gene", "clnid",
]
with open(snakemake.output.pairs, "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda r: (ord38[r["chrom_hg38"]], r["pos_hg38"])):
        out.write("\t".join(str(r[c]) for c in cols) + "\n")

for ref, order, cc, pc in (
    ("hg38", ord38, "chrom_hg38", "pos_hg38"),
    ("chm13", ord13, "chrom_chm13", "pos_chm13"),
):
    tg = sorted({(r[cc], r[pc]) for r in rows}, key=lambda x: (order[x[0]], x[1]))
    with open(getattr(snakemake.output, f"tgt_{ref}"), "w") as out:
        for c, p in tg:
            out.write(f"{c}\t{p}\n")

with open(snakemake.output.log_tsv, "w") as out:
    with open(snakemake.input.prep_counts) as fh:
        out.write(fh.read())
    for k, v in counts.items():
        out.write(f"pair:{k}\t{v}\n")
for k, v in counts.items():
    print(f"{k}\t{v}", file=sys.stderr)
