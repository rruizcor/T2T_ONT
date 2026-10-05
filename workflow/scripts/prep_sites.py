"""Build the hg38 candidate site VCF (sites-only) with ID = hg38 key.

mode == "het"     : biallelic heterozygous SNVs from a germline caller VCF
mode == "clinvar" : ClinVar Pathogenic / Likely_pathogenic SNVs
"""
import bisect
import re
import sys
from collections import Counter, defaultdict

import pysam

sys.stderr = open(snakemake.log[0], "w")

mode = snakemake.params.mode
cfg = snakemake.params.cfg
excl_bed = snakemake.params.exclude_bed

contigs = []
with open(snakemake.input.fai) as fh:
    for line in fh:
        name, length = line.split("\t")[:2]
        contigs.append((name, int(length)))
contig_order = {c: i for i, (c, _) in enumerate(contigs)}


def harmonise(chrom):
    if chrom in contig_order:
        return chrom
    alt = {"MT": "chrM", "M": "chrM"}.get(chrom, "chr" + chrom)
    if alt in contig_order:
        return alt
    if chrom.startswith("chr") and chrom[3:] in contig_order:
        return chrom[3:]
    return None


def load_bed(path):
    iv = defaultdict(list)
    if not path:
        return {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            c = harmonise(f[0]) or f[0]
            iv[c].append((int(f[1]), int(f[2])))
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


bed = load_bed(excl_bed)
fasta = pysam.FastaFile(snakemake.input.fai[: -len(".fai")])
counts = Counter()
records = []


def sanitize(x):
    return re.sub(r"[;=\s,]", "_", str(x))


vcf = pysam.VariantFile(snakemake.input.vcf)
sample0 = list(vcf.header.samples)[0] if mode == "het" and len(vcf.header.samples) else None
if mode == "het" and sample0 is None:
    raise SystemExit("het mode needs a VCF with at least one sample column")

for rec in vcf:
    counts["input"] += 1
    if rec.alts is None or len(rec.alts) != 1:
        counts["not_biallelic"] += 1
        continue
    ref, alt = rec.ref.upper(), rec.alts[0].upper()
    if len(ref) != 1 or len(alt) != 1 or alt not in "ACGT" or ref not in "ACGT":
        counts["not_snv"] += 1
        continue
    chrom = harmonise(rec.chrom)
    if chrom is None:
        counts["contig_not_in_ref"] += 1
        continue
    info = {}

    if mode == "het":
        filt = list(rec.filter.keys())
        if filt and filt != ["PASS"]:
            counts["filter_fail"] += 1
            continue
        if rec.qual is not None and rec.qual < cfg["min_qual"]:
            counts["low_qual"] += 1
            continue
        smp = rec.samples[sample0]
        gt = smp.get("GT")
        if gt is None or None in gt or sorted(gt) != [0, 1]:
            counts["not_het"] += 1
            continue
        dp = smp.get("DP") if "DP" in rec.format else None
        if dp is None:
            dp = rec.info.get("DP")
        if dp is not None and dp < cfg["min_depth"]:
            counts["low_depth"] += 1
            continue
        vaf = None
        if "AF" in rec.format and smp.get("AF") is not None:
            af = smp.get("AF")
            vaf = af[0] if isinstance(af, tuple) else af
        elif "AD" in rec.format and smp.get("AD") is not None:
            ad = smp.get("AD")
            if ad[0] is not None and ad[1] is not None and ad[0] + ad[1] > 0:
                vaf = ad[1] / (ad[0] + ad[1])
        if vaf is not None and not (cfg["vaf_min"] <= vaf <= cfg["vaf_max"]):
            counts["vaf_out_of_range"] += 1
            continue
        info["CAND_VAF"] = f"{vaf:.4f}" if vaf is not None else "."
    else:  # clinvar
        sig = rec.info.get("CLNSIG")
        sig = "|".join(sig) if isinstance(sig, tuple) else (sig or "")
        if not re.search(r"(^|[|/_,])(Pathogenic|Likely_pathogenic)", sig) or "Conflicting" in sig:
            counts["not_pathogenic"] += 1
            continue
        gene = rec.info.get("GENEINFO")
        gene = gene[0] if isinstance(gene, tuple) else gene
        info["CLNSIG"] = sanitize(sig)
        info["GENE"] = sanitize((gene or ".").split(":")[0])
        info["CLNID"] = sanitize(rec.id or ".")

    if in_bed(bed, chrom, rec.pos - 1):
        counts["excluded_bed_hg38"] += 1
        continue
    if fasta.fetch(chrom, rec.pos - 1, rec.pos).upper() != ref:
        counts["ref_mismatch_hg38"] += 1
        continue
    key = f"{chrom}:{rec.pos}:{ref}:{alt}"
    records.append((contig_order[chrom], rec.pos, chrom, key, ref, alt, info))
    counts["kept"] += 1

records.sort()
seen = set()
with open(snakemake.output.vcf, "w") as out:
    out.write("##fileformat=VCFv4.2\n")
    for c, ln in contigs:
        out.write(f"##contig=<ID={c},length={ln}>\n")
    out.write('##INFO=<ID=CAND_VAF,Number=1,Type=String,Description="Caller VAF">\n')
    out.write('##INFO=<ID=CLNSIG,Number=1,Type=String,Description="ClinVar significance">\n')
    out.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">\n')
    out.write('##INFO=<ID=CLNID,Number=1,Type=String,Description="ClinVar ID">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for _, pos, chrom, key, ref, alt, info in records:
        if key in seen:
            continue
        seen.add(key)
        inf = ";".join(f"{k}={v}" for k, v in info.items()) or "."
        out.write(f"{chrom}\t{pos}\t{key}\t{ref}\t{alt}\t.\t.\t{inf}\n")

with open(snakemake.output.counts, "w") as fh:
    fh.write("stage\tn\n")
    for k, v in counts.items():
        fh.write(f"prep:{k}\t{v}\n")
        print(f"{k}\t{v}", file=sys.stderr)
