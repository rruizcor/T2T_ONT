"""Simulate a tiny two-reference test set.

hg38.fa  : random genome (2 contigs)
chm13.fa : hg38 + scattered substitutions (some clustered near het SNPs), small indels
           (shift coordinates -> real chain needed), and "swap" sites where CHM13 carries
           the sample's alt allele.
The sample's haplotypes are built on the CHM13 sequence (CHM13 = "more correct"), with het SNPs.
Long reads (~3% error, indel-rich) are simulated from both haplotypes.
Also writes a hg38->chm13 chain, a ClinVar-like VCF and exclusion BEDs.
"""
import gzip
import os
import random
import subprocess
import sys

out = sys.argv[1] if len(sys.argv) > 1 else "test/data"
os.makedirs(out, exist_ok=True)
rng = random.Random(42)
B = "ACGT"
CONTIGS = {"chr1": 300_000, "chr2": 200_000}
N_HET = {"chr1": 450, "chr2": 300}


def rand_seq(n):
    s = []
    while len(s) < n:
        if rng.random() < 0.004:  # sprinkle homopolymers
            s += [rng.choice(B)] * rng.randint(5, 9)
        else:
            s.append(rng.choice(B))
    return s[:n]


def apply(seq, subs, indels):
    """subs: {pos0: base}; indels: sorted [(pos0, 'ins', seq) | (pos0, 'del', L)]."""
    s = list(seq)
    for p, b in subs.items():
        s[p] = b
    res, last = [], 0
    for p, kind, pay in indels:
        if kind == "ins":
            res += s[last : p + 1] + list(pay)
            last = p + 1
        else:
            res += s[last:p]
            last = p + pay
    res += s[last:]
    return "".join(res)


fa38, fa13, haps, chains, clinvar, het_truth = {}, {}, {}, [], [], []
ex38, ex13 = [], []
for ci, (chrom, L) in enumerate(CONTIGS.items()):
    seq = rand_seq(L)
    fa38[chrom] = "".join(seq)
    # indels in hg38 coords, spaced out
    indels = []
    for p in sorted(rng.sample(range(5_000, L - 5_000, 997), 6)):
        if rng.random() < 0.5:
            indels.append((p, "ins", "".join(rng.choice(B) for _ in range(rng.randint(5, 150)))))
        else:
            indels.append((p, "del", rng.randint(5, 150)))
    blocked = set()
    for p, k, pay in indels:
        blocked.update(range(p - 400, p + 400 + (pay if k == "del" else 0)))
    # het SNPs
    cand = [p for p in range(1_000, L - 1_000) if p not in blocked]
    het = sorted(rng.sample(cand, N_HET[chrom]))
    het = [p for i, p in enumerate(het) if i == 0 or p - het[i - 1] > 60]
    subs13, alt_sub, hap1_override = {}, {}, {}
    near_err = set(rng.sample(het, len(het) // 4))  # 25% of hets get nearby ref differences
    swaps = set(rng.sample([p for p in het if p not in near_err], 12))
    for p in het:
        ref = seq[p]
        alt = rng.choice([b for b in B if b != ref])
        alt_sub[p] = alt
        if p in swaps:  # CHM13 carries the alt; sample is ref/alt
            subs13[p] = alt
            hap1_override[p] = ref
        if p in near_err:
            for _ in range(rng.randint(1, 4)):
                q = p + rng.choice([-1, 1]) * rng.randint(2, 40)
                if q not in het and q not in blocked and 0 < q < L:
                    subs13[q] = rng.choice([b for b in B if b != seq[q]])
    # background differences far from hets
    for q in rng.sample(range(L), 150):
        if q not in blocked and all(abs(q - h) > 160 for h in het):
            subs13[q] = rng.choice([b for b in B if b != seq[q]])
    fa13[chrom] = apply(seq, subs13, indels)
    h1 = dict(subs13)
    h1.update(hap1_override)
    h2 = dict(subs13)
    h2.update(alt_sub)
    haps[chrom] = (apply(seq, h1, indels), apply(seq, h2, indels))
    # chain hg38 (t) -> chm13 (q)
    Lq = len(fa13[chrom])
    lines, t0 = [f"chain 1000 {chrom} {L} + 0 {L} {chrom} {Lq} + 0 {Lq} {ci + 1}"], 0
    for p, k, pay in indels:
        if k == "ins":
            lines.append(f"{p + 1 - t0}\t0\t{len(pay)}")
            t0 = p + 1
        else:
            lines.append(f"{p - t0}\t{pay}\t0")
            t0 = p + pay
    lines.append(f"{L - t0}")
    chains.append("\n".join(lines) + "\n\n")
    # ClinVar-like sites: half at het positions (sample carries them), rest random
    for i, p in enumerate(sorted(rng.sample(het, 40))):
        sig = rng.choice(["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"])
        clinvar.append((chrom.replace("chr", ""), p + 1, f"{ci}{i:04d}", seq[p], alt_sub[p], sig, f"GENE{i % 7}:{i}"))
    for i, p in enumerate(sorted(rng.sample(cand, 40))):
        clinvar.append((chrom.replace("chr", ""), p + 1, f"{ci}9{i:03d}", seq[p],
                        rng.choice([b for b in B if b != seq[p]]), "Benign" if i % 3 == 0 else "Pathogenic", "GENEX:1"))
    het_truth += [(chrom, p + 1, seq[p], alt_sub[p], int(p in swaps), int(p in near_err)) for p in het]
    ex38.append((chrom, 250_000 if L > 250_000 else L - 20_000, L - 10_000))

with open(f"{out}/hg38.fa", "w") as f:
    for c, s in fa38.items():
        f.write(f">{c}\n" + "\n".join(s[i : i + 80] for i in range(0, len(s), 80)) + "\n")
with open(f"{out}/chm13.fa", "w") as f:
    for c, s in fa13.items():
        f.write(f">{c}\n" + "\n".join(s[i : i + 80] for i in range(0, len(s), 80)) + "\n")
with gzip.open(f"{out}/hg38ToChm13.over.chain.gz", "wt") as f:
    f.write("".join(chains))
with open(f"{out}/hg38.exclude.bed", "w") as f:
    for c, s, e in ex38:
        f.write(f"{c}\t{s}\t{e}\n")
open(f"{out}/chm13.exclude.bed", "w").close()
with open(f"{out}/truth_hets.tsv", "w") as f:
    f.write("chrom\tpos\tref\talt\tswap\tnear_ref_diff\n")
    for r in het_truth:
        f.write("\t".join(map(str, r)) + "\n")
with open(f"{out}/clinvar.vcf", "w") as f:
    f.write("##fileformat=VCFv4.1\n")
    for c, l in CONTIGS.items():
        f.write(f"##contig=<ID={c.replace('chr', '')},length={l}>\n")
    f.write('##INFO=<ID=CLNSIG,Number=.,Type=String,Description="x">\n')
    f.write('##INFO=<ID=GENEINFO,Number=1,Type=String,Description="x">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for c, p, i, r, a, sig, g in sorted(clinvar, key=lambda x: (x[0], x[1])):
        f.write(f"{c}\t{p}\t{i}\t{r}\t{a}\t.\t.\tCLNSIG={sig};GENEINFO={g}\n")
subprocess.run(f"bgzip -f {out}/clinvar.vcf && tabix -f -p vcf {out}/clinvar.vcf.gz", shell=True, check=True)


def mutate(s, err=0.03):
    o = []
    for ch in s:
        r = rng.random()
        if r < err / 3:
            o.append(rng.choice(B))
        elif r < 2 * err / 3:
            o.append(ch + rng.choice(B))
        elif r < err:
            continue
        else:
            o.append(ch)
    return "".join(o)


comp = str.maketrans("ACGT", "TGCA")
for sample, depth in (("S1", 45), ("S2", 70)):
    path = f"{out}/{sample}.fastq.gz"
    with gzip.open(path, "wt") as f:
        n = 0
        for chrom, (h1, h2) in haps.items():
            for hap in (h1, h2):
                total = 0
                while total < depth / 2 * len(hap):
                    ln = min(len(hap) - 1, int(rng.lognormvariate(8.7, 0.5)))
                    st = rng.randint(0, len(hap) - ln)
                    s = hap[st : st + ln]
                    if rng.random() < 0.5:
                        s = s.translate(comp)[::-1]
                    s = mutate(s)
                    f.write(f"@{sample}_{n}\n{s}\n+\n{'5' * len(s)}\n")
                    n += 1
                    total += ln
    if sample == "S2":  # also exercise the unaligned-BAM input path
        subprocess.run(f"samtools import -0 {path} -o {out}/{sample}.ubam", shell=True, check=True)
print("ok", {c: len(s) for c, s in fa13.items()})
