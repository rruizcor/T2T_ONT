"""Paired mapping metrics (samtools stats SN section), hg38 vs T2T-CHM13."""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stderr = open(snakemake.log[0], "w")

KEEP = {
    "raw total sequences": "reads_total",
    "reads mapped": "reads_mapped",
    "reads MQ0": "reads_mq0",
    "bases mapped (cigar)": "bases_mapped_cigar",
    "mismatches": "mismatches",
    "error rate": "error_rate",
    "average quality": "average_quality",
}
rows = []
for p in snakemake.input:
    m = re.match(r"(.+)\.(hg38|chm13)\.samtools_stats\.txt$", Path(p).name)
    sample, ref = m.group(1), m.group(2)
    rec = dict(sample=sample, ref=ref)
    with open(p) as fh:
        for line in fh:
            if line.startswith("SN\t"):
                _, key, val = line.rstrip("\n").split("\t")[:3]
                key = key.rstrip(":")
                if key in KEEP:
                    rec[KEEP[key]] = float(val)
    rows.append(rec)
df = pd.DataFrame(rows)
df["pct_mapped"] = 100 * df.reads_mapped / df.reads_total
df["pct_mq0"] = 100 * df.reads_mq0 / df.reads_mapped
df.sort_values(["sample", "ref"]).to_csv(snakemake.output.tsv, sep="\t", index=False)

wide = df.pivot(index="sample", columns="ref")
tests = []
for metric in ["pct_mapped", "pct_mq0", "error_rate", "average_quality", "bases_mapped_cigar"]:
    a, b = wide[(metric, "hg38")], wide[(metric, "chm13")]
    d = (b - a).dropna()
    p = float(stats.wilcoxon(d).pvalue) if len(d) >= 6 and (d != 0).any() else np.nan
    tests.append(
        dict(
            metric=metric, n_samples=len(d), mean_hg38=a.mean(), mean_chm13=b.mean(),
            mean_rel_change_pct=float(100 * ((b - a) / a).mean()),
            n_chm13_higher=int((d > 0).sum()), n_chm13_lower=int((d < 0).sum()),
            p_wilcoxon_paired_two_sided=p,
        )
    )
pd.DataFrame(tests).to_csv(snakemake.output.tests, sep="\t", index=False)
