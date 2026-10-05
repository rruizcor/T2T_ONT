"""VAF / depth at ClinVar P/LP SNVs with alt support, hg38 vs CHM13, stratified by
whether the two references differ within the +/- window around the site."""
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.stderr = open(snakemake.log[0], "w")
OUT = snakemake.params.outdir
MIN_ALT, MIN_COV = snakemake.params.min_alt, snakemake.params.min_cov


def read_counts(path):
    cols = ["chrom", "pos", "ref", "alts", "ad", "dp"]
    try:
        df = pd.read_csv(path, sep="\t", header=None, names=cols, dtype={"chrom": str})
    except pd.errors.EmptyDataError:
        return {}
    out = {}
    for r in df.itertuples(index=False):
        alleles = [r.ref] + str(r.alts).split(",")
        ad = [int(x) if x != "." else 0 for x in str(r.ad).split(",")]
        out[(r.chrom, int(r.pos))] = dict(zip(alleles, ad))
    return out


def wil(x, y, alt):
    d = np.asarray(x) - np.asarray(y)
    d = d[d != 0]
    return float(stats.wilcoxon(d, alternative=alt).pvalue) if len(d) >= 10 else np.nan


pairs = pd.read_csv(snakemake.input.pairs, sep="\t", dtype={"chrom_hg38": str, "chrom_chm13": str})
rows = []
for s in snakemake.params.samples:
    c38 = read_counts(f"{OUT}/counts/{s}.clinvar.hg38.tsv")
    c13 = read_counts(f"{OUT}/counts/{s}.clinvar.chm13.tsv")
    for r in pairs.itertuples(index=False):
        a = c38.get((r.chrom_hg38, int(r.pos_hg38)), {})
        b = c13.get((r.chrom_chm13, int(r.pos_chm13)), {})
        ref38, alt38 = a.get(r.ref_hg38, 0), a.get(r.alt_hg38, 0)
        ref13, alt13 = b.get(r.hg38ref_on_chm13, 0), b.get(r.hg38alt_on_chm13, 0)
        if max(alt38, alt13) < MIN_ALT:
            continue
        rows.append(
            dict(
                sample=s, site_id=r.site_id, gene=r.gene, clnid=r.clnid, clnsig=r.clnsig,
                chrom_chm13=r.chrom_chm13, pos_chm13=r.pos_chm13, swapped=r.swapped,
                prox_diff=r.prox_diff, prox_diff_in_hp=r.prox_diff_in_hp,
                prox_mismatch=r.prox_mismatch, prox_ins=r.prox_ins, prox_del=r.prox_del,
                alt_hg38=alt38, dp_hg38=ref38 + alt38, alt_chm13=alt13, dp_chm13=ref13 + alt13,
            )
        )
df = pd.DataFrame(rows)
if len(df):
    df["vaf_hg38"] = df.alt_hg38 / df.dp_hg38.replace(0, np.nan)
    df["vaf_chm13"] = df.alt_chm13 / df.dp_chm13.replace(0, np.nan)  # same (variant) allele
    df["delta_vaf"] = df.vaf_chm13 - df.vaf_hg38
df.to_csv(snakemake.output.sites, sep="\t", index=False)

tests = []
if len(df):
    ok = df[(df.dp_hg38 >= MIN_COV) & (df.dp_chm13 >= MIN_COV)]
    for lab, g in (("all", ok), ("no_difference", ok[ok.prox_diff == 0]),
                   ("difference", ok[ok.prox_diff == 1]),
                   ("difference_in_homopolymer", ok[ok.prox_diff_in_hp == 1])):
        tests.append(
            dict(
                group=lab, n_sample_sites=len(g),
                mean_vaf_hg38=g.vaf_hg38.mean(), mean_vaf_chm13=g.vaf_chm13.mean(),
                mean_delta_vaf=g.delta_vaf.mean(),
                mean_dp_hg38=g.dp_hg38.mean(), mean_dp_chm13=g.dp_chm13.mean(),
                p_vaf_chm13_gt_hg38=wil(g.vaf_chm13, g.vaf_hg38, "greater"),
                p_dp_chm13_gt_hg38=wil(g.dp_chm13, g.dp_hg38, "greater"),
            )
        )
pd.DataFrame(tests, columns=["group", "n_sample_sites", "mean_vaf_hg38", "mean_vaf_chm13",
                             "mean_delta_vaf", "mean_dp_hg38", "mean_dp_chm13",
                             "p_vaf_chm13_gt_hg38", "p_dp_chm13_gt_hg38"]).to_csv(
    snakemake.output.tests, sep="\t", index=False)
