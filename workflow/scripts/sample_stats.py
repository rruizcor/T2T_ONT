"""Per-sample hg38 vs CHM13 reference-mapping-bias statistics (Cherchi et al. 2026 framework,
adapted for long reads).

VAF definition follows the paper: non-reference reads / (ref + alt) in EACH assembly.
For allele-swapped sites the CHM13 non-reference allele is the hg38 REF, so those sites
are reported separately and excluded from the paired tests.
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.stderr = open(snakemake.log[0], "w")
P = snakemake.params.stats
MIN_COV, THR, HIGH = P["min_cov"], P["delta_vaf"], P["high_cov"]
BINS = list(P["coverage_bins"]) + [np.inf]
C38, C13 = "#2a78d6", "#eb6834"  # categorical slots 1, 2
INK, MUTED, GRID = "#1f1f1e", "#6b6a63", "#e6e5df"


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


def binned_kde(v, lo=0.2, hi=0.8, nbin=1200):
    """Fast binned Gaussian KDE (Scott bandwidth) on [lo, hi]. Returns grid, density."""
    v = np.asarray(v, dtype=float)
    v = v[(v >= lo) & (v <= hi)]
    grid = np.linspace(lo, hi, nbin + 1)
    centers = (grid[:-1] + grid[1:]) / 2
    if len(v) < 50 or np.std(v) == 0:
        return centers, None
    h, _ = np.histogram(v, bins=grid)
    bw = 1.06 * np.std(v) * len(v) ** (-1 / 5)
    sd = bw / (grid[1] - grid[0])
    half = int(np.ceil(4 * sd))
    k = np.exp(-0.5 * (np.arange(-half, half + 1) / sd) ** 2)
    d = np.convolve(h, k / k.sum(), mode="same")
    return centers, d / (d.sum() * (grid[1] - grid[0]))


def kde_peak(v, lo=0.2, hi=0.8):
    """Local maximum of the VAF density closest to 0.5 (paper's localMaxima approach)."""
    grid, d = binned_kde(v, lo, hi)
    if d is None:
        return np.nan
    idx = np.where((d[1:-1] > d[:-2]) & (d[1:-1] >= d[2:]))[0] + 1
    if len(idx) == 0:
        return float(grid[np.argmax(d)])
    return float(grid[idx[np.argmin(np.abs(grid[idx] - 0.5))]])


def wil_greater(x, y):
    d = np.asarray(x) - np.asarray(y)
    d = d[d != 0]
    if len(d) < 10:
        return np.nan
    return float(stats.wilcoxon(d, alternative="greater").pvalue)


pairs = pd.read_csv(snakemake.input.pairs, sep="\t", dtype={"chrom_hg38": str, "chrom_chm13": str})
c38 = read_counts(snakemake.input.c_hg38)
c13 = read_counts(snakemake.input.c_chm13)

rec = []
for r in pairs.itertuples(index=False):
    a = c38.get((r.chrom_hg38, int(r.pos_hg38)), {})
    b = c13.get((r.chrom_chm13, int(r.pos_chm13)), {})
    ref38, alt38 = a.get(r.ref_hg38, 0), a.get(r.alt_hg38, 0)
    # counts of the SAME biological alleles in the CHM13 alignment
    r13_same, a13_same = b.get(r.hg38ref_on_chm13, 0), b.get(r.hg38alt_on_chm13, 0)
    oth38 = sum(a.values()) - ref38 - alt38
    oth13 = sum(b.values()) - r13_same - a13_same
    rec.append((ref38, alt38, oth38, r13_same, a13_same, oth13))
cnt = pd.DataFrame(rec, columns=["ref38", "alt38", "other38", "hg38ref13", "hg38alt13", "other13"])
df = pd.concat([pairs.reset_index(drop=True), cnt], axis=1)
df["dp38"] = df.ref38 + df.alt38
df["dp13"] = df.hg38ref13 + df.hg38alt13
with np.errstate(invalid="ignore", divide="ignore"):
    df["vaf38"] = df.alt38 / df.dp38
    same = df.hg38alt13 / df.dp13
    # CHM13 non-reference VAF: for swapped sites the non-ref allele is the hg38 REF
    df["vaf13"] = np.where(df.swapped == 1, 1 - same, same)
    df["vaf13_same_allele"] = same
df["min_dp"] = df[["dp38", "dp13"]].min(axis=1)
df["delta_vaf"] = df.vaf13 - df.vaf38
df["incoherent"] = (df.delta_vaf.abs() >= THR).astype(int)
df.to_csv(snakemake.output.sites, sep="\t", index=False, compression="gzip")

ok = df[(df.min_dp >= MIN_COV)]
main = ok[ok.swapped == 0]
sw = ok[ok.swapped == 1]
hi = main[main.min_dp >= HIGH]
inc = main[main.incoherent == 1]

summary = dict(
    sample=snakemake.wildcards.sample,
    n_pairs=len(df),
    n_cov_ok=len(ok),
    n_main=len(main),
    n_swapped=len(sw),
    median_depth_hg38=float(main.dp38.median()) if len(main) else np.nan,
    median_depth_chm13=float(main.dp13.median()) if len(main) else np.nan,
    p_depth_chm13_gt_hg38=wil_greater(main.dp13, main.dp38),
    mean_vaf_hg38=float(main.vaf38.mean()),
    mean_vaf_chm13=float(main.vaf13.mean()),
    median_vaf_hg38=float(main.vaf38.median()),
    median_vaf_chm13=float(main.vaf13.median()),
    n_high_cov=len(hi),
    peak_hg38_highcov=kde_peak(hi.vaf38),
    peak_chm13_highcov=kde_peak(hi.vaf13),
    peak_hg38_all=kde_peak(main.vaf38),
    peak_chm13_all=kde_peak(main.vaf13),
    frac_incoherent=float(main.incoherent.mean()) if len(main) else np.nan,
    n_incoherent=len(inc),
    n_inc_chm13_higher=int((inc.delta_vaf > 0).sum()),
    n_inc_hg38_higher=int((inc.delta_vaf < 0).sum()),
    p_wilcoxon_incoherent_chm13_gt_hg38=wil_greater(inc.vaf13, inc.vaf38),
    median_delta_vaf_incoherent=float(inc.delta_vaf.median()) if len(inc) else np.nan,
    swapped_median_vaf_hg38=float(sw.vaf38.median()) if len(sw) else np.nan,
    swapped_median_vaf_chm13=float(sw.vaf13.median()) if len(sw) else np.nan,
)
pd.DataFrame([summary]).to_csv(snakemake.output.summary, sep="\t", index=False)

# coverage-binned
brows = []
main = main.assign(cov_bin=pd.cut(main.min_dp, BINS, right=False))
for b, g in main.groupby("cov_bin", observed=True):
    gi = g[g.incoherent == 1]
    brows.append(
        dict(
            sample=summary["sample"], cov_bin=str(b), n=len(g),
            peak_hg38=kde_peak(g.vaf38), peak_chm13=kde_peak(g.vaf13),
            median_vaf_hg38=g.vaf38.median(), median_vaf_chm13=g.vaf13.median(),
            frac_incoherent=g.incoherent.mean(),
            n_inc_chm13_higher=int((gi.delta_vaf > 0).sum()),
            n_inc_hg38_higher=int((gi.delta_vaf < 0).sum()),
            p_wilcoxon=wil_greater(gi.vaf13, gi.vaf38),
        )
    )
pd.DataFrame(brows).to_csv(snakemake.output.bins, sep="\t", index=False)

# proximal sequence difference between references vs VAF discordance
prow = []
for lab, g in (
    ("no_difference", main[main.prox_diff == 0]),
    ("difference", main[main.prox_diff == 1]),
    ("difference_in_homopolymer", main[main.prox_diff_in_hp == 1]),
):
    prow.append(
        dict(
            sample=summary["sample"], group=lab, n=len(g),
            frac_incoherent=g.incoherent.mean() if len(g) else np.nan,
            median_delta_vaf=g.delta_vaf.median() if len(g) else np.nan,
            median_abs_delta_vaf=g.delta_vaf.abs().median() if len(g) else np.nan,
            p_wilcoxon_chm13_gt_hg38=wil_greater(g.vaf13, g.vaf38),
        )
    )
tab = pd.crosstab(main.prox_diff, main.incoherent).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
fisher_p = float(stats.fisher_exact(tab.values, alternative="greater")[1]) if tab.values.sum() else np.nan
pdf = pd.DataFrame(prow)
pdf["fisher_p_incoherent_enriched_in_diff"] = fisher_p
pdf.to_csv(snakemake.output.prox, sep="\t", index=False)

# ----------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False})
fig, ax = plt.subplots(2, 2, figsize=(10, 8))
src = hi if len(hi) >= 200 else main
grid = np.linspace(0.2, 0.8, 601)
a = ax[0, 0]
for col, c, lab in (("vaf38", C38, "hg38"), ("vaf13", C13, "T2T-CHM13")):
    v = src[col].dropna()
    v = v[(v >= 0.2) & (v <= 0.8)]
    if len(v) > 50:
        g_, d_ = binned_kde(v)
        a.plot(g_, d_, color=c, lw=2, label=lab)
a.axvline(0.5, color=MUTED, lw=1, ls="--")
a.set(xlabel="Het-SNP VAF", ylabel="Density",
      title=f"VAF density ({'depth >= %d' % HIGH if src is hi else 'all sites'}, n={len(src):,})")
a.legend(frameon=False)

a = ax[0, 1]
if len(main):
    a.hexbin(main.vaf38, main.vaf13, gridsize=60, cmap="Blues", mincnt=1, extent=(0, 1, 0, 1))
a.plot([0, 1], [0, 1], color=MUTED, lw=1)
a.set(xlabel="VAF hg38", ylabel="VAF T2T-CHM13", title="Per-site VAF, paired", xlim=(0, 1), ylim=(0, 1))

a = ax[1, 0]
bins_d = np.linspace(-0.1, 0.1, 81)
nz = main[main.delta_vaf != 0]
for g, c, lab in ((nz[nz.prox_diff == 0], MUTED, f"no ref difference within ±window (n={(nz.prox_diff == 0).sum()})"),
                  (nz[nz.prox_diff == 1], C13, f"ref difference within ±window (n={(nz.prox_diff == 1).sum()})")):
    if len(g):
        a.hist(g.delta_vaf.clip(-0.1, 0.1), bins=bins_d, density=True, histtype="step", lw=2, color=c, label=lab)
a.axvline(0, color=MUTED, lw=1, ls="--")
a.set(xlabel="ΔVAF (CHM13 − hg38)", ylabel="Density",
      title=f"Non-zero ΔVAF by proximal ref difference ({len(main) - len(nz):,} sites at 0)")
a.legend(frameon=False, fontsize=8)

a = ax[1, 1]
bdf = pd.DataFrame(brows)
if len(bdf):
    x = np.arange(len(bdf))
    a.plot(x, bdf.peak_hg38 - 0.5, "o-", color=C38, lw=2, ms=6, label="hg38")
    a.plot(x, bdf.peak_chm13 - 0.5, "o-", color=C13, lw=2, ms=6, label="T2T-CHM13")
    labs = []
    for lo_, hi_ in zip(BINS[:-1], BINS[1:]):
        labs.append(f"{int(lo_)}+" if np.isinf(hi_) else f"{int(lo_)}–{int(hi_) - 1}")
    used = {str(iv): l for iv, l in zip(pd.IntervalIndex.from_breaks(BINS, closed="left"), labs)}
    a.set_xticks(x, [used.get(b, b) for b in bdf.cov_bin])
a.axhline(0, color=MUTED, lw=1, ls="--")
a.set(xlabel="Min depth bin", ylabel="VAF peak − 0.5", title="RMB peak shift by depth")
a.legend(frameon=False)
for a in ax.flat:
    a.grid(color=GRID, lw=0.6)
    a.set_axisbelow(True)
fig.suptitle(f"{summary['sample']}: hg38 vs T2T-CHM13 reference mapping bias", color=INK)
fig.tight_layout()
fig.savefig(snakemake.output.png, dpi=150)
