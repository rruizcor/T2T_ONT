"""Cohort-level summary, pooled tests, figure and Markdown report."""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.stderr = open(snakemake.log[0], "w")
P = snakemake.params.stats
MIN_COV, HIGH = P["min_cov"], P["high_cov"]
C38, C13 = "#2a78d6", "#eb6834"
INK, MUTED, GRID = "#1f1f1e", "#6b6a63", "#e6e5df"


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


def bootstrap_peak_p(v38, v13, n=1000, seed=7):
    """P(|peak13-0.5| >= |peak38-0.5|) under paired bootstrap of sites (one-sided:
    small p = CHM13 peak closer to 0.5). Mirrors the paper's bootstrap on 1kGP VCF VAFs."""
    rng = np.random.default_rng(seed)
    v38, v13 = np.asarray(v38), np.asarray(v13)
    if len(v38) < 200:
        return np.nan
    m = len(v38)
    worse = 0
    for _ in range(n):
        i = rng.integers(0, len(v38), m)
        if abs(kde_peak(v13[i]) - 0.5) >= abs(kde_peak(v38[i]) - 0.5):
            worse += 1
    return (worse + 1) / (n + 1)


def fmt(x, nd=4):
    if isinstance(x, (int, np.integer)):
        return f"{x:,}"
    if isinstance(x, (float, np.floating)):
        if np.isnan(x):
            return "NA"
        if x != 0 and (abs(x) < 1e-3 or abs(x) >= 1e5):
            return f"{x:.2e}"
        return f"{x:.{nd}f}"
    return str(x)


def md_table(df, cols=None):
    cols = cols or list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in df[cols].itertuples(index=False):
        out.append("| " + " | ".join(fmt(x) for x in r) + " |")
    return "\n".join(out)


summ = pd.concat([pd.read_csv(p, sep="\t") for p in snakemake.input.summaries], ignore_index=True)
summ.to_csv(snakemake.output.summary, sep="\t", index=False)

frames = []
for p in snakemake.input.sites:
    d = pd.read_csv(p, sep="\t", usecols=["vaf38", "vaf13", "min_dp", "swapped", "incoherent",
                                          "delta_vaf", "prox_diff"])
    d["sample"] = Path(p).name.replace(".sites.tsv.gz", "")
    frames.append(d)
allsites = pd.concat(frames, ignore_index=True)
main = allsites[(allsites.min_dp >= MIN_COV) & (allsites.swapped == 0)].dropna(subset=["vaf38", "vaf13"])
hi = main[main.min_dp >= HIGH]
inc = main[main.incoherent == 1]


def wil(x, y, alt="greater"):
    d = np.asarray(x) - np.asarray(y)
    d = d[d != 0]
    return float(stats.wilcoxon(d, alternative=alt).pvalue) if len(d) >= 10 else np.nan


pooled = [
    dict(test="pooled sites, depth>=min_cov", n=len(main), hg38=kde_peak(main.vaf38),
         chm13=kde_peak(main.vaf13), stat="VAF density peak", p=np.nan),
    dict(test=f"pooled sites, depth>={HIGH}", n=len(hi), hg38=kde_peak(hi.vaf38),
         chm13=kde_peak(hi.vaf13), stat="VAF density peak",
         p=bootstrap_peak_p(hi.vaf38.values, hi.vaf13.values, n=1000)),
    dict(test="incoherent sites (|ΔVAF|>=thr)", n=len(inc), hg38=inc.vaf38.mean(), chm13=inc.vaf13.mean(),
         stat="mean VAF; one-sided Wilcoxon CHM13>hg38", p=wil(inc.vaf13, inc.vaf38)),
    dict(test="fraction incoherent", n=len(main), hg38=np.nan, chm13=np.nan,
         stat=f"{main.incoherent.mean():.4f}" if len(main) else "NA", p=np.nan),
]
# sample-level: is the CHM13 peak closer to 0.5 than hg38?
s = summ.dropna(subset=["peak_hg38_all", "peak_chm13_all"])
dist38 = (s.peak_hg38_all - 0.5).abs()
dist13 = (s.peak_chm13_all - 0.5).abs()
pooled.append(dict(
    test="per-sample |peak-0.5| (all sites)", n=len(s), hg38=dist38.mean(), chm13=dist13.mean(),
    stat="mean distance; one-sided paired Wilcoxon hg38>CHM13",
    p=float(stats.wilcoxon(dist38, dist13, alternative="greater").pvalue)
    if len(s) >= 6 and (dist38 != dist13).any() else np.nan))
pooled = pd.DataFrame(pooled)
pooled.to_csv(snakemake.output.pooled, sep="\t", index=False)

# ----------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                     "axes.spines.right": False})
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
src = hi if len(hi) >= 200 else main
grid = np.linspace(0.2, 0.8, 601)
for col, c, lab in (("vaf38", C38, "hg38"), ("vaf13", C13, "T2T-CHM13")):
    v = src[col]
    v = v[(v >= 0.2) & (v <= 0.8)].values
    if len(v) > 200000:
        v = np.random.default_rng(1).choice(v, 200000, replace=False)
    if len(v) > 50:
        g_, d_ = binned_kde(v)
        ax[0].plot(g_, d_, color=c, lw=2, label=lab)
ax[0].axvline(0.5, color=MUTED, lw=1, ls="--")
ax[0].set(xlabel="Het-SNP VAF", ylabel="Density",
          title=f"Pooled VAF density ({'depth >= %d' % HIGH if src is hi else 'all'}; n={len(src):,})")
ax[0].legend(frameon=False)
if len(s):
    y = np.arange(len(s))
    for i, r in enumerate(s.itertuples()):
        ax[1].plot([r.peak_hg38_all, r.peak_chm13_all], [i, i], color=GRID, lw=2, zorder=1)
    ax[1].scatter(s.peak_hg38_all, y, s=40, color=C38, label="hg38", zorder=2,
                  edgecolor="white", linewidth=1.5)
    ax[1].scatter(s.peak_chm13_all, y, s=40, color=C13, label="T2T-CHM13", zorder=3,
                  edgecolor="white", linewidth=1.5)
    ax[1].set_yticks(y, s["sample"] if len(s) <= 40 else [""] * len(s))
ax[1].axvline(0.5, color=MUTED, lw=1, ls="--")
ax[1].set(xlabel="VAF density peak", title="Per-sample peak (closer to 0.5 = less bias)")
ax[1].legend(frameon=False)
for a in ax:
    a.grid(color=GRID, lw=0.6)
    a.set_axisbelow(True)
fig.tight_layout()
fig.savefig(snakemake.output.png, dpi=150)

# ----------------------------------------------------------------- markdown
filt = []
for p in snakemake.input.filt:
    d = pd.read_csv(p, sep="\t")
    d["sample"] = Path(p).name.replace(".site_filter_counts.tsv", "")
    filt.append(d)
filt = pd.concat(filt).pivot_table(index="stage", columns="sample", values="n", aggfunc="sum").fillna(0).astype(int)
mapping = pd.read_csv(snakemake.input.mapping, sep="\t")
cfg = snakemake.params.config

lines = [
    "# hg38 vs T2T-CHM13 reference mapping bias (ONT)",
    "",
    f"Samples: {len(summ)} · min depth {MIN_COV} in both alignments · high-depth threshold {HIGH} · "
    f"ΔVAF threshold {P['delta_vaf']} · het-site source: `{cfg['het_sites']['source']}` · "
    f"preset `{cfg['alignment']['preset']}`",
    "",
    "## Mapping metrics (paired across samples)",
    "",
    md_table(mapping),
    "",
    "## Reference mapping bias: pooled and sample-level tests",
    "",
    md_table(pooled),
    "",
    "Peak p-value: paired bootstrap of sites; small p means the CHM13 VAF peak sits closer to 0.5 "
    "than the hg38 peak. Allele-swapped sites are excluded from these tests and summarised per sample.",
    "",
    "## Per-sample summary",
    "",
    md_table(summ, ["sample", "n_main", "n_high_cov", "peak_hg38_highcov", "peak_chm13_highcov",
                    "peak_hg38_all", "peak_chm13_all", "frac_incoherent", "n_inc_chm13_higher",
                    "n_inc_hg38_higher", "p_wilcoxon_incoherent_chm13_gt_hg38", "n_swapped"]),
    "",
    "Per-sample coverage-binned and proximal-difference tables: `stats/<sample>.coverage_bins.tsv`, "
    "`stats/<sample>.proximal_diff.tsv`; figures `stats/<sample>.vaf.png`, `report/cohort_vaf.png`.",
    "",
    "## Site filtering",
    "",
    md_table(filt.reset_index()),
    "",
]
if snakemake.input.clinvar:
    cv = pd.read_csv(snakemake.input.clinvar[0], sep="\t")
    lines += ["## ClinVar P/LP sites with alt support", "", md_table(cv) if len(cv) else "No sites passed.", ""]
lines += [
    "## Caveats",
    "",
    "- At ~30× the paper could not resolve RMB peak shifts even in short reads; rely on the "
    f"depth ≥ {HIGH} and coverage-binned results.",
    "- Sites near homopolymers are excluded on both references, which removes some of the regions "
    "where the paper found hg38/CHM13 differences; the proximal-difference table covers the rest.",
    "- Liftover cannot place sites in CHM13-only sequence, so new T2T regions are not assessed.",
    "",
]
Path(snakemake.output.md).write_text("\n".join(lines))
