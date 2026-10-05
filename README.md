# t2t-ont-refbias

Does aligning ONT reads to T2T-CHM13v2.0 instead of hg38 reduce reference mapping bias (RMB)
and shift VAFs at clinically relevant sites? This Snakemake workflow runs the comparison from
Cherchi et al., *iScience* 2026 ([10.1016/j.isci.2026.117283](https://doi.org/10.1016/j.isci.2026.117283))
on long-read data. The paper did the allele-level analysis only on short reads.

## What it does

1. **Align the same reads twice** (minimap2, identical parameters; MM/ML tags kept) or take
   pre-aligned BAM/CRAM pairs (e.g. the 1kGP ONT Vienna CRAMs, which exist for both builds).
2. **Mapping QC**: `samtools stats` per alignment, paired across samples.
3. **Shared het-SNP site set**: per-sample germline het SNVs on hg38 (Clair3, a provided VCF,
   or bcftools for testing), QUAL/DP/VAF filtered, lifted to CHM13 with Picard
   `LiftoverVcf RECOVER_SWAPPED_REF_ALT=true`. Many-to-one lifts are dropped, allele swaps
   flagged, homopolymer-adjacent SNPs (run ≥5 within ±5 bp, ref or alt haplotype, on **both**
   references) and exclusion BEDs (SDs, non-syntenic regions) removed.
4. **Proximal reference comparison**: ±150 bp around every site, hg38 vs CHM13, global alignment
   (edlib). Mismatches, insertions, deletions, and whether a difference falls in a homopolymer.
5. **Allele counting** at the identical sites in both alignments: `bcftools mpileup -B -I`,
   MAPQ/BQ filters, secondary/supplementary/dup removed.
6. **Statistics** (per sample and cohort), following the paper:
   - VAF density peak (closest local max to 0.5), all sites and depth ≥ `high_cov`;
     paired site bootstrap for "CHM13 peak closer to 0.5"
   - ΔVAF = VAF_CHM13 − VAF_hg38; "incoherent" if |ΔVAF| ≥ 0.005; one-sided Wilcoxon
     CHM13 > hg38 at incoherent sites
   - coverage-binned peaks and tests
   - Fisher test: are incoherent sites enriched where the references differ nearby?
   - allele-swapped sites summarised separately
7. **Optional ClinVar module**: P/LP SNVs lifted the same way; VAF/depth in each sample's two
   alignments where the variant has ≥ `min_alt_reads`, stratified by proximal difference.
   Point it at tumor BAMs to replicate the paper's cfDNA analysis on long reads.

Outputs land in `results/report/` (`report.md`, `cohort_summary.tsv`, `cohort_pooled_tests.tsv`,
`mapping_metrics*.tsv`, `cohort_vaf.png`, `clinvar_*.tsv`) and `results/stats/<sample>.*`
(summary, coverage bins, proximal-difference table, per-site table, figure).

## Setup on the cluster

```bash
git clone <this repo> && cd t2t-ont-refbias
mamba create -n snk -c conda-forge -c bioconda snakemake>=8 snakemake-executor-plugin-slurm
conda activate snk
```

References (hg38 via your existing analysis set):

```bash
mkdir -p refs && cd refs
wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz && gunzip hs1.fa.gz   # = chm13v2.0
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz
# Exclusion BEDs: CHM13 SDs + non-syntenic regions from https://github.com/marbl/CHM13
# (chm13v2.0_SD.bed, hg38 non-syntenic BED). Concatenate per build.
```

Check that contig names match across the FASTAs, chain and BAMs (`chr1` vs `1`). Then edit
`config/config.yaml` and `config/samples.tsv`. Two input modes per sample:

| Mode | Columns | Use for |
|---|---|---|
| A: realign | `reads` = Dorado uBAM or FASTQ | your own data; guarantees the reference is the only variable |
| B: pre-aligned | `bam_hg38`, `bam_chm13` | 1kGP ONT CRAMs (confirm both builds used the same reads and aligner) |

`het_vcf` (hg38 coordinates) is required when `het_sites.source: provided`.

## Run

```bash
snakemake -n                                         # dry run
snakemake --profile profiles/slurm                   # set partition/account first
snakemake -c 32 --sw-deployment-method conda         # single node
```

## Test (runs in ~2 min on 2 cores)

```bash
python test/make_test_data.py test/data      # 0.5 Mb toy genome, 45x and 70x simulated reads
snakemake --configfile test/config.yaml -c 2 --sw-deployment-method conda
```

The simulator gives CHM13 known differences from hg38 (clustered substitutions near 25% of het
SNPs, indels that shift coordinates, allele-swap sites), so the test checks liftover, swap
detection (17/17 recovered, 0 false), proximal-difference flagging (134/134 truth sites
flagged), and that incoherent sites are enriched near reference differences (Fisher p ≈ 7e-5).

## Things to decide before the real run

- **Depth.** At 30× the paper saw no RMB peak shift even in short reads. Use deep WGS, pool
  samples, or adaptive-sampling targets; read the coverage-binned tables.
- **Preset.** `lr:hq` needs minimap2 ≥ 2.27 and suits R10.4.1 SUP; use `map-ont` for older data.
- **Het calls on hg38 only.** Sites are defined on hg38 and lifted, which is the paper's design
  but cannot see CHM13-only sequence. A symmetric variant (call on both, intersect) is a small
  change in `rules/sites.smk` if you want it.
- **Accuracy check (not automated).** Lower bias doesn't help if precision drops: run HG002
  through both builds and score with hap.py or `rtg vcfeval` against GIAB truth on each build.

## Layout

```
config/          config.yaml, samples.tsv
workflow/        Snakefile, rules/*.smk, scripts/*.py, envs/*.yaml
profiles/slurm/  Slurm executor profile
test/            simulator + test config
```
