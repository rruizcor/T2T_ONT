rule prep_clinvar_sites:
    """ClinVar pathogenic / likely-pathogenic SNVs (hg38). Chromosome names are
    harmonised to the hg38 FASTA (adds 'chr' if needed)."""
    input:
        vcf=config["clinvar"]["vcf"],
        fai=f"{OUT}/ref/hg38.fa.fai",
    output:
        vcf=f"{OUT}/sites/clinvar.hg38.candidates.vcf",
        counts=f"{OUT}/sites/clinvar.hg38.candidates.counts.tsv",
    params:
        mode="clinvar",
        cfg=config["het_sites"],
        exclude_bed=excl("hg38"),
    log:
        f"{OUT}/logs/prep_sites.clinvar.log",
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/prep_sites.py"


rule clinvar_report:
    input:
        pairs=f"{OUT}/sites/clinvar.pairs.tsv",
        counts=expand(
            f"{OUT}/counts/{{sample}}.clinvar.{{ref}}.tsv", sample=SAMPLES.index, ref=REFS
        ),
    output:
        sites=f"{OUT}/report/clinvar_sites.tsv",
        tests=f"{OUT}/report/clinvar_tests.tsv",
    params:
        samples=list(SAMPLES.index),
        outdir=OUT,
        min_alt=config["clinvar"]["min_alt_reads"],
        min_cov=config["stats"]["min_cov"],
    log:
        f"{OUT}/logs/clinvar_report.log",
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/clinvar_report.py"
