rule sample_stats:
    input:
        pairs=f"{OUT}/sites/{{sample}}.pairs.tsv",
        c_hg38=f"{OUT}/counts/{{sample}}.{{sample}}.hg38.tsv",
        c_chm13=f"{OUT}/counts/{{sample}}.{{sample}}.chm13.tsv",
    output:
        sites=f"{OUT}/stats/{{sample}}.sites.tsv.gz",
        summary=f"{OUT}/stats/{{sample}}.summary.tsv",
        bins=f"{OUT}/stats/{{sample}}.coverage_bins.tsv",
        prox=f"{OUT}/stats/{{sample}}.proximal_diff.tsv",
        png=f"{OUT}/stats/{{sample}}.vaf.png",
    params:
        stats=config["stats"],
    log:
        f"{OUT}/logs/sample_stats.{{sample}}.log",
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/sample_stats.py"


rule cohort_report:
    input:
        summaries=expand(f"{OUT}/stats/{{sample}}.summary.tsv", sample=SAMPLES.index),
        sites=expand(f"{OUT}/stats/{{sample}}.sites.tsv.gz", sample=SAMPLES.index),
        filt=expand(f"{OUT}/sites/{{sample}}.site_filter_counts.tsv", sample=SAMPLES.index),
        mapping=f"{OUT}/report/mapping_metrics_tests.tsv",
        mapping_tsv=f"{OUT}/report/mapping_metrics.tsv",
        clinvar=([f"{OUT}/report/clinvar_tests.tsv"] if RUN_CLINVAR else []),
    output:
        summary=f"{OUT}/report/cohort_summary.tsv",
        pooled=f"{OUT}/report/cohort_pooled_tests.tsv",
        png=f"{OUT}/report/cohort_vaf.png",
        md=f"{OUT}/report/report.md",
    params:
        stats=config["stats"],
        config=config,
    log:
        f"{OUT}/logs/cohort_report.log",
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/cohort_report.py"
