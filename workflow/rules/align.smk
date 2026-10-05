def reads_cmd(sample):
    r = _val(sample, "reads")
    if r.endswith(".bam") or r.endswith(".ubam"):
        # keep modified-base and all other tags; minimap2 -y re-attaches them
        return f"samtools fastq -T '*' {r}"
    return f"zcat -f {r}"


rule align:
    """Same reads, same parameters, only the reference differs."""
    input:
        mmi=f"{OUT}/ref/{{ref}}.mmi",
        reads=lambda wc: _val(wc.sample, "reads"),
    output:
        bam=f"{OUT}/align/{{sample}}.{{ref}}.bam",
        bai=f"{OUT}/align/{{sample}}.{{ref}}.bam.bai",
    params:
        reads_cmd=lambda wc: reads_cmd(wc.sample),
        preset=config["alignment"]["preset"],
        extra=config["alignment"]["extra"],
        rg=lambda wc: f"@RG\\tID:{wc.sample}\\tSM:{wc.sample}\\tPL:ONT",
        sort_threads=config["alignment"].get("sort_threads", 4),
        sort_mem=config["alignment"].get("sort_mem_per_thread", "2G"),
    threads: config["threads"]["align"]
    log:
        f"{OUT}/logs/align.{{sample}}.{{ref}}.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        ( {params.reads_cmd} \
          | minimap2 -a -x {params.preset} {params.extra} -R '{params.rg}' -t {threads} {input.mmi} - \
          | samtools sort -@ {params.sort_threads} -m {params.sort_mem} -o {output.bam} -
          samtools index {output.bam}
        ) 2> {log}
        """


def aln(sample, ref):
    """Path of the alignment used downstream (realigned or user-provided)."""
    if prealigned(sample):
        return _val(sample, f"bam_{ref}")
    return f"{OUT}/align/{sample}.{ref}.bam"


rule mapping_stats:
    input:
        bam=lambda wc: aln(wc.sample, wc.ref),
        fa=f"{OUT}/ref/{{ref}}.fa",
        fai=f"{OUT}/ref/{{ref}}.fa.fai",
    output:
        f"{OUT}/qc/{{sample}}.{{ref}}.samtools_stats.txt",
    threads: config["threads"]["default"]
    log:
        f"{OUT}/logs/mapping_stats.{{sample}}.{{ref}}.log",
    conda:
        "../envs/align.yaml"
    shell:
        "samtools stats -@ {threads} --reference {input.fa} {input.bam} > {output} 2> {log}"


rule mapping_table:
    input:
        expand(f"{OUT}/qc/{{sample}}.{{ref}}.samtools_stats.txt", sample=SAMPLES.index, ref=REFS),
    output:
        tsv=f"{OUT}/report/mapping_metrics.tsv",
        tests=f"{OUT}/report/mapping_metrics_tests.tsv",
    log:
        f"{OUT}/logs/mapping_table.log",
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/mapping_table.py"
