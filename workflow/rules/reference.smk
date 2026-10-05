rule ref_index:
    """faidx + Picard dict (LiftoverVcf needs the target .dict)."""
    input:
        fa=lambda wc: ref_fa(wc.ref),
    output:
        fai=f"{OUT}/ref/{{ref}}.fa.fai",
        dict=f"{OUT}/ref/{{ref}}.dict",
        fa=f"{OUT}/ref/{{ref}}.fa",
    params:
        picard=config["picard"],
    threads: config["threads"]["default"]
    log:
        f"{OUT}/logs/ref_index.{{ref}}.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        ( ln -sf $(readlink -f {input.fa}) {output.fa}
          samtools faidx {output.fa}
          {params.picard} CreateSequenceDictionary R={output.fa} O={output.dict}
        ) > {log} 2>&1
        """


rule ref_mmi:
    input:
        fa=f"{OUT}/ref/{{ref}}.fa",
    output:
        mmi=f"{OUT}/ref/{{ref}}.mmi",
    params:
        preset=config["alignment"]["preset"],
    threads: config["threads"]["default"]
    log:
        f"{OUT}/logs/ref_mmi.{{ref}}.log",
    conda:
        "../envs/align.yaml"
    shell:
        "minimap2 -x {params.preset} -t {threads} -d {output.mmi} {input.fa} > {log} 2>&1"
