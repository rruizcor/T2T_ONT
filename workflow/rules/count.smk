rule count_alleles:
    """Per-allele read counts at the shared site set, identical filters on both references.
    -B: no BAQ (tuned for short reads); -I: SNVs only; secondary/supplementary/dup/QC-fail dropped."""
    input:
        bam=lambda wc: aln(wc.sample, wc.ref),
        fa=f"{OUT}/ref/{{ref}}.fa",
        fai=f"{OUT}/ref/{{ref}}.fa.fai",
        targets=f"{OUT}/sites/{{siteset}}.targets.{{ref}}.tsv",
    output:
        f"{OUT}/counts/{{sample}}.{{siteset}}.{{ref}}.tsv",
    params:
        mq=config["pileup"]["min_mapq"],
        bq=config["pileup"]["min_baseq"],
        maxdp=config["pileup"]["max_depth"],
    log:
        f"{OUT}/logs/count.{{sample}}.{{siteset}}.{{ref}}.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        ( if [ -s {input.targets} ]; then
            bcftools mpileup -f {input.fa} -R {input.targets} \
              -B -I -q {params.mq} -Q {params.bq} -d {params.maxdp} \
              --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY \
              -a FORMAT/AD,FORMAT/DP -Ou {input.bam} \
            | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t[%DP]\n' > {output}
          else
            : > {output}
          fi
        ) 2> {log}
        """
