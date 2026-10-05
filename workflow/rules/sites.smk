SITESETS = list(SAMPLES.index) + (["clinvar"] if RUN_CLINVAR else [])


wildcard_constraints:
    siteset="|".join(map(str, SITESETS)),


# ---------------------------------------------------------------- het calls (hg38)
if config["het_sites"]["source"] == "provided":

    rule het_calls:
        input:
            lambda wc: _val(wc.sample, "het_vcf"),
        output:
            f"{OUT}/sites/{{sample}}.hg38.calls.vcf.gz",
        shell:
            "ln -sf $(readlink -f {input}) {output}"

elif config["het_sites"]["source"] == "clair3":

    rule het_calls:
        input:
            bam=lambda wc: aln(wc.sample, "hg38"),
            fa=f"{OUT}/ref/hg38.fa",
            fai=f"{OUT}/ref/hg38.fa.fai",
        output:
            f"{OUT}/sites/{{sample}}.hg38.calls.vcf.gz",
        params:
            model=config["het_sites"]["clair3_model"],
            platform=config["het_sites"]["clair3_platform"],
            outdir=lambda wc: f"{OUT}/sites/clair3_{wc.sample}",
        threads: config["threads"]["align"]
        log:
            f"{OUT}/logs/clair3.{{sample}}.log",
        conda:
            "../envs/clair3.yaml"
        shell:
            r"""
            run_clair3.sh --bam_fn={input.bam} --ref_fn={input.fa} \
              --threads={threads} --platform={params.platform} \
              --model_path={params.model} --output={params.outdir} > {log} 2>&1
            ln -sf $(readlink -f {params.outdir}/merge_output.vcf.gz) {output}
            """

elif config["het_sites"]["source"] == "bcftools":

    rule het_calls:
        input:
            bam=lambda wc: aln(wc.sample, "hg38"),
            fa=f"{OUT}/ref/hg38.fa",
            fai=f"{OUT}/ref/hg38.fa.fai",
        output:
            f"{OUT}/sites/{{sample}}.hg38.calls.vcf.gz",
        params:
            mq=config["pileup"]["min_mapq"],
            bq=config["pileup"]["min_baseq"],
        log:
            f"{OUT}/logs/bcftools_call.{{sample}}.log",
        conda:
            "../envs/align.yaml"
        shell:
            r"""
            ( bcftools mpileup -B -q {params.mq} -Q {params.bq} -I -a AD,DP -d 10000 \
                -f {input.fa} {input.bam} -Ou \
              | bcftools call -mv -Oz -o {output}
              bcftools index -t {output}
            ) 2> {log}
            """

else:
    raise ValueError("het_sites.source must be provided | clair3 | bcftools")


# ---------------------------------------------------------------- candidates (hg38)
rule prep_het_sites:
    """Biallelic het SNVs passing QUAL/DP/VAF filters; ID := hg38 key for joining."""
    input:
        vcf=f"{OUT}/sites/{{sample}}.hg38.calls.vcf.gz",
        fai=f"{OUT}/ref/hg38.fa.fai",
    output:
        vcf=f"{OUT}/sites/{{sample}}.hg38.candidates.vcf",
        counts=f"{OUT}/sites/{{sample}}.hg38.candidates.counts.tsv",
    params:
        mode="het",
        cfg=config["het_sites"],
        exclude_bed=excl("hg38"),
    log:
        f"{OUT}/logs/prep_sites.{{sample}}.log",
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/prep_sites.py"


# ---------------------------------------------------------------- liftover hg38 -> chm13
rule liftover:
    input:
        vcf=f"{OUT}/sites/{{siteset}}.hg38.candidates.vcf",
        fa=f"{OUT}/ref/chm13.fa",
        dict=f"{OUT}/ref/chm13.dict",
        chain=config["chain_hg38_to_chm13"],
    output:
        vcf=f"{OUT}/sites/{{siteset}}.chm13.lifted.vcf",
        rej=f"{OUT}/sites/{{siteset}}.chm13.rejected.vcf",
    params:
        picard=config["picard"],
    log:
        f"{OUT}/logs/liftover.{{siteset}}.log",
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        {params.picard} LiftoverVcf I={input.vcf} O={output.vcf} REJECT={output.rej} \
          CHAIN={input.chain} R={input.fa} RECOVER_SWAPPED_REF_ALT=true \
          WARN_ON_MISSING_CONTIG=true CREATE_INDEX=false > {log} 2>&1
        """


# ---------------------------------------------------------------- paired site table
rule pair_sites:
    """Join hg38 and lifted sites; flag allele swaps; homopolymer + exclusion filters
    on both references; +/- window sequence comparison between references."""
    input:
        hg38=f"{OUT}/sites/{{siteset}}.hg38.candidates.vcf",
        chm13=f"{OUT}/sites/{{siteset}}.chm13.lifted.vcf",
        rej=f"{OUT}/sites/{{siteset}}.chm13.rejected.vcf",
        prep_counts=f"{OUT}/sites/{{siteset}}.hg38.candidates.counts.tsv",
        fa_hg38=f"{OUT}/ref/hg38.fa",
        fa_chm13=f"{OUT}/ref/chm13.fa",
        fai_hg38=f"{OUT}/ref/hg38.fa.fai",
        fai_chm13=f"{OUT}/ref/chm13.fa.fai",
    output:
        pairs=f"{OUT}/sites/{{siteset}}.pairs.tsv",
        tgt_hg38=f"{OUT}/sites/{{siteset}}.targets.hg38.tsv",
        tgt_chm13=f"{OUT}/sites/{{siteset}}.targets.chm13.tsv",
        log_tsv=f"{OUT}/sites/{{siteset}}.site_filter_counts.tsv",
    params:
        window=config["proximal_window"],
        hp_len=config["het_sites"]["homopolymer_len"],
        hp_flank=config["het_sites"]["homopolymer_flank"],
        exclude_hg38=excl("hg38"),
        exclude_chm13=excl("chm13"),
    log:
        f"{OUT}/logs/pair_sites.{{siteset}}.log",
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/pair_sites.py"
