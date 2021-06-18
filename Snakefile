configfile: "config/config.yml"


rule all:
    input:
        #expand("results/{sample}/{sample}.fastq.gz", sample=config["samples"]),
        #expand("results/{sample}/step_02_qc_reads_with_vsearch/{sample}.ee{ee}minlen{minlen}.fasta", sample=config["samples"], ee=config["pipeline"]["vsearch_filter_maxee"], minlen=config["pipeline"]["vsearch_filter_minlen"]),
        expand("results/{sample}/step_07_combine_tsv/{sample}_trimmed_qcd_frags_interpro_combined.tsv", sample=config["samples"]),


rule bowtie_cleaning:
    input:
        "data/in_small/{base}.fastq.gz"
    output:
        "results/{base}/{base}.fastq.gz"
    params:
        idx="data/bowtie_index/human+phiX",
        image="/groups/bhurwitz/planet-microbe-functional-annotation/singularity/bowtie.simg /bowtie2/bowtie2-2.4.2/bowtie2",
    shell:
        "singularity exec {params.image} -x {params.idx} -U {input} --un-gz {output} -p 4"


rule edit_interproscan:
    output:
        "results/{base}/interproscan.txt"
    shell:
        "sh bash/edit_interproscan.sh {output}"


rule qc_pipeline:
    input:
        fa="results/{base}/{base}.fastq.gz",
    output:
        "results/{base}/step_02_qc_reads_with_vsearch/{base}_trimmed_qcd.fasta"
    params:
        config="config/config.yml",
        outdir="results/{base}/",
    shell:
        """
        python pipeline/qc_pipeline.py -c {params.config} -i {input.fa} -o {params.outdir}
        """


rule check_qc:


rule kraken2:


rule run_pipeline:
    input:
        fa="results/{base}/step_02_qc_reads_with_vsearch/{base}_trimmed_qcd.fasta",
        ips="results/{base}/interproscan.txt",
    output:
        "results/{base}/step_07_combine_tsv/{base}_trimmed_qcd_frags_interpro_combined.tsv"
    params:
        indir="results/{base}/step_02_qc_reads_with_vsearch/",
        config="config/config.yml",
        outdir="results/{base}",
    shell:
        """
        set -u
        python pipeline/pipeline.py -c {params.config} -i {params.indir} -o {params.outdir}
        """
