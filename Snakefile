configfile: "config/config.yml"


rule all:
    input:
        #expand("results/{sample}/{sample}.fastq.gz", sample=config["samples"]),
        #expand("results/{sample}/step_02_qc_reads_with_vsearch/{sample}.ee{ee}minlen{minlen}.fasta", sample=config["samples"], ee=config["pipeline"]["vsearch_filter_maxee"], minlen=config["pipeline"]["vsearch_filter_minlen"]),
        expand("results/{sample}/step_07_combine_tsv/{sample}_trimmed_qcd_frags_interpro_combined.tsv", sample=config["samples"]),
        expand("results/{sample}/bracken/{sample}_profiles.txt", sample=config["samples"]),
        "results/killed_interproscan.txt",


rule bowtie_cleaning:
    input:
        "data/{base}.fastq.gz"
    output:
        "results/{base}/bowtie/{base}.fastq.gz"
    params:
        idx="data/bowtie_index/human+phiX",
        #image="/groups/bhurwitz/planet-microbe-functional-annotation/singularity/bowtie.simg /bowtie2/bowtie2-2.4.2/bowtie2",
        image="singularity/bowtie.simg /bowtie2/bowtie2-2.4.2/bowtie2",
    shell:
        "singularity exec {params.image} -x {params.idx} -U {input} --un-gz {output} -p 4"


rule start_server:
    output:
        "results/interproscan.txt"
    shell:
        """
        set -u
        sh bash/edit_interproscan.sh {output}
        """


rule qc_pipeline:
    input:
        fa="results/{base}/bowtie/{base}.fastq.gz",
    output:
        fa="results/{base}/step_02_qc_reads_with_vsearch/{base}_trimmed_qcd.fasta",
        log="results/{base}/step_02_qc_reads_with_vsearch/log"
    params:
        config="config/config.yml",
        outdir="results/{base}/",
    shell:
        """
        set -u
        python pipeline/qc_pipeline.py -c {params.config} -i {input.fa} -o {params.outdir}
        """


rule check_qc:
    input:
        "results/{base}/step_02_qc_reads_with_vsearch/log"
    output:
        "results/{base}/step_02_qc_reads_with_vsearch/logcheck"
    shell:
        """
        set -u
        bash bash/check_qc.sh {input} {output}
        """

rule kraken2:
    input:
        "results/{base}/step_02_qc_reads_with_vsearch/{base}_trimmed_qcd.fasta"
    output:
        cf="results/{base}/kraken2/{base}_classified.fasta",
        rep="results/{base}/kraken2/{base}_report.tsv"
    params:
        threads={config["pipeline"]["threads"]},
        db="/xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/alise/my_scripts/readbased_metagenomes_snakemake/PBS_pipeline/databases/k2_pluspf_20210127"
    shell:
        """
        bash -c '
            . $HOME/.bashrc
            conda activate kraken2
            kraken2 --db {params.db} --classified-out {output.cf} --report {output.rep} --threads {params.threads} {input}
            '
        """
    

rule bracken:
    input:
        rep="results/{base}/kraken2/{base}_report.tsv",
    output:
        profiles="results/{base}/bracken/{base}_profiles.txt",
        rep="results/{base}/bracken/{base}_report.txt"
    params:
        length=300,
        db="/xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/alise/my_scripts/readbased_metagenomes_snakemake/PBS_pipeline/databases/k2_pluspf_20210127"
    shell:
        """
        bash -c '
            . $HOME/.bashrc
            conda activate braken
            bracken -d {params.db} -i {input.rep} -o {output.profiles} -w {output.rep} -r {params.length} -l S
            '
        """


rule run_pipeline:
    input:
        fa="results/{base}/step_02_qc_reads_with_vsearch/{base}_trimmed_qcd.fasta",
        ips="results/interproscan.txt",
        logcheck="results/{base}/step_02_qc_reads_with_vsearch/logcheck",
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

rule stop_server:
    input:
        expand("results/{sample}/step_07_combine_tsv/{sample}_trimmed_qcd_frags_interpro_combined.tsv", sample=config["samples"])
    output:
        "results/killed_interproscan.txt"
    shell:
        """
        JOB=$(squeue -u $USER | grep "lookup" | head -1 | xargs | cut -d" " -f1)
        scancel $JOB
        touch {output}
        """
