configfile: "config/config.yml"
threads_small=config["snakemake"]["threads_small"]
mem_small=config["snakemake"]["mem_small"]
threads_big=config["snakemake"]["threads_big"]
mem_big=config["snakemake"]["mem_big"]

wildcard_constraints:
    chunk="\d+",

rule all:
    input:
        expand("results/{sample}/bowtie/{sample}.fastq.gz", sample=config["samples"]),
        #expand("results/{sample}/step_02_qc_reads_with_vsearch/{sample}_trimmed_qcd.fasta", sample=config["samples"]),
        expand("results/{sample}/step_07_combine_tsv/{sample}_trimmed_qcd_frags_interpro_combined.tsv", sample=config["samples"]),
        expand("results/{sample}/bracken/{sample}_profiles.txt", sample=config["samples"]),
        "results/killed_interproscan.txt",

bowtie_cmd = 'bowtie2 -x {params.idx} -U {input} --un-gz {output} -p {threads}'
bowtie_output = "results/{base}/bowtie/{base}.fastq.gz"

rule bowtie_fastq:
    input:
        "data/{base}.fastq"
    output:
        bowtie_output
    params:
        idx="data/bowtie_index/human+phiX",
        #image="singularity/bowtie.simg /bowtie2/bowtie2-2.4.2/bowtie2",
    threads: threads_small
    resources:
        mem_mb=mem_small
    shell:
        bowtie_cmd


rule bowtie_fastq_gz:
    input:
        "data/{base}.fastq.gz"
    output:
        bowtie_output
    params:
        idx="data/bowtie_index/human+phiX",
        #image="singularity/bowtie.simg /bowtie2/bowtie2-2.4.2/bowtie2",
    threads: threads_small
    resources:
        mem_mb=mem_small
    shell:
        bowtie_cmd


rule bowtie_fq:
    input:
        "data/{base}.fq"
    output:
        bowtie_output
    params:
        idx="data/bowtie_index/human+phiX",
        #image="singularity/bowtie.simg /bowtie2/bowtie2-2.4.2/bowtie2",
    threads: threads_small
    resources:
        mem_mb=mem_small
    shell:
        bowtie_cmd


rule bowtie_fq_gz:
    input:
        "data/{base}.fq.gz"
    output:
        bowtie_output
    params:
        idx="data/bowtie_index/human+phiX",
        #image="singularity/bowtie.simg /bowtie2/bowtie2-2.4.2/bowtie2",
    threads: threads_small
    resources:
        mem_mb=mem_small
    shell:
        bowtie_cmd


rule start_server:
    output:
        "results/interproscan.txt"
    threads: 1
    resources:
        mem_mb=4000
    shell:
        """
        set -u
        sh bash/edit_interproscan.sh {output}
        """


rule qc_pipeline:
    input:
        #fa="results/{base}.{ext}/bowtie/{base}",
        fa="results/{base}/bowtie/{base}.fastq.gz",
    output:
        #fa="results/{base}.{ext}/step_02_qc_reads_with_vsearch/{base}_trimmed_qcd.fasta",
        #log="results/{base}.{ext}/step_02_qc_reads_with_vsearch/log",
        fa="results/{base}/step_02_qc_reads_with_vsearch/{base}_trimmed_qcd.fasta",
        log="results/{base}/step_02_qc_reads_with_vsearch/log"
    params:
        config="config/config.yml",
        #outdir="results/{base}.{ext}/",
        outdir="results/{base}/"
    threads: threads_small
    resources:
        mem_mb=mem_small
    shell:
        """
        set -u
        python pipeline/qc_pipeline.py -c {params.config} -i {input.fa} -o {params.outdir} -t {threads}
        """

rule check_qc:
    input:
        #"results/{base}.{ext}/step_02_qc_reads_with_vsearch/log"
        "results/{base}/step_02_qc_reads_with_vsearch/log"
    output:
        #"results/{base}.{ext}/step_02_qc_reads_with_vsearch/logcheck"
        "results/{base}/step_02_qc_reads_with_vsearch/logcheck"
    threads: 1
    resources:
        mem_mb=4000
    shell:
        """
        set -u
        bash bash/check_qc.sh {input} {output}
        """

rule kraken2:
    input:
        #"results/{base}.{ext}/step_02_qc_reads_with_vsearch/{base}_trimmed_qcd.fasta"
        "results/{base}/step_02_qc_reads_with_vsearch/{base}_trimmed_qcd.fasta"
    output:
        #cf="results/{base}.{ext}/kraken2/{base}_classified.fasta",
        #rep="results/{base}.{ext}/kraken2/{base}_report.tsv"
        cf="results/{base}/kraken2/{base}_classified.fasta",
        rep="results/{base}/kraken2/{base}_report.tsv"

    params:
        db="tools/k2_pluspf_20230314"
    threads: threads_big
    resources:
        mem_mb=mem_big
    shell:
        """
        bash -c '
            . $HOME/.bashrc
            conda activate kraken2
            kraken2 --db {params.db} --classified-out {output.cf} --report {output.rep} --threads {threads} {input}
            '
        """
    

rule bracken:
    input:
        #rep="results/{base}.{ext}/kraken2/{base}_report.tsv",
        rep="results/{base}/kraken2/{base}_report.tsv",
    output:
        #profiles="results/{base}.{ext}/bracken/{base}_profiles.txt",
        #rep="results/{base}.{ext}/bracken/{base}_report.txt"
        profiles="results/{base}/bracken/{base}_profiles.txt",
        rep="results/{base}/bracken/{base}_report.txt"

    params:
        length=300,
        db="tools/k2_pluspf_20230314"
    threads: threads_big
    resources:
        mem_mb=mem_big
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
        #fa="results/{base}.{ext}/step_02_qc_reads_with_vsearch/{base}_trimmed_qcd.fasta",
        fa="results/{base}/step_02_qc_reads_with_vsearch/{base}_trimmed_qcd.fasta",
        ips="results/interproscan.txt",
        #logcheck="results/{base}.{ext}/step_02_qc_reads_with_vsearch/logcheck",
        logcheck="results/{base}/step_02_qc_reads_with_vsearch/logcheck",
    output:
        #"results/{base}.{ext}/step_07_combine_tsv/{base}_trimmed_qcd_frags_interpro_combined.tsv"
        "results/{base}/step_07_combine_tsv/{base}_trimmed_qcd_frags_interpro_combined.tsv"
        #"results/{base}.{ext}/step_05_chunk_reads/{base}_trimmed_qcd_frags_{chunk}.faa"
        #"results/{base}.{ext}/step_04_get_gene_reads/{base}_trimmed_qcd_frags.faa"
    params:
        #indir="results/{base}.{ext}/step_02_qc_reads_with_vsearch/",
        indir="results/{base}/step_02_qc_reads_with_vsearch/",
        config="config/config.yml",
        #outdir="results/{base}.{ext}",
        outdir="results/{base}",
    threads: threads_small
    resources:
        mem_mb=mem_small
    shell:
        """
        set -u
        python pipeline/pipeline.py -c {params.config} -i {params.indir} -o {params.outdir} -t {threads}
        """


rule stop_server:
    input:
        expand("results/{sample}/step_07_combine_tsv/{sample}_trimmed_qcd_frags_interpro_combined.tsv", sample=config["samples"])
    output:
        "results/killed_interproscan.txt"
    threads: 1
    resources:
        mem_mb=4000
    shell:
        """
        JOB=$(squeue -u $USER | grep "lookup" | head -1 | xargs | cut -d" " -f1)
        scancel $JOB
        touch {output}
        """
