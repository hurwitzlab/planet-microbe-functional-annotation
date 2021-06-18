# planet-microbe-functional-annotation
EBI Functional Annotation Pipeline For Planet Microbe
Matt Miller, Kai Blumberg

Based on the EBI functional pipeline version 4.1 (https://www.ebi.ac.uk/metagenomics/pipelines/4.1)

Pipeline steps: 
1) FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
Quality Control FastQC download 

2) Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)

3) FragGeneScan (https://sourceforge.net/projects/fraggenescan/files/)

ORF prediction 

4) Interproscan (https://github.com/ebi-pf-team/interproscan)
Interproscan on the ORF predictions to get the pfam and interpro to go

Normalize the data after interproscan to get the normalized values or maybe not. 

Flags:
 --goterms to switch on corresponding GO terms to Interpro annotations. 
--applications flag for which analyses we just want to use PFAM

How the pipeline works:
1. Download the sequences (snakemake)
2. Run sh edit_interproscan.sh to start interproscan server
3. Run sh bowtie_cleaning.sh /path/to/input/file /path/to/output/dir (will output unmatched reads to the output dir with the basename of the input file)
4. Run sh run_qc_pipeline.sh /path/to/config /path/to/bowtie/cleaning/output/file /path/to/output/dir
5. Check results of QC (snakemake)
6. Run kraken2 on the qc output (snakemake)
7. Run sh run_pipeline.sh /path/to/config /path/to/step_02_qc_reads_with_vsearch/ /path/to/step_02_qc_reads_with_vsearch/.. (arg 2 is path to the directory not the file in there, can change it if you want)
8. A final combined tsv file with be in step_07_combine_tsv.
