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


