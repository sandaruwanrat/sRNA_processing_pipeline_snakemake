# smallRNA_seq_proccessing_pipeline

## This pipeline process raw smallRNA raw reads and produce bam files


```text
Raw FASTQ files
       ↓
Raw reads QC (FastQC/MultiQC)
       ↓
Adapter & quality trimming (cutadapt/fastp)
       ↓
Bowtie mapping (genome/piRNA clusters)
       ↓
Size selection (15–40 nt)
       ↓
Remove rRNA/tRNA/snRNA/... (sortmerna or bowtie filter)
       ↓
ShortStack
(run locus discovery, annotation, quantification)
       ↓
Final clean BAM files
