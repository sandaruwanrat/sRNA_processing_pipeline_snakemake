# smallRNA_seq_proccessing_pipeline

## This pipeline process raw smallRNA raw reads and produce bam files
Raw FASTQ files
       ↓
Raw reads QC (FastQC/MultiQC)
       ↓
Adapter & quality trimming (fastp)
       ↓
Bowtie mapping (map to genome+ all the transgeens used in the Lab)
       ↓
Size selection (18–30 nt)
       ↓
Remove rRNA/tRNA
       ↓
ShortStack
       ↓
Final BAM files
