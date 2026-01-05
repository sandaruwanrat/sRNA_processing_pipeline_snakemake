
# This pipeline process raw smallRNA raw reads and produce bam files


```text
Raw FASTQ files
       ↓
Raw reads QC (FastQC & MultiQC)
       ↓
Adapter & quality trimming (fastp)
       ↓
Bowtie mapping (genome + transgenes used in the lab)
       ↓
Size selection (18–30 nt)
       ↓
Remove rRNA/tRNA
       ↓
ShortStack
       ↓
Final BAM files
