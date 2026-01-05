# smallRNA_seq_proccessing_pipeline

## This pipeline process raw smallRNA raw reads and produce bam files


flowchart TD
    A[Raw FASTQ files] --> B[Raw reads QC<br>FastQC / MultiQC]
    B --> C[Adapter & quality trimming<br>cutadapt / fastp]
    C --> D[Bowtie mapping<br>genome / piRNA clusters]
    D --> E[Size selection<br>15–40 nt]
    E --> F[Remove rRNA/tRNA/snRNA/...<br>sortmerna or bowtie filter]
    F --> G[ShortStack<br>locus discovery, annotation, quantification]
    G --> H[Final clean BAM files]

    style A fill:#f9f9f9,stroke:#333
    style H fill:#d1e7ff,stroke:#333
