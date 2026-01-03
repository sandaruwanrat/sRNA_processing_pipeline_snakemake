import os
import glob
import pandas as pd
import yaml
import sys


#module load fastp/0.23.4
#snakemake/v8.27.1

#load config file
configfile: "sRNA_R29_config.yaml"

#for DAG and rulegraph print need to comment all the print s
#print("Config keys:", config.keys())


old_name=pd.read_csv(config["filenames"],sep='\t')['old_nm']
new_name=pd.read_csv(config["filenames"],sep='\t')['new_nm']

sample_df=pd.read_csv(config["filenames"],sep='\t',header=0)
old_nm=sample_df.old_nm
new_nm=sample_df.new_nm

rule all:
    input:
        expand(f"{config['rename_dir']}/{{sample}}.fq.gz", sample=sample_df.new_nm),
        expand(f"{config['fastQC_bft_dir']}/{{sample}}_fastqc.html", sample=sample_df.new_nm),
        expand(f"{config['fastQC_bft_dir']}/{{sample}}_fastqc.zip", sample=sample_df.new_nm),
        f"{config['multQC_bft_dir']}/multiqc_report.html",
        expand(f"{config['adpt_trimmed_fp_dir']}/{{sample}}_trimmed_unfiltered.fq.gz", sample=sample_df.new_nm),
        expand(f"{config['adpt_trimmed_fp_dir']}/{{sample}}_untrimmed.fq.gz", sample=sample_df.new_nm),
        expand(f"{config['adpt_trimmed_fp_dir']}/{{sample}}_fastp.html", sample=sample_df.new_nm),
        expand(f"{config['adpt_trimmed_fp_dir']}/{{sample}}_trimmed_unfiltered.fa", sample=sample_df.new_nm),
        expand(f"{config['fastQC_AT_dir']}/{{sample}}_trimmed_unfiltered_fastqc.html", sample=sample_df.new_nm),
        expand(f"{config['fastQC_AT_dir']}/{{sample}}_trimmed_unfiltered_fastqc.zip", sample=sample_df.new_nm),
        f"{config['multQC_AT_dir']}/multiqc_report.html",
        expand(f"{config['mapped_unfiltered_dir']}/{{sample}}.unfiltered.mapped.bam", sample=sample_df.new_nm),
        expand(f"{config['size_filt_rmChloMito_dir']}/{{sample}}.unfiltered.mapped.mito_chloroFree.sizeFilt.fa", sample=sample_df.new_nm),
        expand(f"{config['map_rRNA_tRNA_free_dir']}/{{sample}}.mito_chloroFree.rRNA_tRNA_free.fa", sample=sample_df.new_nm),
        expand(f"{config['sht_dir']}/{{sample}}_ShortStack_out", sample=sample_df.new_nm)
        
        
        




def find_input(wildcards):
    prefix = sample_df[sample_df['new_nm'] == wildcards.sample].old_nm.values[0]
    #example '/indir/path_to_dir/SL2525_BR1*'
    pattern = os.path.join(config["input_dir"], f"{prefix}*")
    matches = glob.glob(pattern)
    if not matches:
        raise ValueError(f"No input file found for prefix '{prefix}' with pattern '{pattern}'")
    return matches[0]


#rename and copy files to new dir
rule rename_files:
    input: 
        srna_in = find_input   
    output:
        f"{config['rename_dir']}/{{sample}}.fq.gz"
    shell:"""
        mkdir -p {config[rename_dir]}
        
        cp {input} {output}
        
        """

rule fastQC_bf_trim:
    input:
        f"{config['rename_dir']}/{{sample}}.fq.gz"
    output:
        f"{config['fastQC_bft_dir']}/{{sample}}_fastqc.html",
        f"{config['fastQC_bft_dir']}/{{sample}}_fastqc.zip"
    shell:"""
        mkdir -p {config[fastQC_bft_dir]}
        fastqc --threads 16 -q -o {config[fastQC_bft_dir]} {input} 
        
        """


#here do not use both files
#f"{config['fastQC_bft_dir']}/{{sample}}_fastqc.html",
#f"{config['fastQC_bft_dir']}/{{sample}}_fastqc.zip"
#multiqc cannot use both files since it's using wildcard and no use
rule multiQC_report:
    input:
        expand(f"{config['fastQC_bft_dir']}/{{sample}}_fastqc.html", sample=sample_df.new_nm)
    output:
        f"{config['multQC_bft_dir']}/multiqc_report.html"
    shell:"""
        mkdir -p {config[multQC_bft_dir]}
        multiqc {config[fastQC_bft_dir]} -o {config[multQC_bft_dir]} --force
        """

# rule run_fastp:
#     input:
#         f"{config['rename_dir']}/{{sample}}.fq.gz"
#     output:
#         trimmed_fq= f"{config['adpt_trimmed_dir']}/{{sample}}_trimmed_unfiltered.fq.gz",
#         untrimmed_fq=f"{config['adpt_trimmed_dir']}/{{sample}}_untrimmed.fq.gz",
#         html = f"{config['adpt_trimmed_dir']}/{{sample}}_fastp.html"
#     log:
#         f"{config['adpt_trimmed_dir']}/{{sample}}_fastp.log"
#     threads: 4
#     shell: """
#             mkdir -p {config[adpt_trimmed_dir]}
#             fastp \
#                 -i {input} \
#                 -o {output.trimmed_fq} \
#                 --adapter_sequence TGGAATTCTCGGGTGCCAAGG \
#                 --qualified_quality_phred 33 \
#                 --length_required 5 \
#                 --failed_out {output.untrimmed_fq} \
#                 --html {output.html} \
#                 --thread {threads} \
#                 > {log} 2>&1
#         """


rule run_fastp:
    input:
        f"{config['rename_dir']}/{{sample}}.fq.gz"
    output:
        trimmed_fq= f"{config['adpt_trimmed_fp_dir']}/{{sample}}_trimmed_unfiltered.fq.gz",
        untrimmed_fq=f"{config['adpt_trimmed_fp_dir']}/{{sample}}_untrimmed.fq.gz",
        html = f"{config['adpt_trimmed_fp_dir']}/{{sample}}_fastp.html"
    log:
        f"{config['adpt_trimmed_fp_dir']}/{{sample}}_fastp.log"
    threads: 4
    shell: """
            module load fastp/0.23.4
            mkdir -p {config[adpt_trimmed_fp_dir]}
            fastp \
                -i {input} \
                -o {output.trimmed_fq} \
                --adapter_sequence TGGAATTCTCGGGTGCCAAGG \
                --qualified_quality_phred 33 \
                --length_required 5 \
                --failed_out {output.untrimmed_fq} \
                --html {output.html} \
                --thread {threads} \
                > {log} 2>&1
        """


rule fq_to_fa:
    input:
        f"{config['adpt_trimmed_fp_dir']}/{{sample}}_trimmed_unfiltered.fq.gz"
    output:
        f"{config['adpt_trimmed_fp_dir']}/{{sample}}_trimmed_unfiltered.fa"
    shell:"""
        seqkit fq2fa {input} -o {output}
        """

rule fastQC_aft_trim:
    input:
        f"{config['adpt_trimmed_fp_dir']}/{{sample}}_trimmed_unfiltered.fq.gz"
    output:
        f"{config['fastQC_AT_dir']}/{{sample}}_trimmed_unfiltered_fastqc.html",
        f"{config['fastQC_AT_dir']}/{{sample}}_trimmed_unfiltered_fastqc.zip"
    shell:"""
        mkdir -p {config[fastQC_AT_dir]}
        fastqc --threads 16 -q -o {config[fastQC_AT_dir]} {input} 
        
        """


rule multiQC_AT_report:
    input:
        expand(f"{config['fastQC_AT_dir']}/{{sample}}_trimmed_unfiltered_fastqc.html", sample=sample_df.new_nm)
    output:
        f"{config['multQC_AT_dir']}/multiqc_report.html"
    shell:"""
        mkdir -p {config[multQC_AT_dir]}
        multiqc {config[fastQC_AT_dir]} -o {config[multQC_AT_dir]} --force
        """

#can add a new rule to check the index and make the index if not present
#bowtie-build <reference_in.fa> <index_name>
#bowtie-build TAIR10_Chr-all_Araport11+current_transgenes_R29.fa TAIR10_Chr-all_Araport11+current_transgenes_R29
rule map_unfiltered_to_TAIRtgs:
    input:
        input_fa=f"{config['adpt_trimmed_fp_dir']}/{{sample}}_trimmed_unfiltered.fa"
    output:
        f"{config['mapped_unfiltered_dir']}/{{sample}}.unfiltered.mapped.bam"
    params:
        gindex=f"{config['INDEX_AT11_plus_TGS']}"
    resources:
        mem_mb=32000 #fails with 16GB
    threads: 4
    shell:"""
        mkdir -p {config[mapped_unfiltered_dir]}
        bowtie -v 0 -p {threads} -fS -x {params.gindex} {input.input_fa} |\
         samtools view -bS -F 4 -@ {threads} - |\
          samtools sort -@ {threads} -  -o {output}


        """



rule size_Chlo_Mito_filt:
    input:
        f"{config['mapped_unfiltered_dir']}/{{sample}}.unfiltered.mapped.bam"
    output:
        f"{config['size_filt_rmChloMito_dir']}/{{sample}}.unfiltered.mapped.mito_chloroFree.sizeFilt.fa"
    resources:
        mem_mb=16000
    shell:"""
        mkdir -p {config[size_filt_rmChloMito_dir]}

        samtools view -h {input} |\
         awk '$1 ~ /^@/ || (length($10) >= 18 && length($10) <= 28 && $3 != "ChrM" && $3 != "ChrC")' |\
          samtools fasta - > {output}

        """


rule remove_tRNA_rRNA:
    input:
        f"{config['size_filt_rmChloMito_dir']}/{{sample}}.unfiltered.mapped.mito_chloroFree.sizeFilt.fa"
    output:
        f"{config['map_rRNA_tRNA_free_dir']}/{{sample}}.mito_chloroFree.rRNA_tRNA_free.fa"
    params:
        tRNArRNA_index=f"{config['INDEX_tRNA_rRNA']}"
    resources:
        mem_mb=16000 
    threads: 4
    shell:"""
        mkdir -p {config[map_rRNA_tRNA_free_dir]}
        
        bowtie -v 0 -p {threads} -fS -x {params.tRNArRNA_index} {input} |\
            samtools fasta -f 4 - >  {output}
        
        """

rule map_to_genome_Shortstack:
    input:
        f"{config['map_rRNA_tRNA_free_dir']}/{{sample}}.mito_chloroFree.rRNA_tRNA_free.fa"
    output:
        directory(f"{config['sht_dir']}/{{sample}}_ShortStack_out")
    params:
        sh_index=f"{config['genome_Shortstack']}"
    resources:
        mem_mb=8200
    threads: 4
    shell:"""
        mkdir -p {config[sht_dir]}

        ShortStack --nohp --mmap f --align_only --sort_mem 8G \
         --bowtie_m all --bowtie_cores {threads} \
          --outdir {output} \
          --readfile {input} \
           --genomefile {params.sh_index}


        """