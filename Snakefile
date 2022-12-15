# Author: Adam Freedman
# Date: 12/15/22
# modified version of workflow by:
    # Author: Maya Bose
    # Date: 5/31/19
    # Run sRNA mapping workflow
    # See README.md for usage instructions


# Get configuration
res_config = yaml.load(open("config/resources.yaml"),Loader=yaml.FullLoader)
configfile: "config/config.yaml"
SAMPLES = config["samples"]
#

# Run workflow
rule all:
    input:
        expand("data/7_fastqs/{sample}_fastqc.zip", sample=SAMPLES),
        expand("data/1_raw/{sample}_fastqc.zip", sample=SAMPLES),
        expand("data/7_fastqs/{sample}_length_profile.txt", sample=SAMPLES)

# Index reference genomes
onstart:
    shell("scripts/index_genomes.sh")

for ext in "fastq fq fastq.gz fq.gz".split():
    rule fastqc:
        input:
            expand("data/1_raw/{{sample}}.{ext}", ext=ext)
        output:
            "data/1_raw/{sample}_fastqc.zip"
        conda:
            "workflow/envs/fastqc.yml"
        threads:
            res_config["fastqc"]["threads"]
        resources:
            mem_mb = res_config['fastqc']['mem_mb']
            time = res_config['fastqc']['time']
        shell:
            '''
            fastqc -o data/1_raw/ -t {threads} {input}
            '''

# Trim reads
for ext in "fastq fq fastq.gz fq.gz".split():
    rule trimgalore:
        input:
            expand("data/1_raw/{{sample}}.{ext}", ext=ext)
        output:
            "data/2_trimmed/{sample}_trimmed.fq.gz"
        conda:
            "workflow/envs/trimgalore.yml"
        threads:
            res_config["trimgalore"]["threads"]
        resources:
            mem_mb = res_config['trimgalore']['mem_mb']
            time = res_config['trimgalore']['time']
        params:
            min_length = config["trimming"]["min_length"],
            max_length = config["trimming"]["max_length"],
            adapter_seq = config["trimming"]["adapter_seq"],
            quality = config["trimming"]["quality"],
        shell:
            '''
            trim_galore \
            --adapter {params.adapter_seq} \
            --gzip \
            --length {params.min_length} \
            --max_length {params.max_length} \
            --output_dir data/2_trimmed/ \
            --quality {params.quality} \
            --fastqc_args "--outdir" \
            {input} 2>> output_logs/2_outlog.txt
            
            '''

# Filter out contaminating highly expressed RNAs
rule filter_rna:
    input:
        "data/2_trimmed/{sample}_trimmed.fq.gz"
    output:
        fqgz = "data/3_ncrna_filtered/{sample}_ncrna_filtered.fq.gz"
    conda:
        "workflow/envs/bowtie.yml"
    threads:
        res_config["filter_rna"]["threads"]
    resources:
        mem_mb = res_config['filter_rna']['mem_mb']
        time = res_config['filter_rna']['time']
    params:
        fq = "data/3_ncrna_filtered/{sample}_ncrna_filtered.fq",
        rna_genome = config["genomes"]["filter_rna"]
    run:
        if {params.rna_genome} == "./genomes/filter_rna/":
            shell("echo No contaminating RNA filter genome provided, \
            skipping this step")
        else:
            shell(
            '''
            bowtie \
            -v 0 \
            -m 50 \
            --best \
            -a \
            --nomaqround \
            --norc \
            --threads {threads} \
            --un {params.fq} \
            {params.rna_genome} \
            {input} 1>> output_logs/3_outlog.txt \

            gzip {params.fq} &&
            
            fastqc -o data/3_ncrna_filtered/ -t 1 {output}
            ''')

# Filter out chloroplast and mitochondrial RNA
for ext in "fastq fq fastq.gz fq.gz".split():
    rule filter_chloroplast_mt_rna:
        input:
            name=expand("data/1_raw/{{sample}}.{ext}" \
                if (str(config["trimming"]["min_length"]).strip() == "" and \
                    str(config["trimming"]["max_length"]).strip() == "" and \
                    str(config["trimming"]["adapter_seq"]).strip() == "" and \
                    str(config["trimming"]["quality"]).strip() == "") \
                    else ("data/2_trimmed/{{sample}}_trimmed.{ext}" \
                    if str(config["genomes"]["filter_rna"]).strip() == "./genomes/filter_rna/" \
                        else "data/3_ncrna_filtered/{{sample}}_ncrna_filtered.{ext}"), ext=ext)
        output:
            fqgz = "data/4_c_m_filtered/{sample}_c_m_filtered.fq.gz"
        conda:
            "workflow/envs/bowtie.yml"
        threads:
            res_config['filter_chloroplast_mt_rna']['threads']
        resources:
            mem_mb = res_config['filter_chloroplast_mt_rna']['mem_mb']
            time = res_config['filter_chloroplast_mt_rna']['time']
        params:
            fq = "data/4_c_m_filtered/{sample}_c_m_filtered.fq",
            c_m_genome = config["genomes"]["chloro_mitochondria"],
            fastqc_threads = 1
        shell:
            '''
            bowtie \
            -v 0 \
            -m 50 \
            --best \
            -a \
            --nomaqround \
            --threads {threads} \
            --un {params.fq} \
            {params.c_m_genome} \
            {input.name} 1>> output_logs/4_outlog.txt \

            gzip {params.fq} &&

            fastqc -o data/4_c_m_filtered/ -t {params.fastqc_threads} {output}
            '''

# Cluster and align reads
rule cluster:
    input:
        expand("data/4_c_m_filtered/{sample}_c_m_filtered.fq.gz",
        sample=SAMPLES)
    output:
        "data/5_clustered/merged.bam"
    conda:
        "workflow/envs/shortstack.yml"
    threads:
        res_config['cluster']['threads']
    resources:
        mem_mb = res_config['cluster']['mem_mb']
        time = res_config['cluster']['time']
    params:
        genome = config["genomes"]["reference_genome"],
        multi_map_handler = config["aligning"]["multi_map_handler"],
        sort_memory = 0.8*res_config['cluster'['mem_mb'],
        nohp = config["aligning"]["no_mirna"],
        mismatches = config["aligning"]["mismatches"]
    shell:
        '''
        if [[ {params.nohp} == "Y" ]]; then
            hp="--nohp"
        else
            hp=""
        fi

        rm -r data/5_clustered && \
        ShortStack \
        --sort_mem {params.sort_memory} \
        --mismatches {params.mismatches}\
        --mmap {params.multi_map_handler} \
        --bowtie_cores {threads} \
        $hp \
        --readfile {input} \
        --genomefile {params.genome}.fasta \
        --outdir data/5_clustered/ 2>> output_logs/5_outlog.txt && \

        mv data/5_clustered/*.bam data/5_clustered/merged.bam && \

        scripts/combine_counts_results.py data/5_clustered/Counts.txt \
        data/5_clustered/Results.txt \
        --output_dir data/5_clustered/
        '''

# Split merged alignments file into multiple BAM files by sample name
rule split_by_sample:
    input:
        "data/5_clustered/merged.bam"
    output:
        expand("data/6_split_by_sample/{sample}_c_m_filtered.bam", sample=SAMPLES)
    conda:
        "workflow/envs/samtools.yml"
    threads:
        res_config['split_by_sample']['threads']
    resources:
        mem_mb = res_config['split_by_sample']['mem_mb']
        time = res_config['split_by_sample']['time']
    shell:
        '''
        mkdir -p data/6_split_by_sample && \

        samtools \
        split \
        -f '%!.bam' \
        {input} 2>> output_logs/6_outlog.txt && \

        mv *.bam data/6_split_by_sample/
        '''

# Extract mapped reads into BAM files
rule convert_1:
    input:
        "data/6_split_by_sample/{sample}_c_m_filtered.bam"
    output:
        temp("data/temp_converted/int1/{sample}_int1.bam")
    conda:
        "workflow/envs/samtools.yml"
    threads:
        res_config["convert_1"]["threads"]
    resources:
        mem_mb = res_config['convert_1']['mem_mb']
        time = res_config['convert_1']['time']
    shell:
        '''
        samtools \
        view \
        -F4 \
        -b \
        -@ {threads} \
        {input} > {output} 2>> Error.txt
        '''

# Convert BAM files to Fastq files
rule convert_2:
    input:
        "data/temp_converted/int1/{sample}_int1.bam"
    output:
        temp("data/temp_converted/{sample}_converted.fq")
    conda:
        "workflow/envs/samtools.yml"
    threads:
        res_config['convert_2']['threads']
    resources:
        mem_mb = res_config['convert_2']['mem_mb']
        time = res_config['convert_2']['time']
    shell:
        '''
        samtools bam2fq -t {input} > {output} 2>> Error.txt
        '''

# Get encoding quality for Fastq files
rule retrieve_encoding_quality:
    input:
        "data/temp_converted/{sample}_converted.fq"
    output:
        "data/7_fastqs/{sample}.fastq.gz"
    script:
        "scripts/match_qual_v2.py"

# Print length profiles of each sample to a log file
rule log_lengths_end:
    input:
        "data/7_fastqs/{sample}.fastq.gz"
    output:
        "data/7_fastqs/{sample}_fastqc.zip"
    conda:
        "workflow/envs/fastqc.yml"
    threads:
        config["threads"]["fastqc_report"]
    shell:
        '''
        fastqc -o data/7_fastqs/ -t {threads} {input} \
        2>> output_logs/7_outlog.txt
        '''

rule size_profiles:
    input:
        "data/7_fastqs/{sample}.fastq.gz"
    output: 
        "data/7_fastqs/{sample}_length_profile.txt"
    shell:
        '''
        scripts/fastq_readlength_profile.py {input} > {output}
        '''

onsuccess:
    shell("rm -r data/temp_converted/")
    print("Workflow finished!")
