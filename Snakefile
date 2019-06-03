# Author: Maya Bose
# Date: 5/31/19
# Run sRNA mapping workflow
# See README.md for usage instructions


# Get configuration
configfile: "config.yaml"
SAMPLES = config["samples"]
#

# Run workflow
rule all:
	input:
		expand("data/8_fastqs/{sample}.fastq",sample=SAMPLES)


# Trim reads
rule trim:
	input:
		"data/1_raw/{sample}.fastq.gz"
	output:
		"data/2_trimmed/{sample}_trimmed.fq.gz"
	params:
		min_length = config["trim"]["min_length"],
		max_length = config["trim"]["max_length"],
		adapter_seq = config["trim"]["adapter_seq"],
		quality = config["trim"]["quality"]
	shell:
		"trim_galore "
		"--adapter {params.adapter_seq} "
		"--gzip "
		"--length {params.min_length} "
		"--max_length {params.max_length} "
                        "--output_dir data/2_trimmed/ "
		"--quality {params.quality} "
		"{input}"


# Filter out junk mRNA
rule filter_rfam:
	input:
		"data/2_trimmed/{sample}_trimmed.fq.gz"
	output:
		"data/3_rfam_filtered/{sample}_rfam_filtered.fq"
	params:
		rfam_genome = config["filter_rfam"]["genome"]
	shell:
			"bowtie "
			"-v 0 "
			"-m 50 "
			"--best "
			"-a "
			"--nomaqround "
			"--norc "
			"--un {output} "
			"{params.rfam_genome} "
			"{input}"	


# Filter out chloroplast and mitochondrial RNA
rule filter_c_m:
	input:
		"data/3_rfam_filtered/{sample}_rfam_filtered.fq"
	output:
		"data/4_c_m_filtered/{sample}_c_m_filtered.fq"
	params:
		c_m_genome = config["filter_c_m"]["genome"]
	shell:
		"bowtie "
		"-v 0 "
		"-m 50 "
		"--best "
		"-a "
		"--nomaqround "
		"--un {output} "
		"{params.c_m_genome} "
		"{input}"


# Cluster and align reads
rule cluster:
	input:
		expand("data/4_c_m_filtered/{sample}_c_m_filtered.fq", sample=SAMPLES)

	output:
		"data/5_clustered/Counts.txt",
		"data/5_clustered/merged_alignments.bam"
	params:
		bowtie_cores = config["cluster"]["bowtie_cores"],
		genome = config["cluster"]["genome"]
	shell:
		"rm -r data/5_clustered && " # Need this line because Snakemake
				     # creates this dir, but ShortStack
				     # won't add files to a dir that already
				     # exists
		"ShortStack "
		"--sort_mem 20G "
		"--mismatches 0 "
		"--mmap u "
		"--bowtie_cores {params.bowtie_cores} "
		"--nohp "
		"--readfile {input} "
		"--genomefile {params.genome} "
		"--outdir data/5_clustered/"

rule split_by_sample:
	input:
		"data/5_clustered/merged_alignments.bam"
	output:
		expand("{sample}_c_m_filtered.bam",sample=SAMPLES)
	shell:
		"samtools split "
  		"-f '%!.bam' "
		"{input}"

rule rename:
	input:
		"{sample}_c_m_filtered.bam"
	output:
		"data/6_split_by_sample/{sample}_aligned.bam"
	shell:
		"mv {input} {output}"


# Extract mapped reads into BAM file
rule convert_1:
	input:
		"data/6_split_by_sample/{sample}_aligned.bam"
	output:
		"data/7_converted/int1/{sample}_int1.bam"
	shell:
		"samtools view -F4 -b {input} > {output}"


# Convert BAM file to Fastq file
rule convert_2:
	input:
		"data/7_converted/int1/{sample}_int1.bam"
	output:
		"data/7_converted/{sample}_converted.fq"
	shell:
		"samtools bam2fq -t {input} > {output}"


# Split Fastq file into multiple Fasta files by Sample name
rule split_fastqs:
	input:
		"data/7_converted/{sample}_converted.fq"
	output:
		"data/8_fastqs/{sample}.fastq"
	script:
		"scripts/match_qual_v2.py"
