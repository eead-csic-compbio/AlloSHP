#!/bin/bash

# sample name
accession=Ae_speltoides_T
# Fastq file (Forward)
FASTQ1=ERR5051909_pass_1.fastq
# Fastq file (Reverse)
FASTQ2=ERR5051909_pass_2.fastq
# Trimmomatic parameters (optional)
crop=250
headcrop=10
minlen=240
# Concatenated reference genome for mapping step
reference=CONCAT_Turartu_AND_Aetauschii_AND_Aespeltoides.vcf2alignment.fna

# Print variables

echo $accession;
echo $FASTQ1;
echo $FASTQ2;
echo $crop;
echo $headcrop;
echo $minlen;
echo $reference;


# Separate PE interleaved file in two files (only if necessary)

perl /SOFT/split_pairs/split_pairs.pl -i $FASTQ -n -1 ${accession}_reads1.fq -2 ${accession}_reads2.fq -p "/" && \

# Trimmomatic (check QC report) (only if necessary)

java -jar /SOFT/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 \
	$FASTQ1 $FASTQ2 \
	${accession}_pass_1_paired.fq ${accession}_pass_1_unpaired.fq \
	${accession}_pass_2_paired.fq ${accession}_pass_2_unpaired.fq \
	ILLUMINACLIP:/SOFT/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True CROP:$crop HEADCROP:$headcrop SLIDINGWINDOW:15:28 MINLEN:$minlen && \


rm -rf *unpaired.fq && \ # (optional)


# QC PE filtered using trimmomatic

mkdir report_filtered && \

/SOFT/FastQC/fastqc ${accession}_pass_1_paired.fq ${accession}_pass_2_paired.fq -o report_filtered/ && \

# WARNING: If the fastq files show the following headers:
   # @SRR4162932.1.1 1 length=151
   # @SRR4162932.1.2 1 length=151
# These headers must be formated as follows:

cat ${accession}_pass_1_paired.fq | cut -d " " -f 1 | sed 's/\./_/1' | sed 's/\./\//1' | gzip > ${accession}_reads1.formatted.fq.gz && \

	rm -rf ./fastq/ERR5051909_pass_1.fastq && \

cat ${accession}_pass_2_paired.fq | cut -d " " -f 1 | sed 's/\./_/1' | sed 's/\./\//1' | gzip > ${accession}_reads2.formatted.fq.gz && \

	rm -rf ./fastq/ERR5051909_pass_2.fastq && \

FASTQ1formatted=${accession}_reads1.formatted.fq.gz
FASTQ2formatted=${accession}_reads2.formatted.fq.gz


gzip ${accession}_pass_1_paired.fq ${accession}_pass_2_paired.fq && \

# Resulting headers:
  # @SRR4162932_1/1
  # @SRR4162932_1/2

# Map paried-end reads (PE) to concatenated reference genome

  # WARNING:  minimap2 suggests using the --split-prefix when there are more than 4G bases in the reference.
  # It is also possible to increase the 4G limit in reference size to avoid using a split reference which produces incorrect MAPQ values

minimap2 -t 10 -I 128G -ax sr $reference $FASTQ1formatted $FASTQ2formatted | \

# Convert minimap2 output (SAM format) to sorted BAM format and keep only the mapped reads (-F 4). This means "exclude unmapped reads"
# Additional filters can be added in this step at the user's discretion (see samtools manual)

samtools view -bS -F 4 | samtools sort -@ 10 -o ${accession}.sort.bam && \

# Remove intermiediate files

rm -rf $FASTQ1formatted $FASTQ2formatted && \

# Compress filtered files (optional)

gzip ${accession}_pass_1_paired.fq ${accession}_pass_2_paired.fq && \

# Compute statistics (optinal but recommended)

samtools flagstat ${accession}.sort.bam > ${accession}.flagstat.txt && \

samtools index -@ 10 -c ${accession}.sort.bam && \

samtools idxstats -@ 10 ${accession}.sort.bam > ${accession}.idxstats.txt && \


## The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped read-segments and # unmapped read-segments.
## It is written to stdout. Note this may count reads multiple times if they are mapped more than once or in multiple fragments.

# Variant calling: BAM to VCF

  # mpileup --> Generate VCF or BCF containing genotype likelihoods for one or multiple alignment (BAM or CRAM) files
  # call --> SNP/indel calling (-O v --> output-type vcf; -m --> alternative model for multiallelic and rare-variant calling designed to overcome known limitations in -c calling model (conflicts with -c))

bcftools mpileup -f $reference -a DP ${accession}.sort.bam | bcftools call -O v -m > ${accession}.vcf && \

# Compress vcf (mandatory using bgzip)

bgzip ${accession}.vcf && \

# Index vcf for downstream analysis

bcftools index ${accession}.vcf.gz
