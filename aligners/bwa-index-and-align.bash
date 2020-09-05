#! /bin/bash

# Create a BWA index using reference.fasta, put it in the subdirectory indices.
# For large references specify the bwtsw algorithm via -a bwtswkk
##mkdir -p indices
##bwa index -p indices/reference reference.fasta 

# Do the BWA alignment, output as SAM.
# Don't forget -a to report all alignments.
##bowtie2 -a -f -U reads.fasta -x indices/reference -S reads_aligned.sam

# Do the BWA alignment, output as sorted BAM, then create an associated BAM index file.
# Don't forget -a to report all alignments.
# Use Samtools to convert the SAM output to sorted BAM.
# Use Samtools to index the sorted BAM.
##bwa align indices/reference reads_1.fq reads_2.fq \
##bwa align indices/reference reads.fq \
bwa mem /scale03/fs0/gpfs0/research/jac/refdata/Zea_mays.AGPv4.fa r1.fq r2.fq \
| samtools view -b \
| samtools sort - -o reads-aligned-sorted.bam
samtools index reads-aligned-sorted.bam reads-aligned-sorted.bai


