#! /bin/bash

# Create a Bowtie2 index using reference.fasta, put it in the subdirectory indices.
mkdir -p indices
bowtie2-build -f reference.fasta indices/reference

# Do the Bowtie2 alignment, output as SAM.
# Don't forget -a to report all alignments.
##bowtie2 -a -f -U reads.fasta -x indices/reference -S reads_aligned.sam

# Do the Bowtie2 alignment, output as sorted BAM, then create an associated BAM index file.
# Don't forget -a to report all alignments.
# Use Samtools to convert the SAM output to sorted BAM.
# Use Samtools to index the sorted BAM.
bowtie2 -a -f  -U reads.fasta -x indices/reference \
| samtools view -bS \
| samtools sort - -o reads-aligned-sorted.bam
samtools index reads-aligned-sorted.bam reads-aligned-sorted.bai


