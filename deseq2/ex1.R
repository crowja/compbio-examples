#! /usr/bin/env Rscript

# Hijacked from https://gist.github.com/stephenturner/f60c1934405c127f09a6,
# thank you Stephen Turner!

# STEP 1: Get the raw count data.
# genename sample01 sample02 sample03 sample04 sample05.
# The first two samples correspond to condition A the final three are condition B.

# Read levels.text, convert the data frame to to a matrix.

countdata <- as.matrix(read.table('levels.txt', header=TRUE, row.names=1))
##head(countdata)

# Set the condition associated with each column
conditions <- factor(c(rep("condition_a", 2), rep("condition_b", 3)))

# STEP 2: DESeq2 analysis
suppressMessages(suppressWarnings(library(DESeq2)))

# Create the coldata data frame for use in instantiating the DESeqDataSet.
coldata <- data.frame(row.names=colnames(countdata), conditions)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~conditions)
dds
