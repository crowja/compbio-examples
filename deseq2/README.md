# deseq2

## Example 1 running the DESeq2 pipeline on a set of count data.

Raw counts for five samples are used. The first two samples correspond to
replicates of Condition A and the remaining three to Condition B. This will
figure in in the DESeq2 analysis.

*   mklevels.py. Creates a fake set of expression counts.
*   run-mklevels.bash. Runs mklevels.
*   levels.txt. Count data produced by mklevels.py.
*   ex1.R. Runs the DESeq2 analysis of levels.txt.
*   run-ex1.bash. Runs ex1.R.
