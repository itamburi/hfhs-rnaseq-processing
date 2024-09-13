# hfhs-rnaseq-processing

This directory demonstrates how two RNA-sequencing batches were batch-corrected using SVA in R to create counts matrices that are reported in our manuscript.

In `data/raw counts` :

-   `counts filtered 10 in min2 rows.csv` is the counts matrix for batch 1 samples

-   `full counts matrix HF 96 files.csv` is the counts matrix for batch 2 samples

The script `sva fitting gene expression and counts.R` can be executed as-is to reproduce the processing steps and outputs that include the final, batch-adjusted counts matrix in `data/processed` named `fitted filtered and batch adjusted counts combined.csv`

Examples of the STAR calls used for alignment are included in our SLURM sub file `star_align_pig.sub` . An example process indicating how bam files were assembled into counts matrices is shown in `counts matrix assembly example.R`
