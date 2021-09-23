## each sample (C1-C4, S1-S4) was aligned individually, shown here is the code for one sample
## the same command was used for all other samples

# C2
/athena/abc/scratch/paz2005/bin/src/cellranger-3.0.2/cellranger count --id=C2 \
--localcores 12 --localmem 128 --transcriptome=cellranger-3.0.2/references/refdata-cellranger-GRCh38-3.0.0/ \
--fastqs=fastq_path/HGJWGBGX9/C22
