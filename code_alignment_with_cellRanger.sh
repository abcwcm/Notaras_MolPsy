## each sample (C1-C4, S1-S4) was aligned individually, shown here is the code for one sample
## the same command was used for all other samples
## some samples (such as C2) were sequenced twice because the first run was done at
## sufficient depth

# C2
/athena/abc/scratch/paz2005/bin/src/cellranger-3.0.2/cellranger count --id=C2 \
--localcores 12 --localmem 128 --transcriptome=cellranger-3.0.2/references/refdata-cellranger-GRCh38-3.0.0/ \
--fastqs=fastq_path/HGJWGBGX9/C22,fastq_path/HYLM5BGX9/C2 ## two fastq paths because this sample was sequenced twice
