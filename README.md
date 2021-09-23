# Notaras et al. (2021) Molecular Psychiatry

Scripts and code details related to Notaras et al., 2021, more specifically to the single-cell RNA-seq data presented there.

Several dozen organoids were dissociated using Accutase and serially filtered to remove debris before being pelleted via slow centrifugation.
Viability and cell numbers were determined prior to 10X Genomics Chromium library preparation, which was carried out in the Tilgner Lab at Weill Cornell Medicine.
Sequencing was was performed at the Genomics Core, Weill Cornell Medicine.

The processing was done in 2019 when Seurat's scTransform method was en vogue.
Raw data (e.g. read counts generated with CellRanger) are available from GEO.
The final, filtered, processed R object including cell labels, clustering, dimensionality reductions etc. can be downloaded [here](https://wcm.box.com/shared/static/a4vfzycxz0j8irg2catitszjzekmg9t9.rda). The object is of the format `SingleCellExperiment`, i.e. you will need the [package of the same name](http://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html).

For questions, don't hesitate to reach out to Friederike Dündar at the [Applied Bioinformatics Core](https://abc.med.cornell.edu/) or by raising an issue here.

![](WCM_MB_LOGO_HZSS1L_CLR_RGB.png)

----------------------------------------------------------

## Alignment, filtering, normalization and batch correction

Initial processing of the raw reads and alignment to GRCh38 was done with the CellRanger pipeline (v. 3.0.2; https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/what-is-cell-ranger).
For details see `code_alignment_with_cellRanger.sh`.

Subsequent analyses were performed in R following the recommendations of Amezquita et al. [Ref: https://www.biorxiv.org/content/10.1101/590562v1] (https://osca.bioconductor.org/) using numerous functions provided in the R packages scater and scran [Refs: https://dx.doi.org/10.12688/f1000research.9501.2, https://dx.doi.org/10.12688/f1000research.9501.2].

For details of the filtering and processing, see `code_processing.md`. 
In brief, the following steps were taken:
Based on the calculation of outliers, we subsequently removed cells with fewer than 317 genes and more than 8% mitochondrial reads.
Genes that were expressed in fewer than 5 cells were removed as well.
We then processed and integrated the different samples using Seurat version 3.1 following the recommendations from the Satija Lab's vignette (https://satijalab.org/seurat/v3.1/integration.html).
More specifically, read counts were first normalized using SCTransform for each sample individually [Ref: https://www.biorxiv.org/content/10.1101/576827v1].

The different samples were then integrated using the top 3,000 most informative genes [Ref: https://doi.org/10.1016/j.cell.2019.05.031] before performing various dimensionality reduction steps including PCA and UMAP [Ref: https://www.nature.com/articles/nbt.4314]. *For details, see `seurat_alignment.R`.*
A shared nearest neighbor graph was constructed using Seurat's FindNeighbors function with default settings (e.g. k = 20) using the first 25 principal components.
Subsequent clustering was performed with Seurat's FindClusters function with the resolution parameter set to 0.2 [Ref: https://link.springer.com/article/10.1140%2Fepjb%2Fe2013-40829-0].
For visualizations and assessments of normalized expression values, the SCTransform-normalized (log-transformed) expression values were used unless noted otherwise.

### Labelling cells

To assign cell type labels to subpopulations of single cells, we used information obtained from (i) genes that were significantly reduced or increased in a given Seurat-identified cluster when compared to the expression in cells belonging to all other clusters ("marker genes") and (ii) automated assignment of cell type labels using SingleR with two different reference data sets [Ref: https://dx.doi.org/10.1038/s41590-018-0276-y].
To identify marker genes of the individual clusters and samples, Seurat's FindMarkers and FindAllMarkers functions were used with default settings as well as scran::FindMarkers(), comparing only cells of the control condition.
For SingleR, we used its in-built reference data set of the human primary cell atlas, which contains RNA-seq from adult human tissues, including neurons and myeloid cells.
In addition, we downloaded the single-cell RNA-seq data from numerous human fetal brain samples that was published by Nowakowski et al. [Ref: https://science.sciencemag.org/content/358/6368/1318.full] and used it as the training data with the refactored Bioconductor version of SingleR (http://www.bioconductor.org/packages/devel/bioc/html/SingleR.html).
For details, see `code_cellLabeling.md`.

### Diffusion maps and pseudotime trajectory analyses

For ordering the cells along pseudotime trajectories, we first constructed diffusion maps using the methods developed by Haghverdi et al. [Ref: https://www.ncbi.nlm.nih.gov/pubmed/27571553] that are implemented in the destiny package [Ref: https://www.ncbi.nlm.nih.gov/pubmed/26668002]. 
We then inferred cell lineages and temporally expressed genes using slingshot [Ref: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4772-0] and following the recommendations of the vignette (https://www.bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html).
In brief, pseudotime values were extracted for each cell along each lineage identified by slingshot.
Then, each gene's normalized expression values were regressed on the pseudotime variable using a general additive model to identify those whose expression patterns most closely correlate with the pseudotime weights.
For details, see `code_pseudotime.md`.

### Identification of differentially expressed genes between control and disease organoids

For determining genes that were differentially expressed between cells that were assigned to the same clusters after Seurat's integration procedure (see above), we used a pseudo-bulk approach [Ref: https://www.ncbi.nlm.nih.gov/pubmed/28334062], summing the reads across all cells belonging to the same type of cells (where "type" would be determined by, for example, the sample and the membership within a given subpopulation of interest).
To identify genes with significantly different expression levels between two types of cells (DEG), we used methods implemented in the edgeR package including dispersion estimation and genewise negative binomial generalized linear models with quasi-likelihood tests [Ref: https://dx.doi.org/10.1093/nar/gks042].
Unless stated otherwise, DEG were based on adjusted p-values lower than 0.01.
For details of the code, see `code_DEG.md`.

### Over-representation of gene sets

To test for over-representation of specific gene sets of KEGG, REACTOME and gene ontology collections in our DEG lists, we used the enrichment functions implemented in clusterProfiler and reactomePA [Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339379/] (enrichKEGG, enrichPathway, enrichGO) with an adjusted p-value cut-off of 0.05.
For details of the code, see `code_DEG.md`.

-----------------

## References

* Amezquita et al.: <https://www.biorxiv.org/content/10.1101/590562v1>
* clusterProfiler: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339379/>
* destiny: <https://www.ncbi.nlm.nih.gov/pubmed/26668002>
* edgeR: <https://dx.doi.org/10.1093/nar/gks042>
* Haghverdi et al: <https://www.ncbi.nlm.nih.gov/pubmed/27571553>
* Nowakowski et al. (human fetal brain data set): <https://science.sciencemag.org/content/358/6368/1318.full>
* pseudo-bulk: <https://www.ncbi.nlm.nih.gov/pubmed/28334062>
* scater: <https://www.ncbi.nlm.nih.gov/pubmed/28088763>
* scTransform: <https://www.biorxiv.org/content/10.1101/576827v1>
* scran: <https://dx.doi.org/10.12688/f1000research.9501.2>
* Seurat integration: <https://doi.org/10.1016/j.cell.2019.05.031>
* SingleR: <https://dx.doi.org/10.1038/s41590-018-0276-y>
* slingshot: <https://www.ncbi.nlm.nih.gov/pubmed/29914354>
* UMAP: <https://www.nature.com/articles/nbt.4314>
