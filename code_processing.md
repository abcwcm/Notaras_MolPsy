## QC and filtering

Cells: min. 10^2.5 genes and max. 8 % mitochondrial reads.
Genes: expr_cutoff <- 10^-2.5


```
library(scater)
library(magrittr)
library(Matrix)
library(ggplot2); theme_set(theme_bw(base_size = 16))
condition_cols <- c("dodgerblue","coral1")
replicate_cols <- paste0("bisque",1:4)
library(patchwork)

## the following two packages are in-house packages, most 
## functions are available in other packages now, such as DropletUtils
library(scABC2) # v. 1.5.0. to provide function for reading in 10X data
library(ABCutilities) # v. 0.3.0

samples <- c("S1","S2","S3","S4","C1","C2","C3")
data_dirs <- paste0("mike/counts/", samples, "/outs/filtered_feature_bc_matrix/")
names(data_dirs) <- samples
sce <- read10XResults.abc(data_dir=data_dirs,
    barcode_file = "barcodes.tsv.gz",
    genes_file = "features.tsv.gz",
    matrix_file = "matrix.mtx.gz",
    expand = TRUE, min_total_cell_counts = 0, min_mean_gene_counts=0) ## these have been filtered by CR anyway

########################## Initial Filtering ###############################
dim(sce)
#[1] 33538 33547
keep_gns <- rowSums(counts(sce)) >0
table(keep_gns)
#keep_gns
#FALSE  TRUE 
# 7003 26535 
sce <- sce[keep_gns,]
dim(sce)
#[1] 26535 33547

## add chromosome information -----------------------------------
#BiocManager::install("EnsDb.Hsapiens.v86")        
library(EnsDb.Hsapiens.v86)                                    
gn_location <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(sce),
                      column = "SEQNAME", keytype = "GENEID")
rowData(sce)$chr <- paste0("chr", gn_location)
summary(gn_location == "MT")
sce <- scater::calculateQCMetrics(sce,  feature_controls = list(mito = which(gn_location == "MT")))

colData(sce)$condition <- ifelse(grepl("^S", colData(sce)$Sample), "Scz", 
                                       ifelse(grepl("^C", colData(sce)$Sample), "Ctrl", NA))
colData(sce)$replicate <- gsub("[SC]", "", colData(sce)$Sample)

########################## PART 2 CELL FILTERING ###############################
high.mito <- lapply(samples, function(x) isOutlier(sce[,sce$Sample == x]$pct_counts_mito, nmads = 5, type = "higher"))
names(high.mito) <- samples

thresh.list <- list( min_log10_total_features_by_counts = 2.5, 
                     max_pct_counts_mito = 8)


## filter cells ---------------------------------------------------------
qc_df <- as.data.frame(colData(sce))
filtered <- list()

keep_cells <- row.names(subset(qc_df,
	log10_total_features_by_counts >= thresh.list$min_log10_total_features_by_counts &
	pct_counts_mito <= thresh.list$max_pct_counts_mito))

if(length(keep_cells) > 0){
  sce <- sce[,keep_cells]
}
dim(sce)
# [1] 26535 31828

## immediately remove no-coverage-genes --------------------------------------
gnsfilt1 <- scABC2::filter_gene_countBased(sce.object = sce, min.rowSums = 0, return_gns= TRUE)

if(length(gnsfilt1$keep_genes) >0){
  sce <- sce[gnsfilt1$keep_genes, ]
}
dim(sce)
# [1] 26471 31828


########################## Calculate cell cycle information #############################

hg.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

cc.noNormAllCells <- list()

for(i in samples){
  set.seed(123)
  message(i)
  p = bpstart(MulticoreParam(12))
  sc.tmp <- sce[, sce$Sample == i]
  cc.noNorm <- scran::cyclone(sc.tmp, pairs=hg.pairs, BPPARAM=p)
  names(cc.noNorm$phases) <- colnames(sc.tmp)
  row.names(cc.noNorm$scores) <- colnames(sc.tmp)
  row.names(cc.noNorm$normalized.scores) <- colnames(sc.tmp)
  cc.noNormAllCells[[i]] <- cc.noNorm
  rm(sc.tmp)
  rm(cc.noNorm)
  gc()
}

phs <- lapply(cc.noNormAllCells, function(x) x$phases) %>% unlist
names(phs) <- gsub(".+\\.","", names(phs))
colData(sce)$cc_phase <- phs[colnames(sce)]

g1phs <- lapply(cc.noNormAllCells, function(x) x$score$G1) %>% unlist
colData(sce)$G1score <- g1phs

g2mphs <- lapply(cc.noNormAllCells, function(x) x$score$G2M) %>% unlist
colData(sce)$G2Mscore <- g2mphs

########################## PART 3 GENE FILTERING ###############################
## make more sensible rownames ---------------------------
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$id, rowData(sce)$symbol)

## update QC measures
sce <- scater::calculateQCMetrics(sce,  feature_controls = list(mito = which(rowData(sce)$chr == "chrMT")))

hk_genes <- unique(c(grep("^mt-", rownames(sce), value=TRUE, ignore.case=TRUE),  # mitochondrial genes
                     grep("^Rp[sl]", rownames(sce), value=TRUE, ignore.case=TRUE))) # ribosomal genes

## calculate drop out rates  per sample
gns_dropouts <- scABC2::calc_dropouts_per_cellGroup(sce, rownames(sce), split_by = "Sample")
gns_dropouts$condition <- gsub("[0-9]","", gns_dropouts$Sample)

cutoff <- -2.5
print(paste("Gene cutoff:", cutoff))

P <- ggplot(data = gns_dropouts,
       aes(x = log10(mean.pct.of.counts),
           y = log10(pct.zeroCov_cells + .1),
           text = paste(gene, condition, sep = "_"))) + 
  geom_point(aes(color = condition), shape = 1, size = .5, alpha = .5) +
  geom_point(data = gns_dropouts[gene %in% hk_genes],
             aes(fill = condition), shape = 22, size = 4, alpha = .8) +
  facet_grid(~condition) + ggtitle("Gene dropout rates -- zoomed in") +
  coord_cartesian(ylim = c(1.9,2.01), xlim = c(-3,-1)) +
  facet_wrap(~Sample) +
  geom_vline(xintercept = cutoff, linetype = "dashed", color = "grey") +
  scale_color_manual(values = condition_cols) +
  scale_fill_manual(values = condition_cols)

expr_cutoff <- 10^-2.5
keep_genes  <- gns_dropouts[mean.pct.of.counts >= expr_cutoff]$gene %>% unique

sce.filt <- sce[keep_genes,]
dim(sce.filt)
#[1] 26276 31828

sce.filt <- scater::calculateQCMetrics(sce.filt, feature_controls = list(mito = which(rowData(sce.filt)$chr == "chrMT")))
save(sce.filt, file = "sce_CellGeneFilt_2019-05.rda")
```

## Normalization & intergration with Seurat

Following recommendations of the [Satija Lab](https://satijalab.org/seurat/v3.1/integration.html) for data normalization and integration.

```
library(scater)
library(Seurat) # v. 3.1.0.
library(ggplot2); theme_set(theme_bw(base_size = 16))

load(file = "sce_CellGeneFilt_2019-05.rda")

## First, setup the Seurat object list, and run SCTransform on each object separately:
srt_list <- list()
for(i in unique(sce.filt$Sample)){
    srt_list[[i]] <- CreateSeuratObject(counts = counts(sce.filt[, sce.filt$Sample == i]), 
        project = "cerOrg", assay = "RNA",
        min.cells = 0, min.features = 0,
        names.field = 1#, will use the first part of the colnames as CELLTYPE
        # meta.data = as.data.frame(colData(sce.filt) )
        )
}

for (i in seq_len(length(srt_list)) ) {
    print(names(srt_list[i]))
    srt_list[[i]] <- SCTransform(srt_list[[i]], verbose = FALSE)
}

## Next, select features for downstream integration,
## and run PrepSCTIntegration, which ensures that all
## necessary Pearson residuals have been calculated.
int_features <- SelectIntegrationFeatures(object.list = srt_list, nfeatures = 3000)

options(future.globals.maxSize = 4000 * 1024^2)
srt_list <- PrepSCTIntegration(object.list = srt_list,
    anchor.features = int_features, 
    verbose = FALSE)

## Next, identify anchors and integrate the datasets.
## Commands are identical to the standard workflow,
## but make sure to set normalization.method = 'SCT':
int_anchors <- FindIntegrationAnchors(object.list = srt_list,
    normalization.method = "SCT", 
    anchor.features = int_features, verbose = FALSE)

srt_integrated <- IntegrateData(anchorset = int_anchors,
    normalization.method = "SCT", 
    verbose = FALSE)

## Now proceed with downstream analysis (i.e. visualization, clustering)
## on the integrated dataset. Commands are identical to the standard workflow,
## but do not run the ScaleData function after integration. 
srt_integrated <- RunPCA(srt_integrated, verbose = FALSE)
srt_integrated <- RunUMAP(srt_integrated, dims = 1:25)

plots <- DimPlot(srt_integrated, group.by = c("orig.ident"), combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x +
    theme(legend.position = "top") +
    guides(color = guide_legend(nrow = 3, byrow = TRUE,
        override.aes = list(size = 3))))

## clustering & marker detection ---------------------------------------
srt_integrated <- FindNeighbors(srt_integrated, dims = 1:25) 

srt_integrated <- FindClusters(srt_integrated, resolution = .3)
markers_seurClsts.3 <- FindAllMarkers(srt_integrated)

## Add Seurat info back to SCE ==========================================
DefaultAssay(srt_integrated) <- "SCT"
dim(srt_integrated)
##[1] 20844 26335

scf <- readRDS("sce_CellGeneFiltWithScTransform_2019-06.rds")

## reduced dimensions -----------------------
reducedDim(scf, "pca_integrated") <- Embeddings(object = srt_integrated,
    reduction = "pca")[colnames(scf),]
reducedDim(scf, "umap_integrated") <- Embeddings(object = srt_integrated,
    reduction = "umap")[colnames(scf),]

## clustering info ----------------------------
cd_srt <-  srt_integrated@meta.data
scf$integClust_res.3 <- cd_srt[colnames(scf),]$integrated_snn_res.0.3

### normalized (scTransformed) data ------------
rnms <- data.table(ori = row.names(scf)) # get original rownames -- Seurat replaces underscores with hyphens...
rnms[,seurat := gsub("_","-", ori)]

## filtering to match the cells of Seurat
normd <- as.data.table(as.matrix(GetAssayData(object = srt_integrated)),
    keep.rownames = TRUE) # log1p(corrected counts, i.e. Pearson residuals)
setnames(normd, "rn", "seurat")
normd <- rnms[normd, on = "seurat"]
scf <- scf[normd$ori, ] ## filtering!

normd.m <- as.matrix(normd[, -c("seurat", "ori"), with=FALSE])
rownames(normd.m) <- normd$ori
normd.m <- normd.m[rownames(scf), colnames(scf)]
assay(scf, "log1p_sctransform") <- normd.m
```
