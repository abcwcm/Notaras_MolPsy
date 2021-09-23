# Differential gene expression between cells of control and disease organoids

Wrapper functions for gene-wise testing of expression differences between two groups of cells:

```
#' Summing counts across cells of the same group
#' @param sce.object
#' @param split_by define the colData accessor that defines individual samples,
#'e.g. "Sample"
#' @return data.frame of gene counts where each row corresponds to a gene and
#' each column to a sample;
#' the value is the sum of counts across all cells of that sample
sum_counts_per_sample <- function(sce.object, split_by = "Sample"){
  
  if(!"counts" %in% assayNames(sce.object)){
    stop("We need the assay 'counts' in the SCE object")}
  n_cells_all <- dim(sce.object)[2]
  
  gns_counts <- as.data.table(as.matrix(counts(sce.object)), keep.rownames = TRUE)
  setnames(gns_counts, "rn", "gene")
  message("Gathering counts in a skinny data.table")
  gns_counts <- melt(gns_counts, id.vars = "gene", variable.name = "cell")
  gns_counts <- gns_counts[value > 0]
  
  ci <- colData(sce.object)[, split_by, drop = FALSE] %>% as.data.frame %>%
    as.data.table(., keep.rownames = TRUE)
  setnames(ci, "rn", "cell")
  gns_counts <- gns_counts[ci, on = "cell"]
  
  message(paste("Summing counts per", split_by))
  countsums <- gns_counts[, sum(value), by = c("gene", split_by)]
  setnames(countsums, split_by, "sample")
  
  ## returning matrix
  cc.df <- dcast(countsums, gene ~ sample, value.var = "V1") %>% as.data.frame
  rownames(cc.df) <- cc.df$gene
  cc.df$gene <- NULL
  cc.df[is.na(cc.df)] <- 0
  
  return(cc.df)
}

#' @title Run edgeR for comparing disease vs. control samples
#' @details The control sample will be the reference set. Colnames of
#' counts_df are expected to start with either C (--> Control) or
#' S (--> disease)
#' @return data.frame of edgeR results
run_edger_ctrl_vs_scz <- function(counts_df, sample_info, min_cpm = 1 , 
  min_samples = 4){
  library(edgeR)
  DGEl <- DGEList(counts = counts_df, group = sample_info)

  keep <- rowSums( cpm(DGEl) >= min_cpm) >= min_samples
  DGEl <- DGEl[keep,]
  DGEl$samples$lib.size <- colSums(DGEl$counts)
  DGEl <- calcNormFactors(DGEl)
  DGEl <- estimateDisp(DGEl, model.matrix(~si.edg) )
  fit <- glmQLFit(DGEl,model.matrix(~si.edg))
  fit.res <- glmQLFTest(fit)
  res <- as.data.frame(topTags(fit.res, n = Inf, sort.by = "none"))
  return(res)
}

```

### Progenitors

Cells of clusters 1,2,3,6,7

```{r de_for_progs, eval=FALSE}
datadir <- "~/Documents/Projects/2019-04_MikeNotaras/data/"

scf.de <- scf[, scf$labels %in% c("Progs.", "Proliferating")]
dim(scf.de)
# [1] 20844  13064
## adding counts back in to calculate the sum of reads across cells of the same sample
cts <- readRDS(file = paste0(datadir, "counts.rds"))
rownames(scf.de) <- rowData(scf.de)$id
assay(scf.de, "counts") <- cts[rownames(scf.de),colnames(scf.de)]
rm(cts);gc()
rownames(scf.de) <- uniquifyFeatureNames(rownames(scf.de), rowData(scf.de)$symbol)

scf.de <- scater::calculateQCMetrics(scf.de)
keep_genes <- rowData(scf.de)$n_cells_by_counts >= 5
scf.de <- scf.de[keep_genes,]
## preparing for edgeR
Progs.df <- sum_counts_per_sample(scf.de, "Sample")
si.edg <- factor(ifelse(grepl("^C", colnames(Progs.df)), "Control",
  ifelse(grepl("^S", colnames(Progs.df)), "Schizophrenia", NA)))
si.edg <- relevel(si.edg, ref = "Control")

de.progs <- run_edger_ctrl_vs_scz(Progs.df, sample_info = si.edg, min_cpm = 1,
  min_samples = 4)
saveRDS(de.progs,
  file = paste0(datadir, "de_res_Progs_2019-09-10.rds"))
rm(scf.de); gc()
```

## Gene set over-representation

We used the tests for over-representation of a given gene set within our list of DEG (`goi`) as implemented in `clusterProfiler` and `reactomePA`. 

```
library(clusterProfiler)
## need ENTREZ IDs
eg <-  clusterProfiler::bitr(rownames(scf), fromType="SYMBOL", 
    toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% as.data.table
setnames(eg, names(eg), c("gene_symbol", "entrez"))

## defining genes of interest
goi <- rownames(subset(de.progs, FDR <= 0.05))
goi <- goi[!grepl("^RP[LS]", goi)]
eg <- eg[!grepl("^RP[LS]", gene_symbol)]


## KEGG =======================
ora_kegg.progs <- enrichKEGG( unique(eg[gene_symbol %in% goi]$entrez), 
                             universe = unique(eg$entrez), 
                             organism = "hsa")

ora_kegg.progs <- setReadable(ora_kegg.progs, 'org.Hs.eg.db', 'ENTREZID')

## REACTOME =======================
ora_reactome.progs <- ReactomePA::enrichPathway( unique(eg[gene_symbol %in% goi]$entrez), 
                             universe = unique(eg$entrez), 
                             organism = "human")

ora_reactome.progs <- setReadable(ora_reactome.progs, 'org.Hs.eg.db', 'ENTREZID')


## GO terms =====================
ora_goBP.progs <- enrichGO(gene  = eg[gene_symbol %in% goi]$entrez,
  universe      = unique(eg$entrez),
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

ora_goMF.progs <- enrichGO(gene = eg[gene_symbol %in% goi]$entrez,
  universe      = unique(eg$entrez),
  OrgDb         = org.Hs.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
```

Visualizations of the ORA results were done with functions of the `clusterProfiler` package:

```
## dotplot
dotplot(ora_goBP.ProliferatingRG) + 
  scale_color_gradientn(colours =rev(c("darksalmon", "firebrick3", "firebrick4")))

## cnetplot
cnetplot(ora_goBP.ProliferatingRG, showCategory = 10, 
         colorEdge = TRUE, foldChange = genes_for_cp) +
  scale_colour_gradient2(name = "log2FC", low = "navyblue", 
    high = "red", mid = "white")

## heatplot
heatplot(ora_goBP.ProliferatingRG,
  foldChange = genes_for_cp, showCategory = 10) +
  scale_fill_gradient2(name = "logFC", 
    low = "blue2", mid = "cornsilk1",
    high = "firebrick4", midpoint = 0) +
    theme(legend.position = "bottom")
```
