# Pseudotime inference

To order the cells along pseudotime trajectories, we first determined diffusion maps using the methods developed by [Haghverdi et al.](https://www.ncbi.nlm.nih.gov/pubmed/27571553) and implemented in the [`destiny` package](https://bioconductor.org/packages/release/bioc/html/destiny.html).

Wrapper functions (originally based off of <https://rdrr.io/bioc/scater/src/R/runDiffusionMap.R>):

```
#' Adding diffusion map coordinates to a SCE object
#' @return SCE object with additional reducedDim entry, see \code{reducedDimNames}
#' @seealso \code{\link{run_slingshot}}
#' @import scater
run_destiny <- function(sce.object, exprs_values = "log1p_sctransform",
  ntop = 1000, seed = 345, use_reddim = NULL, dmname = NULL){
  
  if(is.null(use_reddim)){
    inmat <- scater:::.get_mat_for_reddim(sce.object,
      exprs_values = exprs_values, ntop = ntop, scale=TRUE)
  inmat <- as.matrix(inmat)
  }else{
    inmat <- reducedDim(sce.object, use_reddim, withDimnames=FALSE)
    inmat <- as.matrix(inmat)
  }

  message("Generating diffusion map using destiny")
  set.seed(seed)
  dm <- destiny::DiffusionMap(inmat, distance = 'e')
  message("Finding tip cells.")
  tc <- destiny::find_tips(dm)
  
  message("Calculating pseudotime")
  dptime <- NULL
  tryCatch(
    dptime <- destiny::DPT(dm, tips = tc),
    error = function(e){cat("DPT ERROR :",conditionMessage(e), "\nThere will be no pseudotime output in the SCE object, but the diffusion map should still be there.")}
  )
  
  if(is.null(dmname)){
    if("diffMap" %in% names(colData(sce.object))){
      rdn <- paste0("diffMap_", gsub("[ :]","_", Sys.time() ) )
    }else{ rdn <- "diffMap"}
  }else{ rdn <- dmname}
  
  reducedDim(sce.object, rdn) <- dm@eigenvectors
  
  if(!is.null(dptime)){
    colData(sce.object)[paste(rdn, "dpt1", sep = "_")] <- dptime[tc, ][1,]
    colData(sce.object)[paste(rdn, "dpt2", sep = "_")] <- dptime[tc, ][2,]
    colData(sce.object)[paste(rdn, "dpt3", sep = "_")] <- dptime[tc, ][3,]
  }
  
  return(sce.object)
}


#' Running slingshot
#' @details The clusters identified in this step will be used to determine the global structure of the underlying lineages (that is, their number, when they branch off from one another, and the approximate locations of those branching events). This is different than the typical goal of clustering single-cell data, which is to identify all biologically relevant cell types present in the dataset. For example, when determining global lineage structure, there is no need to distinguish between immature and mature neurons since both cell types will, presumably, fall along the same segment of a lineage. Here, the clustering method implemented in Mclust is used; it assumes that Euclidean distance in a low-dimensional space reflect biological differences between cells and features an automated method for determining the number of clusters based on the Bayesian information criterion (BIC).
#' @return List of 2 \code{\link{SlingshotDataSets}} ((i) lineage and (ii) curves) and a vector of cluster labels
#' @seealso \code{\link{run_destiny}}
#' @import Mclust
#' @import slingshot
run_slingshot <- function(sce.object, rd_type = "diffMap", seed = 123,
  clusterLabels=NULL, ...){
  
  rd <- reducedDim(sce.object, rd_type)[, c(1:2)]
  
  
  if(is.null(clusterLabels)){
    message("Identifying clusters by Gaussian mixture modeling")
    clsres <- mclust::Mclust(rd)
    clusterLabels <- clsres$classification
  }
  
  message("Calculating lineages")
  set.seed(seed)
  lin <- slingshot::getLineages(rd, 
    clusterLabels = clusterLabels,
    ...)
  crvs <- slingshot::getCurves(lin)
  return(list(lineage = lin, curves = crvs, clusters = clusterLabels))
}
```

Using the PCA coordinates of the integrated data sets:


```
library(slingshot); library(scater); library(destiny); library(mclust)
scf <- readRDS(file = "/scratchLocal/frd2007/2019_05_MikeNotaras/2019-09-05/sce_integratedData_noS4_2019-09.rds")

#sc_list.cond <- list()
slng_list.cond_int <- list()

for(CONDITION in unique(as.character(scf$condition))){
  print(CONDITION)
 # sc_list.cond[[CONDITION]] <- scf[, scf$condition == CONDITION]
  
  ## diffusion maps
  sc_list.cond[[CONDITION]] <- run_destiny(sc_list.cond[[CONDITION]],
    exprs_values = "log1p_sctransform", use_reddim = "pca_integrated", seed = 345,
    dmname = "diffMap_integrated")
  
  ## Slingshot with no start/endpoints defined
  slng_list.cond_int[[CONDITION]] <- run_slingshot(sc_list.cond[[CONDITION]],
    rd_type = "diffMap_integrated",seed = 123, clusterLabels = NULL)
}

rdl <- lapply(sc_list.cond, reducedDims)

## add dimension red. coordinates to scf
reducedDim(scf, "diffMap_integrated_perCondition") <- rbind(rdl$Scz$diffMap_integrated,
  rdl$Ctrl$diffMap_integrated)[colnames(scf),]

```

Identifying temporally expressed genes

```
#### temporally expressed genes ===============================================
library(gam)
pt_list.cond_int <- list()
gam_pval_list.cond_int <- list()

for(CONDITION in unique(scf$condition)){
  print(CONDITION)
  
  # matrix of pseudotime values or cells' weights along each lineage
  pt_list.cond_int[[CONDITION]] <- slingPseudotime(slng_list.cond_int[[CONDITION]]$curves)[,1]
  
  # only look at the 500 most variable genes
  Y <- log1p(assay(sc_list.cond[[CONDITION]], "log1p_sctransform"))
  var1K <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:500]
  Y <- Y[var1K,]
  
  # fit a GAM with a loess term for pseudotime
  # We will regress each gene on the pseudotime variable we have generated,
  # using a general additive model (GAM)
  gam_pval_list.cond_int[[CONDITION]] <- apply(Y,1,function(z){
    d <- data.frame(z=z, t = pt_list.cond_int[[CONDITION]])
    tmp <- gam(z ~ lo(pt_list.cond_int[[CONDITION]]), data=d)
    p <- summary(tmp)[4][[1]][1,5]
    return(p)}
    )
}
```
