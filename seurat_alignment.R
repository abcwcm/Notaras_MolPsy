###################################
#### SEURAT ALIGNMENT #############
###################################

library(scater)
library(Seurat) # v. 3.1.0.
library(ggplot2); theme_set(theme_bw(base_size = 16))

load(file = "sce_CellGeneFilt_2019-05.rda")


## First, setup the Seurat object list, and run SCTransform on each object **separately**:
srt_list <- list()
for(i in unique(sce.filt$Sample)){
    srt_list[[i]] <- CreateSeuratObject(counts = counts(sce.filt[, sce.filt$Sample == i]), 
        project = "cerOrg", assay = "RNA",
        min.cells = 0, min.features = 0,
        names.field = 1#, will use the first part of the colnames as CELLTYPE
        # meta.data = as.data.frame(colData(sce.filt) )
        )
}

## remove S4
srt_S4_scTransform <- srt_list$S4
srt_list$S4 <- NULL

for (i in seq_len(length(srt_list)) ) {
    print(names(srt_list[i]))
    srt_list[[i]] <- SCTransform(srt_list[[i]], verbose = FALSE)
}

save(srt_list, srt_S4_scTransform, file = "seurat_list_scTransform.rda")

## Next, select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.
int_features <- SelectIntegrationFeatures(object.list = srt_list, nfeatures = 3000)

options(future.globals.maxSize = 4000 * 1024^2)
srt_list <- PrepSCTIntegration(object.list = srt_list,
    anchor.features = int_features, 
    verbose = FALSE)

## Next, identify anchors and integrate the datasets. Commands are identical to the standard workflow, but make sure to set normalization.method = 'SCT':
int_anchors <- FindIntegrationAnchors(object.list = srt_list,
    normalization.method = "SCT", 
    anchor.features = int_features, verbose = FALSE)

srt_integrated <- IntegrateData(anchorset = int_anchors,
    normalization.method = "SCT", 
    verbose = FALSE)

#Now proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset. Commands are identical to the standard workflow, but do not run the ScaleData function after integration. 
srt_integrated <- RunPCA(srt_integrated, verbose = FALSE)
srt_integrated <- RunUMAP(srt_integrated, dims = 1:25)

plots <- DimPlot(srt_integrated, group.by = c("orig.ident"), combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") +
        guides(color = guide_legend(nrow = 3, 
            byrow = TRUE, override.aes = list(size = 3))))

png("test.png")
CombinePlots(plots)
dev.off()

## clustering & marker detection -------------------------------------------
srt_integrated <- FindNeighbors(srt_integrated, dims = 1:25) 
srt_integrated <- FindClusters(srt_integrated, resolution = .3)
markers_seurClsts.3 <- FindAllMarkers(srt_integrated)
srt_integrated <- FindClusters(srt_integrated, resolution = .2)
markers_seurClsts.2 <- FindAllMarkers(srt_integrated)
srt_integrated <- FindClusters(srt_integrated, resolution = .1)
markers_seurClsts.1 <- FindAllMarkers(srt_integrated)

saveRDS(list(res.3 = markers_seurClsts.3,
    res.2 = markers_seurClsts.2,
    res.1 = markers_seurClsts.1),
    file = "markers_seurat_Sep05.rds")

save(int_anchors, int_features, file = "seurat_integrationAnchors_and_Features.rda")

dim(srt_integrated)
#[1]  3000 26335

srt_integrated@assays
# $RNA
# Assay data with 26276 features for 26335 cells
# First 10 features:
#     ISG15, SDF4, AURKAIP1, MRPL20, SSU72, GNB1, FAAP20, FAM213B, TPRG1L,
# SMIM1 
# 
# $SCT
# Assay data with 20844 features for 26335 cells
# First 10 features:
#     ISG15, SDF4, AURKAIP1, MRPL20, SSU72, GNB1, FAAP20, FAM213B, TPRG1L,
# SMIM1 
# 
# $integrated
# Assay data with 3000 features for 26335 cells
# Top 10 variable features:
#     TTR, HBZ, COL3A1, STMN2, CRABP1, APOA1, DLK1, MEG3, MALAT1, HIST1H4C 
srt_integrated@active.assay
#[1] "integrated"
saveRDS(srt_integrated, file = "seurat_integrated_noS4.rds")

## For visualization etc., they use the COORDINATES of the integrated data set,
## but the expression values of "RNA" 

#########################################################################
## Add Seurat info back to SCE ==========================================
DefaultAssay(srt_integrated) <- "SCT"
dim(srt_integrated)
##[1] 20844 26335

scf <- readRDS("sce_CellGeneFiltWithScTransform_2019-06.rds")
scf.new <- scf[, scf$Sample != "S4"] ## determined by the experimentator to be an outlier/no representative of typical organoids

## reduced dimensions -----------------------
reducedDim(scf.new, "pca_integrated") <- Embeddings(object = srt_integrated,
    reduction = "pca")[colnames(scf.new),]
reducedDim(scf.new, "umap_integrated") <- Embeddings(object = srt_integrated,
    reduction = "umap")[colnames(scf.new),]

## clustering info ----------------------------
cd_srt <-  srt_integrated@meta.data
names(colData(scf.new))[names(colData(scf.new)) == "res.03"] <- "old.clusters_res.03"
names(colData(scf.new))[names(colData(scf.new)) == "res.01"] <- "old.clusters_res.01"
scf.new$integClust_res.1 <- cd_srt[colnames(scf.new),]$integrated_snn_res.0.1
scf.new$integClust_res.3 <- cd_srt[colnames(scf.new),]$integrated_snn_res.0.3

### normalized (scTransformed) data ------------
rnms <- data.table(ori = row.names(scf.new)) # get original rownames -- Seurat replaces underscores with hyphens...
rnms[,seurat := gsub("_","-", ori)]

## filtering to match the cells of Seurat
normd <- as.data.table(as.matrix(GetAssayData(object = srt_integrated)),
    keep.rownames = TRUE) # log1p(corrected counts, i.e. Pearson residuals)
setnames(normd, "rn", "seurat")
normd <- rnms[normd, on = "seurat"]
scf.new <- scf.new[normd$ori, ] ## filtering!

normd.m <- as.matrix(normd[, -c("seurat", "ori"), with=FALSE])
rownames(normd.m) <- normd$ori
normd.m <- normd.m[rownames(scf.new), colnames(scf.new)]
assay(scf.new, "log1p_sctransform") <- normd.m

## save final SCE ------------------------------------------
saveRDS(scf.new, file = "sce_integratedData_2019-09.rds")

rm(normd);rm(normd.m);rm(rnms)
gc()
