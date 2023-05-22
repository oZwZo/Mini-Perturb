library(Seurat)

# read data
hoxb8_2khvg_raw <- "/home/wergillius/Project/Alison_Project/data/HoxB8/TFnet_test_glm_10xmatrices/"
Count.10x <- Read10X(data.dir = hoxb8_2khvg_raw)

meta <- read.csv(paste0(hoxb8_2khvg_raw, "metadata.csv"))
rownames(meta) <- meta$Sample

TFnet_obj <-CreateSeuratObject(counts = Count.10x, meta.data = meta)

# scTransform
correct_TFnet <- SCTransform(TFnet_obj, vars.to.regress = "Plate", verbose = FALSE)


# visualize uncorrected
TFnet_obj <- NormalizeData(TFnet_obj, verbose = FALSE, assay = "RNA")
TFnet_obj <- FindVariableFeatures(TFnet_obj, verbose = FALSE, assay = "RNA")
TFnet_obj <- ScaleData(TFnet_obj, verbose = FALSE)
TFnet_obj <- RunPCA(TFnet_obj,  verbose = FALSE)
TFnet_obj <- RunUMAP(TFnet_obj, dims = 1:50, reduction = "pca")
DimPlot(TFnet_obj, label = TRUE, label.size = 3, group.by = 'Plate')

# visualize corrected
correct_TFnet <- RunPCA(correct_TFnet,  verbose = FALSE)
correct_TFnet <- RunUMAP(correct_TFnet, dims = 1:50, reduction = "pca")
DimPlot(correct_TFnet, label = FALSE, label.size = 3, group.by = 'Plate')


# save result
Matrix::writeMM(correct_TFnet@assays$SCT@data, paste0(hoxb8_2khvg_raw, "scTransform_X.mtx"))

saveRDS(correct_TFnet, "data/HoxB8/TFnext_exon_2khvg_scTransform.rds")
