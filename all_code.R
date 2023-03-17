# Start -------------------------------------------------------------------
work_dir <- "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/"
dir.create(work_dir, recursive = T, showWarnings = FALSE)
setwd(work_dir)

library(SCP)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
library(UCell)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(plotly)
library(openxlsx)

srt_22 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_22_annotation.rds")
srt_h0 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_h0_annotation.rds")
srt_esc <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_esc_annotation.rds")

srt_imp <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_imp_annotation.rds")
srt_epi <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_epi_annotation.rds")
srt_amn <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_amn_annotation.rds")
srt_endo <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_endo_annotation.rds")
srt_meso <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_meso_annotation.rds")
srt_ecto <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_ecto_annotation.rds")

srt_endo_sub1 <- readRDS("srt_endo_sub1.rds")
srt_endo_sub2 <- readRDS("srt_endo_sub2.rds")
srt_endo_sub3 <- readRDS("srt_endo_sub3.rds")
srt_meso_sub1 <- readRDS("srt_meso_sub1.rds")
srt_meso_sub2 <- readRDS("srt_meso_sub2.rds")
srt_ecto_sub <- readRDS("srt_ecto_sub.rds")

srt_CS7 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_E19/human_gastrulation_E19.seurat.rds")
srt_E14 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_3dculture/human_gastrulation_3dculture.seurat.rds")
srt_mGast1 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/mouse_gastrulation_E8.5/mouse_gastrulation_E8.5.corrected.seurat.rds")
srt_mGast2 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/mouse_gastrulation_E8.25/mouse_gastrulation_E8.25.seurat.rds")
srt_mOrg <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/mouse_organogenesis_E13.5_atlas/mouse_organogenesis_E13.5_atlas.seurat.total.de.rds")
srt_hOrg <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_GW4_GW6_organogenesis/human_GW4_GW6_organogenesis.seurat.rds")

srt_mEndo1 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/mouse_gut_endoderm_E8.75/mouse_gut_endoderm_E8.75.seurat.rds")
srt_mEndo2 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/mouse_endoderm_E10.5/mouse_endoderm_E10.5.seurat.rds")
srt_mEndo3 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/mouse_foreguts_endoderm_E9.5/mouse_foreguts_endoderm_E9.5.seurat.rds")
srt_hEndo1 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_endodermal_organ_7W_21W/human_endodermal_organ_7W_21W.seurat.rds")
srt_hEndo2 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gut_6W_11W/human_gut_6W_11W.seurat.rds")

srt_mMeso <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/mouse_mesoderm_E7.75/mouse_mesoderm_E7.75.seurat.rds")

srt_mBrain <- readRDS("/data/tmp/HuangMingQian/scRNA-seq/ref/Gioele_E7_E18_mouse_brain/Gioele_E7_E18_mouse_brain.seurat.de.rds")
srt_hBrain <- readRDS("/data/tmp/HuangMingQian/scRNA-seq/ref/NeMO_GW14_GW25_wholebrain/NeMO_GW14_GW25_wholebrain.raw.seurat.rds")
srt_retina1 <- readRDS("/data/tmp/HuangMingQian/scRNA-seq/ref/Sandra_hESC_RPE/Sandra_hESC_RPE.seurat.rds")
srt_retina2 <- readRDS("/data/tmp/HuangMingQian/scRNA-seq/ref/Hu_retina_5W_24W/Hu_retina_5W_24W.seurat.rds")
srt_retina3 <- readRDS("/data/tmp/HuangMingQian/scRNA-seq/ref/Cameron_retinal_organoid/Cameron_retinal_organoid.seurat.rds")
srt_retina4 <- readRDS("/data/tmp/HuangMingQian/scRNA-seq/ref/Lu_retinal_GW8_PND8/Lu_retinal_GW8_PND8.seurat.rds")
srt_retina5 <- readRDS("/data/tmp/HuangMingQian/scRNA-seq/ref/Lukowski_retinal_adult/Lukowski_retinal_adult.seurat.rds")

srt_teratoma <- readRDS("/data/tmp/HuangMingQian/scRNA-seq/ref/McDonald_teratoma/McDonald_teratoma.seurat.rds")

# Markers -----------------------------------------------------------------
markers <- list(
  "Epiblast" = c("POU5F1", "DNMT3B", "TDGF1"),
  "Primitive streak" = c("TBXT", "EOMES", "FGF4"),
  "Mesoderm" = c("HAND1", "BMP4"),
  "MSC/Fib" = c("THY1", "COL1A1", "COL14A1", "COL15A1"),
  "ExE mesoderm" = c("POSTN", "FRZB", "EZR", "CDH1"), # CDH1(allantois)
  "Limb bud mesenchyme cell" = c("TBX5", "PRRX1", "TBX4"),
  "Endothelial & erythroid cell" = c("MEF2C", "HBE1", "PECAM1", "GATA1"),
  "Endoderm" = c("SOX17", "FOXA2", "HNF4A"),
  "ExE endoderm" = c("AFP", "VTN", "TTR"),
  "Gut" = c("FOXA1", "FAM3B", "VIL1"),
  "Amnion" = c("ISL1", "ABCG2", "TFAP2A"),
  "PGC" = c("PRDM1", "NANOS3", "CD38"),
  "Amnion & PGC" = c("ABCG2", "NANOS3"),
  "Epithelium" = c("KRT17", "WNT6", "KRT4", "KRT18"),
  "Preplacodal ectoderm" = c("SIX4", "DLX3", "SOX3", "FOXI3"),
  "Neural ectoderm" = c("SOX2", "LHX5", "SOX21", "ZIC2"),
  "Radial glial cell" = c("SOX1", "RFX4", "TTYH1", "FABP7", "NES"),
  "Ependymal cell" = c("FOXJ1", "PIFO", "DYNLRB2"),
  "Neuron" = c("DCX", "TUBB3", "TAGLN3", "SLC17A6"),
  "Schwann cell" = c("SOX10", "MPZ", "PLP1"),
  "Sensory neuron" = c("POU4F1", "PIEZO2", "NTRK1"),
  "Retinal progenitor cell" = c("SIX6", "VSX2", "HMX1"),
  "Retinal pigmented epithelium" = c("MLANA", "MITF", "SERPINF1", "RLBP1"),
  "Retinal progenitor cell & RPE" = c("VSX2", "MITF")
)

APS_marker <- c("CER1", "GSC", "EOMES", "MIXL1", "TDGF1", "OTX2") # "LHX1", "NODAL"
PPS_marker <- c("EVX1", "MSGN1", "CDX1", "CDX2", "HOXA1", "HOXB1") # "WNT5A", "WNT8A"


# EDA -----------------------------------------------------
srt_22 <- LoadH5Seurat("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/NGSmodule_SCP_analysis/CellQC/Merge.filtered.h5Seurat")
samples <- c(
  "22_CellLine_injection_10d-3", "22_CellLine_injection_20d-1", "22_CellLine_injection_30d-1",
  "22_CellLine_injection_40d-1", "22_CellLine_injection_50d-1", "22_CellLine_injection_70d-1", "22_CellLine_injection_90d-1"
)
srt_22 <- srt_22[, srt_22$orig.ident %in% samples]
srt_22$Sample <- factor(srt_22$orig.ident, levels = samples)
srt_22$Time <- gsub("(22_CellLine_injection_)|(-.*)", "", srt_22$orig.ident)
srt_22$Time <- factor(srt_22$Time, levels = paste0(c(1, 2, 3, 4, 5, 7, 9), "0d"))

srt_22 <- Integration_SCP(srt_22,
  integration_method = "Uncorrected", cluster_resolution = 5,
  # vars_to_regress = "nCount_RNA", regression_model = "negbinom"
)
for (i in unique(c(c(seq(20, 80, 5))))) {
  srt_22 <- CSS_integrate(srt_22,
    linear_reduction = "Uncorrectedpca", HVF = VariableFeatures(srt_22),
    linear_reduction_dims_use = 1:i, cluster_resolution = 3
  )
  srt_22[[paste0("CSSUMAP2D", i)]] <- srt_22@reductions$CSSUMAP2D
  srt_22[[paste0("CSSUMAP3D", i)]] <- srt_22@reductions$CSSUMAP3D
}
saveRDS(srt_22, "srt_22.rds")

plist <- list()
for (i in sort(unique(c(c(seq(20, 80, 5)))))) {
  plist[[as.character(i)]] <- ExpDimPlot(srt_22, "CDX2", reduction = paste0("CSSUMAP2D", i))
}
p <- plot_grid(plotlist = plist)
ggsave(plot = p, filename = "tmp.png", height = 30, width = 30, limitsize = F)

ClassDimPlot(srt_22, "CSSclusters", reduction = "CSSUMAP2D45")
ClassDimPlot3D(srt_22, "CSSclusters", reduction = "CSSUMAP3D40")

# ncount4000: 26,38,48,64,76,80
# ncount5000: 56,54,44,42,34,24
srt_4000umi <- readRDS("srt_22_filtered4000umi.rds")
srt_5000umi <- readRDS("srt_22_filtered5000umi.rds")
ClassDimPlot(srt_4000umi, "CSSclusters", reduction = "CSSUMAP2D54")
ClassDimPlot3D(srt_4000umi, "CSSclusters", reduction = "CSSUMAP3D54")
ClassDimPlot(srt_5000umi, "CSSclusters", reduction = "CSSUMAP2D54")
ClassDimPlot3D(srt_5000umi, "CSSclusters", reduction = "CSSUMAP3D54")

ClassDimPlot(srt, "CSSclusters",
  reduction = "CSSUMAP2D54",
  cells.highlight = colnames(srt)[!colnames(srt) %in% colnames(srt_5000umi)]
)

out <- proxyC::simil(x = srt_22@assays$RNA@data["TP63", , drop = FALSE], y = srt_22@assays$RNA@data, use_nan = TRUE)
out[is.na(out)] <- 0
sort(out[1, ], decreasing = T) %>% head(20)

srt_hOrg <- readRDS("reference_data/human_GW4_GW6_organogenesis/human_GW4_GW6_organogenesis.seurat.rds")
ClassDimPlot(srt_hOrg, "main_split", label = TRUE)
srt_hOrg <- RunPCAMap(
  srt_query = srt_hOrg, srt_ref = srt_22, ref_pca = "CSSpca",
  ref_group = "CSSclusters", ref_umap = "CSSUMAP2D"
)

srt_22$bg <- NA
p <- ProjectionPlot(
  srt_query = srt_hOrg, srt_ref = srt_22,
  query_group = "main_split", ref_group = "bg",
  query_param = list(palette = "Set1"),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "human org"
  )
)
ClassDimPlot(srt_hOrg, "main_split", label = TRUE, reduction = "ref.umap")

srt_22 <- RunKNNPredict(
  srt_query = srt_22, srt_ref = srt_hOrg,
  query_group = "CSSclusters", ref_group = "main_split"
)
ClassDimPlot(srt_22, "knnpredict_main_split")


srt_22_tmp <- RunKNNPredict(srt_22, bulk_ref = SCP::ref_scHCL)
drop <- table(srt_22_tmp$knnpredict) %>%
  .[. < 10] %>%
  names()
srt_22_tmp$knnpredict[srt_22_tmp$knnpredict %in% drop] <- "unreliable"
ClassDimPlot(srt_22_tmp, group.by = "knnpredict", label = T)

srt_22_tmp <- RunKNNPredict(srt_22, srt_ref = srt_mOrg, ref_group = "Main_cell_type")
ClassDimPlot(srt_22_tmp, group.by = "knnpredict_Main_cell_type", label = T)

srt_22_tmp <- RunKNNPredict(srt_22, srt_ref = srt_hOrg, ref_group = "main_split")
ClassDimPlot(srt_22_tmp, group.by = "knnpredict_main_split", label = T)

srt_22_tmp <- RunKNNPredict(srt_22, srt_ref = srt_mGast, ref_group = "celltype")
ClassDimPlot(srt_22_tmp, group.by = "knnpredict_celltype", label = T)

# Fig2 --------------------------------------------------------------------
dir.create(paste0(work_dir, "/figures/fig2"), recursive = T, showWarnings = FALSE)

## 22-CellLine ---------------------------------------------------------------
### Load the data -----------------------------------------------------
srt_22 <- LoadH5Seurat("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/NGSmodule_SCP_analysis/CellQC/Merge.filtered.h5Seurat")
samples <- c(
  "22_CellLine_injection_10d-3", "22_CellLine_injection_20d-1", "22_CellLine_injection_30d-1",
  "22_CellLine_injection_40d-1", "22_CellLine_injection_50d-1", "22_CellLine_injection_70d-1", "22_CellLine_injection_90d-1"
)
srt_22 <- srt_22[, srt_22$orig.ident %in% samples]
srt_22$Sample <- factor(srt_22$orig.ident, levels = samples)
srt_22$Time <- gsub("(22_CellLine_injection_)|(-.*)", "", srt_22$orig.ident)
srt_22$Time <- factor(srt_22$Time, levels = paste0(c(1, 2, 3, 4, 5, 7, 9), "0d"))

### Integration -------------------------------------------------------------
srt_22 <- Integration_SCP(srt_22, integration_method = "Uncorrected", linear_reduction_dims_use = 1:45, cluster_resolution = 5)
srt_22 <- Integration_SCP(srt_22, integration_method = "CSS", linear_reduction_dims_use = 1:45, cluster_resolution = 5)
# srt_22 <- FindClusters(object = srt_22, resolution = 5, algorithm = 1,graph.name = "CSS_SNN")
# srt_22 <- SrtReorder(srt_22, features = srt_22@misc$CSS_HVF, reorder_by = "seurat_clusters", slot = "data",log = F)
# srt_22[["seurat_clusters"]] <- NULL
# srt_22[["clusters"]] <- Idents(srt_22)

ClassDimPlot(srt_22, "CSSclusters", label = TRUE)
ClassDimPlot3D(srt_22, "CSSclusters")

srt_22[["UMAP"]] <- srt_22[["CSSUMAP2D"]]
srt_22[["UMAP3D"]] <- srt_22[["CSSUMAP3D"]]
srt_22@misc$Default_reduction <- "UMAP"
colnames(srt_22[["UMAP"]]@cell.embeddings) <- c("UMAP_1", "UMAP_2")
colnames(srt_22[["UMAP3D"]]@cell.embeddings) <- c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3")
Key(srt_22[["UMAP"]]) <- "UMAP_"
Key(srt_22[["UMAP3D"]]) <- "UMAP3D_"
srt_22@misc$Default_reduction <- "UMAP"

srt_22 <- RunDEtest(srt_22, group_by = "CSSclusters", test.use = "wilcox", BPPARAM = BiocParallel::MulticoreParam(workers = 30))
srt_22 <- RunDEtest(srt_22, group_by = "CSSclusters", test.use = "roc", BPPARAM = BiocParallel::MulticoreParam(workers = 30))

de_filter <- filter(srt_22@tools$DEtest_CSSclusters$AllMarkers_wilcox, p_val_adj < 0.05 & DE_group_number <= 5)
res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 2)
p <- panel_fix(p, save = "tmp.bar.pdf", height = 1.2, width = 2)
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 5)
p <- panel_fix(p, save = "tmp.word.pdf", height = 3, width = 3)
saveRDS(srt_22, "srt_22.rds")

### paga initiated layout ---------------------------------------------------------------
adata <- srt_to_adata(srt_22)
adata <- RunPAGA(
  adata = adata, group_by = "CSSclusters", linear_reduction = "CSS", nonlinear_reduction = "UMAP",
  n_pcs = ncol(srt_22@reductions$CSS@cell.embeddings), threshold = 0.1,
  embedded_with_PAGA = TRUE, paga_layout = "fa"
)

PAGAUMAP2D <- adata$obsm[["PAGAUMAP2D"]]
colnames(PAGAUMAP2D) <- c("umap_1", "umap_2")
rownames(PAGAUMAP2D) <- adata$obs_names$values
srt_22[["PAGAUMAP2D"]] <- CreateDimReducObject(PAGAUMAP2D, assay = "RNA")
ClassDimPlot(srt_22, group.by = "Time", reduction = "PAGAUMAP2D")

### SCVELO ------------------------------------------------------------------
plt <- reticulate::import("matplotlib")$pyplot
scv <- reticulate::import("scvelo")
adata <- srt_to_adata(srt_22)
adata <- RunSCVELO(
  adata = adata, group_by = "CSSclusters", linear_reduction = "CSS", nonlinear_reduction = "UMAP",
  n_pcs = ncol(srt_22@reductions$CSS@cell.embeddings), n_neighbors = 100
)
ax <- scv$pl$velocity_embedding_stream(adata,
  vkey = "stochastic", basis = "UMAP", title = "RNA Velocity", color = "CSS",
  density = 2, linewidth = 1, arrow_size = 1.5, legend_fontsize = 15, legend_loc = "right_margin",
  save = FALSE, show = FALSE, dpi = 300
)
ax$axis("equal")
plt$show()
plt$savefig("figures/fig2/22_SCVELO.stream.svg", bbox_inches = "tight")

### Annotation --------------------------------------------------------------
ClassDimPlot(srt_22, "CSSclusters", reduction = "CSSUMAP2D", label = T)
# ExpDimPlot(srt_22, c("PDGFRA", "KDR"), calculate_coexp = T) # Cardiomyocytes
# ExpDimPlot(srt_22, c("THY1", "ENG"), calculate_coexp = T)
# ExpDimPlot(srt_22, c("POSTN", "ENG"), calculate_coexp = T)
# epithelium marker, KRT18; epidermal marker, KRT5
# cardiomyocyte marker, TNNT2 ACTC1 NPPA HAND1
# LEF1 presomitic mesoderm
# Single cell analysis of human foetal liver captures the transcriptional profile of hepatobiliary hybrid progenitors
# CDH6, STAT1, CD24, FGFR2, DCDC2 and CTNND2 that we identified are restricted to the ductal plate (DP)
# intestinal stem cell ISCs (LGR5, ASCL2, CD166, and LRIG1), enterocytes (VIL1 and ANPEP),
# and secretory lineage cells (LYZ for Paneth cells, MUC2 for goblet cells, and CHGA for enteroendocrine cells)

# srt_22 <- RunKNNPredict(
#   srt_query = srt_22, query_group = "CSSclusters",
#   srt_ref = srt_22, ref_group = "CSSclusters",
#   query_collapsing = T, ref_collapsing = T,
#   return_full_distance_matrix = TRUE, nn_method = "raw"
# )
# d <- 1 - srt_22@tools$knnpredict$distance_matrix
# d <- as.matrix(d)
# Heatmap(d)

class <- "16,17,18,19,20,21,22,23,26: Epiblast
27: Primitive streak
34,35,36,37,38,39: Mesoderm
53,50,51: ExE mesoderm
46,47,48,49,52,54,55: MSC/Fib
45: Limb bud mesenchyme cell
40: Endothelial & erythroid cell
33,32: Endoderm
30: ExE endoderm
31: Gut
28,44: Amnion
29: PGC
41,42,43: Epithelium
25,24,1,2: Neural ectoderm
3,6: Radial glial cell
7,8: Ependymal cell
15: Neuron
4,5: Retinal progenitor cell
9,10,11,12,13,14: Retinal pigmented epithelium"

class <- readLines(textConnection(class))
srt_22[["Annotation"]] <- srt_22[["CSSclusters"]]
srt_22[["Annotation"]] <- as.character(srt_22[["Annotation", drop = TRUE]])
for (i in 1:length(class)) {
  m <- strsplit(class[i], ": ")[[1]]
  srt_22$Annotation[srt_22$Annotation %in% strsplit(m[1], ",")[[1]]] <- m[[2]]
}
levels <- c(
  "Epiblast", "Primitive streak",
  "Mesoderm", "MSC/Fib", "ExE mesoderm", "Limb bud mesenchyme cell", "Endothelial & erythroid cell",
  "Endoderm", "ExE endoderm", "Gut",
  "Amnion", "PGC", "Epithelium",
  "Neural ectoderm", "Radial glial cell", "Ependymal cell", "Neuron", "Retinal progenitor cell", "Retinal pigmented epithelium"
)
all(levels %in% srt_22$Annotation)
srt_22$Annotation <- factor(srt_22$Annotation, levels = levels)

srt_22$GermLayer <- sapply(srt_22$Annotation, function(x) {
  switch(as.character(x),
    "Epiblast" = "Epiblast",
    "Primitive streak" = "Primitive streak",
    "Mesoderm" = "Mesoderm",
    "ExE mesoderm" = "Mesoderm",
    "MSC/Fib" = "Mesoderm",
    "Limb bud mesenchyme cell" = "Mesoderm",
    "Endothelial & erythroid cell" = "Mesoderm",
    "Endoderm" = "Endoderm",
    "ExE endoderm" = "Endoderm",
    "Gut" = "Endoderm",
    "Amnion" = "Non-neural ectoderm",
    "PGC" = "PGC",
    "Epithelium" = "Non-neural ectoderm",
    "Neural ectoderm" = "Neural ectoderm",
    "Radial glial cell" = "Neural ectoderm",
    "Ependymal cell" = "Neural ectoderm",
    "Neuron" = "Neural ectoderm",
    "Retinal progenitor cell" = "Neural ectoderm",
    "Retinal pigmented epithelium" = "Neural ectoderm"
  )
})
srt_22$GermLayer <- factor(srt_22$GermLayer, levels = c("Epiblast", "Primitive streak", "Mesoderm", "Endoderm", "Neural ectoderm", "Non-neural ectoderm", "PGC"))

ClassDimPlot(srt_22, group.by = "Annotation", reduction = "CSSUMAP2D", label = TRUE, label_insitu = FALSE)
ClassDimPlot(srt_22, c("Annotation", "GermLayer"))

srt_22$Time <- factor(paste0("d", gsub("d", "", srt_22$Time)), levels = c("d10", "d20", "d30", "d40", "d50", "d70", "d90"))
saveRDS(srt_22, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_22_annotation.rds")

### Baisc plot --------------------------------------------------------------------
p <- ClassDimPlot(srt_22, "Annotation", label = TRUE, theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 4, save = "figures/fig2/22_annotation.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_22, "CSSclusters", label = TRUE, label_insitu = TRUE, theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 4, save = "figures/fig2/22_CSSclusters.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_22, "Time", theme_use = "theme_blank", show_stat = TRUE, force = TRUE)
p <- panel_fix(p, height = 4, save = "figures/fig2/22_sample.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_22, "GermLayer", theme_use = "theme_blank", label = TRUE, show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 4, save = "figures/fig2/22_germlayer.umap.pdf", raster = TRUE)

p <- ClassDimPlot(srt_22, "Annotation",
  split.by = "Time", legend.position = "none", bg_color = "grey90", nrow = 2,
  theme_use = "theme_blank", force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig2/22_sample_split.umap.pdf", raster = TRUE)

p <- ClassDimPlot3D(srt_22, "Annotation",
  reduction = "CSSUMAP3D",
  width = 1000, height = 700, axis_labs = paste0("UMAP_", 1:3),
  save = "figures/fig2/22_annotation.umap3D.html"
)

### Annotation stat -------------------------------------------------------------------------
p <- ClassStatPlot(srt_22, stat.by = "Annotation", group.by = "Time", plot_type = "area", aspect.ratio = 0.5)
p <- panel_fix(p, width = 5, save = "figures/fig2/22_annotation.stat.pdf", raster = TRUE)


### Markers -----------------------------------------------------------------
# SLC5A12,LRP2,AQP1,NPHS1 (and LTL): PT and nephron progenitor cells
# "Chondrocytes" = c("FMOD", "IGFL2"),
# "Smooth muscle" = c("ACTA2", "EMILIN1"),

n <- 2
markers_use <- markers[levels(srt_22$Annotation)] %>% lapply(function(x) head(x, 2))
all(srt_22$Annotation %in% names(markers_use))

p <- ExpDimPlot(srt_22,
  features = unlist(markers_use), ncol = 6,
  theme_use = "theme_blank"
)
p <- panel_fix(p, height = 2, save = "figures/fig2/22_markers.umap.pdf", raster = TRUE)

p <- ExpDotPlot(srt_22,
  genes = unlist(markers_use),
  feature_split = unlist(lapply(names(markers_use), function(x) rep(x, length(markers_use[[x]])))),
  cell_split_by = "Annotation"
)
ggsave(plot = p, filename = "figures/fig2/22_markers.dotplot.pdf", width = 10, height = 15)

ht <- GroupHeatmap(srt_22,
  group.by = "Annotation",
  features = unlist(markers_use),
  feature_split = rep(names(markers_use), each = 2),
  show_row_names = TRUE, add_dot = TRUE, add_bg = TRUE,
  feature_split_palette = "Paired", heatmap_palette = "YlOrRd",
  height = 10, width = 8, dot_size = unit(7, "mm")
)
panel_fix(ht$plot, save = "figures_add/22testis.markers.pdf")

### Cell cycle --------------------------------------------------------------
srt_22 <- CellCycleScoring(srt_22,
  s.features = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes
)
srt_22$Phase <- factor(srt_22$Phase, levels = c("G1", "S", "G2M"))
p <- ClassDimPlot(srt_22, "Phase", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig2/22_cellcycle_phase.umap.pdf", raster = TRUE)

p <- ClassStatPlot(srt_22, stat.by = "Phase", group.by = "Annotation", aspect.ratio = 0.5)
p <- panel_fix(p, width = 5, save = "figures/fig2/22_cellcycle.stat.pdf", raster = TRUE)

### PAGA ---------------------------------------------------------------------
srt_22_paga <- RunPAGA(
  srt = srt_22, group_by = "Annotation", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
  n_pcs = ncol(srt_22@reductions$CSS@cell.embeddings), n_neighbors = 100, embedded_with_PAGA = FALSE, return_seurat = TRUE
)
srt_22@misc$paga <- srt_22_paga@misc$paga
p <- ClassDimPlot(srt_22,
  group.by = "Annotation", pt.size = 5, pt.alpha = 0.01,
  label = TRUE, label.size = 3, label_repel = TRUE, label_insitu = TRUE, label_segment_color = "transparent",
  paga = srt_22@misc$paga, paga_edge_threshold = 0.1, paga_edge_color = "black", paga_edge_alpha = 1,
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
)
p <- panel_fix(p, width = 4, save = "figures/fig2/22_paga.pdf", raster = TRUE)

### SCVELO ------------------------------------------------------------------
srt_22_scv <- SCP:::RunSCVELO(
  srt = srt_22, group_by = "Annotation", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
  n_pcs = ncol(srt_22@reductions$CSS@cell.embeddings), n_neighbors = 100, return_seurat = TRUE
)
srt_22[["stochastic_UMAP"]] <- srt_22_scv[["stochastic_CSSUMAP2D"]]
srt_22[["Ms"]] <- srt_22_scv[["Ms"]]
srt_22[["Mu"]] <- srt_22_scv[["Mu"]]
srt_22[["stochastic"]] <- srt_22_scv[["stochastic"]]
srt_22[["variance_stochastic"]] <- srt_22_scv[["variance_stochastic"]]
p <- ClassDimPlot(srt_22,
  group.by = "Annotation", pt.size = 5, pt.alpha = 0.01,
  velocity = "stochastic", velocity_plot_type = "stream",
  velocity_density = 2, velocity_smooth = 1,
  streamline_n = 20, streamline_size = 0.5, streamline_color = "black",
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
)
p <- panel_fix(p, width = 4, save = "figures/fig2/22_scvelo.pdf", raster = TRUE)

### DE analysis ------------------------------------------------------------------
srt_22 <- RunDEtest(srt_22,
  group_by = "CSSclusters", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_22 <- RunDEtest(srt_22,
  group_by = "CSSclusters", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
de_filter <- filter(srt_22@tools$DEtest_CSSclusters$AllMarkers_wilcox, p_val_adj < 0.05 & DE_group_number <= 5)
res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = F, ncol = 3)

srt_22 <- RunDEtest(srt_22,
  group_by = "Annotation", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_22 <- RunDEtest(srt_22,
  group_by = "Annotation", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)

de_filter <- filter(srt_22@tools$DEtest_Annotation$AllMarkers_wilcox, p_val_adj < 0.05)
ribogene <- rownames(srt_22)[grep("^RP[SL]\\d+\\w{0,1}\\d*$", rownames(srt_22))]
de_filter <- de_filter[!de_filter$gene %in% ribogene, ]
write.xlsx(de_filter, file = "figures/fig2/22_allDE.xlsx")

de_top <- de_filter %>%
  group_by(gene) %>%
  top_n(1, avg_log2FC) %>%
  group_by(group1) %>%
  top_n(3, avg_log2FC)
p <- ExpDotPlot(srt_22, features = de_top$gene, feature_split = de_top$group1, cell_split_by = c("Annotation", "Time"))
ggsave(plot = p, filename = "figures/fig2/22_topDE.heatmap.pdf", width = 10, height = 15)

ht <- ExpHeatmap(srt_22,
  features = de_filter$gene, feature_split = de_filter$group1, cell_split_by = "Annotation", use_raster = TRUE,
  anno_terms = TRUE, anno_features = TRUE, topTerm = 4, topWord = 15, db_version = "3.13",
  nlabel = 0, feature_split_palette = "Paired", width = 7, height = 12
)
ggsave(plot = ht$plot, filename = "figures/fig2/22_allDE.heatmap.pdf", width = 20, height = 20, limitsize = FALSE)
saveRDS(ht, "figures/fig2/22_allDE.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig2/22_allDE.enrichment.xlsx")

srt_22 <- RunEnrichment(srt_22, group_by = "Annotation", geneID_exclude = ribogene, db_version = "3.13")
p <- EnrichmentPlot(srt_22, group_by = "Annotation", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "figures/fig2/22_bp_allDE.bar.pdf", raster = TRUE)
p <- EnrichmentPlot(srt_22, group_by = "Annotation", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "figures/fig2/22_bp_allDE.word-term.pdf", raster = TRUE)
p <- EnrichmentPlot(srt_22,
  group_by = "Annotation", plot_type = "wordcloud", word_type = "feature", combine = TRUE, ncol = 6,
  topWord = 30, word_size = c(4, 8), aspect.ratio = 0.3, legend.position = "bottom"
)
p <- panel_fix(p, height = 1.5, save = "figures/fig2/22_bp_allDE.word-gene2.pdf", raster = TRUE)

saveRDS(srt_22, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_22_annotation.rds")

### Cell similarity -----------------------------------------------------
srt_22 <- RunKNNPredict(
  srt_query = srt_22, query_group = "Annotation",
  srt_ref = srt_22, ref_group = "Annotation",
  query_collapsing = T, ref_collapsing = T,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
d <- 1 - srt_22@tools$knnpredict_Annotation$distance_matrix
d <- as.matrix(d)

library(ComplexHeatmap)
Heatmap(d)


### Compare with human E14 --------------------------------------------------
srt_E14 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_3dculture/human_gastrulation_3dculture.seurat.rds")
srt_E14_sub <- srt_E14[, srt_E14$Group %in% c("EPI", "PSA-EPI", "PE")]
srt_E14_sub <- FindVariableFeatures(srt_E14_sub)
srt_E14_sub <- RunKNNMap(
  srt_query = srt_E14_sub, srt_ref = srt_22,
  ref_group = "CSSclusters", ref_umap = "CSSUMAP2D"
)
srt_E14_sub$Group <- factor(srt_E14_sub$Group, levels = c("EPI", "PSA-EPI", "PE"))
srt_E14_sub$cell_type <- srt_E14_sub$Group
srt_22$bg <- NA
p <- ProjectionPlot(
  srt_query = srt_E14_sub, srt_ref = srt_22,
  query_group = "cell_type", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "3D-cultured human pre-gastrulation embryos"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig2/22_E14_projection_celltype.umap.pdf", raster = TRUE)

p <- ProjectionPlot(
  srt_query = srt_E14_sub, srt_ref = srt_22,
  query_group = "Time", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "3D-cultured human pre-gastrulation embryos"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig2/22_E14_projection_time.umap.pdf", raster = TRUE)


# library(epiR)
# df <- table(srt_22$CSSclusters, srt_22$Time)
# df <- as.data.frame(cbind(celltype = rownames(df), df))
# for (i in 2:ncol(df)) {
#   df[, i] <- as.numeric(df[, i])
# }
# df[, "PreGas"] <- table(srt_E14_sub$predicted_ref_group)[rownames(df)]
# df[is.na(df)] <- 0
# df <- df[rowSums(df[, 2:ncol(df)]) != 0, ]
# stat <- epi.psi(df, itno = 99, conf.level = 0.95)
# stat <- stat[stat$v2 == "PreGas", ]
# stat$v1 <- factor(stat$v1, levels = levels(srt_22$Time))
# p <- ggplot(stat, aes(x = v1, y = est, fill = est)) +
#   geom_col(color = "black") +
#   geom_errorbar(aes(ymin = lower, ymax = upper),
#     width = .5,
#     position = position_dodge(.9)
#   ) +
#   scale_fill_gradientn(name = "PSI", colors = palette_scp(stat$est, palette = "material-deep-orange")) +
#   theme_scp() +
#   labs(title = "Similarity with pre-gastrulation embryos", x = "", y = "PSI") +
#   theme(
#     aspect.ratio = 0.8,
#     plot.title = element_text(size = 14, colour = "black", vjust = 1),
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
#   )
# p <- panel_fix(p, height = 2, save = "figures/fig2/22_E14_compare_PSI.stat.pdf")

### Compare with human CS7 --------------------------------------------------------
srt_CS7 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_E19/human_gastrulation_E19.seurat.rds")
srt_CS7 <- RunKNNMap(
  srt_query = srt_CS7, srt_ref = srt_22,
  ref_group = "Annotation", ref_umap = "CSSUMAP2D"
)
levels <- c(
  "Epiblast", "Primitive Streak", "Axial Mesoderm", "Nascent Mesoderm", "Emergent Mesoderm", "Advanced Mesoderm", "ExE Mesoderm", "Endoderm",
  "Non-Neural Ectoderm", "Hemogenic Endothelial Progenitors", "Erythroblasts"
)
srt_CS7$cluster_id <- factor(srt_CS7$cluster_id, levels = levels)
srt_CS7$cell_type <- srt_CS7$cluster_id

srt_22$bg <- NA
p <- ProjectionPlot(
  srt_query = srt_CS7, srt_ref = srt_22,
  query_group = "cell_type", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "CS7 human gastrula"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig2/22_CS7_projection_celltype.umap.pdf", raster = TRUE)

p <- ProjectionPlot(
  srt_query = srt_CS7, srt_ref = srt_22,
  query_group = "spatial", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "CS7 human gastrula"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig2/22_CS7_projection_spatial.umap.pdf", raster = TRUE)

# cepi <- colnames(srt_CS7)[srt_CS7$predicted_ref_group == "Caudal epiblast(CEpi)"]
# cepi <- table(srt_CS7$spatial[cepi])
# epi <- colnames(srt_CS7)[srt_CS7$predicted_ref_group == "Epiblast"]
# epi <- table(srt_CS7$spatial[epi])
# df <- cbind(epi, cepi)
# df <- t(t(df) / colSums(df))
# df <- reshape2::melt(df)
# colnames(df) <- c("Spatial", "CellType", "Proportion")
# p <- ggplot(df, aes(x = CellType, y = Proportion, fill = Spatial)) +
#   geom_col(color = "black") +
#   scale_fill_manual(values = palette_scp(df$Spatial)) +
#   theme_scp() +
#   theme(
#     aspect.ratio = 1,
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
#   )

# srt_22_sub <- srt_22[, colnames(srt_22)[srt_22$Time == "20d"]]
# srt_22_sub <- FindVariableFeatures(srt_22_sub)
# srt_22_sub <- RunKNNPredict(
#   srt_query = srt_22_sub, srt_ref = srt_CS7,
#   query_group = "Annotation", ref_group = "spatial",
#   features = intersect(VariableFeatures(srt_22_sub), VariableFeatures(srt_CS7)),
#   query_collapsing = FALSE, ref_collapsing = FALSE,
#   return_full_distance_matrix = TRUE
# )
# df <- table(srt_22_sub$Annotation, srt_22_sub$knnpredict_spatial)
# df <- df / rowSums(df)
# df <- df[rowSums(df, na.rm = TRUE) != 0, ]
# Heatmap(df)
# d <- as.matrix(1 - srt_22@tools$knnpredict_spatial$distance_matrix)
# ComplexHeatmap::Heatmap(d)
# table(srt_22$Annotation, srt_22$knnpredict_spatial)

library(epiR)
df <- table(srt_22$Annotation, srt_22$Time)
df <- as.data.frame(cbind(celltype = rownames(df), df))
for (i in 2:ncol(df)) {
  df[, i] <- as.numeric(df[, i])
}
df[, "CS7"] <- table(srt_CS7$predicted_ref_group)[rownames(df)]
df[is.na(df)] <- 0
stat <- epi.psi(df, itno = 99, conf.level = 0.95)
stat <- stat[stat$v2 == "CS7", ]
stat$v1 <- factor(stat$v1, levels = levels(srt_22$Time))
p <- ggplot(stat, aes(x = v1, y = est, fill = est)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
    width = .5,
    position = position_dodge(.9)
  ) +
  scale_fill_gradientn(name = "PSI", colors = palette_scp(stat$est, palette = "material-deep-orange")) +
  theme_scp() +
  labs(title = "Similarity with CS7 human gastrula", x = "", y = "PSI") +
  theme(
    aspect.ratio = 0.8,
    plot.title = element_text(size = 14, colour = "black", vjust = 1),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
p <- panel_fix(p, height = 2, save = "figures/fig2/22_CS7_compare_PSI.stat.pdf")

rownames(df) <- df[, 1]
df[, 1] <- NULL
df <- t(df) / colSums(df)
hc <- hclust(as.dist(dist(unclass(df))))
plot(hc)
levels_order <- hc$labels[hc$order]
df <- reshape2::melt(df)
colnames(df) <- c("CellType", "Time", "Proportion")
df$CellType <- factor(df$CellType, levels = levels_order)


# sp_gene <- db$Homo_sapiens$ProteinComplex$TERM2GENE[db$Homo_sapiens$ProteinComplex$TERM2NAME == "ComplexID:351", 2]
# sp_gene_map <- GeneConvert(sp_gene, geneID_from_IDtype = "entrez_id", geneID_to_IDtype = "symbol")
# srt_22[["percent.splice"]] <- PercentageFeatureSet(object = srt_22, features = sp_gene_map$geneID_expand$symbol)
# ExpDimPlot(srt, "percent.splice")

library(openxlsx)
df <- read.xlsx("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-04158-y/MediaObjects/41586_2021_4158_MOESM18_ESM.xlsx", startRow = 2)
srt_22 <- CellClassification(srt_22, features = list(adv_meso = df$Advanced.Mesoderm, exe_meso = df$EXE.Mesoderm), name = "Mesoderm")
ExpDimPlot(srt_22, "Mesoderm_adv_meso")
ExpDimPlot(srt_22, "Mesoderm_exe_meso")

df <- read.xlsx("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-04158-y/MediaObjects/41586_2021_4158_MOESM13_ESM.xlsx", startRow = 2)
srt_22 <- CellClassification(srt_22, features = list(DE1 = df$DE1, DE2 = df$DE2, Hypoblast = df$Hypoblast, YS = df$YS.Endoderm), name = "Endoderm")
ExpDimPlot(srt_22, "Endoderm_DE1")
ExpDimPlot(srt_22, "Endoderm_DE2")
ExpDimPlot(srt_22, "Endoderm_Hypoblast")
ExpDimPlot(srt_22, "Endoderm_YS")
ClassDimPlot(srt_22, "Endoderm_Classification")

df <- read.xlsx("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-04158-y/MediaObjects/41586_2021_4158_MOESM10_ESM.xlsx", startRow = 2)
srt_22 <- CellClassification(srt_22, features = list(Amnion = df$Amnion, NNE = df$NNE), name = "Amn")
ExpDimPlot(srt_22, "Amn_NNE")
ExpDimPlot(srt_22, "Amn_Amnion")
ClassDimPlot(srt_22, "Amn_Classification")

### CS7 two epiblast lineage development -------------------------------------
p <- ExpDimPlot(srt_CS7, "percent.ribo", "ref.umap",
  theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
)
p <- panel_fix(p, height = 2, save = "figures/fig2/22_CS7_ribo.umap.pdf")


p <- ClassDimPlot(srt_CS7, "cluster_id", label = TRUE, theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2") +
  ExpDimPlot(srt_CS7, c("percent.ribo"),
    theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2"
  )
p <- panel_fix(as_gtable(p), height = 3, save = "figures/fig2/22_CS7.newumap.pdf")

p <- ClassDimPlot(srt_CS7, "cluster_id", label = TRUE, reduction = "umap", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2") +
  ExpDimPlot(srt_CS7, c("percent.ribo"),
    reduction = "umap",
    theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2"
  )
p <- panel_fix(p, height = 3, save = "figures/fig2/22_CS7.oldumap.pdf")

adata <- srt_to_adata(srt_CS7)
adata <- RunPAGA(adata, group_by = "cluster_id", nonlinear_reduction = "umap")


library(ggrepel)
out <- proxyC::simil(x = t(srt_CS7[["percent.ribo"]]), y = srt_CS7@assays$RNA@data, use_nan = TRUE)
out[is.na(out)] <- 0
gene <- names(sort(out[1, out[1, ] > 0.9], decreasing = T))
gene_nonribo <- gene[!grepl("^RP[SL]\\d+\\w{0,1}\\d*$", gene)]
df <- data.frame(gene = gene, GeneType = "RP gene")
rownames(df) <- gene
df[gene_nonribo, "GeneType"] <- "non-RP gene"
df <- df %>%
  group_by(GeneType) %>%
  summarise(count = n()) %>%
  as.data.frame()
df$GeneType <- factor(df$GeneType, levels = c("RP gene", "non-RP gene"))
p <- ggplot(df, aes(x = 2, y = count, fill = GeneType)) +
  geom_col(color = "black") +
  geom_text_repel(aes(label = paste0(count, "(", round(count / sum(count) * 100, 1), "%)")),
    size = 4, color = "white", bg.color = "black", fontface = "bold",
    position = position_stack(vjust = 0.5)
  ) +
  xlim(0.5, 2.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = palette_scp(df$GeneType)) +
  theme_blank(add_coord = FALSE)
p <- panel_fix(p, height = 4, save = "figures/fig2/22_CS7_caudalgene.pie.pdf")

res <- RunEnrichment(gene)
p <- EnrichmentPlot(res$enrichment, plot_type = "bar")
p <- panel_fix(p, height = 2, save = "figures/fig2/22_CS7_caudalgene.GO.pdf")

srt_22 <- CellClassification(srt_22, features = list(Caudal = gene))
ExpDimPlot(srt_22, "Caudal_Score")


adata <- srt_to_adata(srt_CS7)
adata <- RunPAGA(
  adata = adata, group_by = "Standardclusters", nonlinear_reduction = "StandardUMAP2D",
  embedded_with_PAGA = T
)
plt <- reticulate::import("matplotlib")$pyplot
sc <- reticulate::import("scanpy")
ax <- sc$pl$paga_compare(adata,
  edges = FALSE,
  threshold = 0.1, max_edge_width = 0.8, basis = "umap", title = "PAGA",
  legend_loc = NULL, legend_fontweight = "bold", legend_fontsize = 10, legend_fontoutline = 1.5,
  frameon = FALSE, save = FALSE, show = FALSE
)
ax[[1]]$axis("none")
ax[[2]]$axis("equal") ## adjust window
plt$show()

### Compare with preimplantation embryo -------------------------------------
# srt_preimp <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_preimplantation_embryo/human_preimplantation_embryo.seurat.rds")
srt_preimp <- readRDS("/data/tmp/HuangMingQian/scRNA-seq/ref/Sophie_preimplantation_embryos/Sophie_preimplantation_embryos.seurat.rds")
# srt_preimp <- srt_preimp[, colnames(srt_preimp)[srt_preimp$CellType == "epiblast"]]
# srt_preimp <- Standard_SCP(srt_preimp)

srt_epiblast <- srt_22[, srt_22$Annotation == "Epiblast"]
srt_epiblast <- RunKNNMap(
  srt_query = srt_epiblast, srt_ref = srt_preimp,
  ref_group = "CellType", ref_umap = "StandardUMAP2D",
  distance_metric = "euclidean"
)
ProjectionPlot(
  srt_query = srt_epiblast, srt_ref = srt_preimp,
  query_group = "Time", ref_group = "CellType"
)
srt <- Integration_SCP(srtList = list(, srt_preimp), integration_method = "CSS")


srt_preimp <- RunKNNMap(
  srt_query = srt_preimp, srt_ref = srt_22,
  ref_group = "Annotation", ref_umap = "CSSUMAP2D"
)
ProjectionPlot(
  srt_query = srt_preimp, srt_ref = srt_22,
  query_group = "Time", ref_group = "bg",
  query_param = list(palette = "Set1"),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "pre-implantation embryo"
  )
)


### Compare with mouse gastrulation -----------------------------------------
srt_mGast <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/mouse_gastrulation_E8.5/mouse_gastrulation_E8.5.corrected.seurat.rds")
celltype_colours <- c(
  "Epiblast" = "#635547",
  "Primitive Streak" = "#DABE99",
  "Caudal epiblast" = "#9e6762",
  "PGC" = "#FACB12",
  "Anterior Primitive Streak" = "#c19f70",
  "Notochord" = "#0F4A9C",
  "Def. endoderm" = "#F397C0",
  "Gut" = "#EF5A9D",
  "Nascent mesoderm" = "#C594BF",
  "Mixed mesoderm" = "#DFCDE4",
  "Intermediate mesoderm" = "#139992",
  "Caudal Mesoderm" = "#3F84AA",
  "Paraxial mesoderm" = "#8DB5CE",
  "Somitic mesoderm" = "#005579",
  "Pharyngeal mesoderm" = "#C9EBFB",
  "Cardiomyocytes" = "#B51D8D",
  "Allantois" = "#532C8A",
  "ExE mesoderm" = "#8870ad",
  "Mesenchyme" = "#cc7818",
  "Haematoendothelial progenitors" = "#FBBE92",
  "Endothelium" = "#ff891c",
  "Blood progenitors 1" = "#f9decf",
  "Blood progenitors 2" = "#c9a997",
  "Erythroid1" = "#C72228",
  "Erythroid2" = "#f79083",
  "Erythroid3" = "#EF4E22",
  "NMP" = "#8EC792",
  "Rostral neurectoderm" = "#65A83E",
  "Caudal neurectoderm" = "#354E23",
  "Neural crest" = "#C3C388",
  "Forebrain/Midbrain/Hindbrain" = "#647a4f",
  "Spinal cord" = "#CDE088",
  "Surface ectoderm" = "#f7f79e",
  "Visceral endoderm" = "#F6BFCB",
  "ExE endoderm" = "#7F6874",
  "ExE ectoderm" = "#989898",
  "Parietal endoderm" = "#1A1A1A"
)
srt_mGast$celltype <- factor(srt_mGast$celltype, levels = names(celltype_colours))
p <- ClassDimPlot(srt_mGast,
  group.by = "celltype", reduction = "umap",
  pal_color = celltype_colours, label = TRUE, theme_use = "theme_blank",
  xlab = "UMAP_1", ylab = "UMAP_2"
)
p <- panel_fix(p, height = 4, save = "figures/fig2/ref_mGas.umap.pdf")

srt_22 <- RunKNNMap(
  srt_query = srt_22, srt_ref = srt_mGast,
  ref_group = "celltype", ref_umap = "umap"
)
srt_mGast <- RunKNNMap(
  srt_query = srt_mGast, srt_ref = srt_22,
  ref_group = "CSSclusters", ref_umap = "CSSUMAP2D"
  # features = intersect(VariableFeatures(srt_22), VariableFeatures(srt_E14_sub))
)
srt_22[["ref.umap.mgas"]] <- srt_22[["ref.umap"]]
srt_22[["predicted.mgas"]] <- srt_22[["predicted_ref_group"]]
srt_22[["predicted.mgas"]] <- factor(srt_22[["predicted.mgas", drop = TRUE]], levels = levels(srt_mGast$celltype))

srt_mGast$bg <- NA
srt_22$bg <- srt_22$Annotation
srt_22$bg[srt_22$bg != "Epiblast"] <- NA
srt_22$bg <- "Our study"
p <- ProjectionPlot(
  srt_query = srt_22, srt_ref = srt_mGast,
  query_group = "bg", ref_group = "bg", pt.size = 0.5,
  query_param = list(pal_color = "red"),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "Mapping to the mouse gastrulation data"
  )
)
p <- panel_fix(p, height = 4, save = "figures/fig2/22_map_to_mGas.umap.pdf")

plist <- list()
for (i in levels(srt_22$Time)) {
  srt_sub <- srt_22[, srt_22$Time == i]
  srt_sub$Time <- i
  plist[[i]] <- ProjectionPlot(
    srt_query = srt_sub, srt_ref = srt_mGast,
    query_group = "Time", ref_group = "bg", pt.size = 0.5,
    query_param = list(pal_color = "red", legend.position = "none"),
    ref_param = list(
      bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
      title = i
    )
  )
}
p <- plot_grid(plotlist = plist, nrow = 2)
p <- panel_fix(p, height = 3, save = "figures/fig2/22_map_to_mGas.TimeSplit.umap.pdf")


library(epiR)
df_base <- table(srt_mGast$celltype, srt_mGast$stage)
df_base <- as.data.frame(cbind(celltype = rownames(df_base), df_base))
for (i in 2:ncol(df_base)) {
  df_base[, i] <- as.numeric(df_base[, i])
}

df <- df_base
df[, "GO_embyro"] <- table(srt_22$predicted.mgas)[rownames(df)]
df[is.na(df)] <- 0
df <- df[rowSums(df[, 2:ncol(df)]) != 0, ]
stat <- epi.psi(df, itno = 99, conf.level = 0.95)
stat <- stat[stat$v2 == "GO_embyro", ]
stat$v1 <- factor(stat$v1, levels = c("E6.5", "E6.75", "E7.0", "E7.25", "E7.5", "E7.75", "E8.0", "E8.25", "E8.5", "mixed_gastrulation"))
p <- ggplot(stat, aes(x = v1, y = est, fill = est)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
    width = .5,
    position = position_dodge(.9)
  ) +
  scale_fill_gradientn(name = "PSI", colors = palette_scp(stat$est, palette = "material-deep-orange")) +
  theme_scp() +
  labs(title = "Similarity with mouse gastrulation embryo", x = "", y = "PSI") +
  theme(
    aspect.ratio = 0.8,
    plot.title = element_text(size = 14, colour = "black", vjust = 1),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
p <- panel_fix(p, height = 2, save = "figures/fig2/22_map_to_mGas.stat.pdf")


statlist <- list()
plist <- list()
for (i in levels(srt_22$Time)) {
  df <- df_base
  df[, "GO_embyro"] <- table(srt_22$predicted.mgas[srt_22$Time == i])[rownames(df)]
  df[is.na(df)] <- 0
  df <- df[rowSums(df[, 2:ncol(df)]) != 0, ]
  stat <- epi.psi(df, itno = 99, conf.level = 0.95)
  stat <- stat[stat$v2 == "GO_embyro", ]
  stat$v1 <- factor(stat$v1, levels = c("E6.5", "E6.75", "E7.0", "E7.25", "E7.5", "E7.75", "E8.0", "E8.25", "E8.5", "mixed_gastrulation"))
  stat$compare <- i
  statlist[[i]] <- stat
  plist[[i]] <- ggplot(stat, aes(x = v1, y = est, fill = est)) +
    geom_col(color = "black") +
    geom_errorbar(aes(ymin = lower, ymax = upper),
      width = .5,
      position = position_dodge(.9)
    ) +
    scale_fill_gradientn(name = "PSI", colors = palette_scp(stat$est, palette = "material-deep-orange")) +
    theme_scp() +
    labs(title = i, x = "", y = "PSI") +
    theme(
      aspect.ratio = 0.8,
      plot.title = element_text(size = 14, colour = "black", vjust = 1),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )
}
p <- plot_grid(plotlist = plist, nrow = 2)
p <- panel_fix(p, height = 2, save = "figures/fig2/22_map_to_mGas.TimeSplit.stat.pdf")

df <- do.call(rbind.data.frame, statlist)
mat <- reshape2::dcast(data = df, v1 ~ compare, value.var = "est")
rownames(mat) <- mat$v1
mat <- as.matrix(mat[, -1])
ht <- Heatmap(mat,
  name = "Estimated PSI", cluster_columns = F, cluster_rows = FALSE,
  split = rownames(mat) == "mixed_gastrulation", row_title = " ",
  col = colorRamp2(seq(0, max(mat), length.out = 3), c("#27408B", "white", "#EE0000")),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = "white", lwd = 1, fill = "white")
    )
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = fill, lwd = 1, fill = alpha(fill, 0.5))
    )
    if (mat[i, j] >= sort(mat[i, ], decreasing = T)[3]) {
      # grid.text("*", x, y, gp = gpar(fontsize = 20))
      grid.text(round(mat[i, j], 2), x, y, gp = gpar(fontsize = 10))
    }
  },
  border = TRUE,
  width = unit(ncol(mat) * 0.8, "cm"),
  height = unit(nrow(mat) * 0.8, "cm"),
  heatmap_legend_param = list(border = TRUE)
)
gTree <- grid.grabExpr({
  draw(ht,
    heatmap_legend_side = "right",
    padding = unit(c(3, 1, 1, 3), "cm")
  ) # bottom, left, top and right
})
p <- plot_grid(gTree)
ggsave(plot = p, filename = "figures/fig2/22_map_to_mGas.TimeSplit.heatmap.pdf", width = 7, height = 10)


saveRDS(srt_22, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_22_annotation.rds")
saveRDS(srt_mGast, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/mouse_gastrulation_E8.5/mouse_gastrulation_E8.5.corrected.seurat.rds")


### Compare with mouse organogenesis -----------------------------------------
srt_mOrg <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/mouse_organogenesis_E13.5_atlas/mouse_organogenesis_E13.5_atlas.seurat.total.de.rds")
srt_22 <- RunKNNPredict(
  srt_query = srt_22, srt_ref = srt_mOrg, features_type = "DE",
  query_group = "CSSclusters", ref_group = "Main_cell_type"
)

ClassDimPlot(srt_22, "knnpredict_Main_cell_type")


# merge 22cellline and h0 (different passage)-----------------------------------------------------------------------------
srt_22cellline1 <- LoadH5Seurat("/ssd/lab/HuangMingQian/scRNAseq/CellCulture/NGSmodule_SCP_analysis/CellQC/22P12.filtered.h5Seurat")
srt_22cellline2 <- LoadH5Seurat("/ssd/lab/HuangMingQian/scRNAseq/CellCulture-project/NGSmodule_SCP_analysis/CellQC/22_CellLine-3.filtered.h5Seurat")
srt_22cellline3 <- LoadH5Seurat("/ssd/lab/HuangMingQian/scRNAseq/CellCulture/NGSmodule_SCP_analysis/CellQC/22P38.filtered.h5Seurat")
srt_h0cellline1 <- LoadH5Seurat("/ssd/lab/HuangMingQian/scRNAseq/CellCulture/NGSmodule_SCP_analysis/CellQC/H0P11.filtered.h5Seurat")
srt_h0cellline2 <- LoadH5Seurat("/ssd/lab/ChenMengQi/scRNAseq/in_vitro_culture/analysis_zh/NGSmodule_SCP_analysis/CellQC/C380h.filtered.h5Seurat")
srt_h0cellline3 <- LoadH5Seurat("/ssd/lab/HuangMingQian/scRNAseq/CellCulture/NGSmodule_SCP_analysis/CellQC/H0P30.filtered.h5Seurat")

srt_22cellline1[["ribo.mito.ratio"]] <- srt_22cellline1$percent.ribo / srt_22cellline1$percent.mito
srt_22cellline1 <- srt_22cellline1[, colnames(srt_22cellline1)[srt_22cellline1$ribo.mito.ratio > 1]]
srt_22cellline2[["ribo.mito.ratio"]] <- srt_22cellline2$percent.ribo / srt_22cellline2$percent.mito
srt_22cellline2 <- srt_22cellline2[, colnames(srt_22cellline2)[srt_22cellline2$ribo.mito.ratio > 1]]
srt_22cellline3[["ribo.mito.ratio"]] <- srt_22cellline3$percent.ribo / srt_22cellline3$percent.mito
srt_22cellline3 <- srt_22cellline3[, colnames(srt_22cellline3)[srt_22cellline3$ribo.mito.ratio > 1]]

srt_h0cellline1[["ribo.mito.ratio"]] <- srt_h0cellline1$percent.ribo / srt_h0cellline1$percent.mito
srt_h0cellline1 <- srt_h0cellline1[, colnames(srt_h0cellline1)[srt_h0cellline1$ribo.mito.ratio > 1]]
srt_h0cellline2[["ribo.mito.ratio"]] <- srt_h0cellline2$percent.ribo / srt_h0cellline2$percent.mito
srt_h0cellline2 <- srt_h0cellline2[, colnames(srt_h0cellline2)[srt_h0cellline2$ribo.mito.ratio > 1]]
srt_h0cellline3[["ribo.mito.ratio"]] <- srt_h0cellline3$percent.ribo / srt_h0cellline3$percent.mito
srt_h0cellline3 <- srt_h0cellline3[, colnames(srt_h0cellline3)[srt_h0cellline3$ribo.mito.ratio > 1]]

srt_22cellline1[["Source"]] <- srt_22cellline2[["Source"]] <- srt_22cellline3[["Source"]] <- "22cellline"
srt_h0cellline1[["Source"]] <- srt_h0cellline2[["Source"]] <- srt_h0cellline3[["Source"]] <- "h0cellline"
srt_22cellline1[["Time"]] <- srt_h0cellline1[["Time"]] <- "Early"
srt_22cellline2[["Time"]] <- srt_h0cellline2[["Time"]] <- "Middle"
srt_22cellline3[["Time"]] <- srt_h0cellline3[["Time"]] <- "Late"

srt_cellline <- Reduce(merge, list(
  srt_22cellline1, srt_22cellline2, srt_22cellline3,
  srt_h0cellline1, srt_h0cellline2, srt_h0cellline3
))
srt_cellline[["Source_Time"]] <- paste0(srt_cellline[["Source", drop = TRUE]], "-", srt_cellline[["Time", drop = TRUE]])
srt_cellline <- srt_cellline[, colnames(srt_cellline)[srt_cellline$nFeature_RNA < 8000]]

srt_cellline <- Integration_SCP(
  srt_cellline,
  batch = "Source_Time",
  integration_method = "Seurat", cluster_resolution = 2
)

ClassDimPlot(srt_cellline, c("Source_Time", "Source", "Time"))
ExpDimPlot(srt_cellline, c("percent.mito", "percent.ribo", "ribo.mito.ratio"))
ExpDimPlot(srt_cellline, c(
  "POU5F1", "NANOG", "SOX2", "DNMT3B", "CDH1", "CDH2", "EOMES", "TBXT",
  "TFAP2A", "ISL1", "SOX17", "FOXA2", "TFAP2C", "NANOS3", "HAND1", "BMP4"
), theme_use = "theme_blank") %>% panel_fix(height = 2)

ClassDimPlot(srt_cellline, "Seuratclusters", label = TRUE)
nameslist <- list(
  "Epiblast(PSA)" = c(1, 2, 3, 4, 11, 12, 13, 14, 15, 16, 17, 18),
  "Epiblast(EMT)" = list(6, 7, 8, 9, 10),
  "Nascent Mesoderm" = list(20),
  "Advanced Mesoderm" = list(21, 23),
  "Endoderm" = list(24, 25, 26, 27),
  "Amnion" = list(5, 22),
  "PGC" = list(19)
)
srt_cellline <- RenameClusters(srt_cellline,
  group.by = "Seuratclusters",
  name = "Annotation",
  nameslist = nameslist
)
srt_cellline$Annotation <- as.character(srt_cellline$Annotation)
srt_cellline$Annotation <- factor(srt_cellline$Annotation, levels = names(nameslist))
srt_cellline$Time <- factor(srt_cellline$Time, levels = c("Early", "Late"))
srt_cellline$TimeAnnotation <- factor(paste0(srt_cellline$Time, "-", srt_cellline$Annotation),
  levels = paste0(levels(srt_cellline$Time), "-", rep(levels(srt_cellline$Annotation), each = 2))
)
srt_cellline$Source <- factor(srt_cellline$Source, levels = c("22cellline", "h0cellline"))
srt_cellline$SourceAnnotation <- factor(paste0(srt_cellline$Source, "-", srt_cellline$Annotation),
  levels = paste0(levels(srt_cellline$Source), "-", rep(levels(srt_cellline$Annotation), each = 2))
)
srt_cellline[["UMAP"]] <- srt_cellline[["SeuratUMAP2D"]]
srt_cellline@misc$Default_reduction <- "UMAP"
ClassDimPlot(srt_cellline, c("Annotation"), label = T, theme_use = "theme_blank") %>%
  panel_fix(height = 3, raster = T, save = "figures_add/cellline.umap.pdf")
ClassDimPlot(srt_cellline, "Time", split.by = "Source", theme_use = "theme_blank") %>%
  panel_fix(height = 3, raster = T, save = "figures_add/cellline.split.pdf")
ExpDimPlot(srt_cellline, c(
  "POU5F1", "NANOG", "SOX2", "DNMT3B", "CDH1", "CDH2", "EOMES", "TBXT",
  "TFAP2A", "ISL1", "SOX17", "FOXA2", "TFAP2C", "NANOS3", "HAND1", "BMP4"
), theme_use = "theme_blank") %>%
  panel_fix(height = 2, raster = T, save = "figures_add/cellline.markers.pdf")


ClassStatPlot(srt_cellline,
  group.by = "Time", stat.by = "Annotation",
  plot_type = "bar", split.by = "Source"
) %>%
  panel_fix(height = 3, width = 3, save = "figures_add/cellline.bar1.pdf")
ClassStatPlot(srt_cellline,
  group.by = "Annotation", stat.by = "Time",
  plot_type = "bar", position = "dodge", split.by = "Source"
) %>%
  panel_fix(height = 3, width = 3, save = "figures_add/cellline.bar2.pdf")

CellCorHeatmap(
  srt_query = srt_cellline[, srt_cellline$Source == "22cellline"],
  srt_ref = srt_cellline[, srt_cellline$Source == "22cellline"],
  query_group = "TimeAnnotation", ref_group = "TimeAnnotation",
  cluster_columns = F, cluster_rows = F, nlabel = 4
) %>% panel_fix(height = 7, width = 7, save = "figures_add/cellline_22.cor.pdf")
CellCorHeatmap(
  srt_query = srt_cellline[, srt_cellline$Source == "h0cellline"],
  srt_ref = srt_cellline[, srt_cellline$Source == "h0cellline"],
  query_group = "TimeAnnotation", ref_group = "TimeAnnotation",
  cluster_columns = F, cluster_rows = F, nlabel = 4
) %>% panel_fix(height = 7, width = 7, save = "figures_add/cellline_h0.cor.pdf")
CellCorHeatmap(
  srt_query = srt_cellline,
  srt_ref = srt_cellline,
  query_group = "SourceAnnotation", ref_group = "SourceAnnotation",
  cluster_columns = F, cluster_rows = F, nlabel = 4
) %>% panel_fix(height = 7, width = 7, save = "figures_add/cellline.cor.pdf")
saveRDS(srt_cellline, "srt_cellline.rds")

srt_E14 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_3dculture/human_gastrulation_3dculture.seurat.rds")
srt_E14$orig.ident <- "3D cultured"
srt_E14_epi <- srt_E14[, srt_E14$Group %in% c("EPI", "PSA-EPI")]
srt_E14_epi$Annotation <- NA
srt_E14_epi$Annotation[srt_E14_epi$Group == "EPI"] <- "EPI(3D cultured)"
srt_E14_epi$Annotation[srt_E14_epi$Group == "PSA-EPI"] <- "PSA-EPI(3D cultured)"
srt_E14_epi$Annotation <- factor(srt_E14_epi$Annotation, levels = c("EPI(3D cultured)", "PSA-EPI(3D cultured)"))
srt_E14_epi$TimeAnnotation <- factor(paste0(srt_E14_epi$Time, "-", srt_E14_epi$Annotation),
  levels = c(
    "D6-EPI(3D cultured)", "D7-EPI(3D cultured)", "D8-EPI(3D cultured)", "D9-EPI(3D cultured)",
    "D10-EPI(3D cultured)", "D12-EPI(3D cultured)", "D14-EPI(3D cultured)", "D14-PSA-EPI(3D cultured)"
  )
)
srt_E14_epi <- RunDEtest(srt_E14_epi, group_by = "TimeAnnotation")
degs <- srt_E14_epi@tools$DEtest_TimeAnnotation$AllMarkers_wilcox[with(srt_E14_epi@tools$DEtest_TimeAnnotation$AllMarkers_wilcox, p_val_adj < 0.05), ]
srt_E14_epi <- Standard_SCP(srt_E14_epi, HVF = degs$gene)
panel_fix(ClassDimPlot(srt_E14_epi, c("Annotation", "Time")), height = 2)

srt <- merge(srt_cellline, srt_E14_epi)
srt$Annotation <- factor(srt$Annotation, levels = c(levels(srt_E14_epi$Annotation), levels(srt_cellline$Annotation)))
srt$Time <- factor(srt$Time, levels = c(levels(srt_E14_epi$Time), levels(srt_cellline$Time)))
srt$TimeAnnotation <- factor(srt$TimeAnnotation, levels = c(levels(srt_E14_epi$TimeAnnotation), levels(srt_cellline$TimeAnnotation)))

srt$TimeAnnotation2 <- as.character(srt$TimeAnnotation)
srt$TimeAnnotation2[which(srt$Source == "22cellline")] <- "22cellline"
srt$TimeAnnotation2[which(srt$Source == "h0cellline")] <- "h0cellline"
srt$TimeAnnotation2 <- factor(srt$TimeAnnotation2, levels = c(levels(srt_E14_epi$TimeAnnotation), "22cellline", "h0cellline"))
ht <- ExpHeatmap(srt,
  features = degs$gene, group.by = "TimeAnnotation2",
  width = 7, height = 5
)
panel_fix(ht$plot, height = 5, width = 10, raster = T, save = "figures_add/cellline_e14.heatmap.pdf")

CellCorHeatmap(
  srt_query = srt, srt_ref = srt, features = degs$gene,
  query_group = "TimeAnnotation2", ref_group = "TimeAnnotation2",
  cluster_columns = FALSE, cluster_rows = FALSE
) %>% panel_fix(height = 7, width = 7, save = "figures_add/cellline_e14.cor.pdf")

# EMT analysis ----------------------------------------------------------------------------------------
# db <- PrepareDB(species = "Homo_sapiens", db = "GO_BP")
# EMT_gene <- filter(db$Homo_sapiens$GO_BP$TERM2GENE, Term == "GO:0001837")[["symbol"]]
# # EMT_gene <- intersect(EMT_gene, VariableFeatures(srt_cellline))
# MET_gene <- filter(db$Homo_sapiens$GO_BP$TERM2GENE, Term == "GO:0060231")[["symbol"]]
# # MET_gene <- intersect(MET_gene, VariableFeatures(srt_cellline))
# srt_cellline <- CellScoring(srt_cellline, features = list(EMT = EMT_gene, MET = MET_gene), name = "Transition")
# ClassDimPlot(srt_cellline, "Transition_classification") %>% panel_fix(height = 3)
# ExpDimPlot(srt_cellline, c("Transition_EMT", "Transition_MET"))
# GroupHeatmap(srt_cellline,
#              features = c(EMT_gene, MET_gene), feature_split = c(rep("EMT", length(EMT_gene)), rep("MET", length(MET_gene))),
#              group.by = "Annotation", show_row_names = TRUE, nlabel = 0, add_dot = F, cluster_rows = T,
# )

CDH_gene <- c("CDH1", "CDH2", "CDH3", "CDH4", "CDH8", "CDH10", "CDH11", "CDH12", "CDH13", "CDH20")
ExpDimPlot(srt_cellline, CDH_gene, theme_use = "theme_blank") %>%
  panel_fix(height = 2, raster = T, save = "figures_add/cellline.cdh.pdf")

db <- PrepareDB(species = "Homo_sapiens", db = "GO_BP", db_version = "3.13")
pEMT_gene <- filter(db$Homo_sapiens$GO_BP$TERM2GENE, Term == "GO:0010718")[["symbol"]]
pEMT_gene <- pEMT_gene[pEMT_gene %in% VariableFeatures(srt_cellline)]
nEMT_gene <- filter(db$Homo_sapiens$GO_BP$TERM2GENE, Term == "GO:0010719")[["symbol"]]
nEMT_gene <- nEMT_gene[nEMT_gene %in% VariableFeatures(srt_cellline)]
srt_cellline <- CellScoring(srt_cellline, features = list(pEMT = pEMT_gene, nEMT = nEMT_gene), name = "Transition")
srt_cellline@meta.data[, "Transition_classification"] <- factor(srt_cellline@meta.data[, "Transition_classification"], levels = c("pEMT", "nEMT"))
srt_cellline$Transition <- srt_cellline$Transition_classification
ClassDimPlot(srt_cellline, "Transition_classification", palette = "Set1", theme_use = "theme_blank") %>%
  panel_fix(height = 3, raster = T, save = "figures_add/emt_class.pdf")
ClassDimPlot(srt_cellline, "Annotation", stat.by = "Transition_classification", theme_use = "theme_blank") %>%
  panel_fix(height = 3, raster = T, save = "figures_add/emt_stat.pdf")
ExpDimPlot(srt_cellline, c("Transition_pEMT", "Transition_nEMT"), theme_use = "theme_blank") %>%
  panel_fix(height = 3, raster = T, save = "figures_add/emt_score.pdf")

ExpHeatmap(srt_cellline,
  group.by = "Annotation",
  features = c(CDH_gene, pEMT_gene, nEMT_gene),
  feature_split = c(rep("CDH", length(CDH_gene)), rep("pEMT", length(pEMT_gene)), rep("nEMT", length(nEMT_gene))),
  nlabel = 0, cluster_rows = T, cluster_columns = T, show_row_names = TRUE,
  width = 7, height = 7
)$plot %>%
  panel_fix(height = 7, width = 9, raster = T, save = "figures_add/emt.expheatmap.pdf")
GroupHeatmap(srt_cellline,
  group.by = "Annotation", heatmap_palette = "YlOrRd",
  features = c(CDH_gene, pEMT_gene, nEMT_gene),
  feature_split = c(rep("CDH", length(CDH_gene)), rep("pEMT", length(pEMT_gene)), rep("nEMT", length(nEMT_gene))),
  # cell_annotation = c("Transition"), cell_annotation_params = list(height = grid::unit(0.5, "in")),
  cluster_rows = T, show_row_names = TRUE, show_column_names = F, nlabel = 0, add_dot = F
)$plot %>%
  panel_fix(height = 7, width = 9, raster = T, save = "figures_add/emt.groupheatmap.pdf")



# Compare with CS7 ------------------------------------------------------------------------------------
srt_CS7 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_E19/human_gastrulation_E19.seurat.rds")
# srt_CS7 <- NormalizeData(srt_CS7)
srt_CS7$orig.ident <- "CS7"
srt_CS7$Time <- "CS7"
srt_CS7$Source <- "CS7"
srt_CS7$Annotation <- factor(srt_CS7$cluster_id,
  levels = c(
    "Epiblast", "Primitive Streak",
    "Nascent Mesoderm", "Emergent Mesoderm", "Advanced Mesoderm",
    "ExE Mesoderm", "Hemogenic Endothelial Progenitors", "Erythroblasts",
    "Axial Mesoderm", "Endoderm", "Non-Neural Ectoderm"
  )
)
srt_CS7$SubAnnotation <- factor(srt_CS7$sub_cluster,
  levels = c(
    "Epiblast", "Primitive Streak", "Nascent Mesoderm", "Emergent Mesoderm", "Advanced Mesoderm",
    "YS Mesoderm", "Hemogenic Endothelium", "Erythro-Myeloid Progenitors", "Myeloid Progenitors", "Blood Progenitors", "Erythroblasts",
    "Axial Mesoderm", "DE(P)", "DE(NP)", "Hypoblast", "YS Endoderm", "Non-Neural Ectoderm", "PGC"
  )
)
srt_CS7$TimeAnnotation <- factor(paste0("CS7-", srt_CS7$Annotation), levels = paste0("CS7-", levels(srt_CS7$Annotation)))
srt_CS7$TimeSubAnnotation <- factor(paste0("CS7-", srt_CS7$SubAnnotation), levels = paste0("CS7-", levels(srt_CS7$SubAnnotation)))
srt_CS7_epi <- srt_CS7[, colnames(srt_CS7)[srt_CS7$Annotation == "Epiblast"]]
srt_CS7 <- RenameClusters(srt_CS7,
  group.by = "Annotation", name = "Annotation_rename",
  nameslist = list(
    "Nascent/Emergent Mesoderm" = c("Nascent Mesoderm", "Emergent Mesoderm"),
    "Advanced/ExE Mesoderm" = c("Advanced Mesoderm", "ExE Mesoderm"),
    "Amnion" = "Non-Neural Ectoderm"
  )
)
srt_CS7$Annotation_rename <- as.character(srt_CS7$Annotation_rename)
srt_CS7$Annotation_rename[srt_CS7$Annotation %in% c("Axial Mesoderm", "Hemogenic Endothelial Progenitors", "Erythroblasts")] <- NA
srt_CS7$Annotation_rename[srt_CS7$SubAnnotation == "PGC"] <- "PGC"
srt_CS7$Annotation_rename <- factor(srt_CS7$Annotation_rename,
  levels = c("Epiblast", "Primitive Streak", "Nascent/Emergent Mesoderm", "Advanced/ExE Mesoderm", "Endoderm", "Amnion", "PGC")
)
ClassDimPlot(srt_CS7, c("Annotation", "Annotation_rename"))

# srt_cellline <- readRDS("srt_cellline.rds")
srt <- Integration_SCP(
  srtList = list(srt_cellline, srt_CS7),
  integration_method = "Seurat"
)
HVF <- VariableFeatures(srt)
srt$Annotation <- factor(srt$Annotation,
  levels = c(
    "Epiblast(PSA)", "Epiblast(EMT)", "Epiblast", "Primitive Streak",
    "Nascent Mesoderm", "Emergent Mesoderm", "Advanced Mesoderm",
    "ExE Mesoderm", "Hemogenic Endothelial Progenitors", "Erythroblasts",
    "Axial Mesoderm", "Endoderm",
    "Amnion", "Non-Neural Ectoderm", "PGC"
  )
)
srt$Annotation_rename <- factor(srt$Annotation_rename, levels = levels(srt_CS7$Annotation_rename))
srt_cellline <- srt_cellline[, srt_cellline$Time != "Middle"]

ClassDimPlot(srt, "Annotation", split.by = "Source", ncol = 1, theme_use = "theme_blank") %>%
  panel_fix(height = 2, raster = T, save = "figures_add/cellline_cs7.umap.pdf")


srt[["CelllLine22"]] <- srt[["CelllLineh0"]] <- srt[["CS7"]] <- NA
cells22 <- colnames(srt_cellline)[srt_cellline$Source == "22cellline"]
cellsh0 <- colnames(srt_cellline)[srt_cellline$Source == "h0cellline"]
srt@meta.data[cells22, "CelllLine22"] <- as.character(srt_cellline$Annotation[cells22])
srt$CelllLine22 <- factor(srt$CelllLine22, levels = levels(srt_cellline$Annotation))
srt@meta.data[cellsh0, "CelllLineh0"] <- as.character(srt_cellline$Annotation[cellsh0])
srt$CelllLineh0 <- factor(srt$CelllLineh0, levels = levels(srt_cellline$Annotation))
srt@meta.data[colnames(srt_CS7), "CS7"] <- as.character(srt_CS7$Annotation_rename)
srt$CS7 <- factor(srt$CS7, levels = levels(srt_CS7$Annotation_rename))

srt_cellline <- RunDEtest(srt_cellline, group_by = "Annotation")
srt_CS7 <- RunDEtest(srt_CS7, group_by = "Annotation_rename")
HVF2 <- Reduce(intersect, list(
  srt_cellline@tools$DEtest_Annotation$AllMarkers_wilcox[srt_cellline@tools$DEtest_Annotation$AllMarkers_wilcox$p_val_adj < 0.05, "gene"],
  srt_CS7@tools$DEtest_Annotation_rename$AllMarkers_wilcox[srt_CS7@tools$DEtest_Annotation_rename$AllMarkers_wilcox$p_val_adj < 0.05, "gene"]
))

GroupHeatmap(srt,
  features = HVF2, group.by = c("CelllLine22", "CelllLineh0", "CS7"),
  group_palette = c("Paired", "Paired", "Paired"),
  cell_annotation = c("CDH1", "CDH2"), slot = "data", lib_normalize = FALSE,
  cluster_rows = TRUE, cluster_features_by = c("CelllLine22", "CelllLineh0", "CS7"),
  show_row_names = FALSE, nlabel = 20,
  height = 8, width = 8, heatmap_palette = "YlOrRd"
)$plot %>%
  panel_fix(height = 7, width = 9, raster = T, save = "figures_add/cellline_cs7.groupheatmap.pdf")

srt$TimeAnnotation2 <- factor(srt$SourceAnnotation,
  levels = c(paste0(c("22cellline-", "h0cellline-"), rep(levels(srt_cellline$Annotation), each = 2)))
)
CellCorHeatmap(
  srt_query = srt[, colnames(srt_cellline)],
  srt_ref = srt[, colnames(srt_CS7)],
  features = HVF2, nlabel = 2, label_by = "row",
  query_group = "TimeAnnotation2", ref_group = "Annotation_rename",
  cluster_columns = F, cluster_rows = F
) %>%
  panel_fix(height = 7, width = 7, raster = T, save = "figures_add/cellline_cs7.cor.pdf")

# Testis-CS7 integration ---------------------------------------------------------------------------------
srt_22$Source <- "22cellline"
srt_22$TimeAnnotation <- paste0("22Testis-", srt_22$Time, "-", srt_22$Annotation)
srt_22$TimeAnnotation <- factor(srt_22$TimeAnnotation,
  levels = paste0("22Testis-", levels(srt_22$Time), "-", rep(levels(srt_22$Annotation), each = nlevels(srt_22$Time)))
)

srt_22testis_d10 <- srt_22[, colnames(srt_22)[srt_22$Time == "d10"]]
srt_22testis_d10$Time <- "22Testis-d10"
srt_22testis_d10$TimeAnnotation <- factor(paste0("22Testis-d10-", srt_22testis_d10$Annotation),
  levels = paste0("22Testis-d10-", levels(srt_22testis_d10$Annotation))
)
srt_22testis_d20 <- srt_22[, colnames(srt_22)[srt_22$Time == "d20"]]
srt_22testis_d20$Time <- "22Testis-d20"
srt_22testis_d20$TimeAnnotation <- factor(paste0("22Testis-d20-", srt_22testis_d20$Annotation),
  levels = paste0("22Testis-d20-", levels(srt_22testis_d20$Annotation))
)
srt_22testis_d30 <- srt_22[, colnames(srt_22)[srt_22$Time == "d30"]]
srt_22testis_d30$Time <- "22Testis-d30"
srt_22testis_d30$TimeAnnotation <- factor(paste0("22Testis-d30-", srt_22testis_d30$Annotation),
  levels = paste0("22Testis-d30-", levels(srt_22testis_d30$Annotation))
)

srt_h0$Source <- "h0cellline"
srt_h0 <- RenameClusters(srt_h0, group.by = "Annotation", nameslist = list("Amnion" = "Amnion & PGC"), name = "Annotation", keep_levels = TRUE)
srt_h0$TimeAnnotation <- paste0("h0Testis-", srt_h0$Time, "-", srt_h0$Annotation)
srt_h0$TimeAnnotation <- factor(srt_h0$TimeAnnotation,
  levels = paste0("h0Testis-", levels(srt_h0$Time), "-", rep(levels(srt_h0$Annotation), each = nlevels(srt_h0$Time)))
)

srt_h0testis_d10 <- srt_h0[, colnames(srt_h0)[srt_h0$Time == "d10"]]
srt_h0testis_d10$Time <- "h0Testis-d10"
srt_h0testis_d10$TimeAnnotation <- factor(paste0("h0Testis-d10-", srt_h0testis_d10$Annotation),
  levels = paste0("h0Testis-d10-", levels(srt_h0testis_d10$Annotation))
)
srt_h0testis_d20 <- srt_h0[, colnames(srt_h0)[srt_h0$Time == "d20"]]
srt_h0testis_d20$Time <- "h0Testis-d20"
srt_h0testis_d20$TimeAnnotation <- factor(paste0("h0Testis-d20-", srt_h0testis_d20$Annotation),
  levels = paste0("h0Testis-d20-", levels(srt_h0testis_d20$Annotation))
)
srt_h0testis_d30 <- srt_h0[, colnames(srt_h0)[srt_h0$Time == "d30"]]
srt_h0testis_d30$Time <- "h0Testis-d30"
srt_h0testis_d30$TimeAnnotation <- factor(paste0("h0Testis-d30-", srt_h0testis_d30$Annotation),
  levels = paste0("h0Testis-d30-", levels(srt_h0testis_d30$Annotation))
)

srt_CS7$Source <- "CS7"
srt_CS7 <- RenameClusters(srt_CS7,
  group.by = "cluster_id",
  nameslist = list(
    "Epiblast" = "Epiblast",
    "Primitive streak" = "Primitive Streak",
    "Mesoderm" = c(
      "Nascent Mesoderm", "Emergent Mesoderm", "Advanced Mesoderm", "ExE Mesoderm",
      "Axial Mesoderm", "Hemogenic Endothelial Progenitors", "Erythroblasts"
    ),
    "Endoderm" = "Endoderm",
    "Non-neural ectoderm" = "Non-Neural Ectoderm"
  ),
  name = "GermLayer"
)
srt_CS7$GermLayer <- as.character(srt_CS7$GermLayer)
srt_CS7$GermLayer[srt_CS7$sub_cluster == "PGC"] <- "PGC"

HVF <- VariableFeatures(srt_CS7)
srt <- Integration_SCP(
  srtList = list(
    srt_22testis_d10,
    srt_22testis_d20,
    srt_22testis_d30,
    srt_h0testis_d10,
    srt_h0testis_d20,
    srt_h0testis_d30,
    srt_CS7
  ),
  HVF = HVF,
  integration_method = "CSS"
)
levels(srt_22$Annotation)
levels(srt_h0$Annotation)
annotation_intersect <- intersect(levels(srt_22$Annotation), levels(srt_h0$Annotation))
annotation_union <- unique(c(levels(srt_22$Annotation), levels(srt_h0$Annotation)))
annotation_diff <- setdiff(annotation_union, annotation_intersect)
srt <- RenameClusters(srt, group.by = "Annotation", nameslist = list("Amnion" = "Amnion & PGC"), name = "Annotation2")
srt$Annotation2 <- factor(srt$Annotation2, levels = unique(c(levels(srt_CS7$Annotation), annotation_union[annotation_union != "Amnion & PGC"])))
srt_testis_cs7 <- srt
saveRDS(srt_testis_cs7, "srt_testis_cs7.rds")

nameslist <- list(
  "Epiblast" = "Epiblast",
  "Primitive streak" = c("Primitive streak", "Primitive Streak"),
  "Nascent/Emergent mesoderm" = c("Mesoderm", "Nascent Mesoderm", "Emergent Mesoderm"),
  "Advanced/ExE mesoderm" = c("Advanced Mesoderm", "ExE Mesoderm", "MSC/Fib", "ExE mesoderm", "Limb bud mesenchyme cell"),
  "Endothelial & erythroid cell" = c("Endothelial & erythroid cell", "Hemogenic Endothelial Progenitors", "Erythroblasts"),
  "Endoderm" = c("Endoderm", "ExE endoderm", "Gut"),
  "Amnion& PGC" = c("Amnion & PGC", "Amnion", "PGC", "Epithelium", "Non-Neural Ectoderm"),
  "Neural ectoderm" = c("Neural ectoderm", "Ependymal cell", "Neuron", "Radial glial cell", "Retinal pigmented epithelium", "Retinal progenitor cell", "Schwann cell", "Sensory neuron"),
  "Axial Mesoderm" = "Axial Mesoderm"
)
srt <- RenameClusters(srt,
  group.by = "Annotation", name = "Annotation_rough",
  nameslist = nameslist
)
srt$Annotation_rough <- factor(srt$Annotation_rough, levels = names(nameslist))
srt$Time <- factor(srt$Time, levels = c("22Testis-d10", "22Testis-d20", "22Testis-d30", "h0Testis-d10", "h0Testis-d20", "h0Testis-d30", "CS7"))
srt$TimeAnnotation_rough <- factor(paste0(srt$Time, "-", srt$Annotation_rough), levels = paste0(levels(srt$Time), "-", rep(names(nameslist), each = nlevels(srt$Time))))
unique(srt$TimeAnnotation_rough)


ClassDimPlot(srt, "Time", cells.highlight = WhichCells(srt, expression = Time == "CS7"), theme_use = "theme_blank") %>%
  panel_fix(height = 3, raster = T, save = "figures_add/testis_cs7.umap.pdf")

ClassDimPlot(srt, "Annotation_rough", split.by = "Time", cells.highlight = T, sizes.highlight = 0.5, theme_use = "theme_blank") %>%
  panel_fix(height = 2, raster = T, save = "figures_add/testis_cs7.split.pdf")

srt$GermLayer[srt$Annotation == "PGC"] <- "PGC"
srt$GermLayer[which(srt$sub_cluster == "PGC")] <- "PGC"
srt$GermLayer <- factor(srt$GermLayer, levels = c("Epiblast", "Primitive streak", "Mesoderm", "Endoderm", "Neural ectoderm", "Non-neural ectoderm", "PGC"))
ClassDimPlot(srt, "GermLayer", theme_use = "theme_blank", label = T) %>%
  panel_fix(height = 3, raster = T, save = "figures_add/testis_cs7.germlayer.pdf")
ClassDimPlot(srt, "GermLayer", split.by = "Time", cells.highlight = T, sizes.highlight = 0.5, theme_use = "theme_blank") %>%
  panel_fix(height = 2, raster = T, save = "figures_add/testis_cs7.split.pdf")

ClassStatPlot(srt, group.by = "Time", stat.by = "GermLayer", plot_type = "bar", stat_type = "count", position = "dodge") %>% panel_fix(height = 3, width = 5)
ClassStatPlot(srt, group.by = "Time", stat.by = "GermLayer", plot_type = "bar", stat_type = "percent", position = "dodge") %>% panel_fix(height = 3, width = 5)
srt$GermLayer2 <- srt$GermLayer
srt$GermLayer2[srt$GermLayer2 == "Epiblast"] <- NA
ClassStatPlot(srt, group.by = "Time", stat.by = "GermLayer2", plot_type = "bar", stat_type = "count", position = "dodge", NA_stat = FALSE) %>% panel_fix(height = 3, width = 5)
ClassStatPlot(srt, group.by = "Time", stat.by = "GermLayer2", plot_type = "bar", stat_type = "percent", position = "dodge", NA_stat = FALSE) %>% panel_fix(height = 3, width = 5)


srt[["CS7"]] <- srt[["Testis_22_d10"]] <- srt[["Testis_h0_d10"]] <-
  srt[["Testis_22_d20"]] <- srt[["Testis_h0_d20"]] <-
  srt[["Testis_22_d30"]] <- srt[["Testis_h0_d30"]] <- NA
cells_CS7 <- colnames(srt_CS7)
srt@meta.data[cells_CS7, "CS7"] <- "CS7"
cells_22_d10 <- colnames(srt_22testis_d10)
srt@meta.data[cells_22_d10, "Testis_22_d10"] <- "Testis_22_d10"
cells_22_d20 <- colnames(srt_22testis_d20)
srt@meta.data[cells_22_d20, "Testis_22_d20"] <- "Testis_22_d20"
cells_22_d30 <- colnames(srt_22testis_d30)
srt@meta.data[cells_22_d30, "Testis_22_d30"] <- "Testis_22_d30"
cells_h0_d10 <- colnames(srt_h0testis_d10)
srt@meta.data[cells_h0_d10, "Testis_h0_d10"] <- "Testis_h0_d10"
cells_h0_d20 <- colnames(srt_h0testis_d20)
srt@meta.data[cells_h0_d20, "Testis_h0_d20"] <- "Testis_h0_d20"
cells_h0_d30 <- colnames(srt_h0testis_d30)
srt@meta.data[cells_h0_d30, "Testis_h0_d30"] <- "Testis_h0_d30"
groups <- c("CS7", "Testis_22_d10", "Testis_22_d20", "Testis_22_d30", "Testis_h0_d10", "Testis_h0_d20", "Testis_h0_d30")
ht <- GroupHeatmap(srt,
  features = VariableFeatures(srt), split.by = "Annotation_rough",
  group.by = groups,
  group_palcolor = as.list(palette_scp(groups, palette = "Set1")),
  cell_split_palette = "Paired", heatmap_palette = "YlOrRd",
  cell_annotation = c("CDH1", "CDH2"), slot = "data", lib_normalize = FALSE,
  cluster_rows = TRUE,
  cluster_features_by = groups,
  show_row_names = FALSE, nlabel = 20,
  height = 8, width = 10
)
panel_fix(ht$plot, height = 7, width = 14, raster = T, save = "figures_add/testis_cs7.heatmap.pdf")

CellCorHeatmap(
  srt_query = srt,
  srt_ref = srt,
  features = VariableFeatures(srt), nlabel = 3, label_by = "column",
  query_group = "TimeAnnotation_rough", ref_group = "TimeAnnotation_rough",
  cluster_columns = F, cluster_rows = F,
  columns = grep("CS7", x = levels(srt$TimeAnnotation_rough), invert = T, value = T),
  rows = grep("CS7", x = levels(srt$TimeAnnotation_rough), value = T),
  # label_cutoff = 0.7,
  column_names_rot = 45
) %>%
  panel_fix(height = 5, width = 20, raster = T, save = "figures_add/testis_cs7.cor.pdf")


# 22CellLine-22Testis-E14-CS7 integration ---------------------------------------------------------------------------------
HVF <- VariableFeatures(srt_CS7)
srt <- Integration_SCP(
  srtList = list(
    srt_E14_epi,
    srt_22cellline,
    srt_22testis_d10,
    srt_22testis_d20,
    srt_22testis_d30,
    srt_CS7
  ),
  HVF = HVF,
  integration_method = "CSS"
)
srt_22cellline_testis_e14_cs7 <- srt
saveRDS(srt_22cellline_testis_e14_cs7, "srt_22cellline_testis_e14_cs7.rds")
panel_fix(ClassDimPlot(srt, "Time",
  cells.highlight = WhichCells(srt, expression = Time == "D14")
), height = 3)
panel_fix(ClassDimPlot(srt, "Annotation", split.by = "orig.ident"), height = 3)
ClassDimPlot3D(srt, "Time", cells.highlight = WhichCells(srt, expression = Time == "D14"))

srt[["cultrued"]] <- srt[["CS7"]] <- srt[["CelllLine22"]] <-
  srt[["TestisD10"]] <- srt[["TestisD30"]] <- srt[["TestisD20"]] <- NA
srt@meta.data[colnames(srt_E14_epi), "cultrued"] <- as.character(srt_E14_epi$TimeAnnotation)
srt$cultrued <- factor(srt$cultrued, levels = levels(srt_E14_epi$TimeAnnotation))
srt@meta.data[colnames(srt_CS7), "CS7"] <- as.character(srt_CS7$Annotation)
srt$CS7 <- factor(srt$CS7, levels = levels(srt_CS7$Annotation))
srt@meta.data[colnames(srt_22cellline), "CelllLine22"] <- as.character(srt_22cellline$Annotation)
srt$CelllLine22 <- factor(srt$CelllLine22, levels = levels(srt_22cellline$Annotation))
srt@meta.data[colnames(srt_22testis_d10), "TestisD10"] <- as.character(srt_22testis_d10$Annotation)
srt$TestisD10 <- factor(srt$TestisD10, levels = levels(srt_22testis_d10$Annotation))
srt@meta.data[colnames(srt_22testis_d20), "TestisD20"] <- as.character(srt_22testis_d20$Annotation)
srt$TestisD20 <- factor(srt$TestisD20, levels = levels(srt_22testis_d20$Annotation))
srt@meta.data[colnames(srt_22testis_d30), "TestisD30"] <- as.character(srt_22testis_d30$Annotation)
srt$TestisD30 <- factor(srt$TestisD30, levels = levels(srt_22testis_d30$Annotation))
srt@meta.data[, "Ratio"] <- srt@assays$RNA@data["CDH2", ] / srt@assays$RNA@data["CDH1", ]

GroupHeatmap(srt,
  features = VariableFeatures(srt),
  group.by = c("cultrued", "CS7", "CelllLine22", "TestisD10", "TestisD20", "TestisD30"),
  group_palette = c("jco", "Set1", "Paired", "jama", "jama", "jama"),
  cell_annotation = c("CDH1", "CDH2"),
  cluster_rows = TRUE,
  cluster_features_by = c("cultrued", "CS7", "CelllLine22", "TestisD10", "TestisD20", "TestisD30"),
  show_row_names = FALSE, nlabel = 20,
  height = 8, width = 8
)
CellCorHeatmap(
  srt_query = srt,
  srt_ref = srt,
  features = HVF, nlabel = 3, label_by = "row",
  query_group = "TimeAnnotation", ref_group = "TimeAnnotation",
  cluster_columns = F, cluster_rows = F,
  columns = c(
    paste0(c("22CellLine-"), rep(levels(srt_22cellline$Annotation), each = 1)),
    paste0(c("Testis-d10", "Testis-d20", "Testis-d30"), rep(levels(srt_22$Annotation), each = 3))
  ),
  rows = c(levels(srt_E14_epi$TimeAnnotation), levels(srt_CS7$TimeAnnotation)),
  label_cutoff = 0.7,
  column_names_rot = 45,
  gird_size = unit(7, "mm")
)

# h0Testis-E14-CS7 integration ---------------------------------------------------------------------------------
srt_h0testis_d10 <- srt_h0[, colnames(srt_h0)[srt_h0$Time == "d10"]]
srt_h0testis_d10$Time <- "H0Testis-d10"
srt_h0testis_d10$TimeAnnotation <- factor(paste0("H0Testis-d10", srt_h0testis_d10$Annotation),
  levels = paste0("H0Testis-d10", levels(srt_h0testis_d10$Annotation))
)
srt_h0testis_d20 <- srt_h0[, colnames(srt_h0)[srt_h0$Time == "d20"]]
srt_h0testis_d20$Time <- "H0Testis-d20"
srt_h0testis_d20$TimeAnnotation <- factor(paste0("H0Testis-d20", srt_h0testis_d20$Annotation),
  levels = paste0("H0Testis-d20", levels(srt_h0testis_d20$Annotation))
)
srt_h0testis_d30 <- srt_h0[, colnames(srt_h0)[srt_h0$Time == "d30"]]
srt_h0testis_d30$Time <- "H0Testis-d30"
srt_h0testis_d30$TimeAnnotation <- factor(paste0("H0Testis-d30", srt_h0testis_d30$Annotation),
  levels = paste0("H0Testis-d30", levels(srt_h0testis_d30$Annotation))
)


HVF <- VariableFeatures(srt_CS7)
srt <- Integration_SCP(
  srtList = list(
    srt_E14_epi,
    srt_h0testis_d10,
    srt_h0testis_d20,
    srt_h0testis_d30,
    srt_CS7
  ),
  HVF = HVF,
  integration_method = "CSS",
  linear_reduction_dims_use = 1:25
)
srt_h0testis_e14_cs7 <- srt
saveRDS(srt_h0testis_e14_cs7, "srt_h0testis_e14_cs7.rds")
panel_fix(ClassDimPlot(srt, "Time",
  cells.highlight = WhichCells(srt, expression = Time == "D14")
), height = 3)
panel_fix(ClassDimPlot(srt, "Annotation", split.by = "orig.ident"), height = 3)
ClassDimPlot3D(srt, "Annotation")

srt[["cultrued"]] <- srt[["CS7"]] <-
  srt[["TestisD10"]] <- srt[["TestisD30"]] <- srt[["TestisD20"]] <- NA
srt@meta.data[colnames(srt_E14_epi), "cultrued"] <- as.character(srt_E14_epi$TimeAnnotation)
srt$cultrued <- factor(srt$cultrued, levels = levels(srt_E14_epi$TimeAnnotation))
srt@meta.data[colnames(srt_CS7), "CS7"] <- as.character(srt_CS7$Annotation)
srt$CS7 <- factor(srt$CS7, levels = levels(srt_CS7$Annotation))
srt@meta.data[colnames(srt_h0testis_d10), "TestisD10"] <- as.character(srt_h0testis_d10$Annotation)
srt$TestisD10 <- factor(srt$TestisD10, levels = levels(srt_h0testis_d10$Annotation))
srt@meta.data[colnames(srt_h0testis_d20), "TestisD20"] <- as.character(srt_h0testis_d20$Annotation)
srt$TestisD20 <- factor(srt$TestisD20, levels = levels(srt_h0testis_d20$Annotation))
srt@meta.data[colnames(srt_h0testis_d30), "TestisD30"] <- as.character(srt_h0testis_d30$Annotation)
srt$TestisD30 <- factor(srt$TestisD30, levels = levels(srt_h0testis_d30$Annotation))
srt@meta.data[, "Ratio"] <- srt@assays$RNA@data["CDH2", ] / srt@assays$RNA@data["CDH1", ]

GroupHeatmap(srt,
  features = VariableFeatures(srt),
  group.by = c("cultrued", "CS7", "TestisD10", "TestisD20", "TestisD30"),
  group_palette = c("jco", "Set1", "Paired", "Paired", "Paired"),
  cell_annotation = c("CDH1", "CDH2"),
  cluster_rows = TRUE,
  cluster_features_by = c("cultrued", "CS7", "TestisD10", "TestisD20", "TestisD30"),
  show_row_names = FALSE, nlabel = 20,
  height = 8, width = 8
)
CellCorHeatmap(
  srt_query = srt,
  srt_ref = srt,
  features = HVF, nlabel = 3, label_by = "column",
  query_group = "TimeAnnotation", ref_group = "TimeAnnotation",
  cluster_columns = F, cluster_rows = F,
  columns = paste0(c("Testis-d10", "Testis-d20", "Testis-d30"), rep(levels(srt_h0$Annotation), each = 3)),
  rows = c(levels(srt_E14_epi$TimeAnnotation), levels(srt_CS7$TimeAnnotation)),
  label_cutoff = 0.7,
  column_names_rot = 45
)

# h0CellLine-h0Testis-E14-CS7 integration ---------------------------------------------------------------------------------
HVF <- VariableFeatures(srt_CS7)
srt <- Integration_SCP(
  srtList = list(
    srt_E14_epi,
    srt_h0cellline,
    srt_testis_d10,
    srt_testis_d20,
    srt_testis_d30,
    srt_CS7
  ),
  HVF = HVF,
  integration_method = "CSS"
)
srt_h0cellline_testis_e14_cs7 <- srt
saveRDS(srt_h0cellline_testis_e14_cs7, "srt_h0cellline_testis_e14_cs7.rds")
panel_fix(ClassDimPlot(srt, "Time"), height = 3)
panel_fix(ClassDimPlot(srt, "Annotation", split.by = "orig.ident"), height = 3)
ClassDimPlot3D(srt, "Time", cells.highlight = WhichCells(srt, expression = Time == "D14"))

srt[["cultrued"]] <- srt[["CS7"]] <- srt[["CelllLineh0"]] <-
  srt[["TestisD10"]] <- srt[["TestisD30"]] <- srt[["TestisD20"]] <- NA
srt@meta.data[colnames(srt_E14_epi), "cultrued"] <- as.character(srt_E14_epi$TimeAnnotation)
srt$cultrued <- factor(srt$cultrued, levels = levels(srt_E14_epi$TimeAnnotation))
srt@meta.data[colnames(srt_CS7), "CS7"] <- as.character(srt_CS7$Annotation)
srt$CS7 <- factor(srt$CS7, levels = levels(srt_CS7$Annotation))
srt@meta.data[colnames(srt_h0cellline), "CelllLineh0"] <- as.character(srt_h0cellline$Annotation)
srt$CelllLineh0 <- factor(srt$CelllLineh0, levels = levels(srt_h0cellline$Annotation))
srt@meta.data[colnames(srt_h0testis_d10), "TestisD10"] <- as.character(srt_h0testis_d10$Annotation)
srt$TestisD10 <- factor(srt$TestisD10, levels = levels(srt_h0testis_d10$Annotation))
srt@meta.data[colnames(srt_h0testis_d20), "TestisD20"] <- as.character(srt_h0testis_d20$Annotation)
srt$TestisD20 <- factor(srt$TestisD20, levels = levels(srt_h0testis_d20$Annotation))
srt@meta.data[colnames(srt_h0testis_d30), "TestisD30"] <- as.character(srt_h0testis_d30$Annotation)
srt$TestisD30 <- factor(srt$TestisD30, levels = levels(srt_h0testis_d30$Annotation))
srt@meta.data[, "Ratio"] <- srt@assays$RNA@data["CDH2", ] / srt@assays$RNA@data["CDH1", ]

GroupHeatmap(srt,
  features = VariableFeatures(srt),
  group.by = c("cultrued", "CS7", "CelllLineh0", "TestisD10", "TestisD20", "TestisD30"),
  group_palette = c("jco", "Set1", "Paired", "jama", "jama", "jama"),
  cell_annotation = c("CDH1", "CDH2"),
  cluster_rows = TRUE,
  cluster_features_by = c("cultrued", "CS7", "CelllLineh0", "TestisD10", "TestisD20", "TestisD30"),
  show_row_names = FALSE, nlabel = 20,
  height = 8, width = 8
)
CellCorHeatmap(
  srt_query = srt,
  srt_ref = srt,
  features = HVF, nlabel = 3, label_by = "column",
  query_group = "TimeAnnotation", ref_group = "TimeAnnotation",
  cluster_columns = F, cluster_rows = F,
  columns = c(
    paste0(c("h0CellLine-"), rep(levels(srt_h0cellline$Annotation), each = 1)),
    paste0(c("Testis-d10", "Testis-d20", "Testis-d30"), rep(levels(srt_h0$Annotation), each = 3))
  ),
  rows = c(levels(srt_E14_epi$TimeAnnotation), levels(srt_CS7$TimeAnnotation)),
  label_cutoff = 0.7,
  column_names_rot = 45,
  gird_size = unit(7, "mm")
)



# E14Epi_CS7Epi_22CellLlineEpi_22TestisEpi ------------------------------------------------------------------------------------
srt_E14 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_3dculture/human_gastrulation_3dculture.seurat.rds")
srt_E14$Project <- "3D cultured"
srt_E14_epi <- srt_E14[, colnames(srt_E14)[srt_E14$Group == "EPI"]]
srt_E14_epi <- srt_E14[, srt_E14$Group %in% c("EPI", "PSA-EPI")]
srt_E14_epi$Annotation <- NA
srt_E14_epi$Annotation[srt_E14_epi$Group == "EPI"] <- "EPI(3D cultured)"
srt_E14_epi$Annotation[srt_E14_epi$Group == "PSA-EPI"] <- "PSA-EPI(3D cultured)"
srt_E14_epi$Annotation <- factor(srt_E14_epi$Annotation, levels = c("EPI(3D cultured)", "PSA-EPI(3D cultured)"))
srt_E14_epi$TimeAnnotation <- factor(paste0(srt_E14_epi$Time, "-", srt_E14_epi$Annotation),
  levels = c(
    "D6-EPI(3D cultured)", "D7-EPI(3D cultured)", "D8-EPI(3D cultured)", "D9-EPI(3D cultured)",
    "D10-EPI(3D cultured)", "D12-EPI(3D cultured)", "D14-EPI(3D cultured)", "D14-PSA-EPI(3D cultured)"
  )
)
srt_E14_epi <- RunDEtest(srt_E14_epi, group_by = "TimeAnnotation")
E14epi_degs <- srt_E14_epi@tools$DEtest_TimeAnnotation$AllMarkers_wilcox[with(srt_E14_epi@tools$DEtest_TimeAnnotation$AllMarkers_wilcox, p_val_adj < 0.05), ]

srt_CS7 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_E19/human_gastrulation_E19.seurat.rds")
srt_CS7 <- NormalizeData(srt_CS7)
srt_CS7$Project <- "CS7"
srt_CS7$Time <- "CS7"
srt_CS7$Annotation <- factor(srt_CS7$cluster_id,
  levels = c(
    "Epiblast", "Primitive Streak", "Axial Mesoderm", "Endoderm",
    "Nascent Mesoderm", "Emergent Mesoderm", "Advanced Mesoderm",
    "ExE Mesoderm", "Hemogenic Endothelial Progenitors", "Erythroblasts",
    "Non-Neural Ectoderm"
  )
)
srt_CS7$SubAnnotation <- factor(srt_CS7$sub_cluster,
  levels = c(
    "Epiblast", "Primitive Streak", "Nascent Mesoderm", "Emergent Mesoderm", "Advanced Mesoderm",
    "YS Mesoderm", "Hemogenic Endothelium", "Erythro-Myeloid Progenitors", "Myeloid Progenitors", "Blood Progenitors", "Erythroblasts",
    "Axial Mesoderm", "DE(P)", "DE(NP)", "Hypoblast", "YS Endoderm", "Non-Neural Ectoderm", "PGC"
  )
)
srt_CS7$TimeAnnotation <- factor(paste0("CS7-", srt_CS7$Annotation), levels = paste0("CS7-", levels(srt_CS7$Annotation)))
srt_CS7$TimeSubAnnotation <- factor(paste0("CS7-", srt_CS7$SubAnnotation), levels = paste0("CS7-", levels(srt_CS7$SubAnnotation)))
srt_CS7 <- RunDEtest(srt_CS7, group_by = "TimeAnnotation")
CS7epi_degs <- srt_CS7@tools$DEtest_TimeAnnotation$AllMarkers_wilcox[with(srt_CS7@tools$DEtest_TimeAnnotation$AllMarkers_wilcox, p_val_adj < 0.05 & group1 %in% c("CS7-Epiblast", "CS7-Primitive Streak")), ]
srt_CS7_epi <- srt_CS7[, colnames(srt_CS7)[srt_CS7$Annotation %in% c("Epiblast", "Primitive Streak")]]

srt_22cellline_epi <- srt_22cellline[, colnames(srt_22cellline)[srt_22cellline$Annotation %in% c("Epiblast(PSA)", "Epiblast(EMT)")]]
srt_22cellline_epi$Project <- "22CellLine"

srt_testis_d10_epi <- srt_22[, colnames(srt_22)[srt_22$Annotation %in% c("Epiblast", "Primitive streak") & srt_22$Time == "d10"]]
srt_testis_d10_epi$Time <- "Testis-d10"
srt_testis_d10_epi$TimeAnnotation <- "Testis-d10-Epiblast"
srt_testis_d20_epi <- srt_22[, colnames(srt_22)[srt_22$Annotation %in% c("Epiblast", "Primitive streak") & srt_22$Time == "d20"]]
srt_testis_d20_epi$Time <- "Testis-d20"
srt_testis_d20_epi$TimeAnnotation <- factor(paste0("Testis-d20-", srt_testis_d20_epi$Annotation), levels = c("Testis-d20-Epiblast", "Testis-d20-Primitive streak"))
srt_testis_d30_epi <- srt_22[, colnames(srt_22)[srt_22$Annotation %in% c("Epiblast", "Primitive streak") & srt_22$Time == "d30"]]
srt_testis_d30_epi$Time <- "Testis-d30"
srt_testis_d30_epi$TimeAnnotation <- factor(paste0("Testis-d30-", srt_testis_d30_epi$Annotation), levels = c("Testis-d30-Epiblast", "Testis-d30-Primitive streak"))
srt_testis_d10_epi$Project <- srt_testis_d20_epi$Project <- srt_testis_d30_epi$Project <- "22Testis"

srt_epi_list <- list(
  srt_E14_epi, srt_CS7_epi, srt_22cellline_epi,
  srt_testis_d10_epi, srt_testis_d20_epi, srt_testis_d30_epi
)
srtMerge <- check_srtList(srt_epi_list, batch = "Project")
srt_compare <- Reduce(merge, srtMerge$srtList)
srt_compare$TimeAnnotation <- factor(srt_compare$TimeAnnotation,
  levels = c(
    levels(srt_E14_epi$TimeAnnotation),
    c("CS7-Epiblast", "CS7-Primitive Streak"),
    c("22CellLine-Epiblast(PSA)", "22CellLine-Epiblast(EMT)"),
    c("Testis-d10-Epiblast", "Testis-d20-Epiblast", "Testis-d30-Epiblast", "Testis-d20-Primitive streak", "Testis-d30-Primitive streak")
  )
)
ht <- GroupHeatmap(srt_compare,
  features = E14epi_degs$gene, group.by = "TimeAnnotation",
  slot = "data", # cluster_rows = T, cluster_columns = T,
  width = 7, height = 7, show_row_names = F,
  show_column_names = T, column_names_side = "bottom", ht_params = list(column_names_rot = 45)
)
ht$plot
CellCorHeatmap(
  srt_query = srt_compare, srt_ref = srt_compare, features = E14epi_degs$gene,
  query_group = "TimeAnnotation", ref_group = "TimeAnnotation",
  cluster_columns = T, cluster_rows = T
)
ht <- GroupHeatmap(srt_compare,
  features = CS7epi_degs$gene, group.by = "TimeAnnotation",
  slot = "data", cluster_rows = T, cluster_columns = T,
  width = 7, height = 7, show_row_names = F,
  show_column_names = T, column_names_side = "bottom", ht_params = list(column_names_rot = 45)
)
ht$plot
CellCorHeatmap(
  srt_query = srt_compare, srt_ref = srt_compare,
  features = CS7epi_degs$gene,
  query_group = "TimeAnnotation", ref_group = "TimeAnnotation",
  cluster_columns = T, cluster_rows = T
)


# DEGs between CS7 epi and d10-testis epi  -------------------------------------------------------------------
candidates_1 <- names(rowSums(srt_CS7@assays$RNA@counts > 0) > 30)
candidates_2 <- names(rowSums(srt_E14@assays$RNA@counts > 0) > 30)
candidates_3 <- names(rowSums(srt_22@assays$RNA@counts > 0) > 30)

srt_compare <- RunDEtest(srt_compare,
  group_by = "TimeAnnotation",
  features = Reduce(intersect, list(candidates_1, candidates_2, candidates_3))
)
epi_degs <- srt_compare@tools$DEtest_TimeAnnotation$AllMarkers_wilcox[with(srt_compare@tools$DEtest_TimeAnnotation$AllMarkers_wilcox, p_val_adj < 0.05), ]
# epi_degs <- epi_degs[epi_degs$group1 %in% c("D14-PSA-EPI(3D cultured)", "CS7-Epiblast", "Testis-d10-Epiblast", "Testis-d20-Epiblast", "Testis-d30-Epiblast"), ]
srt_compare <- AnnotateFeatures(srt_compare, db = "TF")
ht <- GroupHeatmap(srt_compare,
  features = epi_degs$gene, group.by = "TimeAnnotation",
  slot = "data", exp_method = "zscore", cluster_rows = T, n_split = 6,
  anno_terms = T, anno_features = T, feature_annotation = "TF",
  width = 7, height = 7, show_row_names = F,
  show_column_names = T, column_names_side = "bottom", ht_params = list(column_names_rot = 45)
)
ht$plot
ExpStatPlot(srt_compare, "UTF1", group.by = "Time")



# CS7_testisd10_testisd20_testisd30 -------------------------------------------------------------------
srt_CS7 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_E19/human_gastrulation_E19.seurat.rds")
srt_CS7$Time <- "CS7"
srt_CS7$Annotation <- srt_CS7$cluster_id
srt_CS7$TimeAnnotation <- paste0(srt_CS7$Time, "-", srt_CS7$Annotation)

HVF <- intersect(VariableFeatures(srt_CS7), VariableFeatures(srt_22))
srt_d10 <- srt_22[, colnames(srt_22)[srt_22$Time == "d10"]]
srt_d10$TimeAnnotation <- paste0(srt_d10$Time, "-", srt_d10$Annotation)
srt_d20 <- srt_22[, colnames(srt_22)[srt_22$Time == "d20"]]
srt_d20$TimeAnnotation <- paste0(srt_d20$Time, "-", srt_d20$Annotation)
srt_d30 <- srt_22[, colnames(srt_22)[srt_22$Time == "d30"]]
srt_d30$TimeAnnotation <- paste0(srt_d30$Time, "-", srt_d30$Annotation)
srtList <- list(srt_d10, srt_d20, srt_d30)
srt_compare <- Reduce(merge, srtList)
srt_compare <- FindVariableFeatures(srt_compare)
srt_compare <- RunDEtest(srt_compare, group_by = "Annotation")
srt_CS7 <- RunDEtest(srt_CS7, group_by = "Annotation")
srt_22cellline <- RunDEtest(srt_22cellline, group_by = "Annotation")
srt_22cellline$TimeAnnotation <- paste0("22CellLine-", srt_22cellline$Annotation)
degs1 <- srt_compare@tools$DEtest_Annotation$AllMarkers_wilcox
degs2 <- srt_CS7@tools$DEtest_Annotation$AllMarkers_wilcox
degs3 <- srt_22cellline@tools$DEtest_Annotation$AllMarkers_wilcox
features <- Reduce(intersect, list(
  degs1[with(degs1, p_val_adj < 0.05 & avg_log2FC > 1), "gene"],
  degs2[with(degs2, p_val_adj < 0.05 & avg_log2FC > 1), "gene"],
  degs3[with(degs3, p_val_adj < 0.05 & avg_log2FC > 1), "gene"]
))
features <- intersect(VariableFeatures(srt_compare), VariableFeatures(srt_CS7))

srt_Merge <- Reduce(merge, list(srt_22cellline, srt_compare))
srt_Merge <- srt_22cellline
srt_Merge <- RunKNNPredict(
  srt_query = srt_Merge, query_group = "TimeAnnotation",
  srt_ref = srt_CS7, ref_group = "TimeAnnotation",
  query_collapsing = T, ref_collapsing = T,
  features = features,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
d <- 1 - srt_Merge@tools$knnpredict$distance_matrix
d <- as.matrix(t(d))
levels1 <- c(
  c("22CellLine-Epiblast", "22CellLine-Epiblast(EMT)", "22CellLine-Nascent mesoderm", "22CellLine-Advanced Mesoderm", "22CellLine-Endoderm", "22CellLine-Amnion"),
  paste0(
    c("d10-", "d20-", "d30-"),
    rep(c(
      "Epiblast", "Primitive streak", "Mesoderm", "MSC/Fib", "ExE mesoderm", "Limb bud mesenchyme cell", "Endothelial & erythroid cell",
      "Endoderm", "ExE endoderm", "Gut", "Amnion", "PGC", "Epithelium",
      "Neural ectoderm", "Radial glial cell", "Ependymal cell", "Neuron", "Retinal progenitor cell", "Retinal pigmented epithelium"
    ), each = 3)
  )
)
levels2 <- c(
  "CS7-Epiblast", "CS7-Primitive Streak", "CS7-Nascent Mesoderm", "CS7-Emergent Mesoderm", "CS7-Advanced Mesoderm", "CS7-ExE Mesoderm",
  "CS7-Hemogenic Endothelial Progenitors", "CS7-Erythroblasts",
  "CS7-Axial Mesoderm", "CS7-Endoderm",
  "CS7-Non-Neural Ectoderm"
)
d <- d[levels1[levels1 %in% rownames(d)], levels2[levels2 %in% colnames(d)]]
ht <- Heatmap(d,
  name = "Cosine similarity", cluster_columns = F, cluster_rows = F,
  col = colorRamp2(seq(min(d), max(d), length.out = 3), c("#27408B", "white", "#EE0000")),
  cell_fun = function(j, i, x, y, w, h, fill) {
    # grid.rect(x, y,
    #   width = w, height = h,
    #   gp = gpar(col = "white", lwd = 1, fill = "white")
    # )
    # grid.rect(x, y,
    #   width = w, height = h,
    #   gp = gpar(col = fill, lwd = 1, fill = alpha(fill, 0.5))
    # )
    if (d[i, j] >= sort(d[i, ], decreasing = T)[5]) {
      # grid.text("*", x, y, gp = gpar(fontsize = 20))
      # grid.text(round(d[i, j], 2), x, y, gp = gpar(fontsize = 10))
    }
  },
  border = TRUE,
  # width = unit(ncol(d) * 0.8, "cm"),
  # height = unit(nrow(d) * 0.8, "cm"),
)
gTree <- grid.grabExpr({
  draw(ht,
    padding = unit(c(3, 1, 1, 3), "cm")
  ) # bottom, left, top and right
})
p <- plot_grid(gTree)
p

srt <- merge(srt_Merge, srt_CS7)
GroupHeatmap(srt,
  slot = "data", exp_method = "zscore",
  features = features,
  group.by = "TimeAnnotation", nlabel = 0, cluster_rows = T, cluster_columns = T,
  show_row_names = FALSE, show_column_names = TRUE
)






# CS7 -------------------------------------------------------------------------------------------------
srt_CS7 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_E19/human_gastrulation_E19.seurat.rds")
srt_CS7 <- RunKNNMap(
  srt_query = srt_CS7, srt_ref = srt_22,
  ref_group = "Annotation", ref_umap = "CSSUMAP2D"
)
levels <- c(
  "Epiblast", "Primitive Streak", "Axial Mesoderm", "Nascent Mesoderm", "Emergent Mesoderm", "Advanced Mesoderm", "ExE Mesoderm", "Endoderm",
  "Non-Neural Ectoderm", "Hemogenic Endothelial Progenitors", "Erythroblasts"
)
srt_CS7$cluster_id <- factor(srt_CS7$cluster_id, levels = levels)
srt_CS7$cell_type <- srt_CS7$cluster_id
ProjectionPlot(
  srt_query = srt_CS7, srt_ref = srt_22,
  query_group = "cell_type", ref_group = "Annotation"
)



srt_CS7_epi <- srt_CS7[, srt_CS7$cluster_id == "Epiblast"]
srt_CS7_epi[["Time"]] <- "CS7"

srt_22 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_22_annotation.rds")
srt_22_epi <- srt_22[, srt_22$Annotation == "Epiblast"]
srt_22_epi[["Time"]] <- paste0("imp-", srt_22_epi[["Time", drop = TRUE]])

genes <- Reduce(intersect, list(rownames(srt), rownames(srt_CS7_epi), rownames(srt_22_epi)))
srt_epi_all <- Reduce(merge, list(
  srt[genes, ],
  srt_CS7_epi[genes, ],
  srt_22_epi[genes, ]
))

srt_epi_all$Time <- factor(srt_epi_all$Time, levels = c(
  "D6", "D7", "D8", "D9", "D10", "D12", "D14", "22cellline", "h0cellline", "CS7",
  "imp-d10", "imp-d20", "imp-d30", "imp-d40", "imp-d50", "imp-d70", "imp-d90"
))
ht_result <- ExpHeatmap(srt_epi_all,
  features = degs$gene, group.by = "Time", slot = "data", exp_method = "log1p", nlabel = 0, show_row_names = T,
  n_split = 4, anno_terms = T
)
ht_result$plot


srt_epi_all <- RunKNNPredict(
  srt_query = srt_epi_all, query_group = "Time", features = degs$gene,
  srt_ref = srt_epi_all, ref_group = "Time",
  query_collapsing = T, ref_collapsing = T,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
d <- 1 - srt_epi_all@tools$knnpredict$distance_matrix
d <- as.matrix(d)
ht <- Heatmap(d,
  name = "Cosine similarity", cluster_columns = T, cluster_rows = T,
  col = colorRamp2(seq(min(d), max(d), length.out = 3), c("#27408B", "white", "#EE0000")),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = "white", lwd = 1, fill = "white")
    )
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = fill, lwd = 1, fill = alpha(fill, 0.5))
    )
    if (d[i, j] >= sort(d[i, ], decreasing = T)[3]) {
      # grid.text("*", x, y, gp = gpar(fontsize = 20))
      grid.text(round(d[i, j], 2), x, y, gp = gpar(fontsize = 10))
    }
  },
  border = TRUE,
  width = unit(ncol(d) * 0.8, "cm"),
  height = unit(nrow(d) * 0.8, "cm"),
)
gTree <- grid.grabExpr({
  draw(ht,
    padding = unit(c(3, 1, 1, 3), "cm")
  ) # bottom, left, top and right
})
p <- plot_grid(gTree)

genes <- Reduce(intersect, list(rownames(srt_CS7), rownames(srt_22)))
srt_CS7$Annotation <- paste0("CS7-", srt_CS7$cluster_id)
srt_compare <- Reduce(merge, list(
  srt_CS7[genes, ],
  srt_22[genes, ]
))
srt_compare <- RunKNNPredict(
  srt_query = srt_compare, query_group = "Annotation",
  features = intersect(VariableFeatures(srt_CS7), VariableFeatures(srt_22)),
  srt_ref = srt_compare, ref_group = "Annotation",
  query_collapsing = TRUE, ref_collapsing = TRUE,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
d <- 1 - srt_compare@tools$knnpredict$distance_matrix
d <- as.matrix(d)
d <- d[
  setdiff(as.character(unique(srt_22$Annotation)), c("Radial glial cell", "Ependymal cell", "Neuron", "Retinal progenitor cell", "Retinal pigmented epithelium")),
  as.character(unique(srt_CS7$Annotation))
]
ht <- Heatmap(d,
  name = "Cosine similarity", cluster_columns = TRUE, cluster_rows = TRUE,
  col = colorRamp2(seq(min(d), max(d), length.out = 3), c("#27408B", "white", "#EE0000")),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = "white", lwd = 1, fill = "white")
    )
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = fill, lwd = 1, fill = alpha(fill, 0.5))
    )
    if (d[i, j] >= sort(d[i, ], decreasing = T)[3]) {
      # grid.text("*", x, y, gp = gpar(fontsize = 20))
      grid.text(round(d[i, j], 2), x, y, gp = gpar(fontsize = 10))
    }
  },
  border = TRUE,
  width = unit(ncol(d) * 0.8, "cm"),
  height = unit(nrow(d) * 0.8, "cm"),
)
gTree <- grid.grabExpr({
  draw(ht,
    padding = unit(c(3, 1, 1, 3), "cm")
  ) # bottom, left, top and right
})
p <- plot_grid(gTree)



# 22CellLine-H0CellLine-E14PE-E19Endoderm -------------------------------------------------------------
endo_22 <- srt_22cellline[, srt_22cellline$Annotation == "Endoderm"]
endo_22$Annotation <- "22cellline_endo"
endo_22$Time <- "22cellline_endo"
endo_h0 <- srt_h0cellline[, srt_h0cellline$Annotation == "Endoderm"]
endo_h0$Annotation <- "h0cellline_endo"
endo_h0$Time <- "h0cellline_endo"
endo_3d <- srt_E14[, srt_E14$Group == "PE"]
endo_3d$Annotation <- endo_3d$GroupTime
endo_3d$Time <- endo_3d$Time
endo_E19_sub <- srt_CS7[, srt_CS7$cluster_id == "Endoderm"]
endo_E19_sub$Annotation <- endo_E19_sub$sub_cluster
endo_E19_sub$Time <- "CS7"
endo_E19_sub <- NormalizeData(endo_E19_sub)


srt <- Integration_SCP(srtList = list(endo_22, endo_h0, endo_3d, endo_E19_sub), integration_method = "CSS")
ClassDimPlot(srt, c("Annotation", "Time"))
ExpDimPlot(srt, c("AFP", "TTR"))
ExpHeatmap(srt, group.by = "Time", slot = "data", exp_method = "log2fc", cluster_rows = TRUE)

srt <- readRDS("/ssd/lab/wangjiachen/ygg/reference/Cyno/Three_dataset.rds")
srt[["Annotation"]] <- sapply(strsplit(as.character(srt$Cells), split = "_", fixed = T), function(x) x[[1]])
# srt[["Annotation"]] <- factor(srt[["Annotation",drop=TRUE]],levels = c("EmDisc"))
srt[["Stage"]] <- sapply(strsplit(as.character(srt$Cells), split = "_", fixed = T), function(x) x[[2]])
srt[["Stage"]] <- factor(srt[["Stage", drop = TRUE]], levels = c("E50", "CS3", "CS4", "CS5", "CS6", "CS6/7", "CS7"))
ClassDimPlot(srt, c("Cells", "Stage"))

srt <- RunKNNMap(
  srt_query = srt, srt_ref = srt_22cellline, ref_umap = "StandardUMAP2D",
  ref_group = "Annotation"
)
ProjectionPlot(
  srt_query = srt, srt_ref = srt_22cellline,
  query_group = "Stage", ref_group = "Annotation"
)
srt_22cellline <- RunKNNMap(
  srt_query = srt_22cellline, srt_ref = srt, ref_umap = "umap",
  ref_group = "Annotation"
)
ProjectionPlot(
  srt_query = srt_22cellline, srt_ref = srt,
  query_group = "Annotation", ref_group = "Annotation"
)

srt_22cellline <- RunKNNMap(
  srt_query = srt_22cellline, srt_ref = srt_CS7, ref_umap = "umap",
  ref_group = "Annotation"
)
ProjectionPlot(
  srt_query = srt_22cellline, srt_ref = srt_CS7,
  query_group = "Annotation", ref_group = "Annotation"
)


srt_marmoset <- srt[, srt$DataOrder == "1) Marmoset"]
srt_marmoset <- srt_marmoset[, srt_marmoset$Stage != "CS3"]
srt_marmoset <- Integration_SCP(list(srt_marmoset, srt_CS7), integration_method = "CSS")
ClassDimPlot(srt_marmoset, c("Stage", "Annotation"))

srt_marmoset$StageAnnotation <- paste0(srt_marmoset$Stage, "-", srt_marmoset$Annotation)
srt_marmoset <- RunKNNPredict(
  srt_query = srt_marmoset, query_group = "StageAnnotation",
  srt_ref = srt_22cellline, ref_group = "Annotation",
  query_collapsing = TRUE, ref_collapsing = TRUE, features_type = "DE",
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
d <- 1 - srt_marmoset@tools$knnpredict$distance_matrix
d <- as.matrix(d)
# d <- d[
#   setdiff(as.character(unique(srt_22$Annotation)), c("Radial glial cell", "Ependymal cell", "Neuron", "Retinal progenitor cell", "Retinal pigmented epithelium")),
#   as.character(unique(srt_CS7$Annotation))
# ]
ht <- Heatmap(d,
  name = "Cosine similarity", cluster_columns = TRUE, cluster_rows = TRUE,
  col = colorRamp2(seq(min(d), max(d), length.out = 3), c("#27408B", "white", "#EE0000")),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = "white", lwd = 1, fill = "white")
    )
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = fill, lwd = 1, fill = alpha(fill, 0.5))
    )
    if (d[i, j] >= sort(d[i, ], decreasing = T)[3]) {
      # grid.text("*", x, y, gp = gpar(fontsize = 20))
      grid.text(round(d[i, j], 2), x, y, gp = gpar(fontsize = 10))
    }
  },
  border = TRUE,
  width = unit(ncol(d) * 0.8, "cm"),
  height = unit(nrow(d) * 0.8, "cm"),
)
gTree <- grid.grabExpr({
  draw(ht,
    padding = unit(c(3, 1, 1, 3), "cm")
  ) # bottom, left, top and right
})
p <- plot_grid(gTree)

srt_CS7$Stage <- "Human_CS7"
srt_CS7$Batch <- "Human_CS7"
srt_E14$Batch <- "Human_3Dculture"
srt_merge <- Integration_SCP(srtList = list(srt_CS7, srt_E14), integration_method = "CSS", batch = "Batch")
ClassDimPlot(srt_merge, c("Annotation", "Stage"))


srt_merge <- merge(srt_marmoset, srt_22cellline)
GroupHeatmap(srt_marmoset,
  slot = "data",
  features = srt_marmoset@tools$knnpredict_classification$features,
  group.by = "Annotation", nlabel = 0, cluster_rows = T, show_row_names = TRUE, show_column_names = T, add_dot = F
)


srt_sub <- srt[, colnames(srt)[srt$Annotation %in% c("Epi", "ExMes", "VE", "EmDisc", "SYS", "Am", "PGC", "Hyp", "EmDiscPS")]]
ClassDimPlot(srt_sub, c("Annotation", "Stage", "DataOrder"), pt.size = 1, nrow = 1)
srt_sub$Batch <- srt_sub$DataOrder
srt_22cellline$Batch <- "22cellline"
srt_h0cellline$Batch <- "h0cellline"
srt_sub <- FindVariableFeatures(srt_sub)
srt_22cellline <- FindVariableFeatures(srt_22cellline)
srt_h0cellline <- FindVariableFeatures(srt_h0cellline)
HVF <- Reduce(intersect, list(VariableFeatures(srt_22cellline), VariableFeatures(srt_h0cellline), VariableFeatures(srt_sub)))
srt_compare <- Integration_SCP(
  srtList = list(srt_22cellline, srt_h0cellline, srt_sub),
  batch = "Batch", HVF = HVF,
  integration_method = "Harmony"
)
ClassDimPlot(srt_compare, "Annotation", split.by = "Batch")


srt <- RunKNNPredict(
  srt_query = srt, query_group = "Annotation",
  srt_ref = srt, ref_group = "Annotation",
  query_collapsing = TRUE, ref_collapsing = TRUE,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
d <- 1 - srt@tools$knnpredict$distance_matrix
d <- as.matrix(d)
# d <- d[
#   setdiff(as.character(unique(srt_22$Annotation)), c("Radial glial cell", "Ependymal cell", "Neuron", "Retinal progenitor cell", "Retinal pigmented epithelium")),
#   as.character(unique(srt_CS7$Annotation))
# ]
ht <- Heatmap(d,
  name = "Cosine similarity", cluster_columns = TRUE, cluster_rows = TRUE,
  col = colorRamp2(seq(min(d), max(d), length.out = 3), c("#27408B", "white", "#EE0000")),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = "white", lwd = 1, fill = "white")
    )
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = fill, lwd = 1, fill = alpha(fill, 0.5))
    )
    if (d[i, j] >= sort(d[i, ], decreasing = T)[3]) {
      # grid.text("*", x, y, gp = gpar(fontsize = 20))
      grid.text(round(d[i, j], 2), x, y, gp = gpar(fontsize = 10))
    }
  },
  border = TRUE,
  width = unit(ncol(d) * 0.8, "cm"),
  height = unit(nrow(d) * 0.8, "cm"),
)
gTree <- grid.grabExpr({
  draw(ht,
    padding = unit(c(3, 1, 1, 3), "cm")
  ) # bottom, left, top and right
})
p <- plot_grid(gTree)


#############################################################
srt_d0 <- RunPCAMap(
  srt_query = srt_d0, srt_ref = srt_22,
  ref_umap = "UMAP"
)
srt_d0 <- RunKNNMap(
  srt_query = srt_d0, srt_ref = srt_22,
  ref_umap = "UMAP"
)
ProjectionPlot(
  srt_query = srt_d0, srt_ref = srt_22,
  query_group = "Time", ref_group = "Annotation"
)

srt_22_d0 <- merge(srt_d0, srt_22)
srt_22_d0 <- Integration_SCP(srt_22_d0,
  integration_method = "Uncorrected", nHVF = 3000,
  linear_reduction_dims_use = 1:50, cluster_resolution = 5
)

plist <- list()
for (nHVF in c(2000, 3000)) {
  for (dim in seq(25, 80, 5)) {
    srt_22_d0 <- Integration_SCP(srt_22_d0,
      integration_method = "Uncorrected", nHVF = nHVF,
      linear_reduction_dims_use = 1:dim, cluster_resolution = 5
    )
    srt_22_d0[[paste0("UMAP2D", dim)]] <- srt_22_d0@reductions$UncorrectedUMAP2D
    srt_22_d0[[paste0("UMAP3D", dim)]] <- srt_22_d0@reductions$UncorrectedUMAP3D
    plist[[paste0("dim", dim)]] <- ClassDimPlot(srt_22_d0, "Annotation")
  }
}
p <- plot_grid(plotlist = plist)
p <- panel_fix(p, save = "tmp.png", width = 3)

ClassDimPlot(srt_22_d0, "Annotation")
ClassDimPlot3D(srt_22_d0, "Annotation")
saveRDS(srt_22_d0, "srt_22_d0_annotation.rds")


## H0 data --------------------------------------------------------------------
### Load the data -----------------------------------------------------
srt_h0 <- LoadH5Seurat("/ssd/lab/HuangMingQian/scRNAseq/H0CellLine-project/NGSmodule_SCP_analysis/CellQC/Merge.filtered.h5Seurat")
srt_h0[["percent.mito"]] <- PercentageFeatureSet(object = srt_h0, pattern = paste0(c("^MT-", "^Mt-", "^mt-"), collapse = "|"))
srt_h0[["percent.ribo"]] <- PercentageFeatureSet(object = srt_h0, pattern = paste0(c("^RP[SL]\\d+\\w{0,1}\\d*$", "^Rp[sl]\\d+\\w{0,1}\\d*$", "^rp[sl]\\d+\\w{0,1}\\d*$"), collapse = "|"))
srt_h0$Sample <- as.character(srt_h0$orig.ident)
srt_h0$Sample <- sapply(srt_h0$Sample, function(x) {
  switch(x,
    "H0_CellLine_injection_10d-1" = "d10",
    "H0_CellLine_injection_20d-1" = "d20",
    "H0_CellLine_injection_30d-2" = "d30",
    "H0_CellLine_injection_50d-2" = "d50",
    "H0_CellLine_injection_70d-1" = "d70"
  )
})
srt_h0$Sample <- factor(srt_h0$Sample, levels = sort(unique(srt_h0$Sample)))
srt_h0$Time <- srt_h0$Sample

###  EDA ---------------------------------------------------------------------
srt_h0 <- Integration_SCP(srt_h0,
  integration_method = "Uncorrected", cluster_resolution = 5
  # vars_to_regress = "nCount_RNA", regression_model = "negbinom"
)
for (i in unique(c(c(seq(20, 80, 5))))) {
  srt_h0 <- CSS_integrate(srt_h0,
    linear_reduction = "Uncorrectedpca", HVF = VariableFeatures(srt_h0),
    linear_reduction_dims_use = 1:i, cluster_resolution = 3
  )
  srt_h0[[paste0("CSSUMAP2D", i)]] <- srt_h0@reductions$CSSUMAP2D
  srt_h0[[paste0("CSSUMAP3D", i)]] <- srt_h0@reductions$CSSUMAP3D
}
saveRDS(srt_h0, "srt_h0.rds")

srt_h0_tmp <- RunKNNPredict(srt_h0,
  bulk_ref = SCP::ref_scHCL,
  features = intersect(VariableFeatures(srt_h0), rownames(SCP::ref_scHCL))
)
drop <- table(srt_h0_tmp$knnpredict) %>%
  .[. < 30] %>%
  names()
srt_h0_tmp$knnpredict[srt_h0_tmp$knnpredict %in% drop] <- "unreliable"
ClassDimPlot(srt_h0_tmp, group.by = "knnpredict", label = T)

srt_h0_tmp <- RunKNNPredict(srt_h0, srt_ref = srt_mOrg, ref_group = "Main_cell_type")
ClassDimPlot(srt_h0_tmp, group.by = "knnpredict_Main_cell_type", label = T)

srt_h0_tmp <- RunKNNPredict(srt_h0, srt_ref = srt_hOrg, ref_group = "main_split")
ClassDimPlot(srt_h0_tmp, group.by = "knnpredict_main_split", label = T)

srt_h0_tmp <- RunKNNPredict(srt_h0, srt_ref = srt_mGast, ref_group = "celltype")
ClassDimPlot(srt_h0_tmp, group.by = "knnpredict_celltype", label = T)

### Integration -----------------------------------------------------
srt_h0 <- Integration_SCP(srt_h0, integration_method = "CSS", linear_reduction_dims_use = 1:80, cluster_resolution = 7)
# srt_h0 <- FindClusters(object = srt_h0, resolution = 7, algorithm = 1, graph.name = "CSS_SNN")
# srt_h0 <- SrtReorder(srt_h0, features = srt_h0@misc$CSS_HVF, reorder_by = "seurat_clusters", slot = "data")
# srt_h0[["seurat_clusters"]] <- NULL
# srt_h0[["CSSclusters"]] <- Idents(srt_h0)
ClassDimPlot(srt_h0, group.by = "CSSclusters", label = TRUE)

srt_h0[["UMAP"]] <- srt_h0[["CSSUMAP2D"]]
srt_h0[["UMAP3D"]] <- srt_h0[["CSSUMAP3D"]]
srt_h0@misc$Default_reduction <- "UMAP"

colnames(srt_h0[["UMAP"]]@cell.embeddings) <- c("UMAP_1", "UMAP_2")
colnames(srt_h0[["UMAP3D"]]@cell.embeddings) <- c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3")
Key(srt_h0[["UMAP"]]) <- "UMAP_"
Key(srt_h0[["UMAP3D"]]) <- "UMAP3D_"
srt_h0@misc$Default_reduction <- "UMAP"


srt_h0 <- RunKNNPredict(
  srt_query = srt_h0, srt_ref = srt_22,
  ref_group = "CSSclusters",
  query_collapsing = FALSE, ref_collapsing = FALSE
)
srt_h0$knnpredict_CSSclusters <- factor(srt_h0$knnpredict_CSSclusters, levels = levels(srt_22$CSSclusters))
ClassDimPlot(srt_h0, "knnpredict_CSSclusters", label = T)

srt_h0 <- RunDEtest(srt_h0, group_by = "CSSclusters", test.use = "wilcox", BPPARAM = BiocParallel::MulticoreParam(workers = 30))
srt_h0 <- RunDEtest(srt_h0, group_by = "CSSclusters", test.use = "roc", BPPARAM = BiocParallel::MulticoreParam(workers = 30))

de_filter <- filter(srt_h0@tools$DEtest_CSSclusters$AllMarkers_wilcox, p_val_adj < 0.05 & DE_group_number <= 5)
res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 2)
p <- panel_fix(p, height = 1.2, width = 2, save = "tmp_h0.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 5)
p <- panel_fix(p, height = 3, width = 3, save = "tmp_h0.word.pdf")

saveRDS(srt_h0, "srt_h0.rds")

### Annotation --------------------------------------------------------------
ClassDimPlot(srt_h0, "CSSclusters", reduction = "CSSUMAP2D", label = T)

srt_h0 <- RunKNNPredict(
  srt_query = srt_h0, srt_ref = srt_22,
  ref_group = "Annotation",
  query_collapsing = FALSE, ref_collapsing = FALSE
)
srt_h0$knnpredict_Annotation <- factor(srt_h0$knnpredict_Annotation, levels = levels(srt_22$Annotation))
ClassDimPlot(srt_h0, "knnpredict_Annotation", label = T)

class <- "8,9,6,7,12,13,14,15,16: Epiblast
10,11,18: Primitive streak
25,26,53: Mesoderm
54,55: ExE mesoderm
56,57: MSC/Fib
58,59,60,61,62,63,64: Limb bud mesenchyme cell
1,2,34: Endothelial & erythroid cell
19,21,22,23,24,4: Endoderm
3,5: ExE endoderm
30,29,28,27,50: Gut
20: Amnion & PGC
31,32,33: Epithelium
17,38: Neural ectoderm
39,42,43: Radial glial cell
44,45,46,47,48,49: Ependymal cell
40,35,36,37: Schwann cell
52: Neuron
41: Retinal progenitor cell
51: Sensory neuron"

class <- readLines(textConnection(class))
srt_h0[["Annotation"]] <- srt_h0[["CSSclusters"]]
srt_h0[["Annotation"]] <- as.character(srt_h0[["Annotation", drop = TRUE]])
for (i in 1:length(class)) {
  m <- strsplit(class[i], ": ")[[1]]
  srt_h0$Annotation[srt_h0$Annotation %in% strsplit(m[1], ",")[[1]]] <- m[[2]]
}
levels <- c(
  "Epiblast", "Primitive streak",
  "Mesoderm", "MSC/Fib", "ExE mesoderm", "Limb bud mesenchyme cell", "Endothelial & erythroid cell",
  "Endoderm", "ExE endoderm", "Gut",
  "Amnion & PGC", "Epithelium",
  "Neural ectoderm", "Radial glial cell", "Ependymal cell", "Neuron", "Schwann cell", "Sensory neuron", "Retinal progenitor cell"
)
all(levels %in% srt_h0$Annotation)
srt_h0$Annotation <- factor(srt_h0$Annotation, levels = levels)

srt_h0$GermLayer <- sapply(srt_h0$Annotation, function(x) {
  switch(as.character(x),
    "Epiblast" = "Epiblast",
    "Primitive streak" = "Primitive streak",
    "Mesoderm" = "Mesoderm",
    "Mixed mesoderm" = "Mesoderm",
    "ExE mesoderm" = "Mesoderm",
    "MSC/Fib" = "Mesoderm",
    "Limb bud mesenchyme cell" = "Mesoderm",
    "Endothelial & erythroid cell" = "Mesoderm",
    "Endoderm" = "Endoderm",
    "ExE endoderm" = "Endoderm",
    "Gut" = "Endoderm",
    "Amnion" = "Non-neural ectoderm",
    "PGC" = "Non-neural ectoderm",
    "Amnion & PGC" = "Non-neural ectoderm",
    "Epithelium" = "Non-neural ectoderm",
    "Neural ectoderm" = "Neural ectoderm",
    "Radial glial cell" = "Neural ectoderm",
    "Ependymal cell" = "Neural ectoderm",
    "Neuron" = "Neural ectoderm",
    "Schwann cell" = "Neural ectoderm",
    "Sensory neuron" = "Neural ectoderm",
    "Retinal progenitor cell" = "Neural ectoderm",
    "Retinal pigmented epithelium" = "Neural ectoderm"
  )
})
srt_h0$GermLayer <- factor(srt_h0$GermLayer, levels = c("Epiblast", "Primitive streak", "Mesoderm", "Endoderm", "Non-neural ectoderm", "Neural ectoderm"))

ClassDimPlot(srt_h0, group.by = "Annotation", reduction = "CSSUMAP2D", label = TRUE, label_insitu = FALSE)
ClassDimPlot(srt_h0, c("Annotation", "GermLayer"))

srt_h0$Time <- factor(paste0("d", gsub("d", "", srt_h0$Time)), levels = c("d10", "d20", "d30", "d50", "d70"))
# saveRDS(srt_h0, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_h0_annotation.rds")

### Baisc plot --------------------------------------------------------------------
p <- ClassDimPlot(srt_h0, "Annotation", label = TRUE, theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 4, save = "figures/fig2/h0_annotation.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_h0, "CSSclusters", label = TRUE, label_insitu = TRUE, theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 4, save = "figures/fig2/h0_CSSclusters.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_h0, "Time", theme_use = "theme_blank", show_stat = TRUE, force = TRUE)
p <- panel_fix(p, height = 4, save = "figures/fig2/h0_sample.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_h0, "GermLayer", theme_use = "theme_blank", label = T, show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 4, save = "figures/fig2/h0_germlayer.umap.pdf", raster = TRUE)

p <- ClassDimPlot(srt_h0, "Annotation",
  split.by = "Time", legend.position = "none", bg_color = "grey90",
  theme_use = "theme_blank", ncol = 4, force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig2/h0_sample_split.umap.pdf", raster = TRUE)

p <- ClassDimPlot3D(srt_h0, "Annotation",
  reduction = "CSSUMAP3D",
  width = 1000, height = 700, axis_labs = paste0("UMAP_", 1:3),
  save = "figures/fig2/h0_annotation.umap3D.html"
)

### Annotation stat -------------------------------------------------------------------------
p <- ClassStatPlot(srt_h0, stat.by = "Annotation", group.by = "Time", plot_type = "area", aspect.ratio = 0.5)
p <- panel_fix(p, width = 5, save = "figures/fig2/h0_annotation.stat.pdf", raster = TRUE)


### Markers -----------------------------------------------------------------
n <- 2
markers_use <- markers[levels(srt_h0$Annotation)] %>% lapply(function(x) head(x, 2))
all(srt_h0$Annotation %in% names(markers_use))

p <- ExpDimPlot(srt_h0,
  features = unlist(markers_use),
  theme_use = "theme_blank", ncol = 6
)
p <- panel_fix(p, height = 2, save = "figures/fig2/h0_markers.umap.pdf", raster = TRUE)

p <- ExpDotPlot(srt_h0,
  genes = unlist(markers_use),
  feature_split = unlist(lapply(names(markers_use), function(x) rep(x, length(markers_use[[x]])))),
  cell_split_by = "Annotation"
)
ggsave(plot = p, filename = "figures/fig2/h0_markers.dotplot.pdf", width = 10, height = 15)


ht <- GroupHeatmap(srt_h0,
  group.by = "Annotation",
  features = unlist(markers_use),
  feature_split = rep(names(markers_use), each = 2),
  show_row_names = TRUE, add_dot = TRUE, add_bg = TRUE,
  feature_split_palette = "Paired", heatmap_palette = "YlOrRd",
  height = 10, width = 8, dot_size = unit(7, "mm")
)
panel_fix(ht$plot, save = "figures_add/h0testis.markers.pdf")


### Cell cycle --------------------------------------------------------------
srt_h0 <- CellCycleScoring(srt_h0,
  s.features = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes
)
srt_h0$Phase <- factor(srt_h0$Phase, levels = c("G1", "S", "G2M"))
p <- ClassDimPlot(srt_h0, "Phase", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig2/h0_cellcycle_phase.umap.pdf", raster = TRUE)

p <- ClassStatPlot(srt_h0, stat.by = "Phase", group.by = "Annotation", aspect.ratio = 0.5)
p <- panel_fix(p, width = 5, save = "figures/fig2/h0_cellcycle.stat.pdf", raster = TRUE)


### PAGA ---------------------------------------------------------------------
srt_h0_paga <- RunPAGA(
  srt = srt_h0, group_by = "Annotation", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
  n_pcs = ncol(srt_h0@reductions$CSS@cell.embeddings), n_neighbors = 100, embedded_with_PAGA = FALSE, return_seurat = TRUE
)
srt_h0@misc$paga <- srt_h0_paga@misc$paga
p <- ClassDimPlot(srt_h0,
  group.by = "Annotation", pt.size = 5, pt.alpha = 0.01,
  label = TRUE, label.size = 3, label_repel = TRUE, label_insitu = TRUE, label_segment_color = "transparent",
  paga = srt_h0@misc$paga, paga_edge_threshold = 0.1, paga_edge_color = "black", paga_edge_alpha = 1,
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
)
p <- panel_fix(p, width = 4, save = "figures/fig2/h0_paga.pdf", raster = TRUE)

### SCVELO ------------------------------------------------------------------
srt_h0_scv <- SCP:::RunSCVELO(
  srt = srt_h0, group_by = "Annotation", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
  n_pcs = ncol(srt_h0@reductions$CSS@cell.embeddings), n_neighbors = 100, return_seurat = TRUE
)
srt_h0[["stochastic_UMAP"]] <- srt_h0_scv[["stochastic_CSSUMAP2D"]]
srt_h0[["Ms"]] <- srt_h0_scv[["Ms"]]
srt_h0[["Mu"]] <- srt_h0_scv[["Mu"]]
srt_h0[["stochastic"]] <- srt_h0_scv[["stochastic"]]
srt_h0[["variance_stochastic"]] <- srt_h0_scv[["variance_stochastic"]]
p <- ClassDimPlot(srt_h0,
  group.by = "Annotation", pt.size = 5, pt.alpha = 0.01,
  velocity = "stochastic", velocity_plot_type = "stream",
  velocity_density = 2, velocity_smooth = 1,
  streamline_n = 20, streamline_size = 0.5, streamline_color = "black",
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
)
p <- panel_fix(p, width = 4, save = "figures/fig2/h0_scvelo.pdf", raster = TRUE)

### DE analysis ------------------------------------------------------------------
srt_h0 <- RunDEtest(srt_h0,
  group_by = "Annotation", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_h0 <- RunDEtest(srt_h0,
  group_by = "Annotation", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)

de_filter <- filter(srt_h0@tools$DEtest_Annotation$AllMarkers_wilcox, p_val_adj < 0.05)
ribogene <- rownames(srt_h0)[grep("^RP[SL]\\d+\\w{0,1}\\d*$", rownames(srt_h0))]
de_filter <- de_filter[!de_filter$gene %in% ribogene, ]
write.xlsx(de_filter, file = "figures/fig2/h0_allDE.xlsx")

de_top <- de_filter %>%
  group_by(gene) %>%
  top_n(1, avg_log2FC) %>%
  group_by(group1) %>%
  top_n(3, avg_log2FC)
p <- ExpDotPlot(srt_h0, features = de_top$gene, feature_split = de_top$group1, cell_split_by = c("Annotation", "Time"))
ggsave(plot = p, filename = "figures/fig2/h0_topDE.heatmap.pdf", width = 10, height = 15)

ht <- ExpHeatmap(srt_h0,
  features = de_filter$gene, feature_split = de_filter$group1, cell_split_by = "Annotation", use_raster = TRUE,
  anno_terms = TRUE, anno_features = TRUE, topTerm = 4, topWord = 15, db_version = "3.13",
  nlabel = 0, feature_split_palette = "Paired", width = 7, height = 12
)
ggsave(plot = ht$plot, filename = "figures/fig2/h0_allDE.heatmap.pdf", width = 20, height = 20, limitsize = FALSE)
saveRDS(ht, "figures/fig2/h0_allDE.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig2/h0_allDE.enrichment.xlsx")

srt_h0 <- RunEnrichment(srt_h0, group_by = "Annotation", geneID_exclude = ribogene, db_version = "3.13")
p <- EnrichmentPlot(srt_h0, group_by = "Annotation", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "figures/fig2/h0_bp_allDE.bar.pdf", raster = TRUE)
p <- EnrichmentPlot(srt_h0, group_by = "Annotation", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "figures/fig2/h0_bp_allDE.word-term.pdf", raster = TRUE)
p <- EnrichmentPlot(srt_h0,
  group_by = "Annotation", plot_type = "wordcloud", word_type = "feature", combine = TRUE, ncol = 6,
  topWord = 30, word_size = c(4, 8), aspect.ratio = 0.3, legend.position = "bottom"
)
p <- panel_fix(p, height = 1.5, save = "figures/fig2/h0_bp_allDE.word-gene2.pdf", raster = TRUE)

saveRDS(srt_h0, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_h0_annotation.rds")

### Compare with human CS7 --------------------------------------------------------
srt_CS7 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_E19/human_gastrulation_E19.seurat.rds")
srt_CS7[["percent.mito"]] <- PercentageFeatureSet(object = srt_CS7, pattern = paste0(c("^MT-", "^Mt-", "^mt-"), collapse = "|"))
srt_CS7[["percent.ribo"]] <- PercentageFeatureSet(object = srt_CS7, pattern = paste0(c("^RP[SL]\\d+\\w{0,1}\\d*$", "^Rp[sl]\\d+\\w{0,1}\\d*$", "^rp[sl]\\d+\\w{0,1}\\d*$"), collapse = "|"))
srt_CS7 <- RunKNNMap(
  srt_query = srt_CS7, srt_ref = srt_h0,
  ref_group = "CSSclusters", ref_umap = "CSSUMAP2D"
)
levels <- c(
  "Epiblast", "Primitive Streak", "Axial Mesoderm", "Nascent Mesoderm", "Emergent Mesoderm", "Advanced Mesoderm", "ExE Mesoderm", "Endoderm",
  "Non-Neural Ectoderm", "Hemogenic Endothelial Progenitors", "Erythroblasts"
)
srt_CS7$cluster_id <- factor(srt_CS7$cluster_id, levels = levels)
srt_CS7$cell_type <- srt_CS7$cluster_id

srt_h0$bg <- NA
p <- ProjectionPlot(
  srt_query = srt_CS7, srt_ref = srt_h0,
  query_group = "cell_type", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "CS7 human gastrula"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig2/h0_CS7_projection_celltype.umap.pdf", raster = TRUE)

p <- ProjectionPlot(
  srt_query = srt_CS7, srt_ref = srt_h0,
  query_group = "spatial", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "CS7 human gastrula"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig2/h0_CS7_projection_spatial.umap.pdf", raster = TRUE)

### Compare with human E14 --------------------------------------------------
srt_E14 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_3dculture/human_gastrulation_3dculture.seurat.rds")
srt_E14_sub <- srt_E14[, srt_E14$Group %in% c("EPI", "PSA-EPI", "PE")]
srt_E14_sub <- FindVariableFeatures(srt_E14_sub)
srt_E14_sub <- RunKNNMap(
  srt_query = srt_E14_sub, srt_ref = srt_h0,
  ref_group = "CSSclusters", ref_umap = "CSSUMAP2D"
)
srt_E14_sub$Group <- factor(srt_E14_sub$Group, levels = c("EPI", "PSA-EPI", "PE"))
srt_E14_sub$cell_type <- srt_E14_sub$Group
srt_h0$bg <- NA
p <- ProjectionPlot(
  srt_query = srt_E14_sub, srt_ref = srt_h0,
  query_group = "cell_type", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "3D-cultured human pre-gastrulation embryos"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig2/h0_E14_projection_celltype.umap.pdf", raster = TRUE)

p <- ProjectionPlot(
  srt_query = srt_E14_sub, srt_ref = srt_h0,
  query_group = "Time", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "3D-cultured human pre-gastrulation embryos"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig2/h0_E14_projection_time.umap.pdf", raster = TRUE)


## 22ESC data --------------------------------------------------------------------
### Load the data -----------------------------------------------------
srt_esc <- LoadH5Seurat("/ssd/lab/HuangMingQian/scRNAseq/22ESC-project/NGSmodule_SCP_analysis/CellQC/Merge.filtered.h5Seurat")
srt_esc[["percent.mito"]] <- PercentageFeatureSet(object = srt_esc, pattern = paste0(c("^MT-", "^Mt-", "^mt-"), collapse = "|"))
srt_esc[["percent.ribo"]] <- PercentageFeatureSet(object = srt_esc, pattern = paste0(c("^RP[SL]\\d+\\w{0,1}\\d*$", "^Rp[sl]\\d+\\w{0,1}\\d*$", "^rp[sl]\\d+\\w{0,1}\\d*$"), collapse = "|"))
srt_esc$Sample <- as.character(srt_esc$orig.ident)
srt_esc$Sample <- sapply(srt_esc$Sample, function(x) {
  switch(x,
    "22_ESC_injection_10d-1" = "d10-1",
    "22_ESC_injection_10d-2" = "d10-2",
    "22_ESC_injection_20d-1" = "d20-1",
    "22_ESC_injection_20d-2" = "d20-2",
    "22_ESC_injection_30d-1" = "d30-1",
    "22_ESC_injection_30d-2" = "d30-2",
    "22_ESC_injection_30d-3" = "d30-3",
    "22_ESC_injection_40d-1" = "d40-1",
    "22_ESC_injection_50d-1" = "d50-1",
    "22_ESC_injection_50d-2" = "d50-2",
    "22_ESC_injection_70d-1" = "d70-1",
    "22_ESC_injection_70d-2" = "d70-2",
  )
})
srt_esc$Sample <- factor(srt_esc$Sample, levels = sort(unique(srt_esc$Sample)))
srt_esc$Time <- factor(gsub("-.*", "", srt_esc$Sample))

###  EDA ---------------------------------------------------------------------
srt_esc_tmp <- RunKNNPredict(srt_esc,
  bulk_ref = SCP::ref_scHCL,
  features = intersect(VariableFeatures(srt_h0), rownames(SCP::ref_scHCL))
)
drop <- table(srt_esc_tmp$knnpredict) %>%
  .[. < 50] %>%
  names()
srt_esc_tmp$knnpredict[srt_esc_tmp$knnpredict %in% drop] <- "unreliable"
ClassDimPlot(srt_esc_tmp, group.by = "knnpredict", label = T)

srt_esc_tmp <- RunKNNPredict(srt_esc, srt_ref = srt_mOrg, ref_group = "Main_cell_type")
ClassDimPlot(srt_esc_tmp, group.by = "knnpredict_Main_cell_type", label = T)

srt_esc_tmp <- RunKNNPredict(srt_esc, srt_ref = srt_hOrg, ref_group = "main_split")
ClassDimPlot(srt_esc_tmp, group.by = "knnpredict_main_split", label = T)

srt_esc_tmp <- RunKNNPredict(srt_esc, srt_ref = srt_mGast, ref_group = "celltype")
ClassDimPlot(srt_esc_tmp, group.by = "knnpredict_celltype", label = T)


### Integration -----------------------------------------------------
srt_esc <- Integration_SCP(srt_esc, integration_method = "CSS", linear_reduction_dims_use = 1:45, cluster_resolution = 5)
ClassDimPlot(srt_esc, group.by = "CSSclusters", label = TRUE)
ClassDimPlot3D(srt_esc, group.by = "CSSclusters")

srt_esc[["UMAP"]] <- srt_esc[["CSSUMAP2D"]]
srt_esc[["UMAP3D"]] <- srt_esc[["CSSUMAP3D"]]
srt_esc@misc$Default_reduction <- "UMAP"

colnames(srt_esc[["UMAP"]]@cell.embeddings) <- c("UMAP_1", "UMAP_2")
colnames(srt_esc[["UMAP3D"]]@cell.embeddings) <- c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3")
Key(srt_esc[["UMAP"]]) <- "UMAP_"
Key(srt_esc[["UMAP3D"]]) <- "UMAP3D_"
srt_esc@misc$Default_reduction <- "UMAP"

srt_esc <- RunKNNPredict(
  srt_query = srt_esc, srt_ref = srt_22,
  ref_group = "Annotation",
  query_collapsing = FALSE, ref_collapsing = FALSE
)
srt_esc$knnpredict_Annotation_22 <- factor(srt_esc$knnpredict_Annotation, levels = levels(srt_22$Annotation))
ClassDimPlot(srt_esc, "knnpredict_Annotation_22", label = TRUE)

srt_esc <- RunKNNPredict(
  srt_query = srt_esc, srt_ref = srt_h0,
  ref_group = "Annotation",
  query_collapsing = FALSE, ref_collapsing = FALSE
)
srt_esc$knnpredict_Annotation_h0 <- factor(srt_esc$knnpredict_Annotation, levels = levels(srt_h0$Annotation))
ClassDimPlot(srt_esc, "knnpredict_Annotation_h0", label = TRUE)

srt_esc <- RunDEtest(srt_esc, group_by = "CSSclusters", test.use = "wilcox", BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE)
srt_esc <- RunDEtest(srt_esc, group_by = "CSSclusters", test.use = "roc", BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE)

de_filter <- filter(srt_esc@tools$DEtest_CSSclusters$AllMarkers_wilcox, p_val_adj < 0.05 & DE_group_number <= 5)
res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, GO_simplify = FALSE, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 2)
p <- panel_fix(p, height = 1.2, width = 2, save = "tmp_esc.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 5)
p <- panel_fix(p, height = 3, width = 3, save = "tmp_esc.word.pdf")

saveRDS(srt_esc, "srt_esc.rds")

### Annotation --------------------------------------------------------------
ClassDimPlot(srt_esc, "CSSclusters", reduction = "CSSUMAP2D", label = T)

class <- "27,28,29,30,31,32,33: Epiblast
34,35: Primitive streak
41,64: Mesoderm
63,66,67,68: ExE mesoderm
65,69,70: MSC/Fib
58,59,60,61,62,6: Limb bud mesenchyme cell
57: Endothelial & erythroid cell
37,38,39,40,42,43,44,45: Endoderm
46: ExE endoderm
47,48,49,50,51,52: Gut
36: Amnion & PGC
53,54,55,56: Epithelium
1,2,3,4: Neural ectoderm
7,8,9,10,11,12,13,14,15,16,17,18: Radial glial cell
21: Ependymal cell
5: Schwann cell
23,24,25,26: Neuron
22: Sensory neuron
20,19: Retinal progenitor cell & RPE"

class <- readLines(textConnection(class))
srt_esc[["Annotation"]] <- srt_esc[["CSSclusters"]]
srt_esc[["Annotation"]] <- as.character(srt_esc[["Annotation", drop = TRUE]])
for (i in 1:length(class)) {
  m <- strsplit(class[i], ": ")[[1]]
  srt_esc$Annotation[srt_esc$Annotation %in% strsplit(m[1], ",")[[1]]] <- m[[2]]
}
levels <- c(
  "Epiblast", "Primitive streak",
  "Mesoderm", "MSC/Fib", "ExE mesoderm", "Limb bud mesenchyme cell", "Endothelial & erythroid cell",
  "Endoderm", "ExE endoderm", "Gut",
  "Amnion & PGC", "Epithelium",
  "Neural ectoderm", "Radial glial cell", "Ependymal cell", "Neuron", "Schwann cell", "Sensory neuron", "Retinal progenitor cell & RPE"
)
all(levels %in% srt_esc$Annotation)
all(srt_esc$Annotation %in% levels)
srt_esc$Annotation <- factor(srt_esc$Annotation, levels = levels)

srt_esc$GermLayer <- sapply(srt_esc$Annotation, function(x) {
  switch(as.character(x),
    "Epiblast" = "Epiblast",
    "Primitive streak" = "Primitive streak",
    "Mesoderm" = "Mesoderm",
    "Mixed mesoderm" = "Mesoderm",
    "ExE mesoderm" = "Mesoderm",
    "MSC/Fib" = "Mesoderm",
    "Limb bud mesenchyme cell" = "Mesoderm",
    "Endothelial & erythroid cell" = "Mesoderm",
    "Endoderm" = "Endoderm",
    "ExE endoderm" = "Endoderm",
    "Gut" = "Endoderm",
    "Amnion" = "Non-neural ectoderm",
    "PGC" = "Non-neural ectoderm",
    "Amnion & PGC" = "Non-neural ectoderm",
    "Epithelium" = "Non-neural ectoderm",
    "Neural ectoderm" = "Neural ectoderm",
    "Radial glial cell" = "Neural ectoderm",
    "Ependymal cell" = "Neural ectoderm",
    "Neuron" = "Neural ectoderm",
    "Schwann cell" = "Neural ectoderm",
    "Sensory neuron" = "Neural ectoderm",
    "Retinal progenitor cell" = "Neural ectoderm",
    "Retinal pigmented epithelium" = "Neural ectoderm",
    "Retinal progenitor cell & RPE" = "Neural ectoderm"
  )
})
srt_esc$GermLayer <- factor(srt_esc$GermLayer, levels = c("Epiblast", "Primitive streak", "Mesoderm", "Endoderm", "Non-neural ectoderm", "Neural ectoderm"))

ClassDimPlot(srt_esc, group.by = "Annotation", reduction = "CSSUMAP2D", label = TRUE, label_insitu = FALSE)
ClassDimPlot(srt_esc, c("Annotation", "GermLayer"))

srt_esc$Time <- factor(paste0("d", gsub("d", "", srt_esc$Time)), levels = c("d10", "d20", "d30", "d50", "d70"))
saveRDS(srt_esc, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_esc_annotation.rds")

### Baisc plot --------------------------------------------------------------------
p <- ClassDimPlot(srt_esc, "Annotation", label = TRUE, theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 4, save = "figures/fig2/esc_annotation.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_esc, "CSSclusters", label = TRUE, label_insitu = TRUE, theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 4, save = "figures/fig2/esc_CSSclusters.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_esc, "Sample", theme_use = "theme_blank", show_stat = TRUE, force = TRUE)
p <- panel_fix(p, height = 4, save = "figures/fig2/esc_sample.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_esc, "GermLayer", theme_use = "theme_blank", label = TRUE, show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 4, save = "figures/fig2/esc_germlayer.umap.pdf", raster = TRUE)

p <- ClassDimPlot(srt_esc, "Annotation",
  split.by = "Sample", legend.position = "none", bg_color = "grey90",
  theme_use = "theme_blank", ncol = 4, force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig2/esc_sample_split.umap.pdf", raster = TRUE)

p <- ClassDimPlot3D(srt_esc, "Annotation",
  reduction = "CSSUMAP3D",
  width = 1000, height = 700, axis_labs = paste0("UMAP_", 1:3),
  save = "figures/fig2/esc_annotation.umap3D.html"
)

### Annotation stat -------------------------------------------------------------------------
p <- ClassStatPlot(srt_esc, stat.by = "Annotation", group.by = "Sample", plot_type = "area", aspect.ratio = 0.5)
p <- panel_fix(p, width = 5, save = "figures/fig2/esc_annotation.stat.pdf", raster = TRUE)


### Markers -----------------------------------------------------------------
n <- 2
markers_use <- markers[levels(srt_esc$Annotation)] %>% lapply(function(x) head(x, 2))
all(srt_esc$Annotation %in% names(markers_use))

p <- ExpDimPlot(srt_esc,
  features = unlist(markers_use),
  theme_use = "theme_blank", ncol = 6
)
p <- panel_fix(p, height = 2, save = "figures/fig2/esc_markers.umap.pdf", raster = TRUE)

p <- ExpDotPlot(srt_esc,
  genes = unlist(markers_use),
  feature_split = unlist(lapply(names(markers_use), function(x) rep(x, length(markers_use[[x]])))),
  cell_split_by = "Annotation"
)
ggsave(plot = p, filename = "figures/fig2/esc_markers.dotplot.pdf", width = 10, height = 20)

### Cell cycle --------------------------------------------------------------
srt_esc <- CellCycleScoring(srt_esc,
  s.features = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes
)
srt_esc$Phase <- factor(srt_esc$Phase, levels = c("G1", "S", "G2M"))
p <- ClassDimPlot(srt_esc, "Phase", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig2/esc_cellcycle_phase.umap.pdf", raster = TRUE)

p <- ClassStatPlot(srt_esc, stat.by = "Phase", group.by = "Annotation", aspect.ratio = 0.5)
p <- panel_fix(p, width = 5, save = "figures/fig2/esc_cellcycle.stat.pdf", raster = TRUE)


### PAGA ---------------------------------------------------------------------
srt_esc_paga <- RunPAGA(
  srt = srt_esc, group_by = "Annotation", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
  n_pcs = ncol(srt_esc@reductions$CSS@cell.embeddings), n_neighbors = 100, embedded_with_PAGA = FALSE, return_seurat = TRUE
)
srt_esc@misc$paga <- srt_esc_paga@misc$paga
p <- ClassDimPlot(srt_esc,
  group.by = "Annotation", pt.size = 5, pt.alpha = 0.01,
  label = TRUE, label.size = 3, label_repel = TRUE, label_insitu = TRUE, label_segment_color = "transparent",
  paga = srt_esc@misc$paga, paga_edge_threshold = 0.1, paga_edge_color = "black", paga_edge_alpha = 1,
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
)
p <- panel_fix(p, width = 4, save = "figures/fig2/esc_paga.pdf", raster = TRUE)

### SCVELO ------------------------------------------------------------------
srt_esc_scv <- SCP:::RunSCVELO(
  srt = srt_esc, group_by = "Annotation", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
  n_pcs = ncol(srt_esc@reductions$CSS@cell.embeddings), n_neighbors = 100, return_seurat = TRUE
)
srt_esc[["stochastic_UMAP"]] <- srt_esc_scv[["stochastic_CSSUMAP2D"]]
srt_esc[["Ms"]] <- srt_esc_scv[["Ms"]]
srt_esc[["Mu"]] <- srt_esc_scv[["Mu"]]
srt_esc[["stochastic"]] <- srt_esc_scv[["stochastic"]]
srt_esc[["variance_stochastic"]] <- srt_esc_scv[["variance_stochastic"]]
p <- ClassDimPlot(srt_esc,
  group.by = "Annotation", pt.size = 5, pt.alpha = 0.01,
  velocity = "stochastic", velocity_plot_type = "stream",
  velocity_density = 2, velocity_smooth = 1,
  streamline_n = 20, streamline_size = 0.5, streamline_color = "black",
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
)
p <- panel_fix(p, width = 4, save = "figures/fig2/esc_scvelo.pdf", raster = TRUE)


### DE analysis ------------------------------------------------------------------
srt_esc <- RunDEtest(srt_esc,
  group_by = "Annotation", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_esc <- RunDEtest(srt_esc,
  group_by = "Annotation", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)

de_filter <- filter(srt_esc@tools$DEtest_Annotation$AllMarkers_wilcox, p_val_adj < 0.05)
ribogene <- rownames(srt_esc)[grep("^RP[SL]\\d+\\w{0,1}\\d*$", rownames(srt_esc))]
de_filter <- de_filter[!de_filter$gene %in% ribogene, ]
write.xlsx(de_filter, file = "figures/fig2/esc_allDE.xlsx")

de_top <- de_filter %>%
  group_by(gene) %>%
  top_n(1, avg_log2FC) %>%
  group_by(group1) %>%
  top_n(3, avg_log2FC)
p <- ExpDotPlot(srt_esc, features = de_top$gene, feature_split = de_top$group1, cell_split_by = c("Annotation", "Time"))
ggsave(plot = p, filename = "figures/fig2/esc_topDE.heatmap.pdf", width = 10, height = 15)

ht <- ExpHeatmap(srt_esc,
  features = de_filter$gene, feature_split = de_filter$group1, cell_split_by = "Annotation", use_raster = TRUE,
  anno_terms = TRUE, anno_features = TRUE, topTerm = 4, topWord = 15, db_version = "3.13",
  nlabel = 0, feature_split_palette = "Paired", width = 7, height = 12
)
ggsave(plot = ht$plot, filename = "figures/fig2/esc_allDE.heatmap.pdf", width = 20, height = 20, limitsize = FALSE)
saveRDS(ht, "figures/fig2/esc_allDE.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig2/esc_allDE.enrichment.xlsx")

srt_esc <- RunEnrichment(srt_esc, group_by = "Annotation", geneID_exclude = ribogene, db_version = "3.13")
write.xlsx(srt_esc@tools$Enrichment_Annotation_wilcox$enrichment, file = "figures/fig2/esc_allDE.enrichment.xlsx")
p <- EnrichmentPlot(srt_esc, group_by = "Annotation", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "figures/fig2/esc_bp_allDE.bar.pdf", raster = TRUE)
p <- EnrichmentPlot(srt_esc, group_by = "Annotation", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "figures/fig2/esc_bp_allDE.word-term.pdf", raster = TRUE)
p <- EnrichmentPlot(srt_esc,
  group_by = "Annotation", plot_type = "wordcloud", word_type = "feature", combine = TRUE, ncol = 6,
  topWord = 30, word_size = c(4, 8), aspect.ratio = 0.3, legend.position = "bottom"
)
p <- panel_fix(p, height = 1.5, save = "figures/fig2/esc_bp_allDE.word-gene2.pdf", raster = TRUE)

saveRDS(srt_esc, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_esc_annotation.rds")

### Compare with human CS7 --------------------------------------------------------
srt_CS7 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_E19/human_gastrulation_E19.seurat.rds")
srt_CS7[["percent.mito"]] <- PercentageFeatureSet(object = srt_CS7, pattern = paste0(c("^MT-", "^Mt-", "^mt-"), collapse = "|"))
srt_CS7[["percent.ribo"]] <- PercentageFeatureSet(object = srt_CS7, pattern = paste0(c("^RP[SL]\\d+\\w{0,1}\\d*$", "^Rp[sl]\\d+\\w{0,1}\\d*$", "^rp[sl]\\d+\\w{0,1}\\d*$"), collapse = "|"))
srt_CS7 <- RunKNNMap(
  srt_query = srt_CS7, srt_ref = srt_esc,
  ref_group = "CSSclusters", ref_umap = "CSSUMAP2D"
)
levels <- c(
  "Epiblast", "Primitive Streak", "Axial Mesoderm", "Nascent Mesoderm", "Emergent Mesoderm", "Advanced Mesoderm", "ExE Mesoderm", "Endoderm",
  "Non-Neural Ectoderm", "Hemogenic Endothelial Progenitors", "Erythroblasts"
)
srt_CS7$cluster_id <- factor(srt_CS7$cluster_id, levels = levels)
srt_CS7$cell_type <- srt_CS7$cluster_id

srt_esc$bg <- NA
p <- ProjectionPlot(
  srt_query = srt_CS7, srt_ref = srt_esc,
  query_group = "cell_type", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "CS7 human gastrula"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig2/esc_CS7_projection_celltype.umap.pdf", raster = TRUE)

p <- ProjectionPlot(
  srt_query = srt_CS7, srt_ref = srt_esc,
  query_group = "spatial", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "CS7 human gastrula"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig2/esc_CS7_projection_spatial.umap.pdf", raster = TRUE)

### Compare with human E14 --------------------------------------------------
srt_E14 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_3dculture/human_gastrulation_3dculture.seurat.rds")
srt_E14_sub <- srt_E14[, srt_E14$Group %in% c("EPI", "PSA-EPI", "PE")]
srt_E14_sub <- FindVariableFeatures(srt_E14_sub)
srt_E14_sub <- RunKNNMap(
  srt_query = srt_E14_sub, srt_ref = srt_esc,
  ref_group = "CSSclusters", ref_umap = "CSSUMAP2D"
)
srt_E14_sub$Group <- factor(srt_E14_sub$Group, levels = c("EPI", "PSA-EPI", "PE"))
srt_E14_sub$cell_type <- srt_E14_sub$Group
srt_esc$bg <- NA
p <- ProjectionPlot(
  srt_query = srt_E14_sub, srt_ref = srt_esc,
  query_group = "cell_type", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "3D-cultured human pre-gastrulation embryos"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig2/esc_E14_projection_celltype.umap.pdf", raster = TRUE)

p <- ProjectionPlot(
  srt_query = srt_E14_sub, srt_ref = srt_esc,
  query_group = "Time", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "3D-cultured human pre-gastrulation embryos"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig2/esc_E14_projection_time.umap.pdf", raster = TRUE)

## KFGs intersection --------------------------------------------------------
library(ggraph)
library(igraph)
library(ggforce)
library(concaveman)
enrich_22 <- srt_22@tools$Enrichment_Annotation_wilcox$enrichment
enrich_h0 <- srt_h0@tools$Enrichment_Annotation_wilcox$enrichment
enrich_esc <- srt_esc@tools$Enrichment_Annotation_wilcox$enrichment
enrich_list <- list(cellline22 = enrich_22, celllineh0 = enrich_h0)

kfg_list <- list()
for (enrich in names(enrich_list)) {
  df <- enrich_list[[enrich]] %>%
    filter(Enrichment %in% "GO_BP") %>%
    group_by(Enrichment, Groups) %>%
    filter(.data[["p.adjust"]] <= 0.05) %>%
    arrange(desc(.data[["pvalue"]])) %>%
    as.data.frame()
  df_list <- split.data.frame(df, ~ Enrichment + Groups)
  df_list <- df_list[lapply(df_list, nrow) > 0]
  kfg_tmp <- lapply(df_list, function(df) {
    df <- df %>%
      mutate(keyword = strsplit(as.character(.data[["geneID"]]), "/")) %>%
      unnest(cols = "keyword") %>%
      group_by(.data[["keyword"]], Enrichment, Groups) %>%
      summarise(
        keyword = .data[["keyword"]],
        score = sum(-(log10(.data[["p.adjust"]]))),
        count = n(),
        Enrichment = .data[["Enrichment"]],
        Groups = .data[["Groups"]],
        .groups = "keep"
      ) %>%
      distinct() %>%
      mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(60, 40))) %>%
      arrange(desc(.data[["score"]])) %>%
      as.data.frame() %>%
      top_n(100, wt = .data[["score"]])
    return(df)
  })
  kfg <- do.call(rbind.data.frame, kfg_tmp)
  kfg <- split(kfg$keyword, kfg$Groups)
  kfg_list[[enrich]] <- kfg
}
kfg_list <- unlist(kfg_list, recursive = FALSE)

node <- data.frame(
  name = names(kfg_list),
  celltype = gsub(".*\\.", "", names(kfg_list)),
  cellsource = gsub("\\..*", "", names(kfg_list)),
  size = sapply(kfg_list, length)
)

edge <- as.data.frame(t(combn(names(kfg_list), 2)), stringsAsFactors = FALSE)
colnames(edge) <- c("from", "to")
edge[, "size"] <- sapply(seq_len(nrow(edge)), function(x) {
  length(intersect(kfg_list[[edge[x, "from"]]], kfg_list[[edge[x, "to"]]]))
})
edge[, "weight"] <- edge[, "size"] / pmax(node[edge[, "from"], "size"], node[edge[, "to"], "size"])

graph <- graph_from_data_frame(edge[edge$size >= 5, ], directed = FALSE, node)

set.seed(222)
p <- ggraph(graph, layout = "fr", weights = weight) +
  # ggforce::geom_mark_hull(
  #   mapping = aes(x = x, y = y, fill = celltype),
  #   color = "grey",
  #   concavity = 4,
  #   expand = unit(3, "mm"),
  #   alpha = 0.5,
  #   label.fontsize = 12,
  #   label.fill = "grey90",
  #   con.colour = "red",
  #   con.size = 1,
  #   con.cap = 0,
  #   show.legend = FALSE
  # ) +
  geom_edge_fan(aes(edge_width = size), alpha = 1) +
  geom_node_point(aes(size = size, fill = celltype, shape = cellsource), color = "white") +
  scale_fill_manual(name = "Cell type", values = palette_scp(names(markers)), breaks = names(markers), guide = guide_legend(override.aes = list(shape = 21, color = "black", size = 4.5, alpha = 1), order = 1)) +
  scale_shape_manual(name = "Cell source", values = c(21, 23, 24), guide = guide_legend(override.aes = list(color = "black", size = 4, alpha = 1), order = 2)) +
  scale_size_continuous(name = "Number of KFGs", range = c(3, 5), guide = guide_legend(override.aes = list(color = "black", alpha = 1), order = 3)) +
  scale_edge_width(name = "Size of intersection", range = c(0.1, 3), guide = guide_legend(order = 4)) +
  # facet_nodes(~cellsource) +
  labs(caption = "Force-directed graph layout") +
  theme_graph(base_family = "sans", caption_size = 10) +
  theme(aspect.ratio = 1, plot.caption = element_text(hjust = 0))
ggsave(p, filename = "figures/fig2/KFG_intersection.network.pdf", width = 30, height = 15)
ggsave(p, filename = "figures_add/KFG_intersection.network.pdf", width = 10, height = 10)


## KFGs heatmap --------------------------------------------------------
library(ComplexHeatmap)
library(circlize)
kfgs <- unique(unlist(kfg_list))
features <- kfgs
mat_22 <- AverageExpression(srt_22, assays = "RNA", features = features, group.by = "Annotation")[[1]]
mat_h0 <- AverageExpression(srt_h0, assays = "RNA", features = features, group.by = "Annotation")[[1]]
colnames(mat_22) <- paste0("22.", colnames(mat_22))
colnames(mat_h0) <- paste0("h0.", colnames(mat_h0))
mat <- do.call(cbind, list(mat_22, mat_h0))
mat <- t(scale(t(mat)))

celltype <- gsub(".*\\.", "", colnames(mat))
cellsource <- gsub("\\..*", "", colnames(mat))

top_annotation <- HeatmapAnnotation(
  cell_type = celltype,
  cell_source = cellsource,
  col = list(
    cell_type = palette_scp(names(markers)),
    cell_source = palette_scp(c("22", "h0"), palcolor = c("black", "grey50"))
  ),
  annotation_legend_param = list(
    cell_type = list(at = names(markers)),
    cell_source = list(at = c("22", "h0"))
  ),
  gp = gpar(col = "black"),
  simple_anno_size = unit(0.8, "cm"),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  which = "column"
)

b <- ceiling(min(abs(quantile(mat, c(0.01, 0.99), na.rm = TRUE))) * 2) / 2
colors <- colorRamp2(seq(-b, b, length = 100), palette_scp(palette = "RdBu"))
ht <- Heatmap(
  matrix = mat,
  name = "Z-score",
  col = colors,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = top_annotation,
  show_heatmap_legend = TRUE,
  border = TRUE,
  use_raster = FALSE
)
gTree <- grid.grabExpr({
  draw(ht) # bottom, left, top and right
})
p <- plot_grid(gTree)
ggsave(plot = p, filename = "figures/fig2/KFG_expression.ht.pdf", width = 10, height = 5)
ggsave(plot = p, filename = "figures_add/KFG_expression.ht.pdf", width = 10, height = 5)


## DEGs intersection --------------------------------------------------------
de_22 <- srt_22@tools$DEtest_Annotation$AllMarkers_wilcox %>% filter(p_val_adj < 0.05)
de_h0 <- srt_h0@tools$DEtest_Annotation$AllMarkers_wilcox %>% filter(p_val_adj < 0.05)
de_esc <- srt_esc@tools$DEtest_Annotation$AllMarkers_wilcox %>% filter(p_val_adj < 0.05)
de_list <- list(celline22 = de_22, cellineh0 = de_h0)
degs1 <- unique(c(de_22$gene, de_h0$gene))
degs2 <- Reduce(intersect, list(de_22$gene, de_h0$gene))

deg_list <- list()
for (group in names(de_list)) {
  deg <- split(de_list[[group]]$gene, de_list[[group]]$group1)
  deg_list[[group]] <- deg
}
deg_list <- unlist(deg_list, recursive = FALSE)

node <- data.frame(
  name = names(deg_list),
  celltype = gsub(".*\\.", "", names(deg_list)),
  cellsource = gsub("\\..*", "", names(deg_list)),
  size = sapply(deg_list, length)
)

edge <- as.data.frame(t(combn(names(deg_list), 2)), stringsAsFactors = FALSE)
colnames(edge) <- c("from", "to")
edge[, "size"] <- sapply(seq_len(nrow(edge)), function(x) {
  length(intersect(deg_list[[edge[x, "from"]]], deg_list[[edge[x, "to"]]]))
})
edge[, "weight"] <- edge[, "size"] / pmax(node[edge[, "from"], "size"], node[edge[, "to"], "size"])

graph <- graph_from_data_frame(edge[edge$size >= 50, ], directed = FALSE, node)

set.seed(11)
p <- ggraph(graph, layout = "fr", weights = weight) +
  ggforce::geom_mark_hull(
    aes(x, y, group = celltype, fill = celltype, label = celltype),
    color = "grey",
    concavity = 4,
    expand = unit(3, "mm"),
    alpha = 0.5,
    label.fontsize = 12,
    label.fill = "grey90",
    con.colour = "red",
    con.size = 1,
    con.cap = 0,
    show.legend = FALSE
  ) +
  geom_edge_fan(aes(edge_width = size), alpha = 1) +
  geom_node_point(aes(size = size, fill = celltype, shape = cellsource), color = "white") +
  scale_fill_manual(name = "Cell type", values = palette_scp(names(markers)), guide = guide_legend(override.aes = list(shape = 21, color = "black", size = 4.5, alpha = 1), order = 1)) +
  scale_shape_manual(name = "Cell source", values = c(21, 23, 24), guide = guide_legend(override.aes = list(color = "black", size = 4, alpha = 1), order = 2)) +
  scale_size_continuous(name = "Number of KFGs", range = c(3, 5), guide = guide_legend(override.aes = list(color = "black", size = 4, alpha = 1), order = 3)) +
  scale_edge_width(name = "Size of intersection", range = c(0.1, 3), guide = guide_legend(order = 4)) +
  # facet_nodes(~cellsource) +
  labs(caption = "Force-directed graph layout") +
  theme_graph(base_family = "sans", caption_size = 10) +
  theme(aspect.ratio = 1, plot.caption = element_text(hjust = 0))
ggsave(p, filename = "figures/fig2/DEG_intersection.network.pdf", width = 30, height = 15)


## Cell type similarity --------------------------------------------------------------
library(ComplexHeatmap)
features <- Reduce(intersect, list(VariableFeatures(srt_22), VariableFeatures(srt_h0), VariableFeatures(srt_esc)))
# features <- unique(c(VariableFeatures(srt_22), VariableFeatures(srt_h0), VariableFeatures(srt_esc)))
features <- unique(unlist(kfg_list))
mat_22 <- AverageExpression(srt_22, assays = "RNA", features = features, group.by = "Annotation")[[1]]
mat_h0 <- AverageExpression(srt_h0, assays = "RNA", features = features, group.by = "Annotation")[[1]]
colnames(mat_22) <- paste0("22-", colnames(mat_22))
colnames(mat_h0) <- paste0("h0-", colnames(mat_h0))
mat <- do.call(cbind, list(mat_22, mat_h0))
srt <- CreateSeuratObject(counts = mat)
srt$Annotation <- colnames(srt)

srt <- RunKNNPredict(
  srt_query = srt, query_group = "Annotation",
  srt_ref = srt, ref_group = "Annotation", features = features,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
d <- 1 - srt@tools$knnpredict$distance_matrix
d <- as.matrix(d)

markers_use <- markers
markers_use <- markers_use[names(markers_use) %in% c(srt_22$Annotation, srt_h0$Annotation)]

order <- paste0(c("22-", "h0-"), rep(names(markers_use), each = 2))[paste0(c("22-", "h0-"), rep(names(markers_use), each = 2)) %in% colnames(d)]
d <- d[order, order]

bottom_an <- HeatmapAnnotation(
  which = "column",
  cell_type = sapply(strsplit(colnames(d), "-"), function(x) x[[2]]),
  cell_source = sapply(strsplit(colnames(d), "-"), function(x) x[[1]]),
  col = list(
    cell_type = palette_scp(names(markers_use)),
    cell_source = palette_scp(c("22", "h0"), palcolor = c("grey30", "grey90"))
  ),
  annotation_legend_param = list(
    cell_type = list(at = names(markers_use)),
    cell_source = list(at = c("22", "h0"))
  ),
  gp = gpar(col = "black"),
  show_annotation_name = FALSE,
  show_legend = TRUE
)
right_an <- HeatmapAnnotation(
  which = "row",
  cell_type = sapply(strsplit(colnames(d), "-"), function(x) x[[2]]),
  cell_source = sapply(strsplit(colnames(d), "-"), function(x) x[[1]]),
  col = list(
    cell_type = palette_scp(names(markers_use)),
    cell_source = palette_scp(c("22", "h0"), palcolor = c("grey30", "grey90"))
  ),
  gp = gpar(col = "black"),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht <- Heatmap(d,
  name = "Cosine similarity",
  col = colorRamp2(c(min(d), (max(d) + min(d)) / 2, max(d)), c("#27408B", "white", "#EE0000")),
  cluster_columns = F, cluster_rows = F,
  show_column_names = FALSE,
  show_row_names = FALSE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = "white", lwd = 1, fill = "white")
    )
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = fill, lwd = 1, fill = alpha(fill, 0.5))
    )
  },
  border = TRUE,
  bottom_annotation = bottom_an,
  right_annotation = right_an,
  width = unit(ncol(d) * 0.2, "cm"),
  height = unit(nrow(d) * 0.2, "cm"),
)
gTree <- grid.grabExpr({
  draw(ht,
    padding = unit(c(3, 1, 1, 3), "cm")
  ) # bottom, left, top and right
})
p <- plot_grid(gTree)
ggsave(plot = p, filename = "figures/fig2/celltype_similarity.KFGs.pdf", width = 12, height = 12)
ggsave(plot = p, filename = "figures_add/celltype_similarity.KFGs.pdf", width = 12, height = 12)


# Fig3 --------------------------------------------------------------------
dir.create(paste0(work_dir, "/figures/fig3"), recursive = T, showWarnings = FALSE)

## Implantation --------------------------------------------------------------------
srt_d0 <- LoadH5Seurat("../CellCulture-project/NGSmodule_SCP_analysis/CellQC/22_CellLine-3.filtered.h5Seurat")
srt_d0$Time <- "d0"
srt_d0$Annotation <- NA
srt_imp <- merge(srt_d0, srt_22[, srt_22$Time %in% c("d10")])
srt_imp$Time <- factor(srt_imp$Time, levels = c("d0", "d10"))

### Integration  --------------------------------------------------------------------
srt_imp <- CellCycleScoring(srt_imp,
  s.features = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes
)

srt_imp <- Integration_SCP(
  srtMerge = srt_imp, integration_method = "Uncorrected",
  linear_reduction_dims_use = 1:25, cluster_resolution = 2,
  vars_to_regress = c("S.Score", "G2M.Score")
)
srt_imp <- RunDM(srt_imp,
  features = VariableFeatures(srt_imp), slot = "data", dims = NULL,
  k = 15, ndcs = 3, reduction.name = "DM",
)
ClassDimPlot(srt_imp, c("Time", "Annotation", "Uncorrectedclusters", "Phase"), label = T)
ClassDimPlot(srt_imp, c("Annotation", "Uncorrectedclusters", "Phase"), reduction = "DM")
ExpDimPlot(srt_imp, c("POU5F1", "SOX2", "DNMT3B", "TBXT", "EOMES", "MESP1", "SOX17", "ABCG2", "ISL1"))
ExpDimPlot(srt_imp, c("percent.ribo", "percent.mito", "nFeature_RNA", "nCount_RNA"))
ClassDimPlot(srt_imp, "Uncorrectedclusters", label = TRUE)
ClassDimPlot(srt_imp, "Uncorrectedclusters", label = TRUE, reduction = "UncorrectedUMAP3D", dims = c(1, 2))

srt_imp[["UMAP"]] <- srt_imp[["UncorrectedUMAP2D"]]
srt_imp[["UMAP3D"]] <- srt_imp[["UncorrectedUMAP3D"]]
srt_imp@misc$Default_reduction <- "UMAP"
colnames(srt_imp[["UMAP"]]@cell.embeddings) <- c("UMAP_1", "UMAP_2")
colnames(srt_imp[["UMAP3D"]]@cell.embeddings) <- c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3")
Key(srt_imp[["UMAP"]]) <- "UMAP_"
Key(srt_imp[["UMAP3D"]]) <- "UMAP3D_"
srt_imp@misc$Default_reduction <- "UMAP"

srt_imp <- FindClusters(object = srt_imp, resolution = 2, algorithm = 1, graph.name = "Uncorrectedpca_SNN")
srt_imp <- SrtReorder(srt_imp, features = VariableFeatures(srt_imp), reorder_by = "seurat_clusters", slot = "data")
srt_imp[["seurat_clusters"]] <- NULL
srt_imp[["Uncorrectedclusters"]] <- Idents(srt_imp)
ClassDimPlot(srt_imp, c("Uncorrectedclusters"), label = TRUE)

### paga initiated layout ---------------------------------------------------------------
# resolution: 2
plt <- reticulate::import("matplotlib")$pyplot
sc <- reticulate::import("scanpy")
adata <- srt_to_adata(srt_imp)
adata <- RunPAGA(
  adata = adata, group_by = "Uncorrectedclusters",
  linear_reduction = "Uncorrectedpca", nonlinear_reduction = "UncorrectedUMAP2D",
  n_pcs = 25, embedded_with_PAGA = TRUE
)
FA <- adata$obsm[["X_draw_graph_fa"]]
colnames(FA) <- c("FA_1", "FA_2")
rownames(FA) <- adata$obs_names$values
srt_imp[["FA"]] <- CreateDimReducObject(FA, key = "FA_", assay = "RNA")
ClassDimPlot(srt_imp, group.by = "Time", reduction = "FA")

PAGAUMAP2D <- adata$obsm[["PAGAUMAP2D"]]
colnames(PAGAUMAP2D) <- c("PAGAUMAP_1", "PAGAUMAP_2")
rownames(PAGAUMAP2D) <- adata$obs_names$values
srt_imp[["PAGAUMAP2D"]] <- CreateDimReducObject(PAGAUMAP2D, key = "PAGAUMAP_", assay = "RNA")
ClassDimPlot(srt_imp, group.by = "Time", reduction = "PAGAUMAP2D")

srt_imp@misc$Default_reduction <- "FA"

### Monocle2 ----------------------------------------------------------------------------
library(monocle)
expr_matrix <- as(as.matrix(srt_imp@assays$RNA@counts), "sparseMatrix")
p_data <- srt_imp@meta.data
f_data <- data.frame(gene_short_name = row.names(srt_imp), row.names = row.names(srt_imp))
pd <- new("AnnotatedDataFrame", data = p_data)
fd <- new("AnnotatedDataFrame", data = f_data)
cds <- newCellDataSet(expr_matrix,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.1,
  expressionFamily = negbinomial.size()
)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
express_genes <- VariableFeatures(srt_imp)
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)
cds <- reduceDimension(cds,
  max_components = 3,
  method = "DDRTree", norm_method = "log"
)
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = c("Annotation"), cell_size = 0.8)
plot_cell_trajectory(cds, color_by = c("State"), cell_size = 0.8)
cds <- orderCells(cds, root_state = 4)
srt <- as.Seurat(cds)
srt_imp[["DDRTree"]] <- srt[["DDRTree"]]
srt_imp[["Pseudotime"]] <- srt[["Pseudotime"]]
srt_imp[["State"]] <- srt[["State"]]
srt_imp[["State"]] <- sapply(srt_imp[["State", drop = TRUE]], function(x) {
  switch(as.character(x),
    "1" = 6,
    "2" = 5,
    "3" = 4,
    "4" = 1,
    "5" = 2,
    "6" = 3,
    "7" = 2
  )
})
srt_imp[["State"]] <- factor(srt_imp[["State", drop = TRUE]], levels = 1:6)
ClassDimPlot3D(srt_imp, "State", "DDRTree")
ClassDimPlot(srt_imp, "State", "DDRTree")
ExpDimPlot(srt_imp, "Pseudotime", "UncorrectedUMAP2D")
srt_imp@misc$Default_reduction <- "DDRTree"


### FA & DDRTree ------------------------------------------------------------
ClassDimPlot(srt_imp, "State", "DDRTree")
ClassDimPlot(srt_imp, "GermLayer", "DDRTree")

ExpDimPlot(srt_imp, c("POU5F1", "DNMT3B", "SOX2", "ISL1", "TBXT", "NANOS3"), "DDRTree")
srt_imp$SubAnnotation <- paste0("CellLine-", srt_imp$State)
srt_imp$SubAnnotation[!is.na(srt_imp$GermLayer)] <- srt_imp$GermLayer[!is.na(srt_imp$GermLayer)]
srt_imp$SubAnnotation <- factor(srt_imp$SubAnnotation, levels = c(paste0("CellLine-", 1:6), "Epiblast", "PGC", "Primitive streak", "Mesoderm", "Endoderm", "Non-neural ectoderm"))
ClassDimPlot(srt_imp, "SubAnnotation", "DDRTree")
ExpDimPlot(srt_imp, "Pseudotime")


ClassDimPlot(srt_imp, "State", label = TRUE)
srt_imp$TimeState <- paste0(srt_imp$Time, srt_imp$State)
ClassDimPlot(srt_imp, "TimeState", label = T, label_insitu = TRUE)
annotation_list <- list(
  "D0-Epiblast" = c("d04", "d01", "d02", "d03"),
  "D10-Epiblast" = c("d103", "d104", "d101"),
  "D0-ExE mesoderm" = c("d05"),
  "D10-ExE mesoderm" = c("d105"),
  "D0-Primitive endoderm" = c("d06"),
  "D10-Primitive endoderm" = c("d106")
)
annotation <- setNames(rep(names(annotation_list), sapply(annotation_list, length)), nm = unlist(annotation_list))
srt_imp$SubAnnotation <- as.character(srt_imp$SubAnnotation)
srt_imp$SubAnnotation <- annotation[as.character(srt_imp$TimeState)]
srt_imp$SubAnnotation <- factor(srt_imp$SubAnnotation, levels = names(annotation_list))
p <- ClassDimPlot(srt_imp, "SubAnnotation",
  theme_use = "theme_blank", show_stat = FALSE
)
p <- panel_fix(p, height = 4, save = "figures/fig3/imp_annotation.DDRTree.pdf", raster = TRUE)


#### Baisc plot ------------------------------------------------------------
p <- ClassDimPlot(srt_imp, "Time", reduction = "FA", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/imp_sample.fa.pdf", raster = TRUE)
p <- ClassDimPlot(srt_imp, "Time", split.by = "Time", ncol = 1, reduction = "FA", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 2, save = "figures/fig3/imp_sample.fa.pdf", raster = TRUE)
p <- ClassDimPlot(srt_imp, "GermLayer", reduction = "FA", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/imp_germlayer.fa.pdf", raster = TRUE)
p <- ClassDimPlot(srt_imp, "Time", reduction = "DDRTree", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/imp_sample.DDRTree.pdf", raster = TRUE)
p <- ClassDimPlot(srt_imp, "State", reduction = "DDRTree", label = T, label_insitu = T, theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/imp_state.DDRTree.pdf", raster = TRUE)
p <- ClassDimPlot(srt_imp, "SubAnnotation", reduction = "FA", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/imp_annotation.fa.pdf", raster = TRUE)
p <- ClassDimPlot(srt_imp, "SubAnnotation", reduction = "DDRTree", label = TRUE, theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/imp_annotation.DDRTree.pdf", raster = TRUE)
p <- ExpDimPlot(srt_imp, "Pseudotime", palette = "cividis", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/imp_pseudotime.DDRTree.pdf", raster = TRUE)


p <- ExpDimPlot(srt_imp, c("POU5F1", "DNMT3B", "SOX2", "TBXT", "POSTN", "EOMES", "CER1", "SOX17"),
  theme_use = "theme_blank", nrow = 2
)
p <- panel_fix(p, height = 2, save = "figures/fig3/imp_markers.DDRTree.pdf", raster = TRUE)

library(ggridges)
df <- srt_imp[[]]
levels <- df %>%
  group_by(SubAnnotation) %>%
  summarise(m = median(Pseudotime)) %>%
  arrange(desc(m)) %>%
  pull("SubAnnotation") %>%
  as.character()
df$y <- factor(df$SubAnnotation, levels = levels)
p <- ggplot(df, aes(x = Pseudotime, y = y, fill = SubAnnotation)) +
  scale_fill_manual(values = palette_scp(df$SubAnnotation)) +
  geom_density_ridges() +
  labs(y = "") +
  theme_scp(aspect.ratio = 0.5)
p <- panel_fix(p, height = 2, save = "figures/fig3/imp_pseudotime.ridges.pdf", raster = TRUE)

#### Pluripotency markers ----------------------------------------------------
plur_markers <- list(
  core = c("POU5F1", "NANOG", "SOX2"),
  Primed = c("CD24", "SFRP2", "SOX11", "ZIC2", "ZIC3", "FGF2", "SALL2"),
  Naive = c("KLF17", "KLF5", "KLF4", "TFCP2L1", "DNMT3L", "ESRRB", "STAT3", "TBX3")
)
ExpDimPlot(srt_imp, plur_markers$core, nrow = 1)
ExpDimPlot(srt_imp, plur_markers$Primed, nrow = 2)
ExpDimPlot(srt_imp, plur_markers$Naive, nrow = 2)


#### Compare with human E14 --------------------------------------------------
srt_E14 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_3dculture/human_gastrulation_3dculture.seurat.rds")
srt_E14_sub <- srt_E14[, srt_E14$Group %in% c("EPI")]
# srt_E14_sub <- FindVariableFeatures(srt_E14_sub)
srt_imp <- RunKNNPredict(
  srt_query = srt_imp, srt_ref = srt_E14_sub,
  ref_group = "Time", features = degs$gene,
  filter_lowfreq = 30,
  return_full_distance_matrix = TRUE
)
p <- ClassDimPlot(srt_imp, "knnpredict_classification", theme_use = "theme_blank", show_stat = FALSE)
p <- panel_fix(p, height = 3, save = "figures/fig3/imp_knnpredict.DDRTree.pdf", raster = TRUE)

d <- srt_imp@tools$knnpredict_classification$distance_matrix
mat <- as.matrix(1 - d)
Heatmap(mat)

srt_E14_sub <- RunKNNMap(
  srt_query = srt_E14_sub, srt_ref = srt_imp,
  ref_umap = "FA", k = 5,
)
srt_E14_sub$Group <- factor(srt_E14_sub$Groupn, levels = c("EPI", "PSA-EPI", "PE"))
srt_E14_sub$cell_type <- srt_E14_sub$Group
srt_imp$bg <- NA
p <- ProjectionPlot(
  srt_query = srt_E14_sub, srt_ref = srt_imp,
  query_group = "cell_type", ref_group = "bg",
  query_param = list(palette = "Set1"),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "3D-cultured human pre-gastrulation embryos"
  )
)

#### Dynamic analysis --------------------------------------------------------
ExpDimPlot(srt_imp, "Pseudotime")
srt_imp$Lineage_CellLineEpi <- srt_imp@tools$DynamicFeatures_Lineage_CellLineEpi$raw_matrix[, 1]
p <- ExpDimPlot(srt_imp, "Lineage_CellLineEpi", palette = "cividis", theme_use = "theme_blank")
p <- panel_fix(p, height = 2, save = "figures/fig3/imp_pseudotime.Lineage_CellLineEpi.pdf", raster = TRUE)

srt_imp <- RunDynamicFeatures(srt_imp,
  lineages = "Lineage_CellLineEpi",
  n_candidates = 10000
)
saveRDS(srt_imp, "srt_imp_annotation.rds")

colnames(srt_imp@tools$DynamicFeatures_Lineage_CellLineEpi$DynamicFeatures)[5] <- "peaktime"
ht <- DynamicHeatmap(
  srt = srt_imp, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.001, use_fitted = TRUE,
  lineages = c("Lineage_CellLineEpi"),
  cell_annotation = "SubAnnotation",
  n_split = 5, split_method = "kmeans-peaktime",
  use_raster = TRUE, width = 8
)
cluster <- sapply(ht$feature_split, function(x) {
  switch(as.character(x),
    "C1" = "C1",
    "C2" = "C2",
    "C3" = "C3",
    "C4" = "C3",
    "C5" = "C4"
  )
})
cluster <- factor(cluster, levels = paste0("C", 1:4))
cluster <- paste0(cluster, "(", table(cluster)[cluster], ")")
names(cluster) <- names(ht$feature_split)
cluster <- factor(cluster, sort(unique(cluster)))

ht <- DynamicHeatmap(
  srt = srt_imp, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.001, use_fitted = TRUE, feature_split = cluster,
  lineages = c("Lineage_CellLineEpi"),
  cell_annotation = "SubAnnotation",
  n_split = 5, split_method = "kmeans-peaktime",
  anno_terms = TRUE, anno_features = TRUE, padjustCutoff = 1, db_version = "3.13",
  use_raster = TRUE, width = 8
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig3/imp_dynamic.ht.pdf")
saveRDS(ht, "figures/fig3/imp_dynamic.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig3/imp_dynamic.enrichment.xlsx")
srt_imp <- AnnotateFeatures(srt_imp, anno_TF = TRUE, anno_LR = TRUE, anno_db = TRUE, db = c("Chromosome"), overwrite = T)
feature_metadata <- merge(ht$feature_metadata, srt_imp[["RNA"]]@meta.features, by = 0)
feature_metadata <- feature_metadata[order(feature_metadata[["index"]]), ]
write.xlsx(feature_metadata, file = "figures/fig3/imp_dynamic.feature_metadata.xlsx")


srt_imp@tools$Enrichment_CellLineEpi <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_imp@tools$Enrichment_CellLineEpi$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 2)
p <- panel_fix(p, height = 2, width = 2, save = "figures/fig3/imp_dynamic_bp.bar.pdf", raster = TRUE)

p <- DynamicPlot(
  srt = srt_imp,
  lineages = c("Lineage_CellLineEpi"),
  features = c(
    "DNMT3B", "TERF1", "TDGF1", "L1TD1", "FOXH1", "LIN28A", "SOX2", "EOMES", "TBXT"
    # "POU5F1", "SOX2", "NANOG", # Core
    # "CD24", "SFRP2", "ZIC2", "FGF2", "SALL2", # Primed
    # "KLF7", "KLF5", "KLF4", "TFCP2L1", "DNMT3L", "ESRRB", "STAT3", "TBX3" # Naive
  ),
  group.by = "SubAnnotation",
  compare_lineages = TRUE,
  compare_features = FALSE
)
p <- panel_fix(p, height = 1.5, save = "figures/fig3/imp_dynamic_pluripotency.plot.pdf", raster = TRUE)

ExpVlnPlot(srt_E14, "SOX2", group.by = "Time", cells_subset = colnames(srt_E14)[srt_E14$Group == "EPI"])

saveRDS(srt_imp, "srt_imp_annotation.rds")

#### DE analysis -------------------------------------------------------------
srt_imp <- RunDEtest(srt_imp,
  group_by = "SubAnnotation", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_imp <- RunDEtest(srt_imp,
  group_by = "SubAnnotation", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
de_filter <- filter(srt_epi@tools$DEtest_SubAnnotation$AllMarkers_wilcox, p_val_adj < 0.05)
ribogene <- rownames(srt_epi)[grep("^RP[SL]\\d+\\w{0,1}\\d*$", rownames(srt_epi))]
de_filter <- de_filter[!de_filter$gene %in% ribogene, ]

de_top <- de_filter %>%
  group_by(gene) %>%
  top_n(1, avg_log2FC) %>%
  group_by(group1) %>%
  top_n(3, avg_log2FC)
p <- ExpDotPlot(srt_epi, features = de_top$gene, feature_split = de_top$group1, cell_split_by = c("SubAnnotation", "Time"))
ggsave(plot = p, filename = "figures/fig3/epi_topDE.heatmap.pdf", width = 10, height = 15)
p <- ExpHeatmapPlot(srt_epi,
  features = de_filter$gene, feature_split = de_filter$group1, cell_split_by = "SubAnnotation",
  cluster_genes = TRUE, cluster_cells = TRUE
)
ggsave(plot = p, filename = "figures/fig3/epi_allDE.heatmap.png", width = 10, height = 7, limitsize = FALSE)

res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "figures/fig3/epi_allDE_bp.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "figures/fig3/epi_allDE_bp.word.pdf")


### similarity(germ layer level) --------------------------------------------------------------
srt_0d10d$Annotation_split <- paste(srt_0d10d$Annotation, srt_0d10d$Time, sep = "-")
srt <- srt_0d10d[, ]
srt_0d10d <- RunKNNPredict(
  srt_query = srt_0d10d, query_group = "Annotation_split",
  srt_ref = srt_0d10d, ref_group = "Annotation_split",
  query_collapsing = T, ref_collapsing = T,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
d <- 1 - srt_0d10d@tools$knnpredict$distance_matrix
d <- as.matrix(d)

rn <- rep("d0", ncol(d))
rn[grep("d10", rownames(d))] <- "d10"
ht <- Heatmap(d,
  name = "Cosine similarity", cluster_columns = T, cluster_rows = T,
  col = colorRamp2(c(0, 0.5, 1), c("#27408B", "white", "#EE0000")),
  row_split = rn,
  column_split = rn,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = "white", lwd = 1, fill = "white")
    )
    grid.rect(x, y,
      width = w, height = h,
      gp = gpar(col = fill, lwd = 1, fill = alpha(fill, 0.5))
    )
    if (d[i, j] >= sort(d[i, ], decreasing = T)[2]) {
      # grid.text("*", x, y, gp = gpar(fontsize = 20))
      grid.text(round(d[i, j], 2), x, y, gp = gpar(fontsize = 10))
    }
  },
  border = TRUE,
  width = unit(ncol(d) * 0.7, "cm"),
  height = unit(nrow(d) * 0.7, "cm"),
)
gTree <- grid.grabExpr({
  draw(ht,
    padding = unit(c(3, 1, 1, 3), "cm")
  ) # bottom, left, top and right
})
p <- plot_grid(gTree)
ggsave(plot = p, filename = "figures/fig2/d0d10_similarity_layerlevel.pdf", width = 20, height = 20)

srt_0d <- srt_0d10d[, srt_0d10d$Time == "d0"]
srt_10d <- srt_0d10d[, srt_0d10d$Time == "d10"]
srt_0d <- RunKNNPredict(
  srt_query = srt_0d, query_group = "Annotation",
  srt_ref = srt_0d, ref_group = "Annotation",
  query_collapsing = FALSE, ref_collapsing = FALSE,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
sim_0d <- 1 - srt_0d@tools$knnpredict$distance_matrix
sim_0d <- c(as.matrix(sim_0d))

srt_10d <- RunKNNPredict(
  srt_query = srt_10d, query_group = "Annotation",
  srt_ref = srt_10d, ref_group = "Annotation",
  query_collapsing = FALSE, ref_collapsing = FALSE,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
sim_10d <- 1 - srt_10d@tools$knnpredict$distance_matrix
sim_10d <- c(as.matrix(sim_10d))

### SCVELO ------------------------------------------------------------------
adata <- srt_to_adata(srt_imp)
adata <- SCP:::RunSCVELO(
  adata = adata, group_by = "GermLayer", linear_reduction = "Uncorrectedpca", nonlinear_reduction = "FA",
  n_pcs = 25, n_neighbors = 100
)
plt <- reticulate::import("matplotlib")$pyplot
scv <- reticulate::import("scvelo")
ax <- scv$pl$velocity_embedding_stream(adata,
  vkey = "stochastic", basis = "FA", title = "RNA Velocity", color = "GermLayer",
  density = 3, linewidth = 1, arrow_size = 1.5, legend_fontsize = 15, legend_loc = "right_margin",
  save = FALSE, show = FALSE, dpi = 300
)
# ax$axis("equal")
plt$show()
plt$savefig("figures/fig3/imp.SCVELO.svg", bbox_inches = "tight")

## ESC compare ------------------------------------------------------------------
srt_E14 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_3dculture/human_gastrulation_3dculture.seurat2.rds")
srt_E14_sub <- srt_E14[, srt_E14$Group %in% c("EPI", "PSA-EPI", "PE")]
srt_E14_sub$Group <- factor(srt_E14_sub$Group)
srt_esccellline <- LoadH5Seurat("../CellCulture-project/NGSmodule_SCP_analysis/CellQC/22_ESC.filtered.h5Seurat")
srt_esccellline$Time <- "ESC"
srt_pgclc <- LoadH5Seurat("../CellCulture-project/NGSmodule_SCP_analysis/CellQC/22_PGCLC_d6.filtered.h5Seurat")
srt_d0 <- LoadH5Seurat("../CellCulture-project/NGSmodule_SCP_analysis/CellQC/22_CellLine-3.filtered.h5Seurat")
srt_d0$Time <- "D0"
srt_compare1 <- Reduce(merge, list(srt_esccellline, srt_d0))
srt_compare1 <- srt_compare[Reduce(intersect, list(rownames(srt_esccellline), rownames(srt_d0))), ]
srt_compare2 <- Reduce(merge, list(srt_esccellline, srt_pgclc, srt_d0))
srt_compare2 <- srt_compare[Reduce(intersect, list(rownames(srt_esccellline), rownames(srt_pgclc), rownames(srt_d0))), ]
srt_compare3 <- Reduce(merge, list(srt_pgclc, srt_d0))
srt_compare3 <- srt_compare[Reduce(intersect, list(rownames(srt_pgclc), rownames(srt_d0))), ]

srt_compare <- Integration_SCP(srt_compare, batch = "Time", integration_method = "Uncorrected")
srt_compare <- Integration_SCP(srt_compare, batch = "Time", integration_method = "Harmony")
srt_compare <- Integration_SCP(srt_compare, batch = "Time", integration_method = "CSS")
srt_compare <- Integration_SCP(srt_compare,
  batch = "Time", integration_method = "Scanorama",
  linear_reduction_dims_use = 1:80
)
ClassDimPlot(srt_compare, "Time", label = T)

srt_compare <- srt_compare[, !srt_compare$Scanoramaclusters %in% c(11, 12, 13, 14, 15, 4, 5)]
srt_compare <- Integration_SCP(srt_compare, batch = "Time", integration_method = "Scanorama")

srt_compare <- srt_compare[, !srt_compare$Time %in% c("ESC")]
srt_compare <- Integration_SCP(srt_compare, batch = "Time", integration_method = "CSS")
ClassDimPlot(srt_compare, "Time", reduction = "CSSUMAP", label = T)

ExpDimPlot(srt_compare, c("ITGA6", "EPCAM"),
  reduction = "ScanoramaUMAP",
  calculate_coexp = T, add_density = T
)

srt_E14_sub <- RunKNNMap(
  srt_query = srt_E14_sub, srt_ref = srt_compare,
  ref_umap = "ScanoramaUMAP2D"
)
ProjectionPlot(
  srt_query = srt_E14_sub, srt_ref = srt_compare,
  query_group = "Group", ref_group = "Time", ref_reduction = "UncorrectedUMAP2D"
)

### Annotation --------------------------------------------------------------
ClassDimPlot(srt_compare, "Annotation")
ClassDimPlot(srt_compare, "CSSclusters", label = TRUE)

srt_compare <- RunKNNMap(
  srt_query = srt_compare, srt_ref = srt_esc,
  ref_group = "CSSclusters", ref_umap = "UMAP"
)
ProjectionPlot(
  srt_query = srt_compare, srt_ref = srt_esc,
  query_group = "CSSclusters", ref_group = "Annotation"
)
ClassDimPlot(srt_compare, group.by = "CSSclusters", reduction = "ref.umap", label = T)
srt_compare <- RunKNNMap(
  srt_query = srt_compare, srt_ref = srt_h0,
  ref_group = "CSSclusters", ref_umap = "UMAP"
)
ProjectionPlot(
  srt_query = srt_compare, srt_ref = srt_h0,
  query_group = "CSSclusters", ref_group = "Annotation"
)
ClassDimPlot(srt_compare, group.by = "CSSclusters", reduction = "ref.umap", label = T)


srt_compare <- RunKNNPredict(
  srt_query = srt_compare, query_group = "CSSclusters",
  srt_ref = srt_CS7, ref_group = "sub_cluster",
  query_collapsing = T, ref_collapsing = T,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
ClassDimPlot(srt_compare, "knnpredict_classification")

srt_compare <- RunKNNPredict(
  srt_query = srt_compare, bulk_ref = SCP::ref_scHCL,
  filter_lowfreq = 10, prefix = "HCL"
)
ClassDimPlot(srt_compare, group.by = "HCL_classification", label = T)

srt_compare <- RunKNNPredict(
  srt_query = srt_compare, srt_ref = srt_CS7,
  ref_group = "sub_cluster",
  filter_lowfreq = 10, prefix = "E19"
)
ClassDimPlot(srt_compare, group.by = "E19_classification", label = T)

srt_compare <- RunKNNPredict(
  srt_query = srt_compare, srt_ref = srt_mGast1,
  ref_group = "celltype",
  filter_lowfreq = 10, prefix = "mGast1"
)
ClassDimPlot(srt_compare, group.by = "mGast1_classification", label = T)

srt_compare <- RunKNNPredict(
  srt_query = srt_compare, srt_ref = srt_mGast2,
  ref_group = "cellType",
  filter_lowfreq = 10, prefix = "mGast2"
)
ClassDimPlot(srt_compare, group.by = "mGast2_classification", label = T)

srt_compare <- RunKNNPredict(
  srt_query = srt_compare, srt_ref = srt_mOrg,
  ref_group = "Sub_trajectory_name",
  filter_lowfreq = 10, prefix = "mOrg"
)
ClassDimPlot(srt_compare, group.by = "mOrg_classification", label = T)

srt_compare <- RunKNNPredict(
  srt_query = srt_compare, srt_ref = srt_hOrg,
  ref_group = "annotation",
  filter_lowfreq = 10, prefix = "hOrg"
)
ClassDimPlot(srt_compare, group.by = "hOrg_classification", label = T)

srt_compare <- RunKNNPredict(
  srt_query = srt_compare, srt_ref = srt_mEndo1,
  ref_group = "CellType",
  filter_lowfreq = 10, prefix = "mEndo1"
)
ClassDimPlot(srt_compare, group.by = "mEndo1_classification", label = T)

srt_compare <- RunKNNPredict(
  srt_query = srt_compare, srt_ref = srt_mEndo2,
  ref_group = "celltype",
  filter_lowfreq = 10, prefix = "mEndo2"
)
ClassDimPlot(srt_compare, group.by = "mEndo2_classification", label = T)

srt_compare <- RunKNNPredict(
  srt_query = srt_compare, srt_ref = srt_mEndo3,
  ref_group = "LineageAnnotations",
  filter_lowfreq = 10, prefix = "mEndo3"
)
ClassDimPlot(srt_compare, group.by = "mEndo3_classification", label = T)

srt_compare <- RunKNNPredict(
  srt_query = srt_compare, srt_ref = srt_hEndo1,
  ref_group = "Cell_type",
  filter_lowfreq = 10, prefix = "hEndo1"
)
ClassDimPlot(srt_compare, group.by = "hEndo1_classification", label = T)

srt_compare <- RunKNNPredict(
  srt_query = srt_compare, srt_ref = srt_hEndo2,
  ref_group = "cell_name_detailed",
  filter_lowfreq = 10, prefix = "hEndo2"
)
ClassDimPlot(srt_compare, group.by = "hEndo2_classification", label = T)

srt_compare <- RunKNNPredict(
  srt_query = srt_compare, srt_ref = srt_mMeso,
  ref_group = "celltype",
  filter_lowfreq = 10, prefix = "mMeso"
)
ClassDimPlot(srt_compare, group.by = "mMeso_classification", label = T)
saveRDS(srt_compare, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_compare_annotation.rds")


### paga initiated layout ---------------------------------------------------------------
adata <- srt_to_adata(srt_compare)
adata <- RunPAGA(
  adata = adata, group_by = "Scanoramaclusters", linear_reduction = "Scanoramapca", nonlinear_reduction = "ScanoramaUMAP2D",
  n_pcs = 16, threshold = 0.1,
  embedded_with_PAGA = TRUE, paga_layout = "fa"
)
FA <- adata$obsm[["X_draw_graph_fa"]]
colnames(FA) <- c("FA_1", "FA_2")
rownames(FA) <- adata$obs_names$values
srt_compare[["FA"]] <- CreateDimReducObject(FA, key = "FA_", assay = "RNA")
ClassDimPlot(srt_compare, group.by = "Time", reduction = "FA")
srt_compare@misc$Default_reduction <- "FA"

PAGAUMAP2D <- adata$obsm[["PAGAUMAP2D"]]
colnames(PAGAUMAP2D) <- c("UMAP_1", "UMAP_2")
rownames(PAGAUMAP2D) <- adata$obs_names$values
srt_compare[["PAGAUMAP2D"]] <- CreateDimReducObject(PAGAUMAP2D, key = "UMAP_", assay = "RNA")
ClassDimPlot(srt_compare, group.by = "Time", reduction = "PAGAUMAP2D")
srt_compare@misc$Default_reduction <- "PAGAUMAP2D"


## Epiblast --------------------------------------------------------------------
srt_epi <- srt_22[, srt_22$CSSclusters %in% c(16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 44, 25, 34)]
srt_epi <- srt_epi[, !srt_epi$Time %in% c("d70", "d90")]
srt_epi$Annotation <- factor(srt_epi$Annotation, levels = levels(srt_22$Annotation))

### Integration  --------------------------------------------------------------------
srt_epi <- CellCycleScoring(srt_epi,
  s.features = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes
)

srt_epi <- Integration_SCP(
  srtMerge = srt_epi, integration_method = "Uncorrected",
  linear_reduction_dims_use = 1:25, cluster_resolution = 2,
  vars_to_regress = c("S.Score", "G2M.Score", "percent.mito")
)
srt_epi <- RunDM(srt_epi,
  features = VariableFeatures(srt_epi), slot = "data", dims = NULL,
  k = 15, ndcs = 3, reduction.name = "DM"
)
ClassDimPlot(srt_epi, c("Time", "Annotation", "Uncorrectedclusters", "Phase"), label = T)
ClassDimPlot(srt_epi, c("Annotation", "Uncorrectedclusters", "Phase"), reduction = "DM")
ExpDimPlot(srt_epi, c("percent.ribo", "percent.mito", "nFeature_RNA", "nCount_RNA"))
ClassDimPlot(srt_epi, "Uncorrectedclusters", label = TRUE)
ClassDimPlot(srt_epi, "Uncorrectedclusters", label = TRUE, reduction = "UncorrectedUMAP3D", dims = c(1, 2))

srt_epi[["UMAP"]] <- srt_epi[["UncorrectedUMAP2D"]]
srt_epi[["UMAP3D"]] <- srt_epi[["UncorrectedUMAP3D"]]
srt_epi@misc$Default_reduction <- "UMAP"
colnames(srt_epi[["UMAP"]]@cell.embeddings) <- c("UMAP_1", "UMAP_2")
colnames(srt_epi[["UMAP3D"]]@cell.embeddings) <- c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3")
Key(srt_epi[["UMAP"]]) <- "UMAP_"
Key(srt_epi[["UMAP3D"]]) <- "UMAP3D_"
srt_epi@misc$Default_reduction <- "UMAP"

srt_epi <- FindClusters(object = srt_epi, resolution = 1.5, algorithm = 1, graph.name = "Uncorrectedpca_SNN")
srt_epi <- SrtReorder(srt_epi, features = VariableFeatures(srt_epi), reorder_by = "seurat_clusters", slot = "data")
srt_epi[["seurat_clusters"]] <- NULL
srt_epi[["Uncorrectedclusters"]] <- Idents(srt_epi)
ClassDimPlot(srt_epi, c("Uncorrectedclusters"), label = TRUE)

### paga initiated layout ---------------------------------------------------------------
# resolution: 1.5
plt <- reticulate::import("matplotlib")$pyplot
sc <- reticulate::import("scanpy")
adata <- srt_to_adata(srt_epi)
adata <- RunPAGA(
  adata = adata, group_by = "Uncorrectedclusters",
  linear_reduction = "Uncorrectedpca", nonlinear_reduction = "UncorrectedUMAP2D",
  n_pcs = 25, embedded_with_PAGA = TRUE
)
FA <- adata$obsm[["X_draw_graph_fa"]]
colnames(FA) <- c("FA_1", "FA_2")
rownames(FA) <- adata$obs_names$values
srt_epi[["FA"]] <- CreateDimReducObject(FA, key = "FA_", assay = "RNA")
ClassDimPlot(srt_epi, group.by = "State", reduction = "FA")
ClassDimPlot(srt_epi, group.by = "Uncorrectedclusters", reduction = "FA")

PAGAUMAP2D <- adata$obsm[["PAGAUMAP2D"]]
colnames(PAGAUMAP2D) <- c("PAGAUMAP_1", "PAGAUMAP_2")
rownames(PAGAUMAP2D) <- adata$obs_names$values
srt_epi[["PAGAUMAP2D"]] <- CreateDimReducObject(PAGAUMAP2D, key = "PAGAUMAP_", assay = "RNA")
ClassDimPlot(srt_epi, group.by = "State", reduction = "PAGAUMAP2D")


### Monocle2 ----------------------------------------------------------------------------
library(monocle)
expr_matrix <- as(as.matrix(srt_epi@assays$RNA@counts), "sparseMatrix")
p_data <- srt_epi@meta.data
f_data <- data.frame(gene_short_name = row.names(srt_epi), row.names = row.names(srt_epi))
pd <- new("AnnotatedDataFrame", data = p_data)
fd <- new("AnnotatedDataFrame", data = f_data)
cds <- newCellDataSet(expr_matrix,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.1,
  expressionFamily = negbinomial.size()
)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
express_genes <- VariableFeatures(srt_epi)
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)
cds <- clusterCells(cds, num_clusters = 2)
cds <- reduceDimension(cds,
  max_components = 3,
  method = "DDRTree", norm_method = "log"
)
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = c("State"), cell_size = 0.8)
plot_cell_trajectory(cds, color_by = c("Cluster"), cell_size = 0.8)
plot_cell_trajectory(cds, color_by = c("Pseudotime"), cell_size = 0.8)
cds <- orderCells(cds, root_state = 2)
srt <- as.Seurat(cds)
srt_epi[["DDRTree"]] <- srt[["DDRTree"]]
srt_epi[["Cluster"]] <- srt[["Cluster"]]
srt_epi[["Pseudotime"]] <- srt[["Pseudotime"]]
srt_epi[["State"]] <- srt[["State"]]

srt_epi@tools$Monocle2 <- list(
  reducedDimS = cds@reducedDimS,
  reducedDimW = cds@reducedDimW,
  reducedDimK = cds@reducedDimK
)
srt_epi[["State"]] <- sapply(srt_epi[["State", drop = TRUE]], function(x) {
  switch(as.character(x),
    "1" = 2,
    "2" = 1,
    "3" = 3,
    "4" = 4,
    "5" = 5
  )
})
ClassDimPlot(srt_epi, "State", "DDRTree")
ClassDimPlot(srt_epi, "State", "UncorrectedUMAP2D")
ExpDimPlot(srt_epi, "Pseudotime", "UncorrectedUMAP2D")
DynamicScatter(srt_epi, group.by = "SubAnnotation", features = "SOX2", lineages = "Pseudotime")

ClassDimPlot(srt_epi, "Time", "DDRTree", split.by = "Time")
srt_epi@misc$Default_reduction <- "DDRTree"


### FA & DDRTree ------------------------------------------------------------
ClassDimPlot(srt_epi, "State", "DDRTree")
ClassDimPlot(srt_epi, "Annotation", "DDRTree")

ExpDimPlot(srt_epi, c("POU5F1", "DNMT3B", "SOX2", "ISL1", "TBXT"), "DDRTree")
srt_epi$SubAnnotation <- paste0("Epiblast-", srt_epi$State)
srt_epi$SubAnnotation[srt_epi$Annotation != "Epiblast"] <- as.character(srt_epi$Annotation[srt_epi$Annotation != "Epiblast"])
srt_epi$SubAnnotation <- factor(srt_epi$SubAnnotation, levels = c(paste0("Epiblast-", 1:5), "Primitive streak", "Mesoderm", "Amnion", "Neural ectoderm"))
ClassDimPlot(srt_epi, "SubAnnotation", "DDRTree")

srt_epi$Lineage_EpiEct <- srt_epi$Pseudotime
srt_epi$Lineage_EpiEct[!srt_epi$SubAnnotation %in% c("Epiblast-1", "Epiblast-2")] <- NA
srt_epi$Lineage_EpiBi <- srt_epi$Pseudotime
srt_epi$Lineage_EpiBi[!srt_epi$SubAnnotation %in% c("Epiblast-1", "Epiblast-3")] <- NA
srt_epi$Lineage_BiAmn <- srt_epi$Pseudotime
srt_epi$Lineage_BiAmn[!srt_epi$SubAnnotation %in% c("Epiblast-3", "Epiblast-5")] <- NA
srt_epi$Lineage_BiPS <- srt_epi$Pseudotime
srt_epi$Lineage_BiPS[!srt_epi$SubAnnotation %in% c("Epiblast-3", "Epiblast-4")] <- NA

srt_epi$Lineage_EpiEct[which(srt_epi$Lineage_EpiEct > max(srt_epi$Lineage_EpiBi, na.rm = T))] <- NA
ExpDimPlot(srt_epi, "Lineage_EpiEct")
ExpDimPlot(srt_epi, c("Lineage_EpiEct", "Lineage_EpiBi"))


#### Baisc plot ------------------------------------------------------------
srt_22_tmp <- srt_22
srt_22_tmp$GermLayer[!colnames(srt_22) %in% colnames(srt_epi)] <- NA
srt_22_tmp$GermLayer <- factor(srt_22_tmp$GermLayer, levels = levels(srt_22$GermLayer))
p <- ClassDimPlot(srt_22_tmp,
  group.by = "GermLayer", theme_use = "theme_blank",
  cells.highlight = colnames(srt_epi), show_stat = FALSE
)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_highlight.umap.pdf", raster = TRUE)

p <- ClassDimPlot(srt_epi, "Annotation", reduction = "FA", label = TRUE, theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_annotation.fa.pdf", raster = TRUE)
p <- ClassDimPlot(srt_epi, "Time", reduction = "FA", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_time.fa.pdf", raster = TRUE)
p <- ClassDimPlot(srt_epi, "Annotation", reduction = "DDRTree", label = TRUE, theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_annotation.DDRTree.pdf", raster = TRUE)
p <- ClassDimPlot(srt_epi, "Time", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_sample.DDRTree.pdf", raster = TRUE)
p <- ClassDimPlot(srt_epi, "SubAnnotation", label = TRUE, show_stat = FALSE, theme_use = "theme_blank", force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_subannotation.DDRTree.pdf", raster = TRUE)
p <- ExpDimPlot(srt_epi, "Pseudotime", palette = "cividis", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_pseudotime.DDRTree.pdf", raster = TRUE)

p <- ExpDimPlot(srt_epi, c("Lineage_EpiEct", "Lineage_EpiBi", "Lineage_BiAmn", "Lineage_BiPS"),
  palette = "cividis", byrow = FALSE, theme_use = "theme_blank", force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig3/epi_lineage.DDRTree.pdf", raster = TRUE)
p <- ExpDimPlot(srt_epi, c("POU5F1", "DNMT3B", "MYC", "ISL1", "TBXT", "SOX2"),
  theme_use = "theme_blank", nrow = 2
)
p <- panel_fix(p, height = 2, save = "figures/fig3/epi_markers.DDRTree.pdf", raster = TRUE)

library(ggridges)
df <- srt_epi[[]]
levels <- df %>%
  group_by(SubAnnotation) %>%
  summarise(m = median(Pseudotime)) %>%
  arrange(desc(m)) %>%
  pull("SubAnnotation") %>%
  as.character()
df$y <- factor(df$SubAnnotation, levels = levels)
p <- ggplot(df, aes(x = Pseudotime, y = y, fill = SubAnnotation)) +
  scale_fill_manual(values = palette_scp(df$SubAnnotation)) +
  geom_density_ridges() +
  labs(y = "") +
  theme_scp(aspect.ratio = 0.5)
p <- panel_fix(p, height = 2, save = "figures/fig3/epi_pseudotime.ridges.pdf", raster = TRUE)


#### Compare with human E14 --------------------------------------------------
srt_E14 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_3dculture/human_gastrulation_3dculture.seurat.rds")
srt_E14_sub <- srt_E14[, srt_E14$Group %in% c("EPI")]
srt_epi <- RunKNNPredict(
  srt_query = srt_epi, srt_ref = srt_E14_sub,
  query_group = "Annotation", ref_group = "Time",
  return_full_distance_matrix = TRUE
)
ClassDimPlot(srt_epi, "knnpredict_classification")

srt_E14_sub <- RunKNNMap(
  srt_query = srt_E14_sub, srt_ref = srt_epi,
  ref_group = "SubAnnotation", ref_umap = "FA"
)
srt_E14_sub$Group <- factor(srt_E14_sub$Group, levels = c("EPI"))
srt_E14_sub$cell_type <- srt_E14_sub$Group
srt_epi$bg <- NA
p <- ProjectionPlot(
  srt_query = srt_E14_sub, srt_ref = srt_epi,
  query_group = "Time", ref_group = "bg",
  query_param = list(palette = "Set1"),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "FA_1", ylab = "FA_2",
    title = "3D-cultured human pre-gastrulation embryos"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_E14.fa.pdf")

#### Compare with human CS7 --------------------------------------------------------
srt_CS7 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_E19/human_gastrulation_E19.seurat.rds")
srt_CS7[["percent.mito"]] <- PercentageFeatureSet(object = srt_CS7, pattern = paste0(c("^MT-", "^Mt-", "^mt-"), collapse = "|"))
srt_CS7[["percent.ribo"]] <- PercentageFeatureSet(object = srt_CS7, pattern = paste0(c("^RP[SL]\\d+\\w{0,1}\\d*$", "^Rp[sl]\\d+\\w{0,1}\\d*$", "^rp[sl]\\d+\\w{0,1}\\d*$"), collapse = "|"))
srt_CS7_sub <- srt_CS7[, srt_CS7$cluster_id %in% c("Epiblast")]
srt_CS7_sub <- RunKNNMap(
  srt_query = srt_CS7_sub, srt_ref = srt_epi,
  ref_group = "SubAnnotation", ref_umap = "FA"
)
srt_CS7_sub$cell_type <- srt_CS7_sub$cluster_id

srt_E14_sub$Epiblast <- "Xiang et al."
srt_CS7_sub$Epiblast <- "RCV Tyser et al."
srt <- merge(srt_CS7_sub, srt_E14_sub)
srt[["ref.umap"]] <- CreateDimReducObject(embeddings = rbind(
  srt_E14_sub[["ref.umap"]]@cell.embeddings,
  srt_CS7_sub[["ref.umap"]]@cell.embeddings
))

srt_epi$bg <- NA
p <- ProjectionPlot(
  srt_query = srt, srt_ref = srt_epi,
  query_group = "Epiblast", ref_group = "bg",
  query_reduction = "ref.umap", ref_reduction = "FA",
  query_param = list(palette = "Set1"),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank",
    title = ""
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_merge.fa.pdf")


#### Dynamic analysis --------------------------------------------------------
srt_epi <- RunDynamicFeatures(srt_epi,
  lineages = c("Lineage_EpiEct", "Lineage_EpiBi", "Lineage_BiAmn", "Lineage_BiPS"),
  n_candidates = 10000
)
saveRDS(srt_epi, "srt_epi_annotation.rds")

colnames(srt_epi@tools$DynamicFeatures_Lineage_EpiEct$DynamicFeatures)[5] <- "peaktime"
colnames(srt_epi@tools$DynamicFeatures_Lineage_EpiBi$DynamicFeatures)[5] <- "peaktime"
ht <- DynamicHeatmap(
  srt = srt_epi, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE, split_method = "mfuzz",
  lineages = c("Lineage_EpiEct", "Lineage_EpiBi"),
  cell_annotation = "SubAnnotation",
  n_split = 4, reverse_ht = 1, use_raster = TRUE, width = 8
)
cluster <- sapply(ht$feature_split, function(x) {
  switch(as.character(x),
    "C1" = "C2",
    "C2" = "C1",
    "C3" = "C3",
    "C4" = "C4"
  )
})
cluster <- factor(cluster, levels = paste0("C", 1:4))
cluster <- paste0(cluster, "(", table(cluster)[cluster], ")")
names(cluster) <- names(ht$feature_split)
cluster <- factor(cluster, sort(unique(cluster)))
ht <- DynamicHeatmap(
  srt = srt_epi, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE, split_method = "mfuzz", feature_split = cluster,
  lineages = c("Lineage_EpiEct", "Lineage_EpiBi"),
  cell_annotation = "SubAnnotation",
  anno_terms = TRUE, anno_features = TRUE, padjustCutoff = 1, db_version = "3.13",
  n_split = 4, reverse_ht = 1, use_raster = TRUE, width = 8
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig3/epi_dynamic1.ht.pdf")
saveRDS(ht, "figures/fig3/epi_dynamic1.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig3/epi_dynamic1.enrichment.xlsx")
srt_epi <- AnnotateFeatures(srt_epi, anno_TF = TRUE, anno_LR = TRUE, anno_db = TRUE, db = c("Chromosome"), overwrite = TRUE)
feature_metadata <- merge(ht$feature_metadata, srt_epi[["RNA"]]@meta.features, by = 0)
feature_metadata <- feature_metadata[order(feature_metadata[["index"]]), ]
write.xlsx(feature_metadata, file = "figures/fig3/epi_dynamic1.feature_metadata.xlsx")

srt_epi@tools$Enrichment_dynamic1 <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_epi@tools$Enrichment_dynamic1$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 2)
p <- panel_fix(p, height = 2, width = 2, save = "figures/fig3/epi_dynamic1_bp.bar.pdf", raster = TRUE)
a <- srt_epi@tools$Enrichment_dynamic1$enrichment

p <- DynamicPlot(
  srt = srt_epi,
  lineages = c("Lineage_EpiEct", "Lineage_EpiBi"),
  features = c("POU5F1", "DNMT3B", "MYCN"),
  group.by = "SubAnnotation",
  compare_lineages = T, add_point = TRUE,
  compare_features = FALSE
)

colnames(srt_epi@tools$DynamicFeatures_Lineage_BiAmn$DynamicFeatures)[5] <- "peaktime"
colnames(srt_epi@tools$DynamicFeatures_Lineage_BiPS$DynamicFeatures)[5] <- "peaktime"
ht <- DynamicHeatmap(
  srt = srt_epi, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 1, use_fitted = TRUE,
  lineages = c("Lineage_BiAmn", "Lineage_BiPS"),
  cell_annotation = "SubAnnotation",
  anno_terms = TRUE, anno_features = TRUE, db_version = "3.13",
  n_split = 4, reverse_ht = 1, use_raster = TRUE, width = 8
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig3/epi_dynamic2.ht.pdf")
saveRDS(ht, "figures/fig3/epi_dynamic2.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig3/epi_dynamic2.enrichment.xlsx")
feature_metadata <- merge(ht$feature_metadata, srt_epi[["RNA"]]@meta.features, by = 0)
feature_metadata <- feature_metadata[order(feature_metadata[["index"]]), ]
write.xlsx(feature_metadata, file = "figures/fig3/epi_dynamic2.feature_metadata.xlsx")

srt_epi@tools$Enrichment_dynamic2 <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_epi@tools$Enrichment_dynamic2$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 2)
p <- panel_fix(p, height = 2, width = 2, save = "figures/fig3/epi_dynamic2_bp.bar.pdf")
b <- srt_epi@tools$Enrichment_dynamic2$enrichment

# glycolytic process
gene1 <- a[a$Description == "glycolytic process" & a$Groups == "C4(232)", "geneID"]
gene1 <- unlist(strsplit(gene1, "/"))
gene2 <- b[b$Description == "oxidative phosphorylation" & b$Groups == "C2", "geneID"]
gene2 <- unlist(strsplit(gene2, "/"))
p <- DynamicPlot(
  srt = srt_epi,
  lineages = c("Lineage_EpiEct", "Lineage_EpiBi"),
  features = gene1, line_pal_color = palette_scp(n = 4, palette = "Dark2")[1:2],
  group.by = "SubAnnotation",
  compare_lineages = TRUE,
  compare_features = FALSE,
  ncol = 3
)
p <- panel_fix(p, height = 1.5, save = "figures/fig3/epi_dynamicplot.glycolytic.pdf", raster = TRUE)

p <- DynamicPlot(
  srt = srt_epi,
  lineages = c("Lineage_BiAmn", "Lineage_BiPS"),
  features = gene2, line_pal_color = palette_scp(n = 4, palette = "Dark2")[3:4],
  group.by = "SubAnnotation",
  compare_lineages = TRUE,
  compare_features = FALSE,
  ncol = 3
)
p <- panel_fix(p, height = 1.5, save = "figures/fig3/epi_dynamicplot.oxidative_phos.pdf", raster = TRUE)
saveRDS(srt_epi, "srt_epi_annotation.rds")

#### DE analysis -------------------------------------------------------------
srt_epi <- RunDEtest(srt_epi,
  group_by = "SubAnnotation", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_epi <- RunDEtest(srt_epi,
  group_by = "SubAnnotation", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
de_filter <- filter(srt_epi@tools$DEtest_SubAnnotation$AllMarkers_wilcox, p_val_adj < 0.05)
ribogene <- rownames(srt_epi)[grep("^RP[SL]\\d+\\w{0,1}\\d*$", rownames(srt_epi))]
de_filter <- de_filter[!de_filter$gene %in% ribogene, ]

de_top <- de_filter %>%
  group_by(gene) %>%
  top_n(1, avg_log2FC) %>%
  group_by(group1) %>%
  top_n(3, avg_log2FC)
p <- ExpDotPlot(srt_epi, features = de_top$gene, feature_split = de_top$group1, cell_split_by = c("SubAnnotation", "Time"))
ggsave(plot = p, filename = "figures/fig3/epi_topDE.heatmap.pdf", width = 10, height = 15)
p <- ExpHeatmapPlot(srt_epi,
  features = de_filter$gene, feature_split = de_filter$group1, cell_split_by = "SubAnnotation",
  cluster_genes = TRUE, cluster_cells = TRUE
)
ggsave(plot = p, filename = "figures/fig3/epi_allDE.heatmap.png", width = 10, height = 7, limitsize = FALSE)

res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "figures/fig3/epi_allDE_bp.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "figures/fig3/epi_allDE_bp.word.pdf")


### UMAP ------------------------------------------------------------
#### Annotation --------------------------------------------------------------
# resolution=1.8
ClassDimPlot(srt_epi, c("Uncorrectedclusters"), label = TRUE)
ExpDimPlot(srt = srt_epi, "SOX2")
# oral-pharyngeal ectoderm/preplacodal ectoderm : WNT9A HES1/SOX3 DLX2 DLX3 SIX4 GBX2 FOXI3 FOXD3 ZIC1 PITX3 CXCR4 SOX3

class <- "2,3,4,5,6: Epiblast-1
1,11: Epiblast-2
7,9,10,12,13,14: Epiblast-3
15: Preplacodal ectoderm
16,17: Neural ectoderm
19,20: Amniotic ectoderm
8,18: Primitive streak & Mesoderm"
class <- readLines(textConnection(class))
srt_epi[["SubAnnotation"]] <- srt_epi[["Uncorrectedclusters"]]
srt_epi[["SubAnnotation"]] <- as.character(srt_epi[["SubAnnotation", drop = TRUE]])
for (i in 1:length(class)) {
  m <- strsplit(class[i], ": ")[[1]]
  srt_epi$SubAnnotation[srt_epi$SubAnnotation %in% strsplit(m[1], ",")[[1]]] <- m[[2]]
}
levels <- c(
  "Epiblast-1", "Epiblast-2",
  "Epiblast-3",
  "Preplacodal ectoderm", "Neural ectoderm",
  "Amniotic ectoderm", "Primitive streak & Mesoderm"
)
all(levels %in% srt_epi$SubAnnotation)
srt_epi$SubAnnotation <- factor(srt_epi$SubAnnotation, levels = levels)
ClassDimPlot(srt_epi, "SubAnnotation")
ClassDimPlot3D(srt_epi, "SubAnnotation")

#### Baisc plot --------------------------------------------------------------------
srt_22$bg <- srt_22$Annotation
srt_22$bg[setdiff(colnames(srt_22), colnames(srt_epi))] <- NA
p <- ClassDimPlot(srt_22,
  group.by = "bg",
  label = FALSE, label_insitu = FALSE, theme_use = "theme_blank",
  cells.highlight = colnames(srt_epi)
)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_highlight.umap.pdf")

p <- ClassDimPlot(srt_epi, "Annotation", label = TRUE, theme_use = "theme_blank", force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_annotation.umap.pdf")
p <- ClassDimPlot(srt_epi, "SubAnnotation", label = TRUE, theme_use = "theme_blank", force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_subannotation.umap.pdf")
p <- ClassDimPlot(srt_epi, "Uncorrectedclusters", label = TRUE, theme_use = "theme_blank", force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_clusters.umap.pdf")
p <- ClassDimPlot(srt_epi, "Time", theme_use = "theme_blank", force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_sample.umap.pdf")
p <- ExpDimPlot(srt_epi, c("SIX4", "DLX3", "FOXI3"), theme_use = "theme_blank", nrow = 1)
p <- panel_fix(p, height = 2, save = "figures/fig3/epi_PPEmarker.umap.pdf")

p <- ClassDimPlot(srt_epi, "Annotation",
  split.by = "Time", legend.position = "none", bg_color = "grey90", nrow = 2,
  theme_use = "theme_blank", force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig3/epi_sample_split.umap.pdf")

#### Cell cycle --------------------------------------------------------------
p <- ClassDimPlot(srt_epi, "Phase", theme_use = "theme_blank", force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_cellcycle_phase.umap.pdf")

df <- table(srt_epi$Phase, srt_epi$SubAnnotation)
df <- t(t(df) / colSums(df))
df <- reshape2::melt(df)
colnames(df) <- c("Phase", "CellType", "Proportion")
p <- ggplot(df, aes(x = CellType, y = Proportion, fill = Phase)) +
  geom_col(color = "black") +
  labs(x = "") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = palette_scp(df$Phase)) +
  theme_scp() +
  theme(
    aspect.ratio = 0.5,
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
p <- panel_fix(p, height = 2, save = "figures/fig3/epi_cellcycle.stat.pdf")


#### PAGA ---------------------------------------------------------------------
adata <- srt_to_adata(srt_epi)
adata <- RunPAGA(
  adata = adata, group_by = "SubAnnotation", linear_reduction = "Uncorrectedpca", nonlinear_reduction = "UncorrectedUMAP2D",
  n_pcs = 25, embedded_with_PAGA = FALSE
)
plt <- reticulate::import("matplotlib")$pyplot
sc <- reticulate::import("scanpy")
ax <- sc$pl$paga_compare(adata,
  edges = FALSE,
  threshold = 0.1, min_edge_width = 0.8, basis = "UncorrectedUMAP2D", title = "PAGA",
  legend_loc = NULL, legend_fontweight = "bold", legend_fontsize = 10, legend_fontoutline = 1.5,
  frameon = FALSE, save = FALSE, show = T
)
ax[[1]]$axis("none")
ax[[2]]$axis("equal") ## adjust window
plt$show()
plt$savefig("figures/fig3/epi.PAGA.svg", bbox_inches = "tight")

#### SCVELO ------------------------------------------------------------------
adata <- srt_to_adata(srt_epi)
adata <- SCP:::RunSCVELO(
  adata = adata, group_by = "SubAnnotation", linear_reduction = "Uncorrectedpca", nonlinear_reduction = "UncorrectedUMAP2D",
  n_pcs = 25, n_neighbors = 30
)
plt <- reticulate::import("matplotlib")$pyplot
scv <- reticulate::import("scvelo")
ax <- scv$pl$velocity_embedding_stream(adata,
  vkey = "stochastic", basis = "UncorrectedUMAP2D", title = "RNA Velocity", color = "SubAnnotation",
  density = 3, linewidth = 1, arrow_size = 1.5, legend_fontsize = 15, legend_loc = "right_margin",
  save = FALSE, show = FALSE, dpi = 300
)
ax$axis("equal")
plt$show()
plt$savefig("figures/fig3/epi.SCVELO.svg", bbox_inches = "tight")

#### DE analysis -------------------------------------------------------------
srt_epi <- RunDEtest(srt_epi,
  group_by = "Uncorrectedclusters", test.use = "wilcox", fc.threshold = 1.5,
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
de_filter <- filter(srt_epi@tools$DEtest_Uncorrectedclusters$AllMarkers_wilcox, p_val_adj < 0.05)
annotation_list <- list(
  "Epiblast-1" = c(2, 3, 4, 5, 6),
  "Epiblast-2" = c(1, 11),
  "Epiblast-3" = c(7, 9, 10, 12, 13, 14),
  "Preplacodal ectoderm" = c(15),
  "Neural ectoderm" = c(16, 17),
  "Amniotic ectoderm" = c(19, 20),
  "Primitive streak & Mesoderm" = c(8, 18)
)
annotation <- setNames(rep(names(annotation_list), sapply(annotation_list, length)), nm = unlist(annotation_list))
de_filter$group1 <- annotation[as.character(de_filter$group1)]
de_filter$group1 <- factor(de_filter$group1, levels = names(annotation_list))

srt_epi <- RunDEtest(srt_epi,
  group_by = "SubAnnotation", test.use = "wilcox", fc.threshold = 1.2,
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
de_filter <- filter(srt_epi@tools$DEtest_SubAnnotation$AllMarkers_wilcox, p_val_adj < 0.05)
ribogene <- rownames(srt_epi)[grep("^RP[SL]\\d+\\w{0,1}\\d*$", rownames(srt_epi))]
de_filter <- de_filter[!de_filter$gene %in% ribogene, ]

# srt_epi <- RunDEtest(srt_epi,
#                      group_by = "SubAnnotation", test.use = "roc", fc.threshold = 1.2,
#                      BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
# )
# de_roc <- filter(srt_epi@tools$DEtest_SubAnnotation$AllMarkers_roc, power > 0.7)
# de_filter <- de_filter[de_filter$gene %in% de_roc$gene, ]

de_top <- de_filter %>%
  group_by(gene) %>%
  top_n(1, avg_log2FC) %>%
  group_by(group1) %>%
  top_n(3, avg_log2FC)
p <- ExpDotPlot(srt_epi, features = de_top$gene, feature_split = de_top$group1, cell_split_by = c("SubAnnotation", "Time"))
ggsave(plot = p, filename = "figures/fig3/epi_topDE.heatmap.pdf", width = 10, height = 15)
p <- ExpHeatmapPlot(srt_epi,
  features = de_filter$gene, feature_split = de_filter$group1, cell_split_by = "SubAnnotation",
  cluster_genes = FALSE, cluster_cells = TRUE
)
ggsave(plot = p, filename = "figures/fig3/epi_allDE.heatmap.png", width = 15, height = 10, limitsize = FALSE)

res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "figures/fig3/epi_allDE_bp.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "figures/fig3/epi_allDE_bp.word.pdf")

#### Compare with human E14 --------------------------------------------------
srt_E14 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_3dculture/human_gastrulation_3dculture.seurat.rds")
srt_E14_sub <- srt_E14[, srt_E14$Group %in% c("EPI")]
srt_epi <- RunKNNPredict(
  srt_query = srt_epi, srt_ref = srt_E14_sub,
  query_group = "Uncorrectedclusters", ref_group = "Time",
  return_full_distance_matrix = TRUE
)
ClassDimPlot(srt_epi, "knnpredict_Time")
srt_epi$knnpredict_Time <- factor(srt_epi$knnpredict_Time, levels = levels(srt_E14_sub$Time))
d <- 1 - srt_E14_sub@tools$knnpredict_Uncorrectedclusters$distance_matrix
d <- t(scale(t(d)))
Heatmap(as.matrix(d), cluster_columns = F)

srt_E14_sub <- RunKNNMap(
  srt_query = srt_E14_sub, srt_ref = srt_epi,
  ref_group = "SubAnnotation", ref_umap = "UMAP"
)
srt_E14_sub$Group <- factor(srt_E14_sub$Group, levels = c("EPI"))
srt_E14_sub$cell_type <- srt_E14_sub$Group
srt_epi$bg <- NA
p <- ProjectionPlot(
  srt_query = srt_E14_sub, srt_ref = srt_epi,
  query_group = "Time", ref_group = "bg",
  query_param = list(palette = "Set1"),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "3D-cultured human pre-gastrulation embryos"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_E14.umap.pdf")

#### Compare with human CS7 --------------------------------------------------------
srt_CS7 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_E19/human_gastrulation_E19.seurat.rds")
srt_CS7[["percent.mito"]] <- PercentageFeatureSet(object = srt_CS7, pattern = paste0(c("^MT-", "^Mt-", "^mt-"), collapse = "|"))
srt_CS7[["percent.ribo"]] <- PercentageFeatureSet(object = srt_CS7, pattern = paste0(c("^RP[SL]\\d+\\w{0,1}\\d*$", "^Rp[sl]\\d+\\w{0,1}\\d*$", "^rp[sl]\\d+\\w{0,1}\\d*$"), collapse = "|"))
srt_CS7_sub <- srt_CS7[, srt_CS7$cluster_id %in% c("Epiblast")]
srt_CS7_sub <- RunKNNMap(
  srt_query = srt_CS7_sub, srt_ref = srt_epi,
  ref_group = "Uncorrectedclusters", ref_umap = "UncorrectedUMAP2D"
)
srt_CS7_sub$cell_type <- srt_CS7_sub$cluster_id

srt_E14_sub$Epiblast <- "Xiang et al."
srt_CS7_sub$Epiblast <- "RCV Tyser et al."
srt <- merge(srt_CS7_sub, srt_E14_sub)
srt[["ref.umap"]] <- CreateDimReducObject(embeddings = rbind(
  srt_E14_sub[["ref.umap"]]@cell.embeddings,
  srt_CS7_sub[["ref.umap"]]@cell.embeddings
))

srt_epi$bg <- NA
p <- ProjectionPlot(
  srt_query = srt, srt_ref = srt_epi,
  query_group = "Epiblast", ref_group = "bg",
  query_reduction = "ref.umap", ref_reduction = "UMAP",
  query_param = list(palette = "Set1"),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank",
    title = ""
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig3/epi_merge.umap.pdf")



### hc --------------------------------------------------------------
data.avg <- AverageExpression(
  object = srt_epi, features = VariableFeatures(srt_epi),
  slot = "data", assays = "RNA", group.by = "Uncorrectedclusters", verbose = FALSE
)[[1]]
data.avg <- log1p(data.avg)
mat <- t(x = data.avg[VariableFeatures(srt_epi), ])
mat <- as(mat, "dgCMatrix")
d <- dist(as(mat, "dgCMatrix"), method = "euclidean")
data.dist <- as.dist(d)
hc <- hclust(d = data.dist, method = "complete")
dd <- as.dendrogram(hc)
plot(hc)


dm1 <- srt_epi@reductions$UncorrectedDM2D@cell.embeddings[, 1]
dm.results <- srt_epi@reductions$UncorrectedDM2D@misc$dm.results
dpt <- destiny::DPT(dm.results, tips = which.min(dm1))
plot(dpt)
plot(dpt, col_by = "branch")

gene.relevance <- destiny::gene_relevance(
  coords = dm.results@eigenvectors,
  exprs = t(srt_epi@assays$RNA@data)
)
srt_epi$dpt <- dpt[["dpt"]]
ExpDimPlot(srt_epi, c("dpt"), reduction = "CSSDM2D", palette = "Reds")

df <- srt_epi@meta.data
library(ggforce)
ggplot(df, aes(x = Pseudotime, fill = SubAnnotation)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = palette_scp(df$CSSclusters)) +
  scale_fill_brewer(palette = "Paired") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  facet_zoom(xlim = c(0, 3)) +
  theme_classic()

library(ggridges)
levels <- df %>%
  group_by(SubAnnotation) %>%
  summarise(m = median(Pseudotime)) %>%
  arrange(desc(m)) %>%
  pull("SubAnnotation") %>%
  as.character()
df$CSSclusters2 <- factor(df$SubAnnotation, levels = levels)
ggplot(df, aes(x = Pseudotime, y = CSSclusters2, fill = SubAnnotation)) +
  scale_fill_manual(values = palette_scp(df$SubAnnotation)) +
  geom_density_ridges() +
  theme_scp(aspect.ratio = 0.8)

ClassDimPlot(srt_epi, c("Time", "Annotation", "Uncorrectedclusters"), reduction = "DM")
ExpDimPlot(srt_epi22, c("TBXT", "EOMES", "ISL1", "ABCG2", "SOX2", "POU5F1", "DPPA3", "DNMT3B"))


# Fig4 --------------------------------------------------------------------
dir.create(paste0(work_dir, "/figures/fig4"), recursive = T, showWarnings = FALSE)

## Primitive streak, nascent mesoderm and endoderm --------------------------------------------------------------------
epi <- colnames(srt_epi)[srt_epi$SubAnnotation %in% paste0("Epiblast-", 4)]
cells <- c(epi, colnames(srt_22)[srt_22$CSSclusters %in% c(27, 34, 33, 32, 30, 31)])
srt_endo <- srt_22[, cells]

# cells <- colnames(srt_22)[srt_22$CSSclusters %in% c(33, 32, 30, 31)]
# srt_endo_only <- srt_22[, cells]
# srt_endo_only <- Integration_SCP(
#   srtMerge = srt_endo_only, integration_method = "Uncorrected",
#   linear_reduction_dims_use = 1:30, cluster_resolution = 3
# )
# srt_endo_only$SubAnnotation <- srt_endo$SubAnnotation[colnames(srt_endo_only)]
# ClassDimPlot(srt_endo_only, "SubAnnotation", label = T)
# ClassDimPlot(srt_endo_only, "Uncorrectedclusters", label = T)
# ClassDimPlot(srt_endo, "SubAnnotation",
#   cells.highlight = WhichCells(srt_endo_only, expression = Uncorrectedclusters == 9)
# )
# ExpDimPlot(srt_endo_only, "CDX2")


### Integration  --------------------------------------------------------------------
srt_endo <- CellCycleScoring(srt_endo,
  s.features = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes
)

plist <- list()
pagalist <- list()
progress <- 1
total <- length(c(2000, 3000)) * length(seq(15, 50, 5)) * length(seq(0.5, 5, 0.5))
for (nhvf in c(2000, 3000)) {
  for (i in seq(15, 50, 5)) {
    for (j in seq(0.5, 5, 0.5)) {
      cat("progress:", progress, "/", total, "\n")
      progress <- progress + 1
      srt_endo <- CSS_integrate(
        srtMerge = srt_endo, nHVF = nhvf,
        linear_reduction_dims_use = 1:i, cluster_resolution = j,
        CSS_param = list(cluster_resolution = j)
      )
      srt_endo[[paste0("CSSUMAP2Dnhvf", nhvf, "dim", i, "res", j)]] <- srt_endo@reductions$CSSUMAP2D
      srt_endo[[paste0("CSSUMAP3Dnhvf", nhvf, "dim", i, "res", j)]] <- srt_endo@reductions$CSSUMAP3D
      plist[[paste0("CSSUMAP2Dnhvf", nhvf, "dim", i, "res", j)]] <- ClassDimPlot(srt_endo, "Annotation", title = paste("hvf:", nhvf, "dim", i, "res", j))

      adata <- srt_to_adata(srt_endo)
      adata <- RunPAGA(
        adata = adata, group_by = "CSSclusters", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
        n_pcs = ncol(srt_endo@reductions$CSS@cell.embeddings), threshold = 0.1,
        embedded_with_PAGA = TRUE, paga_layout = "fa"
      )
      PAGAUMAP2D <- adata$obsm[["PAGAUMAP2D"]]
      colnames(PAGAUMAP2D) <- c("UMAP_1", "UMAP_2")
      rownames(PAGAUMAP2D) <- adata$obs_names$values
      srt_endo[["PAGAUMAP2D"]] <- CreateDimReducObject(PAGAUMAP2D, key = "UMAP_", assay = "RNA")
      ClassDimPlot(srt_endo, group.by = "Time", reduction = "PAGAUMAP2D")
      pagalist[[paste0("PAGAUMAP2Dnhvf", nhvf, "dim", i, "res", j)]] <- ClassDimPlot(srt_endo, "Annotation", reduction = "PAGAUMAP2D", title = paste("hvf:", nhvf, "dim", i, "res", j))
    }
  }
}
p <- plot_grid(plotlist = plist)
p <- panel_fix(p, save = "tmp_endo1.png", width = 3)
p <- plot_grid(plotlist = pagalist)
p <- panel_fix(p, save = "tmp_endo2.png", width = 3)

srt_endo <- Integration_SCP(
  srtMerge = srt_endo, integration_method = "Uncorrected",
  linear_reduction_dims_use = 1:25, cluster_resolution = 5
)
srt_endo <- Integration_SCP(
  srtMerge = srt_endo, integration_method = "CSS",
  linear_reduction_dims_use = 1:25, cluster_resolution = 5,
  CSS_param = list(cluster_resolution = 5)
)
ClassDimPlot(srt_endo, c("Annotation", "CSSclusters", "Phase", "Time"), label = T) # 25-5
ExpDimPlot(srt_endo, c("percent.ribo", "percent.mito", "nFeature_RNA", "nCount_RNA"))

# srt_endo <- FindClusters(object = srt_endo, resolution = 10, algorithm = 1, graph.name = "CSS_SNN")
# srt_endo <- SrtReorder(srt_endo, features = VariableFeatures(srt_endo), reorder_by = "seurat_clusters", slot = "data")
# srt_endo[["seurat_clusters"]] <- NULL
# srt_endo[["CSSclusters"]] <- Idents(srt_endo)
# ClassDimPlot(srt_endo, c("CSSclusters"), label = TRUE)

### paga initiated layout ---------------------------------------------------------------
adata <- srt_to_adata(srt_endo)
adata <- RunPAGA(
  adata = adata, group_by = "CSSclusters", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
  n_pcs = ncol(srt_endo@reductions$CSS@cell.embeddings), threshold = 0.1,
  embedded_with_PAGA = TRUE, paga_layout = "fa"
)
PAGAUMAP2D <- adata$obsm[["PAGAUMAP2D"]]
colnames(PAGAUMAP2D) <- c("UMAP_1", "UMAP_2")
rownames(PAGAUMAP2D) <- adata$obs_names$values
srt_endo[["PAGAUMAP2D"]] <- CreateDimReducObject(PAGAUMAP2D, key = "UMAP_", assay = "RNA")
ClassDimPlot(srt_endo, group.by = "Time", reduction = "PAGAUMAP2D")
srt_endo[["UMAP"]] <- srt_endo[["PAGAUMAP2D"]]
srt_endo@misc$Default_reduction <- "UMAP"

### Annotation --------------------------------------------------------------
ClassDimPlot(srt_endo, "Annotation")
ClassDimPlot(srt_endo, "CSSclusters", label = TRUE)
ExpDimPlot(srt_endo, c("POU5F1", "TBXT", "EOMES", "NANOS3", "TFAP2C", "TFAP2A", "ABCG2", "ISL1", "TP63"))

srt_endo1 <- RunKNNMap(
  srt_query = srt_endo, srt_ref = srt_esc,
  ref_group = "CSSclusters", ref_umap = "UMAP"
)
ProjectionPlot(
  srt_query = srt_endo1, srt_ref = srt_esc,
  query_group = "CSSclusters", ref_group = "Annotation"
)
ClassDimPlot(srt_endo1, group.by = "CSSclusters", reduction = "ref.umap", label = T)
srt_endo2 <- RunKNNMap(
  srt_query = srt_endo, srt_ref = srt_h0,
  ref_group = "CSSclusters", ref_umap = "UMAP"
)
ProjectionPlot(
  srt_query = srt_endo2, srt_ref = srt_h0,
  query_group = "CSSclusters", ref_group = "Annotation"
)
ClassDimPlot(srt_endo2, group.by = "CSSclusters", reduction = "ref.umap", label = T)
ClassDimPlot(srt_endo2,
  group.by = "CSSclusters", reduction = "ref.umap", label = T,
  cells.highlight = WhichCells(srt_endo2, expression = CSSclusters %in% c(6, 7))
)

srt_endo <- RunKNNPredict(
  srt_query = srt_endo, query_group = "CSSclusters",
  srt_ref = srt_endo, ref_group = "CSSclusters",
  query_collapsing = T, ref_collapsing = T,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
d <- 1 - srt_endo@tools$knnpredict_classification$distance_matrix
d <- as.matrix(d)
Heatmap(d)

srt_endo <- RunKNNPredict(
  srt_query = srt_meso, query_group = "CSSclusters",
  srt_ref = srt_CS7, ref_group = "sub_cluster",
  query_collapsing = T, ref_collapsing = T,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
ClassDimPlot(srt_endo, "knnpredict_sub_cluster")

srt_endo <- RunKNNPredict(
  srt_query = srt_endo, bulk_ref = SCP::ref_scHCL,
  filter_lowfreq = 10, prefix = "HCL"
)
ClassDimPlot(srt_endo, group.by = "HCL_classification", label = T)

srt_endo <- RunKNNPredict(
  srt_query = srt_endo, srt_ref = srt_CS7,
  ref_group = "sub_cluster",
  filter_lowfreq = 10, prefix = "E19"
)
ClassDimPlot(srt_endo, group.by = "E19_classification", label = T)

srt_endo <- RunKNNPredict(
  srt_query = srt_endo, srt_ref = srt_mGast,
  ref_group = "celltype",
  filter_lowfreq = 10, prefix = "mGast"
)
ClassDimPlot(srt_endo, group.by = "mGast_classification", label = T)

srt_endo <- RunKNNPredict(
  srt_query = srt_endo, srt_ref = srt_mOrg,
  ref_group = "Sub_trajectory_name",
  filter_lowfreq = 10, prefix = "mOrg"
)
ClassDimPlot(srt_endo, group.by = "mOrg_classification", label = T)

srt_endo <- RunKNNPredict(
  srt_query = srt_endo, srt_ref = srt_mEndo1,
  ref_group = "CellType",
  filter_lowfreq = 10, prefix = "mEndo1"
)
ClassDimPlot(srt_endo, group.by = "mEndo1_classification", label = T)

srt_endo <- RunKNNPredict(
  srt_query = srt_endo, srt_ref = srt_mEndo2,
  ref_group = "celltype",
  filter_lowfreq = 10, prefix = "mEndo2"
)
ClassDimPlot(srt_endo, group.by = "mEndo2_classification", label = T)

srt_endo <- RunKNNPredict(
  srt_query = srt_endo, srt_ref = srt_mEndo3,
  ref_group = "LineageAnnotations_DESM",
  filter_lowfreq = 10, prefix = "mEndo3"
)
ClassDimPlot(srt_endo, group.by = "mEndo3_classification", label = T)

srt_endo <- RunKNNPredict(
  srt_query = srt_endo, srt_ref = srt_hEndo1,
  ref_group = "Cell_type",
  filter_lowfreq = 10, prefix = "hEndo1"
)
ClassDimPlot(srt_endo, group.by = "hEndo1_classification", label = T)

srt_endo <- RunKNNPredict(
  srt_query = srt_endo, srt_ref = srt_hEndo2,
  ref_group = "cell_name_detailed",
  filter_lowfreq = 10, prefix = "hEndo2"
)
ClassDimPlot(srt_endo, group.by = "hEndo2_classification", label = T)

# KRT4 KRT13 KRT17
# dermis: ACTA2 COL1A1
# epidermis: KRT5 KRT10

annotation_list <- list(
  "Epiblast" = c(16, 18),
  "Primitive streak" = c(17, 19, 11, 12),
  "Nascent mesoderm" = c(13, 14, 15),
  "Primitive endoderm" = c(27, 28, 29),
  "Definitive endoderm" = c(26, 25, 24, 22),
  "Gut endoderm-1" = c(21, 20, 23),
  "Gut endoderm-2" = c(1, 2),
  "ExE endoderm" = c(3, 4),
  "Mid/Hindgut" = c(5, 8, 9),
  "Foregut" = c(6, 7, 10)
)
annotation <- setNames(rep(names(annotation_list), sapply(annotation_list, length)), nm = unlist(annotation_list))
srt_endo$SubAnnotation <- annotation[as.character(srt_endo$CSSclusters)]
srt_endo$SubAnnotation <- factor(srt_endo$SubAnnotation, levels = names(annotation_list))
ClassDimPlot(srt_endo, "SubAnnotation", label = T, cells.highlight = WhichCells(srt_endo, expression = SubAnnotation == "Primitive streak"))
srt_endo$SubAnnotation[srt_endo$SubAnnotation == "Primitive streak" & srt_endo$Annotation == "Gut"] <- "Foregut"


# ExpDimPlot(srt_endo, c(APS_marker, PPS_marker), "UMAP", nrow = 2)
srt_endo <- CellClassification(srt_endo,
  features = list(Anterior = APS_marker, Posterior = PPS_marker),
  prefix = "AP"
)
ClassDimPlot(srt_endo, "AP_classification")
srt_endo$AP_classification[is.na(srt_endo$AP_classification)] <- "Posterior"
srt_endo$SubAnnotation <- as.character(srt_endo$SubAnnotation)
srt_endo$SubAnnotation[srt_endo$SubAnnotation == "Primitive streak"] <- srt_endo$AP_classification[srt_endo$SubAnnotation == "Primitive streak"]
srt_endo$SubAnnotation[srt_endo$SubAnnotation == "Anterior"] <- "Anterior primitive streak"
srt_endo$SubAnnotation[srt_endo$SubAnnotation == "Posterior"] <- "Posterior primitive streak"
srt_endo$SubAnnotation <- factor(srt_endo$SubAnnotation,
  levels = c(
    "Epiblast",
    "Anterior primitive streak",
    "Posterior primitive streak",
    "Nascent mesoderm",
    "Primitive endoderm",
    "Definitive endoderm",
    "Gut endoderm-1",
    "Gut endoderm-2",
    "ExE endoderm",
    "Mid/Hindgut",
    "Foregut"
  )
)
ClassDimPlot(srt_endo, "SubAnnotation", label = T)
srt_endo_sub1 <- srt_endo[, srt_endo$SubAnnotation %in% c(
  "Epiblast",
  "Anterior primitive streak",
  "Posterior primitive streak",
  "Nascent mesoderm",
  "Definitive endoderm"
)]
srt_endo_sub1 <- RunSlingshot(srt_endo_sub1,
  group.by = "SubAnnotation", reduction = "UMAP", prefix = "EpiDE",
  start = c("Epiblast"),
  end = c("Nascent mesoderm", "Definitive endoderm")
)
srt_endo_sub1$Lineage_EpiDE <- srt_endo_sub1$EpiDE_Lineage1
srt_endo_sub1$Lineage_EpiNM <- srt_endo_sub1$EpiDE_Lineage2
ExpDimPlot(srt_endo_sub1, c("Lineage_EpiDE", "Lineage_EpiNM"), palette = "cividis")

saveRDS(srt_endo, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_endo_annotation.rds")


### Baisc plot ------------------------------------------------------------
srt_22_tmp <- srt_22
srt_22_tmp$GermLayer[!colnames(srt_22) %in% colnames(srt_endo)] <- NA
srt_22_tmp$GermLayer <- factor(srt_22_tmp$GermLayer, levels = levels(srt_22$GermLayer))
p <- ClassDimPlot(srt_22_tmp,
  group.by = "GermLayer", theme_use = "theme_blank",
  cells.highlight = colnames(srt_endo), show_stat = FALSE
)
p <- panel_fix(p, height = 3, save = "figures/fig4/endo_highlight.umap.pdf", raster = TRUE)

p <- ClassDimPlot(srt_endo, "SubAnnotation", label = TRUE, theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig4/endo_annotation.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_endo, "Time", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig4/endo_time.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_endo, "Phase", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig4/endo_cellcycle.umap.pdf", raster = TRUE)

markers_endo <- list(
  "Epiblast" = c("POU5F1", "DNMT3B"),
  "Primitive streak" = c("TBXT", "FGF4"),
  "Nascent mesoderm" = c("MESP1", "MIXL1"),
  "Primitive endoderm" = c("CXCR4", "GATA6"),
  "Definitive endoderm" = c("SOX17", "FOXA3"),
  "Gut endoderm" = c("KITLG", "HNF4A"),
  "ExE endoderm" = c("AFP", "TTR"),
  "Mid/Hindgut" = c("CDX2", "CDX1"),
  "Foregut" = c("CD47", "SOX2")
)
p <- ExpDimPlot(srt_endo, unlist(markers_endo),
  theme_use = "theme_blank", ncol = 4
)
p <- panel_fix(p, height = 2, save = "figures/fig4/endo_markers.umap.pdf", raster = TRUE)

ht <- GroupHeatmap(srt_endo,
  group.by = "SubAnnotation", heatmap_palette = "YlOrRd",
  features = unlist(markers_endo),
  feature_split = rep(names(markers_endo), each = 2),
  show_row_names = TRUE, add_dot = TRUE, add_bg = TRUE,
  height = 6, width = 6, dot_size = unit(7, "mm")
)
panel_fix(ht$plot, save = "figures_add/endoderm.markers.pdf")

p <- ExpDimPlot(srt_endo, c("DMKN", "FN1", "KITLG", "HNF4A"),
  theme_use = "theme_blank", ncol = 4
)
p <- panel_fix(p, height = 2, save = "figures/fig4/endo_mix_markers.umap.pdf", raster = TRUE)

p <- ExpDimPlot(srt_endo, c("CER1", "LEFTY1", "LEFTY2", "NODAL", "FZD5"),
  theme_use = "theme_blank", ncol = 4, calculate_coexp = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig4/endo_ave_markers.umap.pdf", raster = TRUE)



p <- ExpDotPlot(srt_endo,
  genes = unlist(markers_endo),
  feature_split = unlist(lapply(names(markers_endo), function(x) rep(x, length(markers_endo[[x]])))),
  cell_split_by = "SubAnnotation"
)
ggsave(plot = p, filename = "figures/fig2/22_markers.dotplot.pdf", width = 10, height = 15)

p <- ExpDimPlot(srt_endo, c(APS_marker, PPS_marker), theme_use = "theme_blank")
p <- panel_fix(p, height = 2, save = "figures/fig4/endo_ap.umap.pdf", raster = TRUE)

p <- ExpDimPlot(srt_endo, c("AP_Anterior", "AP_Posterior"), theme_use = "theme_blank")
p <- panel_fix(p, height = 2, save = "figures/fig4/endo_ap_score.umap.pdf", raster = TRUE)

### PAGA ---------------------------------------------------------------------
srt_endo_paga <- RunPAGA(
  srt = srt_endo, group_by = "SubAnnotation", linear_reduction = "CSS", nonlinear_reduction = "UMAP",
  n_pcs = ncol(srt_endo@reductions$CSS@cell.embeddings), n_neighbors = 100, embedded_with_PAGA = FALSE, return_seurat = TRUE
)
srt_endo@misc$paga <- srt_endo_paga@misc$paga
p <- ClassDimPlot(srt_endo,
  group.by = "SubAnnotation", pt.size = 5, pt.alpha = 0.1,
  label = TRUE, label.size = 3, label_repel = TRUE, label_insitu = TRUE, label_segment_color = "transparent",
  paga = srt_endo@misc$paga, paga_edge_threshold = 0.1, paga_edge_color = "black", paga_edge_alpha = 1,
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
)
p <- panel_fix(p, width = 4, save = "figures/fig4/endo_paga.pdf", raster = TRUE)

### SCVELO ------------------------------------------------------------------
srt_endo_scv <- SCP:::RunSCVELO(
  srt = srt_endo, group_by = "SubAnnotation", linear_reduction = "CSS", nonlinear_reduction = "UMAP",
  n_pcs = ncol(srt_endo@reductions$CSS@cell.embeddings), n_neighbors = 100, return_seurat = TRUE
)
srt_endo[["stochastic_UMAP"]] <- srt_endo_scv[["stochastic_UMAP"]]
srt_endo[["Ms"]] <- srt_endo_scv[["Ms"]]
srt_endo[["Mu"]] <- srt_endo_scv[["Mu"]]
srt_endo[["stochastic"]] <- srt_endo_scv[["stochastic"]]
srt_endo[["variance_stochastic"]] <- srt_endo_scv[["variance_stochastic"]]
p <- ClassDimPlot(srt_endo,
  group.by = "SubAnnotation", pt.size = 5, pt.alpha = 0.1,
  velocity = "stochastic", velocity_plot_type = "stream",
  velocity_density = 2, velocity_smooth = 1,
  streamline_n = 15, streamline_size = 0.5, streamline_color = "black",
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
)
p <- panel_fix(p, width = 4, save = "figures/fig4/endo_scvelo.pdf", raster = TRUE)

saveRDS(srt_endo, "srt_endo_annotation.rds")


### Compare with human CS7 --------------------------------------------------------
srt_CS7 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_E19/human_gastrulation_E19.seurat.rds")
srt_CS7[["percent.mito"]] <- PercentageFeatureSet(object = srt_CS7, pattern = paste0(c("^MT-", "^Mt-", "^mt-"), collapse = "|"))
srt_CS7[["percent.ribo"]] <- PercentageFeatureSet(object = srt_CS7, pattern = paste0(c("^RP[SL]\\d+\\w{0,1}\\d*$", "^Rp[sl]\\d+\\w{0,1}\\d*$", "^rp[sl]\\d+\\w{0,1}\\d*$"), collapse = "|"))
levels <- c("Epiblast", "Primitive Streak", "Endoderm", "Nascent Mesoderm")
srt_CS7_sub <- srt_CS7[, srt_CS7$cluster_id %in% levels]
srt_CS7_sub <- FindVariableFeatures(srt_CS7_sub)
srt_CS7_sub <- RunKNNMap(
  srt_query = srt_CS7_sub, srt_ref = srt_endo,
  ref_group = "CSSclusters", ref_umap = "UMAP"
)
srt_CS7_sub$cluster_id <- factor(srt_CS7_sub$cluster_id, levels = levels)
srt_CS7_sub$cell_type <- srt_CS7_sub$cluster_id

srt_endo$bg <- NA
p <- ProjectionPlot(
  srt_query = srt_CS7_sub, srt_ref = srt_endo,
  query_group = "cell_type", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "CS7 human gastrula"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig4/endo_CS7_projection_celltype.umap.pdf", raster = TRUE)

p <- ProjectionPlot(
  srt_query = srt_CS7_sub, srt_ref = srt_endo,
  query_group = "spatial", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "CS7 human gastrula"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig4/endo_CS7_projection_spatial.umap.pdf", raster = TRUE)

### Slingshot ---------------------------------------------------------------
ClassDimPlot(srt_endo, "SubAnnotation", "UMAP", label = T)

#### srt_endo_sub1 -----------------------------------------------------------
srt_endo_sub1 <- srt_endo[, srt_endo$SubAnnotation %in% c(
  "Epiblast",
  "Anterior primitive streak",
  "Posterior primitive streak",
  "Nascent mesoderm",
  "Definitive endoderm"
)]
srt_endo_sub1 <- RunSlingshot(srt_endo_sub1,
  group.by = "SubAnnotation", reduction = "UMAP", prefix = "EpiDE",
  start = c("Epiblast"),
  end = c("Nascent mesoderm", "Definitive endoderm")
)
srt_endo_sub1$Lineage_EpiDE <- srt_endo_sub1$EpiDE_Lineage1
srt_endo_sub1$Lineage_EpiNM <- srt_endo_sub1$EpiDE_Lineage2
ExpDimPlot(srt_endo_sub1, c("Lineage_EpiDE", "Lineage_EpiNM"), palette = "cividis")

srt_endo_sub1 <- RunDynamicFeatures(srt_endo_sub1,
  lineages = c("Lineage_EpiDE", "Lineage_EpiNM"),
  n_candidates = 10000
)
colnames(srt_endo_sub1@tools$DynamicFeatures_Lineage_EpiDE$DynamicFeatures)[5] <- "peaktime"
colnames(srt_endo_sub1@tools$DynamicFeatures_Lineage_EpiNM$DynamicFeatures)[5] <- "peaktime"

srt_endo_sub1 <- AnnotateFeatures(srt_endo_sub1, anno_TF = T)
ht <- DynamicHeatmap(
  srt = srt_endo_sub1, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_EpiDE", "Lineage_EpiNM"),
  cell_annotation = "SubAnnotation",
  feature_annotation = c("TF", "TF_cofactors"), feature_annotation_palcolor = list(c("black", "transparent")),
  anno_terms = TRUE, anno_features = TRUE, db_version = "3.13",
  n_split = 4, reverse_ht = 1, use_raster = TRUE, width = 8
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig4/endo_lineageNMDE.ht.pdf")
saveRDS(ht, "figures/fig4/endo_lineageNMDE.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig4/endo_lineageNMDE.enrichment.xlsx")
write.xlsx(ht$feature_metadata, file = "figures/fig4/endo_lineageNMDE.feature_metadata.xlsx")

p <- ExpDimPlot(srt_endo_sub1,
  features = c("Lineage_EpiDE", "Lineage_EpiNM"), ncol = 1,
  palette = "cividis", byrow = FALSE, theme_use = "theme_blank", force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig4/endo_lineageNMDE.pseudotime.umap.pdf", raster = TRUE)

srt_endo_sub1@tools$Enrichment_NMDE <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(
  res = srt_endo_sub1@tools$Enrichment_NMDE$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE,
  ncol = 2
)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig4/endo_lineageNMDE_bp.bar.pdf", raster = TRUE)


ht <- DynamicHeatmap(
  srt = srt_endo_sub1, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_EpiDE"),
  cell_annotation = "SubAnnotation",
  n_split = NULL, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig4/endo_lineageEpiDE.ht.pdf")

ht <- DynamicHeatmap(
  srt = srt_endo_sub1, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_EpiNM"),
  cell_annotation = "SubAnnotation",
  n_split = NULL, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig4/endo_lineageEpiNM.ht.pdf")

saveRDS(srt_endo_sub1, "srt_endo_sub1.rds")

#### srt_endo_sub2 -----------------------------------------------------------
srt_endo_sub2 <- srt_endo[, srt_endo$SubAnnotation %in% c(
  "Primitive endoderm",
  "Definitive endoderm",
  "Gut endoderm-1"
)]
srt_endo_sub2 <- RunSlingshot(srt_endo_sub2,
  group.by = "SubAnnotation", reduction = "UMAP", prefix = "PE",
  start = c("Gut endoderm-1"),
  end = c("Primitive endoderm", "Definitive endoderm"),
  reverse = TRUE
)
ExpDimPlot(srt_endo_sub2, paste0("PE_Lineage", 1:2), palette = "cividis")
srt_endo_sub2$Lineage_PE <- srt_endo_sub2$PE_Lineage1
srt_endo_sub2$Lineage_DE <- srt_endo_sub2$PE_Lineage2
ExpDimPlot(srt_endo_sub2, c("Lineage_PE", "Lineage_DE"), palette = "cividis")
srt_endo_sub2 <- RunDynamicFeatures(srt_endo_sub2,
  lineages = c("Lineage_PE", "Lineage_DE"),
  n_candidates = 10000
)

colnames(srt_endo_sub2@tools$DynamicFeatures_Lineage_PE$DynamicFeatures)[5] <- "peaktime"
colnames(srt_endo_sub2@tools$DynamicFeatures_Lineage_DE$DynamicFeatures)[5] <- "peaktime"
ht <- DynamicHeatmap(
  srt = srt_endo_sub2, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_PE", "Lineage_DE"),
  cell_annotation = "SubAnnotation",
  n_split = 6, reverse_ht = 2, use_raster = TRUE, width = 8
)
cluster <- sapply(ht$feature_split, function(x) {
  switch(as.character(x),
    "C1" = "C3",
    "C2" = "C1",
    "C3" = "C2",
    "C4" = "C4",
    "C5" = "C5",
    "C6" = "C6"
  )
})
cluster <- factor(cluster, levels = paste0("C", 1:6))
cluster <- paste0(cluster, "(", table(cluster)[cluster], ")")
names(cluster) <- names(ht$feature_split)
cluster <- factor(cluster, sort(unique(cluster)))
ht <- DynamicHeatmap(
  srt = srt_endo_sub2, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE, feature_split = cluster,
  lineages = c("Lineage_PE", "Lineage_DE"),
  cell_annotation = "SubAnnotation",
  anno_terms = TRUE, anno_features = TRUE, db_version = "3.13",
  n_split = 6, reverse_ht = 2, use_raster = TRUE, width = 8
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig4/endo_lineagePEDE.ht.pdf")
saveRDS(ht, "figures/fig4/endo_lineagePEDE.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig4/endo_lineagePEDE.enrichment.xlsx")

p <- ExpDimPlot(srt_endo_sub2,
  features = c("Lineage_PE", "Lineage_DE"), ncol = 1,
  palette = "cividis", byrow = FALSE, theme_use = "theme_blank", force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig4/endo_lineagePEDE.pseudotime.umap.pdf", raster = TRUE)
srt_endo_sub2@tools$Enrichment_PEDE <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(
  res = srt_endo_sub2@tools$Enrichment_PEDE$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE,
  ncol = 3
)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig4/endo_lineagePEDE_bp.bar.pdf", raster = TRUE)

ht <- DynamicHeatmap(
  srt = srt_endo_sub2, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_PE"),
  cell_annotation = "SubAnnotation",
  n_split = NULL, use_raster = TRUE, width = 5
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig4/endo_lineagePE.ht.pdf")
ht <- DynamicHeatmap(
  srt = srt_endo_sub2, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_DE"),
  cell_annotation = "SubAnnotation",
  n_split = NULL, use_raster = TRUE, width = 5
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig4/endo_lineageDE.ht.pdf")

saveRDS(srt_endo_sub2, "srt_endo_sub2.rds")

#### srt_endo_sub3 -----------------------------------------------------------
srt_endo_sub3 <- srt_endo[, srt_endo$SubAnnotation %in% c(
  "Gut endoderm-1",
  "Gut endoderm-2",
  "ExE endoderm",
  "Mid/Hindgut",
  "Foregut"
)]
srt_endo_sub3 <- RunSlingshot(srt_endo_sub3,
  group.by = "SubAnnotation", reduction = "UMAP", prefix = "DE",
  start = c("Gut endoderm-1"),
  end = c("ExE endoderm", "Mid/Hindgut", "Foregut")
)
ExpDimPlot(srt_endo_sub3, paste0("DE_Lineage", 1:3), palette = "cividis")
srt_endo_sub3$Lineage_MHG <- srt_endo_sub3$DE_Lineage1
srt_endo_sub3$Lineage_ExE <- srt_endo_sub3$DE_Lineage2
srt_endo_sub3$Lineage_FG <- srt_endo_sub3$DE_Lineage3
ExpDimPlot(srt_endo_sub3, c("Lineage_FG", "Lineage_MHG", "Lineage_ExE"), palette = "cividis")
srt_endo_sub3 <- RunSlingshot(srt_endo_sub3,
  group.by = "SubAnnotation", reduction = "UMAP", prefix = "possible",
  start = c("Gut endoderm-1")
)
ExpDimPlot(srt_endo_sub3, c("possible_Lineage1"), palette = "cividis")
srt_endo_sub3 <- RunDynamicFeatures(srt_endo_sub3,
  lineages = c("Lineage_FG", "Lineage_MHG", "Lineage_ExE"),
  n_candidates = 10000
)
ht <- DynamicHeatmap(
  srt = srt_endo_sub3, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("possible_Lineage1"),
  cell_annotation = "SubAnnotation",
  n_split = NULL, reverse_ht = NULL, use_raster = TRUE, width = 8
)
p <- ht$plot

colnames(srt_endo_sub3@tools$DynamicFeatures_Lineage_FG$DynamicFeatures)[5] <- "peaktime"
colnames(srt_endo_sub3@tools$DynamicFeatures_Lineage_MHG$DynamicFeatures)[5] <- "peaktime"
colnames(srt_endo_sub3@tools$DynamicFeatures_Lineage_ExE$DynamicFeatures)[5] <- "peaktime"
ht <- DynamicHeatmap(
  srt = srt_endo_sub3, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_FG", "Lineage_MHG", "Lineage_ExE"),
  cell_annotation = "SubAnnotation",
  anno_terms = TRUE, anno_features = TRUE, db_version = "3.13",
  n_split = 6, reverse_ht = NULL, use_raster = TRUE, width = 10
)
p <- ht$plot
ggsave(p, height = 12, width = 25, filename = "figures/fig4/endo_lineageGut.ht.pdf")
saveRDS(ht, "figures/fig4/endo_lineageGut.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig4/endo_lineageGut.enrichment.xlsx")


p <- ExpDimPlot(srt_endo_sub3,
  features = c("Lineage_FG", "Lineage_MHG", "Lineage_ExE"), ncol = 3,
  palette = "cividis", byrow = FALSE, theme_use = "theme_blank", force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig4/endo_lineageGut.pseudotime.umap.pdf", raster = TRUE)
srt_endo_sub3@tools$Enrichment_Gut <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(
  res = srt_endo_sub3@tools$Enrichment_Gut$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE,
  ncol = 3
)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig4/endo_lineageGut_bp.bar.pdf", raster = TRUE)


ht <- DynamicHeatmap(
  srt = srt_endo_sub3, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_FG"), split_method = "kmeans-peaktime",
  cell_annotation = "SubAnnotation",
  n_split = 6, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig4/endo_lineageFG.ht.pdf")
srt_endo_sub3@tools$Enrichment_FG <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(
  res = srt_endo_sub3@tools$Enrichment_FG$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE,
  ncol = 2
)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig4/endo_lineageFG_bp.bar.pdf", raster = TRUE)

ht <- DynamicHeatmap(
  srt = srt_endo_sub3, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_MHG"), split_method = "kmeans-peaktime",
  cell_annotation = "SubAnnotation",
  n_split = 6, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig4/endo_lineageMHG.ht.pdf")
srt_endo_sub3@tools$Enrichment_MHG <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(
  res = srt_endo_sub3@tools$Enrichment_MHG$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE,
  ncol = 2
)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig4/endo_lineageMHG_bp.bar.pdf", raster = TRUE)

ht <- DynamicHeatmap(
  srt = srt_endo_sub3, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_ExE"), split_method = "kmeans-peaktime",
  cell_annotation = "SubAnnotation",
  n_split = 6, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig4/endo_lineageExE.ht.pdf")
srt_endo_sub3@tools$Enrichment_ExE <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(
  res = srt_endo_sub3@tools$Enrichment_ExE$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE,
  ncol = 2
)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig4/endo_lineageExE_bp.bar.pdf", raster = TRUE)

# ht <- DynamicHeatmap(
#   srt = srt_endo_sub3, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
#   cell_density = 0.01, use_fitted = TRUE,
#   lineages = c("Lineage_MHG","Lineage_ExE"),
#   cell_annotation = "SubAnnotation",reverse_ht = 1,
#   n_split = 4, use_raster = TRUE, width = 7
# )
# p <- ht$plot
# ggsave(p, height = 12, width = 12, filename = "figures/fig4/endo_lineageMHGExE.ht.pdf")
# srt_endo_sub3@tools$Enrichment_MHGExE <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
# p <- EnrichmentPlot(
#   res = srt_endo_sub3@tools$Enrichment_MHGExE$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE,
#   ncol = 2
# )
# p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig4/endo_lineageMHGExE_bp.bar.pdf", raster = TRUE)

saveRDS(srt_endo_sub3, "srt_endo_sub3.rds")


#### all lieages ############
srt_endo[["Lineage_EpiDE"]] <- srt_endo_sub1[["Lineage_EpiDE"]]
srt_endo[["Lineage_EpiNM"]] <- srt_endo_sub1[["Lineage_EpiNM"]]
srt_endo[["Lineage_PE"]] <- srt_endo_sub2[["Lineage_PE"]]
srt_endo[["Lineage_DE"]] <- srt_endo_sub2[["Lineage_DE"]]
srt_endo[["Lineage_FG"]] <- srt_endo_sub3[["Lineage_FG"]]
srt_endo[["Lineage_MHG"]] <- srt_endo_sub3[["Lineage_MHG"]]
srt_endo[["Lineage_ExE"]] <- srt_endo_sub3[["Lineage_ExE"]]

srt_endo2 <- srt_endo
srt_endo2$SubAnnotation[setdiff(colnames(srt_endo2), colnames(srt_endo_sub1))] <- NA
p <- ClassDimPlot(srt_endo2, "SubAnnotation",
  cells.highlight = colnames(srt_endo_sub1),
  lineages = c("Lineage_EpiDE", "Lineage_EpiNM"),
  show_stat = FALSE, theme_use = "theme_blank"
)
p <- panel_fix(p, height = 3, save = "figures/fig4/endo_lineageNMDE.umap.pdf", raster = TRUE)

srt_endo2 <- srt_endo
srt_endo2$SubAnnotation[setdiff(colnames(srt_endo2), colnames(srt_endo_sub2))] <- NA
p <- ClassDimPlot(srt_endo2, "SubAnnotation",
  cells.highlight = colnames(srt_endo_sub2),
  lineages = c("Lineage_PE", "Lineage_DE"),
  show_stat = FALSE, theme_use = "theme_blank"
)
p <- panel_fix(p, height = 3, save = "figures/fig4/endo_lineagePEDE.umap.pdf", raster = TRUE)

srt_endo2 <- srt_endo
srt_endo2$SubAnnotation[setdiff(colnames(srt_endo2), colnames(srt_endo_sub3))] <- NA
p <- ClassDimPlot(srt_endo2, "SubAnnotation",
  cells.highlight = colnames(srt_endo_sub3),
  lineages = c("Lineage_FG", "Lineage_MHG", "Lineage_ExE"),
  show_stat = FALSE, theme_use = "theme_blank"
)
p <- panel_fix(p, height = 3, save = "figures/fig4/endo_lineageGut.umap.pdf", raster = TRUE)

p <- ClassDimPlot(srt_endo2, "SubAnnotation",
  lineages = c(
    "Lineage_EpiDE", "Lineage_EpiNM",
    "Lineage_PE", "Lineage_DE",
    "Lineage_FG", "Lineage_MHG", "Lineage_ExE"
  ),
  show_stat = FALSE, theme_use = "theme_blank",
  lineages_trim = c(0.1, 0.9)
)
panel_fix(p, height = 3, raster = T, save = "figures_add/endoderm.lineages.pdf")

### DE analysis ------------------------------------------------------------------
srt_endo <- RunDEtest(srt_endo,
  group_by = "CSSclusters", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_endo <- RunDEtest(srt_endo,
  group_by = "CSSclusters", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
de_filter <- filter(srt_endo@tools$DEtest_CSSclusters$AllMarkers_wilcox, p_val_adj < 0.05)
res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "tmp_endo_bp.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "tmp_endo_bp.word.pdf")

srt_endo <- RunDEtest(srt_endo,
  group_by = "SubAnnotation", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_endo <- RunDEtest(srt_endo,
  group_by = "SubAnnotation", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
de_filter <- filter(srt_endo@tools$DEtest_SubAnnotation$AllMarkers_wilcox, p_val_adj < 0.05)
res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "tmp_endo_bp2.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "tmp_endo_bp2.word.pdf")

## Mesoderm --------------------------------------------------------------------
epi <- colnames(srt_epi)[srt_epi$SubAnnotation %in% paste0("Epiblast-", 4)]
cells <- c(epi, colnames(srt_22)[srt_22$CSSclusters %in% c(27, 34:40, 45:55)])
srt_meso <- srt_22[, cells]

### Integration  --------------------------------------------------------------------
srt_meso <- CellCycleScoring(srt_meso,
  s.features = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes
)

plist <- list()
pagalist <- list()
progress <- 1
total <- length(c(2000, 3000)) * length(seq(15, 50, 5)) * length(seq(0.5, 5, 0.5))
for (nhvf in c(2000)) {
  for (i in seq(15, 50, 5)) {
    for (j in seq(0.5, 5, 0.5)) {
      cat("progress:", progress, "/", total, "\n")
      progress <- progress + 1
      srt_meso <- CSS_integrate(
        srtMerge = srt_meso, nHVF = nhvf,
        linear_reduction_dims_use = 1:i, cluster_resolution = j,
        CSS_param = list(cluster_resolution = j)
      )
      srt_meso[[paste0("CSSUMAP2Dnhvf", nhvf, "dim", i, "res", j)]] <- srt_meso@reductions$CSSUMAP2D
      srt_meso[[paste0("CSSUMAP3Dnhvf", nhvf, "dim", i, "res", j)]] <- srt_meso@reductions$CSSUMAP3D
      plist[[paste0("CSSUMAP2Dnhvf", nhvf, "dim", i, "res", j)]] <- ClassDimPlot(srt_meso, "Annotation", title = paste("hvf:", nhvf, "dim", i, "res", j))

      adata <- srt_to_adata(srt_meso)
      adata <- RunPAGA(
        adata = adata, group_by = "CSSclusters", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
        n_pcs = ncol(srt_meso@reductions$CSS@cell.embeddings), threshold = 0.1,
        embedded_with_PAGA = TRUE, paga_layout = "fa"
      )
      PAGAUMAP2D <- adata$obsm[["PAGAUMAP2D"]]
      colnames(PAGAUMAP2D) <- c("UMAP_1", "UMAP_2")
      rownames(PAGAUMAP2D) <- adata$obs_names$values
      srt_meso[["PAGAUMAP2D"]] <- CreateDimReducObject(PAGAUMAP2D, key = "UMAP_", assay = "RNA")
      ClassDimPlot(srt_meso, group.by = "Time", reduction = "PAGAUMAP2D")
      pagalist[[paste0("PAGAUMAP2Dnhvf", nhvf, "dim", i, "res", j)]] <- ClassDimPlot(srt_meso, "Annotation", reduction = "PAGAUMAP2D", title = paste("hvf:", nhvf, "dim", i, "res", j))
    }
  }
}
p <- plot_grid(plotlist = plist)
p <- panel_fix(p, save = "tmp_meso3.png", width = 3)
p <- plot_grid(plotlist = pagalist)
p <- panel_fix(p, save = "tmp_meso4.png", width = 3)

srt_meso <- Integration_SCP(
  srtMerge = srt_meso, integration_method = "Uncorrected",
  linear_reduction_dims_use = 1:40, cluster_resolution = 5
)
srt_meso <- Integration_SCP(
  srtMerge = srt_meso, integration_method = "CSS",
  linear_reduction_dims_use = 1:40, cluster_resolution = 5,
  CSS_param = list(cluster_resolution = 5)
)
ClassDimPlot(srt_meso, c("Annotation", "CSSclusters", "Phase", "Time"), label = T) # 25-5
ExpDimPlot(srt_meso, c("percent.ribo", "percent.mito", "nFeature_RNA", "nCount_RNA"))

# srt_meso <- FindClusters(object = srt_meso, resolution = 5, algorithm = 1, graph.name = "CSS_SNN")
# srt_meso <- SrtReorder(srt_meso, features = VariableFeatures(srt_meso), reorder_by = "seurat_clusters", slot = "data")
# srt_meso[["seurat_clusters"]] <- NULL
# srt_meso[["CSSclusters"]] <- Idents(srt_meso)
# ClassDimPlot(srt_meso, c("CSSclusters"), label = TRUE)

### paga initiated layout ---------------------------------------------------------------
adata <- srt_to_adata(srt_meso)
adata <- RunPAGA(
  adata = adata, group_by = "CSSclusters", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
  n_pcs = ncol(srt_meso@reductions$CSS@cell.embeddings), threshold = 0.1,
  embedded_with_PAGA = TRUE, paga_layout = "fa"
)
PAGAUMAP2D <- adata$obsm[["PAGAUMAP2D"]]
colnames(PAGAUMAP2D) <- c("UMAP_1", "UMAP_2")
rownames(PAGAUMAP2D) <- adata$obs_names$values
srt_meso[["PAGAUMAP2D"]] <- CreateDimReducObject(PAGAUMAP2D, key = "UMAP_", assay = "RNA")
ClassDimPlot(srt_meso, group.by = "Time", reduction = "PAGAUMAP2D")
srt_meso[["UMAP"]] <- srt_meso[["PAGAUMAP2D"]]
srt_meso@misc$Default_reduction <- "UMAP"


### Annotation --------------------------------------------------------------
ClassDimPlot(srt_meso, "Annotation")
ClassDimPlot(srt_meso, "CSSclusters", label = TRUE)

srt_meso1 <- RunKNNMap(
  srt_query = srt_meso, srt_ref = srt_esc,
  ref_group = "CSSclusters", ref_umap = "UMAP"
)
ProjectionPlot(
  srt_query = srt_meso1, srt_ref = srt_esc,
  query_group = "CSSclusters", ref_group = "Annotation"
)
ClassDimPlot(srt_meso1, group.by = "CSSclusters", reduction = "ref.umap", label = T)
srt_meso2 <- RunKNNMap(
  srt_query = srt_meso, srt_ref = srt_h0,
  ref_group = "CSSclusters", ref_umap = "UMAP"
)
ProjectionPlot(
  srt_query = srt_meso2, srt_ref = srt_h0,
  query_group = "CSSclusters", ref_group = "Annotation"
)
ClassDimPlot(srt_meso2, group.by = "CSSclusters", reduction = "ref.umap", label = T)
ClassDimPlot(srt_meso2,
  group.by = "CSSclusters", reduction = "ref.umap", label = T,
  cells.highlight = WhichCells(srt_meso2, expression = CSSclusters %in% c(6, 7))
)

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, query_group = "CSSclusters",
  srt_ref = srt_meso, ref_group = "CSSclusters",
  query_collapsing = T, ref_collapsing = T,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
d <- 1 - srt_meso@tools$knnpredict_classification$distance_matrix
d <- as.matrix(d)
Heatmap(d)

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, query_group = "CSSclusters",
  srt_ref = srt_CS7, ref_group = "sub_cluster",
  query_collapsing = T, ref_collapsing = T,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
ClassDimPlot(srt_meso, "knnpredict_classification")

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, bulk_ref = SCP::ref_scHCL,
  filter_lowfreq = 10, prefix = "HCL"
)
ClassDimPlot(srt_meso, group.by = "HCL_classification", label = T)

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, srt_ref = srt_CS7,
  ref_group = "sub_cluster",
  filter_lowfreq = 10, prefix = "E19"
)
ClassDimPlot(srt_meso, group.by = "E19_classification", label = T)

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, srt_ref = srt_mGast1,
  ref_group = "celltype",
  filter_lowfreq = 10, prefix = "mGast1"
)
ClassDimPlot(srt_meso, group.by = "mGast1_classification", label = T)

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, srt_ref = srt_mGast2,
  ref_group = "cellType",
  filter_lowfreq = 10, prefix = "mGast2"
)
ClassDimPlot(srt_meso, group.by = "mGast2_classification", label = T)

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, srt_ref = srt_mOrg,
  ref_group = "Sub_trajectory_name",
  filter_lowfreq = 10, prefix = "mOrg"
)
ClassDimPlot(srt_meso, group.by = "mOrg_classification", label = T)

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, srt_ref = srt_hOrg,
  ref_group = "annotation",
  filter_lowfreq = 10, prefix = "hOrg"
)
ClassDimPlot(srt_meso, group.by = "hOrg_classification", label = T)

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, srt_ref = srt_mEndo1,
  ref_group = "CellType",
  filter_lowfreq = 10, prefix = "mEndo1"
)
ClassDimPlot(srt_meso, group.by = "mEndo1_classification", label = T)

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, srt_ref = srt_mEndo2,
  ref_group = "celltype",
  filter_lowfreq = 10, prefix = "mEndo2"
)
ClassDimPlot(srt_meso, group.by = "mEndo2_classification", label = T)

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, srt_ref = srt_mEndo3,
  ref_group = "LineageAnnotations",
  filter_lowfreq = 10, prefix = "mEndo3"
)
ClassDimPlot(srt_meso, group.by = "mEndo3_classification", label = T)

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, srt_ref = srt_hEndo1,
  ref_group = "Cell_type",
  filter_lowfreq = 10, prefix = "hEndo1"
)
ClassDimPlot(srt_meso, group.by = "hEndo1_classification", label = T)

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, srt_ref = srt_hEndo2,
  ref_group = "cell_name_detailed",
  filter_lowfreq = 10, prefix = "hEndo2"
)
ClassDimPlot(srt_meso, group.by = "hEndo2_classification", label = T)

srt_meso <- RunKNNPredict(
  srt_query = srt_meso, srt_ref = srt_mMeso,
  ref_group = "celltype",
  filter_lowfreq = 10, prefix = "mMeso"
)
ClassDimPlot(srt_meso, group.by = "mMeso_classification", label = T)
saveRDS(srt_meso, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_meso_annotation.rds")


# Axial mesoderm: CHRD NOTO
# Paraial mesoderm: TBX6 SFRP2
# LPM: BMP4 PDGFRA GATA6
# Cardiac progenitor: ISL1 MAB21L2 HOPX
# Contractile: MYL7 TNNI1 TNNT2 MYH10

# Paraxial mesoderm: "FOXC1", "FOXC2", "TBX6", "CDH2"
# Presomitic mesoderm: "CYP26A1", "FGF8", "TBX6", "DKK1", "FGF17"
# Intermediate mesoderm: "PAX2", "PAX8", "OSR1", "SALL1", "HOXA9"
# Lateral plate mesoderm: "FOXF1", "PDGFRA", "BMP4"
# Cardiomyocytes: "HAND1", "MYL4", "MYL7", "MYH9", "MYH10", "TNNI1", "TNNT2", "TBX20"
# Splanchnic mesoderm: "FOXF1"
# pro-hematopoietic transcription factors: GATA2 LMO2 RUNX1 TAL1 ETV2 CD34

annotation_list <- list(
  "Epiblast" = c(4, 5),
  "Primitive streak" = c(2),
  "Nascent mesoderm" = c(3),
  "Paraxial mesoderm" = c(14, 15, 16),
  "Lateral plate mesoderm" = c(1, 12, 6, 7, 10, 25),
  "Cardiomyocytes" = c(13, 8, 9, 26),
  "Endothelial & erythroid cell" = c(17, 18),
  "MSC/Fib" = c(11, 22, 24, 27, 28, 29, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43),
  "Limb bud mesenchyme cell" = c(44, 45, 46, 47),
  "ExE mesoderm" = c(19, 20, 21, 23, 30, 31, 32, 33)
)
annotation <- setNames(rep(names(annotation_list), sapply(annotation_list, length)), nm = unlist(annotation_list))
srt_meso$SubAnnotation <- annotation[as.character(srt_meso$CSSclusters)]
ClassDimPlot(srt_meso, "SubAnnotation", label = T, cells.highlight = WhichCells(srt_meso, expression = SubAnnotation == "Primitive streak"))
cells <- colnames(srt_endo)[srt_endo$SubAnnotation %in% c("Anterior primitive streak", "Posterior primitive streak")]
srt_meso$SubAnnotation[cells] <- "Primitive streak"
cells <- colnames(srt_endo)[srt_endo$SubAnnotation %in% c("Nascent mesoderm")]
srt_meso$SubAnnotation[cells] <- "Nascent mesoderm"
ClassDimPlot(srt_meso, "SubAnnotation", label = T)

srt_meso$SubAnnotation <- factor(srt_meso$SubAnnotation, levels = names(annotation_list))

# ExpDimPlot(srt_meso, c(APS_marker, PPS_marker), "UMAP", nrow = 2)
# srt_meso_sub <- srt_meso[, srt_meso$SubAnnotation %in% c("Nascent mesoderm")]
# srt_meso_sub <- CellClassification(srt_meso_sub,
#   features = list(Anterior = APS_marker, Posterior = PPS_marker),
#   prefix = "AP"
# )
# ClassDimPlot(srt_meso_sub, "AP_classification")
# srt_meso_sub$AP_classification[is.na(srt_meso_sub$AP_classification)] <- "Posterior"
# srt_meso$AP_classification <- srt_meso_sub$AP_classification
# srt_meso$SubAnnotation <- as.character(srt_meso$SubAnnotation)
# srt_meso$SubAnnotation[srt_meso$SubAnnotation == "Nascent mesoderm"] <- srt_meso$AP_classification[srt_meso$SubAnnotation == "Nascent mesoderm"]
# srt_meso$SubAnnotation[srt_meso$SubAnnotation == "Anterior"] <- "Anterior nascent mesoderm"
# srt_meso$SubAnnotation[srt_meso$SubAnnotation == "Posterior"] <- "Posterior nascent mesoderm"
# srt_meso$SubAnnotation <- factor(srt_meso$SubAnnotation,
#   levels = c(
#     "Epiblast",
#     "Primitive streak",
#     "Anterior primitive streak",
#     "Posterior primitive streak",
#     "Anterior nascent mesoderm",
#     "Posterior nascent mesoderm",
#     "Nascent mesoderm",
#     "Paraxial mesoderm",
#     "Lateral plate mesoderm",
#     "Cardiomyocytes",
#     "Endothelial & erythroid cell",
#     "MSC/Fib",
#     "Limb bud mesenchyme cell",
#     "ExE mesoderm"
#   )
# )
# ClassDimPlot(srt_meso, "SubAnnotation", label = T)

srt_meso <- RunSlingshot(srt_meso,
  group.by = "SubAnnotation", reduction = "UMAP", prefix = "Meso",
  start = c("Epiblast"),
  end = c(
    "Paraxial mesoderm",
    "Cardiomyocytes",
    "Endothelial & erythroid cell",
    "Limb bud mesenchyme cell",
    "ExE mesoderm"
  )
)
ExpDimPlot(srt_meso, c("Meso_Lineage3"))

saveRDS(srt_meso, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_meso_annotation.rds")


### Baisc plot ------------------------------------------------------------
srt_22_tmp <- srt_22
srt_22_tmp$GermLayer[!colnames(srt_22) %in% colnames(srt_meso)] <- NA
srt_22_tmp$GermLayer <- factor(srt_22_tmp$GermLayer, levels = levels(srt_22$GermLayer))
p <- ClassDimPlot(srt_22_tmp,
  group.by = "GermLayer", theme_use = "theme_blank",
  cells.highlight = colnames(srt_meso), show_stat = FALSE
)
p <- panel_fix(p, height = 3, save = "figures/fig4/meso_highlight.umap.pdf", raster = TRUE)

p <- ClassDimPlot(srt_meso, "SubAnnotation", label = TRUE, show_stat = FALSE, theme_use = "theme_blank", force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig4/meso_annotation.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_meso, "Time", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig4/meso_time.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_meso, "Phase", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig4/meso_cellcycle.umap.pdf", raster = TRUE)

markers_meso <- list(
  "Epiblast" = c("POU5F1", "DNMT3B"),
  "Primitive streak" = c("TBXT", "EOMES"),
  "Nascent mesoderm" = c("MESP1", "MIXL1"),
  "Paraxial mesoderm" = c("DKK1", "FOXC1"),
  "Lateral plate mesoderm" = c("FOXF1", "PDGFRA"),
  "Cardiomyocytes" = c("HAND1", "MYL7"),
  "Endothelial & erythroid cell" = c("MEF2C", "HBE1"),
  "MSC/Fib" = c("THY1", "COL1A1"),
  "Limb bud mesenchyme cell" = c("TBX5", "PRRX1"),
  "ExE mesoderm" = c("POSTN", "FRZB")
)
p <- ExpDimPlot(srt_meso, unlist(markers_meso),
  theme_use = "theme_blank", ncol = 4
)
p <- panel_fix(p, height = 2, save = "figures/fig4/meso_markers.umap.pdf", raster = TRUE)

ht <- GroupHeatmap(srt_meso,
  group.by = "SubAnnotation", heatmap_palette = "YlOrRd",
  features = unlist(markers_meso),
  feature_split = rep(names(markers_meso), each = 2),
  show_row_names = TRUE, add_dot = TRUE, add_bg = TRUE,
  height = 6, width = 6, dot_size = unit(7, "mm")
)
panel_fix(ht$plot, save = "figures_add/mesoderm.markers.pdf")


p <- ExpDimPlot(srt_meso, c("PAX2", "PAX8", "OSR1", "SALL1"), theme_use = "theme_blank", calculate_coexp = TRUE)
p <- ExpDimPlot(srt_meso, c("PAX2", "OSR1", "SALL1"), theme_use = "theme_blank", calculate_coexp = TRUE)
p <- panel_fix(p, height = 2, save = "figures/fig4/meso_IM.marker.pdf", raster = TRUE)

p <- ExpDimPlot(srt_meso, c(APS_marker, PPS_marker), theme_use = "theme_blank")
p <- panel_fix(p, height = 2, save = "figures/fig4/meso_ap.umap.pdf", raster = TRUE)

srt_meso <- CellScoring(srt_meso,
  features = list("Anterior" = APS_marker, "Posterior" = PPS_marker),
  name = "AP"
)
p <- ExpDimPlot(srt_meso, c("AP_Anterior", "AP_Posterior"), theme_use = "theme_blank")
p <- panel_fix(p, height = 2, save = "figures/fig4/meso_apscore.umap.pdf", raster = TRUE)


### PAGA ---------------------------------------------------------------------
srt_meso_paga <- RunPAGA(
  srt = srt_meso, group_by = "SubAnnotation", linear_reduction = "CSS", nonlinear_reduction = "UMAP",
  n_pcs = ncol(srt_meso@reductions$CSS@cell.embeddings), n_neighbors = 100, embedded_with_PAGA = FALSE, return_seurat = TRUE
)
srt_meso@misc$paga <- srt_meso_paga@misc$paga
p <- ClassDimPlot(srt_meso,
  group.by = "SubAnnotation", pt.size = 5, pt.alpha = 0.05,
  label = TRUE, label.size = 3, label_repel = TRUE, label_insitu = TRUE, label_segment_color = "transparent",
  paga = srt_meso@misc$paga, paga_edge_threshold = 0.1, paga_edge_color = "black", paga_edge_alpha = 1,
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
)
p <- panel_fix(p, width = 4, save = "figures/fig4/meso_paga.pdf", raster = TRUE)

### SCVELO ------------------------------------------------------------------
srt_meso_scv <- SCP:::RunSCVELO(
  srt = srt_meso, group_by = "SubAnnotation", linear_reduction = "CSS", nonlinear_reduction = "UMAP",
  n_pcs = ncol(srt_meso@reductions$CSS@cell.embeddings), n_neighbors = 100, return_seurat = TRUE
)
srt_meso[["stochastic_UMAP"]] <- srt_meso_scv[["stochastic_UMAP"]]
srt_meso[["Ms"]] <- srt_meso_scv[["Ms"]]
srt_meso[["Mu"]] <- srt_meso_scv[["Mu"]]
srt_meso[["stochastic"]] <- srt_meso_scv[["stochastic"]]
srt_meso[["variance_stochastic"]] <- srt_meso_scv[["variance_stochastic"]]
p <- ClassDimPlot(srt_meso,
  group.by = "SubAnnotation", pt.size = 5, pt.alpha = 0.1,
  velocity = "stochastic", velocity_plot_type = "stream",
  velocity_density = 2, velocity_smooth = 1,
  streamline_n = 15, streamline_size = 0.5, streamline_color = "black",
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
)
p <- panel_fix(p, width = 4, save = "figures/fig4/meso_scvelo.pdf", raster = TRUE)


### Compare with human CS7 --------------------------------------------------------
srt_CS7 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_E19/human_gastrulation_E19.seurat.rds")
srt_CS7[["percent.mito"]] <- PercentageFeatureSet(object = srt_CS7, pattern = paste0(c("^MT-", "^Mt-", "^mt-"), collapse = "|"))
srt_CS7[["percent.ribo"]] <- PercentageFeatureSet(object = srt_CS7, pattern = paste0(c("^RP[SL]\\d+\\w{0,1}\\d*$", "^Rp[sl]\\d+\\w{0,1}\\d*$", "^rp[sl]\\d+\\w{0,1}\\d*$"), collapse = "|"))
levels <- c(
  "Epiblast", "Primitive Streak", "Axial Mesoderm", "Nascent Mesoderm", "Emergent Mesoderm", "Advanced Mesoderm", "YS Mesoderm",
  "Hemogenic Endothelium", "Blood Progenitors"
  # "Hemogenic Endothelium", "Erythroblasts", "Blood Progenitors", "Erythro-Myeloid Progenitors", "Myeloid Progenitors"
)
srt_CS7_sub <- srt_CS7[, srt_CS7$sub_cluster %in% levels]
srt_CS7_sub <- FindVariableFeatures(srt_CS7_sub)
srt_CS7_sub <- RunKNNMap(
  srt_query = srt_CS7_sub, srt_ref = srt_meso,
  ref_group = "CSSclusters", ref_umap = "UMAP"
)

srt_CS7_sub$sub_cluster <- factor(srt_CS7_sub$sub_cluster, levels = levels)
srt_CS7_sub$cell_type <- srt_CS7_sub$sub_cluster

srt_meso$bg <- NA
p <- ProjectionPlot(
  srt_query = srt_CS7_sub, srt_ref = srt_meso,
  query_group = "cell_type", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "CS7 human gastrula"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig4/meso_CS7_celltype.umap.pdf", raster = TRUE)

p <- ProjectionPlot(
  srt_query = srt_CS7_sub, srt_ref = srt_meso,
  query_group = "spatial", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "CS7 human gastrula"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig4/meso_CS7_spatial.umap.pdf", raster = TRUE)

### Slingshot ---------------------------------------------------------------
ClassDimPlot(srt_meso, "SubAnnotation", "UMAP", label = T)
srt_meso <- RunSlingshot(srt_meso,
  group.by = "SubAnnotation", reduction = "UMAP", prefix = "Meso",
  start = c("Epiblast"),
  end = c(
    "Paraxial mesoderm", "Cardiomyocytes", "Endothelial & erythroid cell",
    "Limb bud mesenchyme cell", "ExE mesoderm"
  )
)
srt_meso$Lineage_EpiPM <- srt_meso$Meso_Lineage5
srt_meso$Lineage_EpiHE <- srt_meso$Meso_Lineage4
srt_meso$Lineage_EpiCM <- srt_meso$Meso_Lineage3
srt_meso$Lineage_EpiLB <- srt_meso$Meso_Lineage1
srt_meso$Lineage_EpiExEM <- srt_meso$Meso_Lineage2
p <- ClassDimPlot(srt_meso, "SubAnnotation",
  lineages = c("Lineage_EpiPM", "Lineage_EpiHE", "Lineage_EpiCM", "Lineage_EpiLB", "Lineage_EpiExEM"),
  theme_use = "theme_blank"
)
p <- panel_fix(p, height = 3, save = "figures/fig4/meso_slingshot.umap.pdf")

srt_meso <- RunDynamicFeatures(srt_meso,
  lineages = c("Lineage_EpiPM", "Lineage_EpiHE", "Lineage_EpiCM", "Lineage_EpiLB", "Lineage_EpiExEM"),
  n_candidates = 10000
)
ht <- DynamicHeatmap(
  srt = srt_meso, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Meso_Lineage1"),
  cell_annotation = "SubAnnotation",
  n_split = 6, reverse_ht = 1, use_raster = TRUE, width = 10
)
ht$plot
saveRDS(srt_meso, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_meso_annotation.rds")

#### srt_meso_sub1 -----------------------------------------------------------
srt_meso_sub1 <- srt_meso[, srt_meso$SubAnnotation %in% c(
  "Epiblast",
  "Primitive streak",
  "Nascent mesoderm",
  "Paraxial mesoderm",
  "Lateral plate mesoderm",
  "Cardiomyocytes",
  "Endothelial & erythroid cell"
)]
srt_meso_sub1$Lineage_EpiLM <- srt_meso_sub1$Lineage_EpiExEM
ExpDimPlot(srt_meso_sub1, c("Lineage_EpiPM", "Lineage_EpiLM", "Lineage_EpiCM", "Lineage_EpiHE"), palette = "cividis")

# srt_meso_sub1 <- RunSlingshot(srt_meso_sub1,
#   group.by = "SubAnnotation", reduction = "UMAP", prefix = "Meso",
#   start = c("Epiblast"),
#   end = c(
#     "Paraxial mesoderm", "Lateral plate mesoderm", "Cardiomyocytes", "Endothelial & erythroid cell"
#   )
# )
# srt_meso_sub1$Lineage_EpiPM <- srt_meso$Meso_Lineage3
# srt_meso_sub1$Lineage_EpiHE <- srt_meso$Meso_Lineage4
# srt_meso_sub1$Lineage_EpiCM <- srt_meso$Meso_Lineage2
# srt_meso_sub1$Lineage_EpiLB <- srt_meso$Meso_Lineage1

srt_meso_sub1 <- RunDynamicFeatures(srt_meso_sub1,
  lineages = c("Lineage_EpiPM", "Lineage_EpiLM", "Lineage_EpiCM", "Lineage_EpiHE"),
  n_candidates = 10000
)
saveRDS(srt_meso_sub1, "srt_meso_sub1.rds")

colnames(srt_meso_sub1@tools$DynamicFeatures_Lineage_EpiPM$DynamicFeatures)[5] <- "peaktime"
colnames(srt_meso_sub1@tools$DynamicFeatures_Lineage_EpiLM$DynamicFeatures)[5] <- "peaktime"
ht <- DynamicHeatmap(
  srt = srt_meso_sub1, r.sq = 0.1, dev.expl = 0.1,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_EpiPM", "Lineage_EpiLM"),
  cell_annotation = "SubAnnotation",
  anno_terms = TRUE, anno_features = TRUE, db_version = "3.13",
  n_split = 4, reverse_ht = 1, use_raster = TRUE, width = 8
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig4/meso_lineagePMLM.ht.pdf")
saveRDS(ht, "figures/fig4/meso_lineagePMLM.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig4/meso_lineagePMLM.enrichment.xlsx")

p <- ExpDimPlot(srt_meso_sub1,
  features = c("Lineage_EpiPM", "Lineage_EpiLM"), ncol = 1,
  palette = "cividis", byrow = FALSE, theme_use = "theme_blank", force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig4/meso_lineagePMLM.pseudotime.umap.pdf", raster = TRUE)
srt_meso_sub1@tools$Enrichment_PMLM <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_meso_sub1@tools$Enrichment_PMLM$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 2)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig4/meso_lineagePMLM_bp.bar.pdf", raster = TRUE)

ht <- DynamicHeatmap(
  srt = srt_meso_sub1, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_EpiPM"),
  cell_annotation = "SubAnnotation",
  n_split = NULL, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig4/meso_lineageEpiPM.ht.pdf")
ht <- DynamicHeatmap(
  srt = srt_meso_sub1, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_EpiLM"),
  cell_annotation = "SubAnnotation",
  n_split = NULL, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig4/meso_lineageEpiLM.ht.pdf")

#############################
colnames(srt_meso_sub1@tools$DynamicFeatures_Lineage_EpiCM$DynamicFeatures)[5] <- "peaktime"
colnames(srt_meso_sub1@tools$DynamicFeatures_Lineage_EpiHE$DynamicFeatures)[5] <- "peaktime"
ht <- DynamicHeatmap(
  srt = srt_meso_sub1, r.sq = 0.1, dev.expl = 0.1,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_EpiCM", "Lineage_EpiHE"),
  cell_annotation = "SubAnnotation",
  anno_terms = TRUE, anno_features = TRUE, db_version = "3.13",
  n_split = 4, reverse_ht = 1, use_raster = TRUE, width = 8
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig4/meso_lineageCMHE.ht.pdf")
saveRDS(ht, "figures/fig4/meso_lineageCMHE.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig4/meso_lineageCMHE.enrichment.xlsx")

p <- ExpDimPlot(srt_meso_sub1,
  features = c("Lineage_EpiCM", "Lineage_EpiHE"), ncol = 1,
  palette = "cividis", byrow = FALSE, theme_use = "theme_blank", force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig4/meso_lineageCMHE.pseudotime.umap.pdf", raster = TRUE)
srt_meso_sub1@tools$Enrichment_CMHE <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_meso_sub1@tools$Enrichment_CMHE$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 2)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig4/meso_lineageCMHE_bp.bar.pdf", raster = TRUE)

ht <- DynamicHeatmap(
  srt = srt_meso_sub1, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_EpiCM"),
  cell_annotation = "SubAnnotation",
  n_split = NULL, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig4/meso_lineageEpiCM.ht.pdf")
ht <- DynamicHeatmap(
  srt = srt_meso_sub1, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_EpiHE"),
  cell_annotation = "SubAnnotation",
  n_split = NULL, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig4/meso_lineageEpiHE.ht.pdf")


#### srt_meso_sub2 -----------------------------------------------------------
srt_meso_sub2 <- srt_meso[, srt_meso$SubAnnotation %in% c(
  "MSC/Fib",
  "Limb bud mesenchyme cell",
  "ExE mesoderm"
)]
srt_meso_sub2 <- srt_meso_sub2[, !srt_meso_sub2$CSSclusters %in% c(1, 2, 3, 4)]
srt_meso_sub2$Lineage_LMLB <- srt_meso_sub2$Lineage_EpiLB
srt_meso_sub2$Lineage_LMExEM <- srt_meso_sub2$Lineage_EpiExEM
ExpDimPlot(srt_meso_sub2, c("Lineage_LMLB", "Lineage_LMExEM"), palette = "cividis")

# srt_meso_sub2 <- RunSlingshot(srt_meso_sub2,
#                          group.by = "SubAnnotation", reduction = "UMAP", prefix = "Meso",
#                          start = c("MSC/Fib"),
#                          end = c(
#                            "Limb bud mesenchyme cell", "ExE mesoderm"
#                          )
# )
# srt_meso_sub2$Lineage_EpiLB <- srt_meso$Meso_Lineage1
# srt_meso_sub2$Lineage_EpiExEM <- srt_meso$Meso_Lineage2

srt_meso_sub2 <- RunDynamicFeatures(srt_meso_sub2,
  lineages = c("Lineage_LMLB", "Lineage_LMExEM"),
  n_candidates = 10000
)
saveRDS(srt_meso_sub2, "srt_meso_sub2.rds")

colnames(srt_meso_sub2@tools$DynamicFeatures_Lineage_LMLB$DynamicFeatures)[5] <- "peaktime"
colnames(srt_meso_sub2@tools$DynamicFeatures_Lineage_LMExEM$DynamicFeatures)[5] <- "peaktime"
ht <- DynamicHeatmap(
  srt = srt_meso_sub2, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_LMLB", "Lineage_LMExEM"),
  cell_annotation = "SubAnnotation",
  n_split = 8, reverse_ht = 1, use_raster = TRUE, width = 10
)
cluster <- sapply(ht$feature_split, function(x) {
  switch(as.character(x),
    "C1" = "C1",
    "C2" = "C2",
    "C3" = "C2",
    "C4" = "C2",
    "C6" = "C3",
    "C5" = "C4",
    "C7" = "C5",
    "C8" = "C5",
  )
})
cluster <- factor(cluster, levels = paste0("C", 1:5))
cluster <- paste0(cluster, "(", table(cluster)[cluster], ")")
names(cluster) <- names(ht$feature_split)
cluster <- factor(cluster, sort(unique(cluster)))
ht <- DynamicHeatmap(
  srt = srt_meso_sub2, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE, feature_split = cluster,
  lineages = c("Lineage_LMLB", "Lineage_LMExEM"),
  cell_annotation = "SubAnnotation",
  anno_terms = TRUE, anno_features = TRUE, db_version = "3.13",
  n_split = 8, reverse_ht = 1, use_raster = TRUE, width = 8
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig4/meso_lineageLBExEM.ht.pdf")
saveRDS(ht, "figures/fig4/meso_lineageLBExEM.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig4/meso_lineageLBExEM.enrichment.xlsx")

p <- ExpDimPlot(srt_meso_sub2,
  features = c("Lineage_LMLB", "Lineage_LMExEM"), ncol = 1,
  palette = "cividis", byrow = FALSE, theme_use = "theme_blank", force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig4/meso_lineageLBExEM.pseudotime.umap.pdf", raster = TRUE)
srt_meso_sub2@tools$Enrichment_LBExEM <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_meso_sub2@tools$Enrichment_LBExEM$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 2)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig4/meso_lineageLBExEM_bp.bar.pdf", raster = TRUE)

ht <- DynamicHeatmap(
  srt = srt_meso_sub2, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_LMLB"), split_method = "kmeans-peaktime",
  cell_annotation = "SubAnnotation",
  anno_terms = TRUE, anno_features = TRUE, db_version = "3.13",
  n_split = 6, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig4/meso_lineageLMLB.ht.pdf")
srt_meso_sub2@tools$Enrichment_LMLB <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_meso_sub2@tools$Enrichment_LMLB$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig4/meso_lineageLMLB_bp.bar.pdf", raster = TRUE)

ht <- DynamicHeatmap(
  srt = srt_meso_sub2, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_LMExEM"), split_method = "kmeans-peaktime",
  cell_annotation = "SubAnnotation",
  anno_terms = TRUE, anno_features = TRUE, db_version = "3.13",
  n_split = 6, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig4/meso_lineageLMExEM.ht.pdf")
srt_meso_sub2@tools$Enrichment_LMExEM <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_meso_sub2@tools$Enrichment_LMExEM$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig4/meso_lineageLMExEM_bp.bar.pdf", raster = TRUE)

saveRDS(srt_meso_sub2, "srt_meso_sub2.rds")

#### all lieages ############
srt_meso[["Lineage_EpiPM"]] <- srt_meso_sub1[["Lineage_EpiPM"]]
srt_meso[["Lineage_EpiLM"]] <- srt_meso_sub1[["Lineage_EpiLM"]]
srt_meso[["Lineage_EpiCM"]] <- srt_meso_sub1[["Lineage_EpiCM"]]
srt_meso[["Lineage_EpiHE"]] <- srt_meso_sub1[["Lineage_EpiHE"]]
srt_meso[["Lineage_LMLB"]] <- srt_meso_sub2[["Lineage_LMLB"]]
srt_meso[["Lineage_LMExEM"]] <- srt_meso_sub2[["Lineage_LMExEM"]]

srt_meso2 <- srt_meso
cell_highlight <- colnames(srt_meso_sub1)[!is.na(srt_meso_sub1[["Lineage_EpiPM", drop = TRUE]]) | !is.na(srt_meso_sub1[["Lineage_EpiLM", drop = TRUE]])]
srt_meso2$SubAnnotation[setdiff(colnames(srt_meso2), cell_highlight)] <- NA
p <- ClassDimPlot(srt_meso2, "SubAnnotation",
  cells.highlight = cell_highlight,
  lineages = c("Lineage_EpiPM", "Lineage_EpiLM"),
  show_stat = FALSE, theme_use = "theme_blank"
)
p <- panel_fix(p, height = 3, save = "figures/fig4/meso_lineagePMLM.umap.pdf", raster = TRUE)

srt_meso2 <- srt_meso
cell_highlight <- colnames(srt_meso_sub1)[!is.na(srt_meso_sub1[["Lineage_EpiCM", drop = TRUE]]) | !is.na(srt_meso_sub1[["Lineage_EpiHE", drop = TRUE]])]
srt_meso2$SubAnnotation[setdiff(colnames(srt_meso2), cell_highlight)] <- NA
p <- ClassDimPlot(srt_meso2, "SubAnnotation",
  cells.highlight = cell_highlight,
  lineages = c("Lineage_EpiCM", "Lineage_EpiHE"),
  show_stat = FALSE, theme_use = "theme_blank"
)
p <- panel_fix(p, height = 3, save = "figures/fig4/meso_lineageCMHE.umap.pdf", raster = TRUE)

srt_meso2 <- srt_meso
cell_highlight <- colnames(srt_meso_sub2)[!is.na(srt_meso_sub2[["Lineage_LMLB", drop = TRUE]]) | !is.na(srt_meso_sub2[["Lineage_LMExEM", drop = TRUE]])]
srt_meso2$SubAnnotation[setdiff(colnames(srt_meso2), cell_highlight)] <- NA
p <- ClassDimPlot(srt_meso2, "SubAnnotation",
  cells.highlight = cell_highlight, lineages_trim = c(0.05, 0.95),
  lineages = c("Lineage_LMLB", "Lineage_LMExEM"),
  show_stat = FALSE, theme_use = "theme_blank"
)
p <- panel_fix(p, height = 3, save = "figures/fig4/meso_lineageLBExEM.umap.pdf", raster = TRUE)



### DE analysis ------------------------------------------------------------------
srt_meso <- RunDEtest(srt_meso,
  group_by = "CSSclusters", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_meso <- RunDEtest(srt_meso,
  group_by = "CSSclusters", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
de_filter <- filter(srt_meso@tools$DEtest_CSSclusters$AllMarkers_wilcox, p_val_adj < 0.05)
res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "tmp_meso_bp.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "tmp_meso_bp.word.pdf")

srt_meso <- RunDEtest(srt_meso,
  group_by = "SubAnnotation", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_meso <- RunDEtest(srt_meso,
  group_by = "SubAnnotation", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
de_filter <- filter(srt_meso@tools$DEtest_SubAnnotation$AllMarkers_wilcox, p_val_adj < 0.05)
res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "tmp_meso_bp2.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "tmp_meso_bp2.word.pdf")

# Fig5 --------------------------------------------------------------------
dir.create(paste0(work_dir, "/figures/fig5"), recursive = T, showWarnings = FALSE)

## Amnion and PGC --------------------------------------------------------------------
epi <- colnames(srt_epi)[srt_epi$SubAnnotation %in% paste0("Epiblast-", 3:5)]
cells <- c(epi, colnames(srt_22)[srt_22$CSSclusters %in% c(27, 28, 29, 41, 42, 43, 44)])
srt_amn <- srt_22[, cells]

### Integration  --------------------------------------------------------------------
srt_amn <- CellCycleScoring(srt_amn,
  s.features = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes
)

plist <- list()
pagalist <- list()
progress <- 1
total <- length(c(2000, 3000)) * length(seq(15, 50, 5)) * length(seq(0.5, 5, 0.5))
for (nhvf in c(2000, 3000)) {
  for (i in seq(15, 50, 5)) {
    for (j in seq(1, 2, 0.5)) {
      cat("progress:", progress, "/", total, "\n")
      progress <- progress + 1
      srt_amn <- CSS_integrate(
        srtMerge = srt_amn, nHVF = nhvf,
        linear_reduction_dims_use = 1:i, cluster_resolution = j,
        CSS_param = list(cluster_resolution = j)
      )
      srt_amn[[paste0("CSSUMAP2Dnhvf", nhvf, "dim", i, "res", j)]] <- srt_amn@reductions$CSSUMAP2D
      srt_amn[[paste0("CSSUMAP3Dnhvf", nhvf, "dim", i, "res", j)]] <- srt_amn@reductions$CSSUMAP3D
      plist[[paste0("CSSUMAP2Dnhvf", nhvf, "dim", i, "res", j)]] <- ClassDimPlot(srt_amn, "Annotation", title = paste("hvf:", nhvf, "dim", i, "res", j))

      adata <- srt_to_adata(srt_amn)
      adata <- RunPAGA(
        adata = adata, group_by = "CSSclusters", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
        n_pcs = ncol(srt_amn@reductions$CSS@cell.embeddings), threshold = 0.1,
        embedded_with_PAGA = TRUE, paga_layout = "fa"
      )
      PAGAUMAP2D <- adata$obsm[["PAGAUMAP2D"]]
      colnames(PAGAUMAP2D) <- c("UMAP_1", "UMAP_2")
      rownames(PAGAUMAP2D) <- adata$obs_names$values
      srt_amn[["PAGAUMAP2D"]] <- CreateDimReducObject(PAGAUMAP2D, key = "UMAP_", assay = "RNA")
      ClassDimPlot(srt_amn, group.by = "Time", reduction = "PAGAUMAP2D")
      pagalist[[paste0("PAGAUMAP2Dnhvf", nhvf, "dim", i, "res", j)]] <- ClassDimPlot(srt_amn, "Annotation", reduction = "PAGAUMAP2D", title = paste("hvf:", nhvf, "dim", i, "res", j))
    }
  }
}
p <- plot_grid(plotlist = plist)
p <- panel_fix(p, save = "tmp_amn1.png", width = 3)
p <- plot_grid(plotlist = pagalist)
p <- panel_fix(p, save = "tmp_amn2.png", width = 3)

# 20:3 25:3 30:3
ClassDimPlot3D(srt_amn, group.by = "SubAnnotation", "CSSUMAP3Ddim20res3")

srt_amn <- Integration_SCP(
  srtMerge = srt_amn, integration_method = "Uncorrected",
  linear_reduction_dims_use = 1:50, cluster_resolution = 2
)
srt_amn <- Integration_SCP(
  srtMerge = srt_amn, integration_method = "CSS",
  linear_reduction_dims_use = 1:50, cluster_resolution = 2,
  CSS_param = list(cluster_resolution = 2)
)
ClassDimPlot(srt_amn, "CSSclusters", label = TRUE) # 50-2
srt_amn$SubAnnotation <- srt_epi$SubAnnotation
ExpDimPlot(srt_amn, c("POU5F1", "TBXT", "EOMES", "TFAP2A", "TFAP2C", "NANOS3", "ISL1", "ABCG2", "CDH3"))

srt_amn <- RunDM(srt_amn,
  features = VariableFeatures(srt_amn), slot = "data", dims = NULL,
  k = 15, ndcs = 3, reduction.name = "DM"
)
ClassDimPlot(srt_amn, c("Time", "Annotation", "CSSclusters", "Phase"),
  label = T,
  cells.highlight = WhichCells(srt_amn, expression = Annotation == "PGC")
)
ClassDimPlot(srt_amn, c("CSSclusters"), label = T)
ExpDimPlot(srt_amn, c("percent.ribo", "percent.mito", "nFeature_RNA", "nCount_RNA"))

srt_amn[["UMAP"]] <- srt_amn[["CSSUMAP2D"]]
srt_amn[["UMAP3D"]] <- srt_amn[["CSSUMAP3D"]]
srt_amn@misc$Default_reduction <- "UMAP"
colnames(srt_amn[["UMAP"]]@cell.embeddings) <- c("UMAP_1", "UMAP_2")
colnames(srt_amn[["UMAP3D"]]@cell.embeddings) <- c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3")
Key(srt_amn[["UMAP"]]) <- "UMAP_"
Key(srt_amn[["UMAP3D"]]) <- "UMAP3D_"
srt_amn@misc$Default_reduction <- "UMAP"


srt_amn_sub <- Standard_SCP(srt_amn[, srt_amn$CSSclusters %in% c(17)],
  cluster_resolution = 2, linear_reduction_dims = 50
)
srt_amn$subclusters1 <- srt_amn_sub$Standardclusters
ClassDimPlot(srt_amn, "subclusters1",
  pt.size = 1,
  cells.highlight = WhichCells(srt_amn, expression = subclusters1 %in% c(6, 7))
)

srt_amn_sub <- Standard_SCP(srt_amn[, srt_amn$CSSclusters %in% c(18)],
  cluster_resolution = 2, linear_reduction_dims = 50
)
srt_amn$subclusters2 <- srt_amn_sub$Standardclusters
ClassDimPlot(srt_amn, "subclusters2",
  pt.size = 1,
  cells.highlight = WhichCells(srt_amn, expression = subclusters2 == 3)
)

# srt_amn <- FindClusters(object = srt_amn, resolution = 2, algorithm = 1, graph.name = "CSS_SNN")
# srt_amn <- SrtReorder(srt_amn, features = VariableFeatures(srt_amn), reorder_by = "seurat_clusters", slot = "data")
# srt_amn[["seurat_clusters"]] <- NULL
# srt_amn[["CSSclusters"]] <- Idents(srt_amn)
# ClassDimPlot(srt_amn, c("CSSclusters"), label = TRUE)




### paga initiated layout ---------------------------------------------------------------
adata <- srt_to_adata(srt_amn)
adata <- RunPAGA(
  adata = adata, group_by = "CSSclusters", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
  n_pcs = ncol(srt_amn@reductions$CSS@cell.embeddings), threshold = 0.1,
  embedded_with_PAGA = TRUE, paga_layout = "fa"
)
FA <- adata$obsm[["X_draw_graph_fa"]]
colnames(FA) <- c("FA_1", "FA_2")
rownames(FA) <- adata$obs_names$values
srt_amn[["FA"]] <- CreateDimReducObject(FA, key = "FA_", assay = "RNA")
ClassDimPlot(srt_amn, group.by = "Annotation", reduction = "FA")
srt_amn@misc$Default_reduction <- "FA"

PAGAUMAP2D <- adata$obsm[["PAGAUMAP2D"]]
colnames(PAGAUMAP2D) <- c("UMAP_1", "UMAP_2")
rownames(PAGAUMAP2D) <- adata$obs_names$values
srt_amn[["PAGAUMAP2D"]] <- CreateDimReducObject(PAGAUMAP2D, key = "UMAP_", assay = "RNA")
ClassDimPlot(srt_amn, group.by = "Annotation", reduction = "PAGAUMAP2D")
srt_amn@misc$Default_reduction <- "PAGAUMAP2D"

### Annotation --------------------------------------------------------------
ClassDimPlot(srt_amn, "Annotation")
ClassDimPlot(srt_amn, "CSSclusters", label = TRUE)
ExpDimPlot(srt_amn, c("POU5F1", "TBXT", "EOMES", "NANOS3", "TFAP2C", "TFAP2A", "ABCG2", "ISL1", "TP63"))

srt_amn <- RunKNNPredict(
  srt_query = srt_amn, query_group = "CSSclusters",
  srt_ref = srt_amn, ref_group = "CSSclusters",
  query_collapsing = T, ref_collapsing = T,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
d <- 1 - srt_amn@tools$knnpredict_CSSclusters$distance_matrix
d <- as.matrix(d)
Heatmap(d)

# KRT4 KRT13 KRT17
# dermis: ACTA2 COL1A1
# epidermis: KRT5 KRT10


annotation_list <- list(
  "Epiblast" = c(12, 13),
  "Preplacodal ectoderm" = c(15),
  "Primitive streak" = c(14),
  "PGC" = c(16),
  "Amnion-1" = c(17),
  "Amnion-2" = c(18),
  "Amnion-3" = c(10, 11),
  "Amnion-4" = c(6, 7, 8, 9),
  "Surface ectoderm-1" = c(19),
  "Surface ectoderm-2" = c(2, 3, 4, 5),
  "Epidermis" = c(1)
)
annotation <- setNames(rep(names(annotation_list), sapply(annotation_list, length)), nm = unlist(annotation_list))
srt_amn$SubAnnotation <- annotation[as.character(srt_amn$CSSclusters)]
srt_amn$SubAnnotation <- factor(srt_amn$SubAnnotation, levels = names(annotation_list))
ClassDimPlot(srt_amn, "SubAnnotation", label = T)

srt_amn$SubAnnotation[WhichCells(srt_amn, expression = subclusters1 %in% c(6, 7))] <- "PGC"
srt_amn$SubAnnotation[WhichCells(srt_amn, expression = subclusters2 %in% c(3))] <- "Surface ectoderm-1"
ClassDimPlot(srt_amn, "SubAnnotation", label = T)

srt_amn$Epiblast <- NA
srt_amn$Epiblast[epi] <- as.character(srt_epi$SubAnnotation[epi])
ClassDimPlot(srt_amn, "Epiblast")

### Baisc plot ------------------------------------------------------------
srt_22_tmp <- srt_22
srt_22_tmp$GermLayer[!colnames(srt_22) %in% colnames(srt_amn)] <- NA
srt_22_tmp$GermLayer <- factor(srt_22_tmp$GermLayer, levels = levels(srt_22$GermLayer))
p <- ClassDimPlot(srt_22_tmp,
  group.by = "GermLayer", theme_use = "theme_blank",
  cells.highlight = colnames(srt_amn), show_stat = FALSE
)
p <- panel_fix(p, height = 3, save = "figures/fig5/amn_rawpos.umap.pdf", raster = TRUE)

p <- ClassDimPlot(srt_amn, "SubAnnotation", label = TRUE, theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig5/amn_annotation.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_amn, "Time", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig5/amn_time.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_amn, "Phase", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig5/amn_cellcycle.umap.pdf", raster = TRUE)
p <- ExpDimPlot(srt_amn, "nFeature_RNA", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 2, save = "figures/fig5/amn_nfeature.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_amn, "Epiblast", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
p <- panel_fix(p, height = 3, save = "figures/fig5/amn_epiblast.umap.pdf", raster = TRUE)

srt_amn$Epiblast[srt_amn$Epiblast == "Epiblast-3"] <- "BDA-Epi"
srt_amn$Epiblast[srt_amn$Epiblast == "Epiblast-4"] <- "PSA-Epi"
srt_amn$Epiblast[srt_amn$Epiblast == "Epiblast-5"] <- "nNEA-Epi"

markers_amn <- list(
  "Epiblast" = c("POU5F1", "DNMT3B"),
  # "Preplacodal ectoderm" = c("SIX4", "DLX3"),
  "Primitive streak" = c("TBXT", "EOMES"),
  "PGC" = c("NANOS3", "TFAP2C"),
  "Amnion" = c("TFAP2A", "ISL1"), # "IGFBP3", "KCNMA1"
  "Surface ectoderm" = c("FOXG1", "TP63"), # "PDGFC", "TFAP2B"
  "Epidermis" = c("KRT4", "KRT13")
)
p <- ExpDimPlot(srt_amn, unlist(markers_amn),
  theme_use = "theme_blank", ncol = 6
)
p <- panel_fix(p, height = 2, save = "figures/fig5/amn_markers.umap.pdf", raster = TRUE)

nameslist <- list(
  "Amnion & Epidermis-precursor" = c("Amnion-1", "Amnion-2"),
  "Amnion-1" = "Amnion-3",
  "Amnion-2" = "Amnion-4"
)
srt_amn <- RenameClusters(srt_amn, group.by = "SubAnnotation", nameslist = nameslist, name = "SubAnnotation2", keep_levels = T)
ClassDimPlot(srt_amn, "SubAnnotation2")
srt_amn$SubAnnotation3 <- as.character(srt_amn$SubAnnotation2)
epis <- colnames(srt_amn)[srt_amn$Epiblast %in% c("BDA-Epi", "PSA-Epi", "nNEA-Epi")]
srt_amn$SubAnnotation3[epis] <- srt_amn$Epiblast[epis]
srt_amn$SubAnnotation3 <- factor(srt_amn$SubAnnotation3, c("BDA-Epi", "PSA-Epi", "nNEA-Epi", setdiff(levels(srt_amn$SubAnnotation2), c("Epiblast", "Preplacodal ectoderm")), "Unassigned"))
srt_amn$SubAnnotation3[is.na(srt_amn$SubAnnotation3)] <- "Unassigned"
ClassDimPlot(srt_amn, "SubAnnotation3")
palcolor <- palette_scp(srt_amn$SubAnnotation3)
palcolor["Unassigned"] <- "grey"
ClassDimPlot(srt_amn, "SubAnnotation3",
  palcolor = palcolor,
  theme_use = "theme_blank"
) %>% panel_fix(height = 3, save = "figures_add/amn.new_anno.pdf")

ClassStatPlot(srt_amn, c("Epiblast"), group.by = "SubAnnotation2") %>%
  panel_fix(height = 2, width = 3, save = "figures_add/amn_epi.stat.pdf")
ClassDimPlot(srt_amn, "Epiblast", split.by = "Epiblast", cells.highlight = T, nrow = 1, theme_use = "theme_blank") %>%
  panel_fix(height = 2, save = "figures_add/amn_epi.split.pdf")

srt_amn$SubAnnotation4 <- srt_amn$SubAnnotation3
srt_amn$SubAnnotation4[srt_amn$SubAnnotation4 == "Unassigned"] <- NA
ht <- GroupHeatmap(srt_amn,
  group.by = "SubAnnotation4", heatmap_palette = "YlOrRd",
  features = unlist(markers_amn),
  feature_split = rep(names(markers_amn), each = 2),
  show_row_names = TRUE, add_dot = TRUE, add_bg = TRUE,
  height = 4.4, width = 6, dot_size = unit(7, "mm")
)
panel_fix(ht$plot, save = "figures_add/non_ectoderm.markers.pdf")

p <- ExpDimPlot(srt_amn, c("CDH1", "CDH2", "CDH3", "CDH4"),
  theme_use = "theme_blank", ncol = 4
)
p <- panel_fix(p, height = 2, save = "figures/fig5/amn_CDH.umap.pdf", raster = TRUE)


p <- ExpDimPlot(srt_amn, c("COL1A1", "COL3A1", "COL5A1", "FN1", "LAMB1", "HSPG2"),
  theme_use = "theme_blank", ncol = 3
)
p <- panel_fix(p, height = 2, save = "figures/fig5/amn_basement_markers.umap.pdf", raster = TRUE)

srt_amn$Epiblast_type <- srt_epi$SubAnnotation[colnames(srt_amn)]
srt_amn_sub <- srt_amn[, intersect(epi, WhichCells(srt_amn, expression = SubAnnotation %in% c("Epiblast", "Preplacodal ectoderm", "Primitive streak")))]
srt_amn_sub$Epiblast_type <- factor(srt_amn_sub$Epiblast_type, levels = unique(srt_amn_sub$Epiblast_type))
srt_amn_sub$SubAnnotation <- factor(srt_amn_sub$SubAnnotation, levels = unique(srt_amn_sub$SubAnnotation))
p <- ClassStatPlot(srt_amn_sub, stat.by = "Epiblast_type", group.by = "SubAnnotation", aspect.ratio = 0.5)
p <- panel_fix(p, height = 2, save = "figures/fig5/amn_epiblast.bar.pdf", raster = TRUE)

### PAGA ---------------------------------------------------------------------
srt_amn_paga <- RunPAGA(
  srt = srt_amn, group_by = "SubAnnotation2", linear_reduction = "CSS", nonlinear_reduction = "UMAP",
  n_pcs = ncol(srt_amn@reductions$CSS@cell.embeddings), embedded_with_PAGA = FALSE, return_seurat = TRUE
)
srt_amn@misc$paga <- srt_amn_paga@misc$paga
ClassDimPlot(srt_amn,
  group.by = "SubAnnotation2", pt.size = 5, pt.alpha = 0.05,
  label = TRUE, label.size = 3, label_repel = TRUE, label_insitu = TRUE, label_segment_color = "transparent",
  paga = srt_amn@misc$paga, paga_edge_threshold = 0.1, paga_edge_color = "black", paga_edge_alpha = 1,
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
) %>%
  panel_fix(width = 4, save = "figures_add/amn_paga.pdf", raster = TRUE)

srt_amn_paga <- srt_amn[, !is.na(srt_amn$SubAnnotation4)]
srt_amn_paga$SubAnnotation4 <- factor(srt_amn_paga$SubAnnotation4, levels = setdiff(levels(srt_amn$SubAnnotation4), "Unassigned"))
srt_amn_paga <- RunPAGA(
  srt = srt_amn_paga, group_by = "SubAnnotation4", linear_reduction = "CSS", nonlinear_reduction = "UMAP",
  n_pcs = ncol(srt_amn_paga@reductions$CSS@cell.embeddings), embedded_with_PAGA = FALSE, return_seurat = TRUE
)
ClassDimPlot(srt_amn_paga,
  group.by = "SubAnnotation4", pt.size = 5, pt.alpha = 0.05,
  label = TRUE, label.size = 3, label_repel = TRUE, label_insitu = TRUE, label_segment_color = "transparent",
  paga = srt_amn_paga@misc$paga, paga_edge_threshold = 0.1, paga_edge_color = "black", paga_edge_alpha = 1,
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
) %>%
  panel_fix(width = 4, save = "figures_add/amn_paga2.pdf", raster = TRUE)

### SCVELO ------------------------------------------------------------------
srt_amn_scv <- SCP:::RunSCVELO(
  srt = srt_amn, group_by = "SubAnnotation", linear_reduction = "CSS", nonlinear_reduction = "UMAP",
  n_pcs = ncol(srt_amn@reductions$CSS@cell.embeddings), n_neighbors = 100, return_seurat = TRUE
)
srt_amn[["stochastic_UMAP"]] <- srt_amn_scv[["stochastic_UMAP"]]
srt_amn[["Ms"]] <- srt_amn_scv[["Ms"]]
srt_amn[["Mu"]] <- srt_amn_scv[["Mu"]]
srt_amn[["stochastic"]] <- srt_amn_scv[["stochastic"]]
srt_amn[["variance_stochastic"]] <- srt_amn_scv[["variance_stochastic"]]
p <- ClassDimPlot(srt_amn,
  group.by = "SubAnnotation3", palcolor = palcolor, pt.size = 5, pt.alpha = 0.1,
  velocity = "stochastic", velocity_plot_type = "stream",
  velocity_density = 2, velocity_smooth = 1,
  streamline_n = 15, streamline_size = 0.5, streamline_color = "black",
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
)
p <- panel_fix(p, width = 4, save = "figures_add/amn_scvelo.pdf", raster = TRUE)

### Compare with human CS7 --------------------------------------------------------
srt_CS7 <- readRDS("/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/reference_data/human_gastrulation_E19/human_gastrulation_E19.seurat.rds")
srt_CS7_sub <- srt_CS7[, srt_CS7$sub_cluster %in% c("Epiblast", "Primitive Streak", "Non-Neural Ectoderm", "PGC")]
srt_CS7_sub <- FindVariableFeatures(srt_CS7_sub)
srt_CS7_sub <- RunKNNMap(
  srt_query = srt_CS7_sub, srt_ref = srt_amn,
  ref_group = "SubAnnotation", ref_umap = "CSSUMAP2D"
)
levels <- c("Epiblast", "Primitive Streak", "Non-Neural Ectoderm", "PGC")
srt_CS7_sub$sub_cluster <- factor(srt_CS7_sub$sub_cluster, levels = levels)
srt_CS7_sub$cell_type <- srt_CS7_sub$sub_cluster

srt_amn$bg <- NA
p <- ProjectionPlot(
  srt_query = srt_CS7_sub, srt_ref = srt_amn,
  query_group = "cell_type", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "CS7 human gastrula"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig5/amn_CS7_celltype.umap.pdf", raster = TRUE)

p <- ProjectionPlot(
  srt_query = srt_CS7_sub, srt_ref = srt_amn,
  query_group = "spatial", ref_group = "bg",
  query_param = list(palette = "Set1", show_stat = FALSE),
  ref_param = list(
    bg_color = "grey80", legend.position = "none", theme_use = "theme_blank", xlab = "UMAP_1", ylab = "UMAP_2",
    title = "CS7 human gastrula"
  )
)
p <- panel_fix(p, height = 3, save = "figures/fig5/amn_CS7_spatial.umap.pdf", raster = TRUE)

### Slingshot ---------------------------------------------------------------
ClassDimPlot(srt_amn, "CSSclusters", "UMAP", label = T)
srt_amn <- RunSlingshot(srt_amn,
  group.by = "CSSclusters", reduction = "UMAP3D",
  start = "12", end = c("14", "16", "1"), prefix = "Amn"
)
ExpDimPlot(srt_amn, paste0("Amn_Lineage", 1:7))

ClassDimPlot(srt_amn, "SubAnnotation", "UMAP", label = T)
ClassDimPlot(srt_amn, "SubAnnotation2", "UMAP", label = T)
srt_amn <- RunSlingshot(srt_amn,
  group.by = "SubAnnotation", reduction = "UMAP", prefix = "Amn",
  start = "Epiblast", end = c("Preplacodal ectoderm", "Primitive streak", "PGC", "Epidermis", "Amnion-4")
)
ExpDimPlot(srt_amn, paste0("Amn_Lineage", 1:5), palette = "cividis")

srt_amn$Amn_Lineage1_sub1 <- srt_amn$Amn_Lineage1
srt_amn$Amn_Lineage1_sub1[which(srt_amn$Amn_Lineage1_sub1 > max(srt_amn$Amn_Lineage3, na.rm = T) + 1)] <- NA
ExpDimPlot(srt_amn, "Amn_Lineage1_sub1")
srt_amn$Amn_Lineage1_sub2 <- srt_amn$Amn_Lineage1
srt_amn$Amn_Lineage1_sub2[which(srt_amn$Amn_Lineage1_sub2 < max(srt_amn$Amn_Lineage3, na.rm = T) - 2)] <- NA
ExpDimPlot(srt_amn, "Amn_Lineage1_sub2")
srt_amn$Amn_Lineage2_sub2 <- srt_amn$Amn_Lineage2
srt_amn$Amn_Lineage2_sub2[which(srt_amn$Amn_Lineage2_sub2 < max(srt_amn$Amn_Lineage3, na.rm = T) - 2)] <- NA
ExpDimPlot(srt_amn, "Amn_Lineage2_sub2")

srt_amn <- RunDynamicFeatures(srt_amn,
  lineages = c(paste0("Amn_Lineage", 1:5), "Amn_Lineage1_sub1", "Amn_Lineage1_sub2", "Amn_Lineage2_sub2"),
  n_candidates = 10000
)
saveRDS(srt_amn, "srt_amn_annotation.rds")

srt_amn$Lineage_EpiAmn <- srt_amn$Amn_Lineage1_sub1
srt_amn$Lineage_EpiPGC <- srt_amn$Amn_Lineage3
srt_amn@tools$DynamicFeatures_Lineage_EpiAmn <- srt_amn@tools$DynamicFeatures_Amn_Lineage1_sub1
srt_amn@tools$DynamicFeatures_Lineage_EpiPGC <- srt_amn@tools$DynamicFeatures_Amn_Lineage3
colnames(srt_amn@tools$DynamicFeatures_Lineage_EpiAmn$DynamicFeatures)[5] <- "peaktime"
colnames(srt_amn@tools$DynamicFeatures_Lineage_EpiPGC$DynamicFeatures)[5] <- "peaktime"
ht <- DynamicHeatmap(
  srt = srt_amn, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 1, use_fitted = TRUE,
  lineages = c("Lineage_EpiAmn", "Lineage_EpiPGC"),
  cell_annotation = "SubAnnotation",
  n_split = 6, reverse_ht = 1, use_raster = TRUE, width = 8
)
cluster <- sapply(ht$feature_split, function(x) {
  switch(as.character(x),
    "C1" = "C1",
    "C2" = "C2",
    "C3" = "C3",
    "C4" = "C4",
    "C5" = "C6",
    "C6" = "C5"
  )
})
cluster <- factor(cluster, levels = levels(ht$feature_split))
# cluster <- paste0(cluster, "(", table(cluster)[cluster], ")")
# names(cluster) <- names(ht$feature_split)
ht <- DynamicHeatmap(
  srt = srt_amn, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 1, use_fitted = TRUE, feature_split = cluster,
  lineages = c("Lineage_EpiAmn", "Lineage_EpiPGC"),
  cell_annotation = c("SubAnnotation", "Epiblast"),
  anno_terms = TRUE, anno_features = TRUE, db_version = "3.13",
  n_split = 6, reverse_ht = 1,
  height = 7, width = 5, use_raster = TRUE
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig5/amn_lineage1.ht.pdf")
ggsave(p, height = 12, width = 20, filename = "figures_add/amn_lineage1.ht.pdf")
saveRDS(ht, "figures/fig5/amn_lineage1.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig5/amn_lineage1.enrichment.xlsx")
write.xlsx(ht$feature_metadata, file = "figures/fig5/amn_lineage1.feature_metadata.xlsx")


ht <- DynamicHeatmap(
  srt = srt_amn, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 1, use_fitted = TRUE,
  lineages = c("Lineage_EpiPGC"),
  cell_annotation = "SubAnnotation",
  anno_terms = TRUE, anno_features = TRUE, db_version = "3.13",
  n_split = 6, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig5/amn_pgc.ht.pdf")

p <- ExpDimPlot(srt_amn,
  features = c("Lineage_EpiAmn", "Lineage_EpiPGC"), ncol = 1,
  palette = "cividis", byrow = FALSE, theme_use = "theme_blank", force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig5/amn_lineage1.umap.pdf", raster = TRUE)
srt_amn@tools$Enrichment_lineage1 <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_amn@tools$Enrichment_lineage1$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig5/amn_lineage1_bp.bar.pdf", raster = TRUE)


srt_amn$Lineage_Sfc <- srt_amn$Amn_Lineage1_sub2
srt_amn$Lineage_Amn <- srt_amn$Amn_Lineage2_sub2
srt_amn@tools$DynamicFeatures_Lineage_Sfc <- srt_amn@tools$DynamicFeatures_Amn_Lineage1_sub2
srt_amn@tools$DynamicFeatures_Lineage_Amn <- srt_amn@tools$DynamicFeatures_Amn_Lineage2_sub2

ht <- DynamicHeatmap(
  srt = srt_amn, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.05, use_fitted = TRUE,
  lineages = c("Lineage_Amn"), split_method = "kmeans-peaktime",
  cell_annotation = "SubAnnotation",
  n_split = 6, use_raster = TRUE, width = 8
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig5/amn_lineage2.ht.pdf")

p <- ExpDimPlot(srt_amn,
  features = c("Lineage_Amn"), ncol = 1,
  palette = "cividis", byrow = FALSE, theme_use = "theme_blank", force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig5/amn_lineage2.umap.pdf", raster = TRUE)
srt_amn@tools$Enrichment_lineage2 <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_amn@tools$Enrichment_lineage2$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig5/amn_lineage2_bp.bar.pdf", raster = TRUE)

ht <- DynamicHeatmap(
  srt = srt_amn, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.05, use_fitted = TRUE,
  lineages = c("Lineage_Sfc"), split_method = "kmeans-peaktime",
  cell_annotation = "SubAnnotation",
  n_split = 6, use_raster = TRUE, width = 8
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig5/amn_lineage3.ht.pdf")
p <- ExpDimPlot(srt_amn,
  features = c("Lineage_Sfc"), ncol = 1,
  palette = "cividis", byrow = FALSE, theme_use = "theme_blank", force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig5/amn_lineage3.umap.pdf", raster = TRUE)
srt_amn@tools$Enrichment_lineage3 <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_amn@tools$Enrichment_lineage3$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig5/amn_lineage3_bp.bar.pdf", raster = TRUE)

colnames(srt_amn@tools$DynamicFeatures_Lineage_Amn$DynamicFeatures)[5] <- "peaktime"
colnames(srt_amn@tools$DynamicFeatures_Lineage_Sfc$DynamicFeatures)[5] <- "peaktime"
ht <- DynamicHeatmap(
  srt = srt_amn, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.1, use_fitted = TRUE,
  lineages = c("Lineage_Amn", "Lineage_Sfc"), split_method = "mfuzz",
  cell_annotation = "SubAnnotation",
  # anno_terms = TRUE,  anno_features = TRUE, padjustCutoff = 1, db_version = "3.13",
  n_split = 8, reverse_ht = 1, use_raster = TRUE, width = 8
)
cluster <- sapply(ht$feature_split, function(x) {
  switch(as.character(x),
    "C1" = "C1",
    "C2" = "C2",
    "C4" = "C3",
    "C6" = "C3",
    "C3" = "C4",
    "C5" = "C4",
    "C7" = "C5",
    "C8" = "C5"
  )
})
cluster <- factor(cluster, levels = levels(ht$feature_split))
cluster <- paste0(cluster, "(", table(cluster)[cluster], ")")
names(cluster) <- names(ht$feature_split)
ht <- DynamicHeatmap(
  srt = srt_amn, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.1, use_fitted = TRUE, feature_split = cluster,
  lineages = c("Lineage_Amn", "Lineage_Sfc"), split_method = "mfuzz",
  cell_annotation = "SubAnnotation",
  feature_annotation = c("TF", "TF_cofactors"), feature_annotation_palcolor = list(c("black", "transparent")),
  anno_terms = TRUE, anno_features = TRUE, padjustCutoff = 1, db_version = "3.13",
  n_split = 8, reverse_ht = 1, use_raster = TRUE, width = 8
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig5/amn_AmnSfc.ht.pdf")
saveRDS(ht, "figures/fig5/amn_AmnSfc.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig5/amn_AmnSfc.enrichment.xlsx")
write.xlsx(ht$feature_metadata, file = "figures/fig5/amn_AmnSfc.feature_metadata.xlsx")

srt_amn@tools$Enrichment_AmnSfc <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_amn@tools$Enrichment_AmnSfc$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig5/amn_Amn&Sfc_5c_bp.bar.pdf", raster = TRUE)

saveRDS(srt_amn, "srt_amn_annotation.rds")

#### all lieages ############
p <- ClassDimPlot(srt_amn, "SubAnnotation3",
  palcolor = palcolor,
  lineages = paste0("Amn_Lineage", 1:5), lineages_trim = c(0.05, 0.95),
  show_stat = FALSE, theme_use = "theme_blank"
)
p <- panel_fix(p, height = 3, save = "figures_add/amn_lineages.umap.pdf", raster = TRUE)


### DE analysis ------------------------------------------------------------------
srt_amn <- RunDEtest(srt_amn,
  group_by = "CSSclusters", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_amn <- RunDEtest(srt_amn,
  group_by = "CSSclusters", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
de_filter <- filter(srt_amn@tools$DEtest_CSSclusters$AllMarkers_wilcox, p_val_adj < 0.05)
res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "tmp_amn_bp.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "tmp_amn_bp.word.pdf")


srt_amn <- RunDEtest(srt_amn,
  group_by = "SubAnnotation", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_amn <- RunDEtest(srt_amn,
  group_by = "SubAnnotation", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
de_filter <- filter(srt_amn@tools$DEtest_SubAnnotation$AllMarkers_wilcox, p_val_adj < 0.05)

ht <- ExpHeatmap(srt_amn,
  group.by = "SubAnnotation", features = de_filter$gene, feature_split = de_filter$group1,
  cluster_rows = T, cluster_row_slices = T, cluster_columns = T, cluster_column_slices = T
)
ggsave(plot = p, filename = "figures/fig5/amn_allDE.heatmap.png", width = 15, height = 10, limitsize = FALSE)

res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "figures/fig5/amn_allDE_bp.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "figures/fig5/amn_allDE_bp.word.pdf")

saveRDS(srt_amn, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_amn_annotation.rds")

### PGC source ------------------------------------------------------------------------------------------
srt_amn$Epiblast_PGC <- srt_amn$Epiblast
srt_amn$Epiblast_PGC[srt_amn$SubAnnotation == "PGC"] <- "PGC"
srt_amn$PSA_PGC <- srt_amn$Epiblast_PGC
srt_amn$PSA_PGC[!srt_amn$PSA_PGC %in% c("PSA-Epi", "PGC")] <- NA
srt_amn$nNEA_PGC <- srt_amn$Epiblast_PGC
srt_amn$nNEA_PGC[!srt_amn$nNEA_PGC %in% c("nNEA-Epi", "PGC")] <- NA
srt_amn <- RunSlingshot(srt_amn,
  group.by = "PSA_PGC", reduction = "UMAP", prefix = "PSA_PGC",
  start = "PSA-Epi", end = c("PGC"),
)
srt_amn <- RunSlingshot(srt_amn,
  group.by = "nNEA_PGC", reduction = "UMAP", prefix = "nNEA_PGC",
  start = "nNEA-Epi", end = c("PGC")
)
ClassDimPlot(srt_amn, "Epiblast_PGC",
  lineages = c("PSA_PGC_Lineage1", "nNEA_PGC_Lineage1"),
  lineages_trim = c(0.2, 0.95), theme_use = "theme_blank"
) %>%
  panel_fix(height = 3)
srt_amn <- RunDynamicFeatures(srt_amn, lineages = c("PSA_PGC_Lineage1", "nNEA_PGC_Lineage1"), n_candidates = 10000)
saveRDS(srt_amn, "srt_amn_pgc_source.rds")
ht1 <- DynamicHeatmap(
  srt = srt_amn,
  lineages = c("PSA_PGC_Lineage1"),
  cell_annotation = "Epiblast_PGC",
  use_fitted = T, height = 5, width = 4
)
ht1$plot

ht2 <- DynamicHeatmap(
  srt = srt_amn,
  lineages = c("nNEA_PGC_Lineage1"),
  cell_annotation = "Epiblast_PGC",
  use_fitted = T, height = 5, width = 4
)
ht2$plot

ht3 <- DynamicHeatmap(
  srt = srt_amn, dev.expl = 0.3, r.sq = 0.3,
  lineages = c("PSA_PGC_Lineage1", "nNEA_PGC_Lineage1"),
  cell_annotation = "Epiblast_PGC",
  n_split = 4,
  reverse_ht = "PSA_PGC_Lineage1",
  use_fitted = T, height = 5, width = 6
)
ht3$plot

library(reticulate)
plt <- import("matplotlib.pyplot")
np <- import("numpy")
pd <- import("pandas")
SCP::check_Python("wot")
wot <- import("wot")

############# example #################
setwd("/home/zhanghao/wot_example/")
# input paths
FULL_DS_PATH <- "data/ExprMatrix.h5ad"
CELL_DAYS_PATH <- "data/cell_days.txt"
VAR_DS_PATH <- "data/ExprMatrix.var.genes.h5ad"
TMAP_PATH <- "tmaps/serum"
CELL_SETS_PATH <- "data/major_cell_sets.gmt"
COORDS_PATH <- "data/fle_coords.txt"

tmap_model <- wot$tmap$TransportMapModel$from_directory(TMAP_PATH)
cell_sets <- wot$io$read_sets(CELL_SETS_PATH, as_dict = TRUE)
populations <- tmap_model$population_from_cell_sets(cell_sets, at_time = 12)
trajectory_ds <- tmap_model$trajectories(populations)
#######################################
srt_amn$time_num <- as.numeric(srt_amn$Time)
adata <- srt_to_adata(srt_amn)
ot_model <- wot$ot$OTModel(adata, epsilon = 0.05, lambda1 = 1, lambda2 = 50, day_field = "time_num")
tmap_annotated <- ot_model$compute_transport_map(1, 2)
ot_model_strict <- wot$ot$OTModel(adata, epsilon = 0.05, lambda1 = 3, lambda2 = 50, day_field = "time_num")
tmap_anno_strict <- ot_model_strict$compute_transport_map(1, 2)
ot_model <- wot$ot$OTModel(adata, epsilon = 0.05, lambda1 = 1, lambda2 = 50, growth_iters = as.integer(3), day_field = "time_num")
ot_model$compute_all_transport_maps(tmap_out = "tmaps/tmap_out")

tmap_model <- wot$tmap$TransportMapModel$from_directory("tmaps/tmap_out")
cell_sets <- split(colnames(srt_amn), srt_amn$SubAnnotation3)
populations <- tmap_model$population_from_cell_sets(cell_sets, at_time = max(srt_amn$time_num))
trajectory_ds <- tmap_model$trajectories(populations)
fates_ds <- tmap_model$fates(populations)

start_populations <- tmap_model$population_from_cell_sets(cell_sets, at_time = 3)
end_populations <- tmap_model$population_from_cell_sets(cell_sets, at_time = max(srt_amn$time_num))
transition_table <- tmap_model$transition_table(start_populations, end_populations)

srt_trajectory <- adata_to_srt(trajectory_ds)
srt_trajectory[["UMAP"]] <- srt_amn[["UMAP"]]
srt_trajectory[["SubAnnotation3"]] <- srt_amn[["SubAnnotation3"]]
srt_fate <- adata_to_srt(fates_ds)
srt_fate[["UMAP"]] <- srt_amn[["UMAP"]]
srt_fate[["SubAnnotation3"]] <- srt_amn[["SubAnnotation3"]]
p1 <- ExpDimPlot(srt_fate, features = "PGC", title = "PGC_fate", theme_use = "theme_blank") %>% panel_fix(height = 2)
p2 <- ExpDimPlot(srt_trajectory, features = "PGC", title = "PGC_ancestors", theme_use = "theme_blank") %>% panel_fix(height = 2)
panel_fix(p1 + p2, height = 2, width = 2, save = "figures_add/PGC_wot.umap.pdf")
p1 <- ExpStatPlot(srt_fate, cells = colnames(srt_fate)[srt_fate$SubAnnotation3 != "PGC"], group.by = "SubAnnotation3", features = "PGC", title = "PGC_fate", y.max = 0.2) %>% panel_fix(height = 2, width = 3)
p2 <- ExpStatPlot(srt_trajectory, cells = colnames(srt_trajectory)[srt_trajectory$SubAnnotation3 != "PGC"], group.by = "SubAnnotation3", features = "PGC", title = "PGC_ancestors", y.max = 0.1) %>% panel_fix(height = 2, width = 3)
p1 + p2
panel_fix(p1 + p2, height = 2, width = 3, save = "figures_add/PGC_wot.stat.pdf")

srt_trajectory$ancestors <- rownames(srt_trajectory@assays$RNA@counts)[apply(srt_trajectory@assays$RNA@counts, 2, which.max)]
srt_fate$fate <- rownames(srt_fate@assays$RNA@counts)[apply(srt_fate@assays$RNA@counts, 2, which.max)]
p1 <- ClassStatPlot(srt_fate, stat.by = "fate", group.by = "SubAnnotation3", title = "fate") %>% panel_fix(height = 2, width = 3)
p2 <- ClassStatPlot(srt_trajectory, stat.by = "SubAnnotation3", group.by = "ancestors", title = "ancestors") %>% panel_fix(height = 2, width = 3)
panel_fix(p1 / p2, height = 2, width = 3, save = "figures_add/PGC_wot.bar.pdf")



srt_transition_table <- adata_to_srt(transition_table)
ht <- ExpHeatmap(srt_transition_table, features = rownames(srt_transition_table), exp_method = "raw", show_column_names = T)
ht$plot





ClassDimPlot(srt_amn, c("SubAnnotation", "Epiblast"), theme_use = "theme_blank", nrow = 2, force = TRUE)
srt_amn2 <- srt_amn[, srt_amn$SubAnnotation == "PGC" | !is.na(srt_amn$Epiblast)]
srt_amn2 <- Integration_SCP(srt_amn2,
  batch = "Time", integration_method = "Uncorrected",
  nonlinear_reduction = c("umap", "tsne", "dm", "phate", "pacmap", "trimap", "largevis")
)
srt_amn2$Epiblast[srt_amn2$SubAnnotation == "PGC"] <- "PGC"
plist <- list()
for (i in c("umap", "tsne", "dm", "phate", "pacmap", "trimap", "largevis")) {
  plist[[i]] <- ClassDimPlot(srt_amn2, "Epiblast", pt.size = 1, reduction = i)
}
p <- plot_grid(plotlist = plist)
panel_fix(p, height = 2)

srt_amn2 <- RunMonocle2(srt = srt_amn2, annotation = "Epiblast", max_components = 2)
ClassDimPlot(srt_amn2, "Epiblast", reduction = "DDRTree", pt.size = 1) %>% panel_fix(height = 3)

srt_amn3 <- srt_amn2[, srt_amn2$Epiblast != "Epiblast-3"]
srt_amn3 <- RunMonocle2(srt = srt_amn3, annotation = "Epiblast", max_components = 2)
ClassDimPlot(srt_amn3, "Epiblast", reduction = "DDRTree", pt.size = 1) %>% panel_fix(height = 3)

ht <- GroupHeatmap(srt_amn2,
  features = VariableFeatures(srt_amn2), group.by = "Epiblast", n_split = 5,
  slot = "data", cluster_rows = T, cluster_columns = T,
  width = 7, height = 7, show_row_names = F,
  nlabel = 20, anno_terms = T, anno_features = T,
  show_column_names = T, column_names_side = "bottom"
)
ht$plot
CellCorHeatmap(
  srt_query = srt_amn2, srt_ref = srt_amn2,
  features = VariableFeatures(srt_amn2),
  query_group = "Epiblast", ref_group = "Epiblast",
  cluster_columns = T, cluster_rows = T
)

srt_amn2 <- RunDEtest(srt_amn2,
  group_by = "Epiblast", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
PGC_markers <- filter(srt_amn2@tools$DEtest_Epiblast$AllMarkers_wilcox, p_val_adj < 0.05 & group1 == "PGC") %>% pull("gene")
ht <- ExpHeatmap(srt_amn2,
  group.by = "Epiblast", features = PGC_markers,
  # anno_terms = T,anno_keys = T,anno_features = T,
  n_split = 4, cluster_rows = T, cluster_column_slices = T, cluster_columns = T,
  height = 7, width = 7
)
ht$plot
Amnion_1_markers <- filter(srt_amn@tools$DEtest_SubAnnotation$AllMarkers_wilcox, p_val_adj < 0.05 & group1 == "Amnion-1") %>% pull("gene")
ht <- ExpHeatmap(srt_amn,
  group.by = "SubAnnotation", features = Amnion_1_markers,
  # anno_terms = T,anno_keys = T,anno_features = T,
  n_split = 4, cluster_rows = T, cluster_column_slices = T, cluster_columns = T,
  height = 7, width = 7
)
ht$plot
PS_markers <- filter(srt_amn@tools$DEtest_SubAnnotation$AllMarkers_wilcox, p_val_adj < 0.05 & group1 == "Primitive streak") %>% pull("gene")
ht <- ExpHeatmap(srt_amn,
  group.by = "SubAnnotation", features = PS_markers,
  # anno_terms = T,anno_keys = T,anno_features = T,
  n_split = 4, cluster_rows = T, cluster_column_slices = T, cluster_columns = T,
  height = 7, width = 7
)
ht$plot

srt_amn <- RunDEtest(srt_amn,
  group_by = "SubAnnotation", test.use = "wilcox",
  group1 = "PGC", group2 = "Amnion-1", only.pos = FALSE,
  BPPARAM = BiocParallel::MulticoreParam(workers = 30)
)
PGC_up <- filter(srt_amn@tools$DEtest_custom$AllMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 0) %>% pull("gene")
PGC_down <- filter(srt_amn@tools$DEtest_custom$AllMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC < 0) %>% pull("gene")
ht <- ExpHeatmap(srt_amn,
  group.by = "SubAnnotation", features = c(PGC_up, PGC_down),
  feature_split = c(rep(paste0("PGC_up", "(", length(PGC_up), ")"), length(PGC_up)), rep(paste0("PGC_down", "(", length(PGC_down), ")"), length(PGC_down))),
  # anno_terms = T,anno_keys = T,anno_features = T,
  n_split = 4, cluster_rows = T, cluster_column_slices = T, cluster_columns = T,
  height = 7, width = 7
)
ht$plot

srt_amn <- RunDEtest(srt_amn,
  group_by = "SubAnnotation", test.use = "wilcox",
  group1 = "PGC", group2 = "Primitive streak", only.pos = FALSE,
  BPPARAM = BiocParallel::MulticoreParam(workers = 30)
)
PGC_up <- filter(srt_amn@tools$DEtest_custom$AllMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC > 0) %>% pull("gene")
PGC_down <- filter(srt_amn@tools$DEtest_custom$AllMarkers_wilcox, p_val_adj < 0.05 & avg_log2FC < 0) %>% pull("gene")
ht <- ExpHeatmap(srt_amn,
  group.by = "SubAnnotation", features = c(PGC_up, PGC_down),
  feature_split = c(rep(paste0("PGC_up", "(", length(PGC_up), ")"), length(PGC_up)), rep(paste0("PGC_down", "(", length(PGC_down), ")"), length(PGC_down))),
  # anno_terms = T,anno_keys = T,anno_features = T,
  n_split = 4, cluster_rows = T, cluster_column_slices = T, cluster_columns = T,
  height = 7, width = 7
)
ht$plot




## Ectoderm --------------------------------------------------------------------
epi <- colnames(srt_epi)[srt_epi$SubAnnotation %in% paste0("Epiblast-", 2)]
cells <- c(epi, colnames(srt_22)[srt_22$CSSclusters %in% c(25, 24, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)])
srt_ecto <- srt_22[, cells]

### Integration  --------------------------------------------------------------------
srt_ecto <- CellCycleScoring(srt_ecto,
  s.features = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes
)

plist <- list()
pagalist <- list()
progress <- 1
total <- length(c(2000, 3000)) * length(seq(15, 50, 5)) * length(seq(0.5, 5, 0.5))
for (nhvf in c(2000, 3000)) {
  for (i in seq(15, 50, 5)) {
    for (j in seq(0.5, 5, 0.5)) {
      cat("progress:", progress, "/", total, "\n")
      progress <- progress + 1
      srt_ecto <- CSS_integrate(
        srtMerge = srt_ecto, nHVF = nhvf,
        linear_reduction_dims_use = 1:i, cluster_resolution = j,
        CSS_param = list(cluster_resolution = j)
      )
      srt_ecto[[paste0("CSSUMAP2Dnhvf", nhvf, "dim", i, "res", j)]] <- srt_ecto@reductions$CSSUMAP2D
      srt_ecto[[paste0("CSSUMAP3Dnhvf", nhvf, "dim", i, "res", j)]] <- srt_ecto@reductions$CSSUMAP3D
      plist[[paste0("CSSUMAP2Dnhvf", nhvf, "dim", i, "res", j)]] <- ClassDimPlot(srt_ecto, "Annotation", title = paste("hvf:", nhvf, "dim", i, "res", j))

      adata <- srt_to_adata(srt_ecto)
      adata <- RunPAGA(
        adata = adata, group_by = "CSSclusters", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
        n_pcs = ncol(srt_ecto@reductions$CSS@cell.embeddings), threshold = 0.1,
        embedded_with_PAGA = TRUE, paga_layout = "fa"
      )
      PAGAUMAP2D <- adata$obsm[["PAGAUMAP2D"]]
      colnames(PAGAUMAP2D) <- c("UMAP_1", "UMAP_2")
      rownames(PAGAUMAP2D) <- adata$obs_names$values
      srt_ecto[["PAGAUMAP2D"]] <- CreateDimReducObject(PAGAUMAP2D, key = "UMAP_", assay = "RNA")
      ClassDimPlot(srt_ecto, group.by = "Time", reduction = "PAGAUMAP2D")
      pagalist[[paste0("PAGAUMAP2Dnhvf", nhvf, "dim", i, "res", j)]] <- ClassDimPlot(srt_ecto, "Annotation", reduction = "PAGAUMAP2D", title = paste("hvf:", nhvf, "dim", i, "res", j))
    }
  }
}
p <- plot_grid(plotlist = plist)
p <- panel_fix(p, save = "tmp_ect1.png", width = 3)
p <- plot_grid(plotlist = pagalist)
p <- panel_fix(p, save = "tmp_ect2.png", width = 3)

srt_ecto <- Integration_SCP(
  srtMerge = srt_ecto, integration_method = "Uncorrected",
  linear_reduction_dims_use = 1:35, cluster_resolution = 5
)
srt_ecto <- Integration_SCP(
  srtMerge = srt_ecto, integration_method = "CSS",
  linear_reduction_dims_use = 1:35, cluster_resolution = 5,
  CSS_param = list(cluster_resolution = 5)
)
# 50-5 potential energy
ClassDimPlot(srt_ecto, c("Annotation", "CSSclusters", "Phase", "Time"), label = T) # 25-5
ExpDimPlot(srt_ecto, c("percent.ribo", "percent.mito", "nFeature_RNA", "nCount_RNA"))

srt_ecto <- FindClusters(object = srt_ecto, resolution = 8, algorithm = 1, graph.name = "CSS_SNN")
srt_ecto <- SrtReorder(srt_ecto, features = VariableFeatures(srt_ecto), reorder_by = "seurat_clusters", slot = "data")
srt_ecto[["seurat_clusters"]] <- NULL
srt_ecto[["CSSclusters"]] <- Idents(srt_ecto)
ClassDimPlot(srt_ecto, c("CSSclusters"), label = TRUE)

### paga initiated layout ---------------------------------------------------------------
adata <- srt_to_adata(srt_ecto)
adata <- RunPAGA(
  adata = adata, group_by = "CSSclusters", linear_reduction = "CSS", nonlinear_reduction = "CSSUMAP2D",
  n_pcs = ncol(srt_ecto@reductions$CSS@cell.embeddings), threshold = 0.1,
  embedded_with_PAGA = TRUE, paga_layout = "fa"
)
PAGAUMAP2D <- adata$obsm[["PAGAUMAP2D"]]
colnames(PAGAUMAP2D) <- c("UMAP_1", "UMAP_2")
rownames(PAGAUMAP2D) <- adata$obs_names$values
srt_ecto[["PAGAUMAP2D"]] <- CreateDimReducObject(PAGAUMAP2D, key = "UMAP_", assay = "RNA")
ClassDimPlot(srt_ecto, group.by = "Time", reduction = "PAGAUMAP2D")
srt_ecto[["UMAP"]] <- srt_ecto[["PAGAUMAP2D"]]
srt_ecto@misc$Default_reduction <- "UMAP"


### Annotation --------------------------------------------------------------
ClassDimPlot(srt_ecto, "Annotation")
ClassDimPlot(srt_ecto, "CSSclusters", label = TRUE)

srt_ecto <- RunKNNMap(
  srt_query = srt_ecto, srt_ref = srt_esc,
  ref_group = "CSSclusters", ref_umap = "UMAP"
)
ProjectionPlot(
  srt_query = srt_ecto, srt_ref = srt_esc,
  query_group = "CSSclusters", ref_group = "Annotation"
)
ClassDimPlot(srt_ecto1, group.by = "CSSclusters", reduction = "ref.umap", label = T)
srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, srt_ref = srt_esc,
  ref_group = "Annotation", prefix = "esc"
)
ClassDimPlot(srt_ecto, "esc_classification", label = TRUE)

srt_ecto <- RunKNNMap(
  srt_query = srt_ecto, srt_ref = srt_h0,
  ref_group = "CSSclusters", ref_umap = "UMAP"
)
ProjectionPlot(
  srt_query = srt_ecto, srt_ref = srt_h0,
  query_group = "CSSclusters", ref_group = "Annotation"
)
ClassDimPlot(srt_ecto, group.by = "CSSclusters", reduction = "ref.umap", label = T)
srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, srt_ref = srt_h0,
  ref_group = "Annotation", prefix = "h0"
)
ClassDimPlot(srt_ecto, "h0_classification", label = TRUE)

srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, query_group = "CSSclusters",
  srt_ref = srt_ecto, ref_group = "CSSclusters",
  query_collapsing = T, ref_collapsing = T,
  return_full_distance_matrix = TRUE, nn_method = "raw"
)
d <- 1 - srt_ecto@tools$knnpredict_classification$distance_matrix
d <- as.matrix(d)
Heatmap(d)

srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, bulk_ref = SCP::ref_scHCL,
  filter_lowfreq = 10, prefix = "HCL"
)
ClassDimPlot(srt_ecto, group.by = "HCL_classification", label = T)

srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, srt_ref = srt_mOrg,
  ref_group = "Sub_trajectory_name",
  filter_lowfreq = 10, prefix = "mOrg1"
)
ClassDimPlot(srt_ecto, group.by = "mOrg1_classification", label = T)

srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, srt_ref = srt_mOrg,
  ref_group = "Main_cell_type",
  filter_lowfreq = 10, prefix = "mOrg2"
)
ClassDimPlot(srt_ecto, group.by = "mOrg2_classification", label = T)

srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, srt_ref = srt_hOrg,
  ref_group = "annotation",
  filter_lowfreq = 10, prefix = "hOrg"
)
ClassDimPlot(srt_ecto, group.by = "hOrg_classification", label = T)

srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, srt_ref = srt_mBrain,
  ref_group = "Subclass",
  filter_lowfreq = 10, prefix = "mBrain"
)
ClassDimPlot(srt_ecto, group.by = "mBrain_classification", label = T)

srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, srt_ref = srt_hBrain,
  ref_group = "cell.type",
  filter_lowfreq = 10, prefix = "hBrain"
)
ClassDimPlot(srt_ecto, group.by = "hBrain_classification", label = T)

srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, srt_ref = srt_retina1,
  ref_group = "CELLTYPE1",
  filter_lowfreq = 10, prefix = "retina1"
)
ClassDimPlot(srt_ecto, group.by = "retina1_classification", label = T)

srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, srt_ref = srt_retina2,
  ref_group = "type",
  filter_lowfreq = 10, prefix = "retina2"
)
ClassDimPlot(srt_ecto, group.by = "retina2_classification", label = T)

srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, srt_ref = srt_retina3,
  ref_group = "cell_type",
  filter_lowfreq = 10, prefix = "retina3"
)
ClassDimPlot(srt_ecto, group.by = "retina3_classification", label = T)

srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, srt_ref = srt_retina4,
  ref_group = "umap2_CellType",
  filter_lowfreq = 10, prefix = "retina4"
)
ClassDimPlot(srt_ecto, group.by = "retina4_classification", label = T)

srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, srt_ref = srt_retina5,
  ref_group = "CellType",
  filter_lowfreq = 10, prefix = "retina5"
)
ClassDimPlot(srt_ecto, group.by = "retina5_classification", label = T)

srt_ecto <- RunKNNPredict(
  srt_query = srt_ecto, srt_ref = srt_teratoma,
  ref_group = "celltype",
  filter_lowfreq = 10, prefix = "teratoma"
)
ClassDimPlot(srt_ecto, group.by = "teratoma_classification", label = T)

saveRDS(srt_ecto, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_ecto_annotation.rds")


# Axial ectoderm: CHRD NOTO
# Paraial ectoderm: TBX6 SFRP2
# LPM: BMP4 PDGFRA GATA6
# Cardiac progenitor: ISL1 MAB21L2 HOPX
# Contractile: MYL7 TNNI1 TNNT2 MYH10

# Radial glial cell: "SOX1","RFX4","TTYH1","FABP7","NES"
# Ependymal cell: "FOXJ1","PIFO","DYNLRB2"
# Neuron: "DCX","TUBB3","TAGLN3","SLC17A6"
# Schwann cell: "SOX10","MPZ","PLP1"
# Retinal progenitor cell: "SIX6", "VSX2", "HMX1"
# Retinal pigmented epithelium: "MLANA","MITF","SERPINF1","RLBP1"
# Dividing: "MKI67", "TOP2A", "CDK1"
# Preplacodal ectoderm: "SIX4","DLX3","SOX3","FOXI3"
# ventral:"PAX6", "FOXP1"
# dorsal: "PAX2", "PAX7","FABP7"
# Astrocytes: c("SLC1A3", "FABP7", "FGFR3","CD44", "NFIA","NFIB" ,"SLIT1", "PLPP3")
# Oligodendrocyte: c("CSPG4","OLIG1","OLIG2","PDGFRA","SOX9")

c("PAX6", "FOXP1") # ventral
c("PAX2", "PAX7") # dorsal FABP7
c(
  "POU5F1", "NANOG", # Epiblast
  "CDH1", "NOTCH1", # Neuroepithelial cell
  "CDH2", "BMP7", # Radial glial cell
  "SOX2", "NES", "VIM", # Neuroepithelial cell & Radial glial cell
  "PAX6", # Retina & & Radial glial cell
  "PDGFRA", "TMEM132D", # Oligodendrocyte
  "SLC1A3", "FABP7", "SLIT1", # Astrocytes
  "RAX", "SIX6", "VSX2", # Retinal progenitor cells
  "MITF", "SERPINF1", # Retinal pigment epithelial cell
  "FOXJ1", "PIFO", # Ependymal cell
  "MKI67" # Prolifrating cell
)



annotation_list <- list(
  "Epiblast" = c(12, 13, 14, 15, 17, 18, 19, 20),
  "Preplacodal ectoderm" = c(16),
  "Neuroepithelial cell-1" = c(4, 5, 1),
  "Neuroepithelial cell-2" = c(2, 3),
  "Neuroepithelial cell-3" = c(7, 8, 9),
  "Neuroepithelial cell-4" = c(6, 10),
  "Radial glial cell(dorsal)" = c(11),
  "Radial glial cell(ventral)" = c(42, 41),
  "Ependymal cell" = c(40, 43),
  "Astrocytes" = c(44),
  "Neuron" = c(45),
  "Retinal progenitor cell" = c(21, 22, 23, 24),
  "Retinal pigmented epithelium" = c(25:39)
)
annotation <- setNames(rep(names(annotation_list), sapply(annotation_list, length)), nm = unlist(annotation_list))
srt_ecto$SubAnnotation <- annotation[as.character(srt_ecto$CSSclusters)]
srt_ecto$SubAnnotation <- factor(srt_ecto$SubAnnotation, levels = names(annotation_list))
ClassDimPlot(srt_ecto, "SubAnnotation", label = T)

ExpDimPlot(srt_ecto, c("SIX6", "VSX2", "HMX1"))
ExpDimPlot(srt_ecto, c("SIX4", "DLX3", "SOX3", "FOXI3"))
ExpDimPlot(srt_ecto, c("MLANA", "MITF", "SERPINF1", "RLBP1"))
ExpDimPlot(srt_ecto, c("SOX10", "MPZ", "PLP1"))

saveRDS(srt_ecto, "/ssd/lab/HuangMingQian/scRNAseq/22CellLine-project/srt_ecto_annotation.rds")

### Baisc plot ------------------------------------------------------------
srt_22_tmp <- srt_22
srt_22_tmp$GermLayer[!colnames(srt_22) %in% colnames(srt_ecto)] <- NA
srt_22_tmp$GermLayer <- factor(srt_22_tmp$GermLayer, levels = levels(srt_22$GermLayer))
p <- ClassDimPlot(srt_22_tmp,
  group.by = "GermLayer", theme_use = "theme_blank",
  cells.highlight = colnames(srt_ecto), show_stat = FALSE
)
p <- panel_fix(p, height = 3, save = "figures/fig5/ecto_highlight.umap.pdf", raster = TRUE)

srt_ecto$SubAnnotation[srt_ecto$SubAnnotation == "Preplacodal ectoderm"] <- "Epiblast"
srt_ecto$SubAnnotation <- factor(srt_ecto$SubAnnotation, levels = setdiff(levels(srt_ecto$SubAnnotation), "Preplacodal ectoderm"))
p <- ClassDimPlot(srt_ecto, "SubAnnotation", label = TRUE, theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
panel_fix(p, height = 3, save = "figures_add/ecto_annotation.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_ecto, "Time", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
panel_fix(p, height = 3, save = "figures/fig5/ecto_time.umap.pdf", raster = TRUE)
p <- ClassDimPlot(srt_ecto, "Phase", theme_use = "theme_blank", show_stat = FALSE, force = TRUE)
panel_fix(p, height = 3, save = "figures/fig5/ecto_cellcycle.umap.pdf", raster = TRUE)

markers_ecto <- list(
  "Epiblast" = c("POU5F1", "DNMT3B"),
  # "Preplacodal ectoderm" = c("SIX4", "DLX3"),
  "Neuroepithelial cell" = c("SOX2", "CDH2"),
  "Radial glial cell" = c("NES", "TTYH1"),
  "Ependymal cell" = c("FOXJ1", "CFAP157"),
  "Astrocytes" = c("SLC1A3", "FABP7"),
  "Neuron" = c("DCX", "TAGLN3"),
  "Retinal progenitor cell" = c("SIX6", "VSX2"),
  "Retinal pigmented epithelium" = c("MITF", "MLANA")
  # "Ventral/Middle" = c("FOXP1", "ARX", "PAX6"), # PAX6 middle
  # "Dorsal" = c("PAX7", "PAX3", "BMP5")
)
p <- ExpDimPlot(srt_ecto, unlist(markers_ecto),
  theme_use = "theme_blank", ncol = 6
)
p <- panel_fix(p, height = 2, save = "figures/fig5/ecto_markers.umap.pdf", raster = TRUE)

ht <- GroupHeatmap(srt_ecto,
  group.by = "SubAnnotation", heatmap_palette = "YlOrRd",
  features = unlist(markers_ecto),
  feature_split = rep(names(markers_ecto), each = 2),
  show_row_names = TRUE, add_dot = TRUE, add_bg = TRUE,
  height = 5.5, width = 7, dot_size = unit(7, "mm")
)
panel_fix(ht$plot, save = "figures_add/ectoderm.markers.pdf")


p <- ExpDimPlot(srt_ecto, c("MKI67", "TOP2A"), theme_use = "theme_blank")
p <- panel_fix(p, height = 2, save = "figures/fig5/ecto_prolif.umap.pdf", raster = TRUE)

p <- ExpStatPlot(srt_ecto, c("CDH1", "CDH2"), group.by = "SubAnnotation")
panel_fix(p, height = 2, save = "figures_add/ecto_CDH.vln.pdf", raster = TRUE)


p <- ExpDimPlot(srt_ecto, c("SOX10", "MPZ", "PLP1"),
  calculate_coexp = TRUE, ncol = 4,
  theme_use = "theme_blank"
)
p <- panel_fix(p, height = 2, save = "figures/fig5/ecto_schwann.umap.pdf", raster = TRUE)

markers_vd <- list(
  "Ventral/Middle" = c("FOXP1", "ARX", "PAX6"), # PAX6 middle
  "Dorsal" = c("PAX7", "PAX3", "BMP5")
)
p <- ExpDimPlot(srt_ecto, features = unlist(markers_vd), theme_use = "theme_blank")
p <- panel_fix(p, height = 2, save = "figures/fig5/ecto_vd.umap.pdf", raster = TRUE)
p <- ExpVlnPlot(srt_ecto,
  features = unlist(markers_vd), group.by = "SubAnnotation", legend.position = "none", xlab = "",
  cells_subset = WhichCells(srt_ecto, expression = SubAnnotation %in% c("Radial glial cell(dorsal)", "Radial glial cell(ventral)"))
)
p <- panel_fix(p, height = 1, save = "figures/fig5/ecto_vd.vln.pdf", raster = TRUE)

### PAGA ---------------------------------------------------------------------
srt_ecto_paga <- RunPAGA(
  srt = srt_ecto, group_by = "SubAnnotation", linear_reduction = "CSS", nonlinear_reduction = "UMAP",
  n_pcs = ncol(srt_ecto@reductions$CSS@cell.embeddings), embedded_with_PAGA = FALSE, return_seurat = TRUE
)
srt_ecto@misc$paga <- srt_ecto_paga@misc$paga
p <- ClassDimPlot(srt_ecto,
  group.by = "SubAnnotation", pt.size = 5, pt.alpha = 0.05,
  label = TRUE, label.size = 3, label_repel = TRUE, label_insitu = TRUE, label_segment_color = "transparent",
  paga = srt_ecto@misc$paga, paga_edge_threshold = 0.1, paga_edge_color = "black", paga_edge_alpha = 1,
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
)
panel_fix(p, width = 4, save = "figures_add/ecto_paga.pdf", raster = TRUE)

### SCVELO ------------------------------------------------------------------
srt_ecto_scv <- SCP:::RunSCVELO(
  srt = srt_ecto, group_by = "SubAnnotation", linear_reduction = "CSS", nonlinear_reduction = "UMAP",
  n_pcs = ncol(srt_ecto@reductions$CSS@cell.embeddings), n_neighbors = 100, return_seurat = TRUE
)
srt_ecto[["stochastic_UMAP"]] <- srt_ecto_scv[["stochastic_UMAP"]]
srt_ecto[["Ms"]] <- srt_ecto_scv[["Ms"]]
srt_ecto[["Mu"]] <- srt_ecto_scv[["Mu"]]
srt_ecto[["stochastic"]] <- srt_ecto_scv[["stochastic"]]
srt_ecto[["variance_stochastic"]] <- srt_ecto_scv[["variance_stochastic"]]
p <- ClassDimPlot(srt_ecto,
  group.by = "SubAnnotation", pt.size = 5, pt.alpha = 0.1,
  velocity = "stochastic", velocity_plot_type = "stream",
  velocity_density = 2, velocity_smooth = 1,
  streamline_n = 15, streamline_size = 0.5, streamline_color = "black",
  show_stat = FALSE, legend.position = "none", theme_use = "theme_blank"
)
panel_fix(p, width = 4, save = "figures_add/ecto_scvelo.pdf", raster = TRUE)

saveRDS(srt_ecto, "srt_ecto_annotation.rds")

### Slingshot ---------------------------------------------------------------
ClassDimPlot(srt_ecto, "SubAnnotation", label = T)
srt_ecto_sub <- srt_ecto[, !srt_ecto$SubAnnotation %in% c("Preplacodal ectoderm", "Neuron")]
srt_ecto_sub <- RunSlingshot(srt_ecto_sub,
  group.by = "SubAnnotation", reduction = "UMAP", prefix = "ecto",
  start = c("Epiblast"),
  end = c(
    "Retinal pigmented epithelium",
    "Ependymal cell", "Astrocytes"
  )
)
srt_ecto_sub@tools$Slingshot_SubAnnotation_UMAP@metadata$plot
ExpDimPlot(srt_ecto_sub, paste0("ecto_Lineage", 1:3))
srt_ecto_sub$Lineage_EpiEP <- srt_ecto_sub$ecto_Lineage1
srt_ecto_sub$Lineage_EpiRPE <- srt_ecto_sub$ecto_Lineage2
srt_ecto_sub$Lineage_EpiAst <- srt_ecto_sub$ecto_Lineage3
srt_ecto_sub <- RunDynamicFeatures(srt_ecto_sub,
  lineages = c("Lineage_EpiEP", "Lineage_EpiRPE", "Lineage_EpiAst"),
  n_candidates = 10000
)
saveRDS(srt_ecto_sub, "srt_ecto_sub.rds")

colnames(srt_ecto_sub@tools$DynamicFeatures_Lineage_EpiEP$DynamicFeatures)[5] <- "peaktime"
colnames(srt_ecto_sub@tools$DynamicFeatures_Lineage_EpiRPE$DynamicFeatures)[5] <- "peaktime"
ht <- DynamicHeatmap(
  srt = srt_ecto_sub, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_EpiEP", "Lineage_EpiRPE"),
  cell_annotation = "SubAnnotation", nlabel = 30,
  # anno_terms = TRUE,  anno_features = TRUE, padjustCutoff = 1, db_version = "3.13",
  n_split = 8, reverse_ht = 1, use_raster = TRUE
)
cluster <- sapply(ht$feature_split, function(x) {
  switch(as.character(x),
    "C1" = "C1",
    "C2" = "C2",
    "C3" = "C2",
    "C4" = "C3",
    "C5" = "C4",
    "C7" = "C5",
    "C6" = "C6",
    "C8" = "C6"
  )
})
cluster <- factor(cluster, levels = levels(ht$feature_split))
cluster <- paste0(cluster, "(", table(cluster)[cluster], ")")
names(cluster) <- names(ht$feature_split)
ht <- DynamicHeatmap(
  srt = srt_ecto_sub, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE, feature_split = cluster,
  lineages = c("Lineage_EpiEP", "Lineage_EpiRPE"),
  cell_annotation = "SubAnnotation", nlabel = 30,
  anno_terms = TRUE, anno_features = TRUE, padjustCutoff = 1, db_version = "3.13",
  n_split = 8, reverse_ht = 1, use_raster = TRUE
)
p <- ht$plot
ggsave(p, height = 12, width = 20, filename = "figures/fig5/ecto_lineageEPRPE.ht.pdf")
saveRDS(ht, "figures/fig5/ecto_lineageEPRPE.heatmap.rds")
write.xlsx(ht$enrichment$enrichment, file = "figures/fig5/ecto_lineageEPRPE.enrichment.xlsx")


p <- ExpDimPlot(srt_ecto_sub,
  features = c("Lineage_EpiEP", "Lineage_EpiRPE"), ncol = 1,
  palette = "cividis", byrow = FALSE, theme_use = "theme_blank", force = TRUE
)
p <- panel_fix(p, height = 2, save = "figures/fig5/ecto_lineageEPRPE.umap.pdf", raster = TRUE)
srt_ecto_sub@tools$Enrichment_EPRPE <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_ecto_sub@tools$Enrichment_EPRPE$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig5/ecto_lineageEPRPE_bp.bar.pdf", raster = TRUE)

ht <- DynamicHeatmap(
  srt = srt_ecto_sub, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_EpiEP"), split_method = "kmeans-peaktime",
  cell_annotation = "SubAnnotation",
  n_split = 6, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig5/ecto_lineageEpiEP.ht.pdf")
srt_ecto_sub@tools$Enrichment_EpiEP <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_ecto_sub@tools$Enrichment_EpiEP$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig5/ecto_lineageEpiEP_bp.bar.pdf", raster = TRUE)

ht <- DynamicHeatmap(
  srt = srt_ecto_sub, r.sq = 0.1, dev.expl = 0.1, min_expcells = 20,
  cell_density = 0.01, use_fitted = TRUE,
  lineages = c("Lineage_EpiRPE"), split_method = "kmeans-peaktime",
  cell_annotation = "SubAnnotation",
  n_split = 6, use_raster = TRUE, width = 7
)
p <- ht$plot
ggsave(p, height = 12, width = 12, filename = "figures/fig5/ecto_lineageEpiRPE.ht.pdf")
srt_ecto_sub@tools$Enrichment_EpiRPE <- RunEnrichment(geneID = names(ht$feature_split), geneID_groups = ht$feature_split, db_version = "3.13")
p <- EnrichmentPlot(res = srt_ecto_sub@tools$Enrichment_EpiRPE$enrichment, db = "GO_BP", plot_type = "bar", padjustCutoff = 1, combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 2, width = 2, margin = 0, save = "figures/fig5/ecto_lineageEpiRPE_bp.bar.pdf", raster = TRUE)

saveRDS(srt_ecto_sub, "srt_ecto_sub.rds")

#### all lieages ############
srt_ecto[["Lineage_EpiEP"]] <- srt_ecto_sub[["Lineage_EpiEP"]]
srt_ecto[["Lineage_EpiRPE"]] <- srt_ecto_sub[["Lineage_EpiRPE"]]
p <- ClassDimPlot(srt_ecto, "SubAnnotation",
  lineages = c("Lineage_EpiEP", "Lineage_EpiRPE"),
  show_stat = FALSE, theme_use = "theme_blank"
)
panel_fix(p, height = 3, save = "figures_add/ecto_lineages.umap.pdf", raster = TRUE)

### DE analysis ------------------------------------------------------------------
srt_ecto <- RunDEtest(srt_ecto,
  group_by = "CSSclusters", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_ecto <- RunDEtest(srt_ecto,
  group_by = "CSSclusters", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
de_filter <- filter(srt_ecto@tools$DEtest_CSSclusters$AllMarkers_wilcox, p_val_adj < 0.05)
res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "tmp_ecto_bp.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "tmp_ecto_bp.word.pdf")

srt_ecto <- RunDEtest(srt_ecto,
  group_by = "SubAnnotation", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
srt_ecto <- RunDEtest(srt_ecto,
  group_by = "SubAnnotation", test.use = "roc",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)
de_filter <- filter(srt_ecto@tools$DEtest_SubAnnotation$AllMarkers_wilcox, p_val_adj < 0.05)
res <- RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "tmp_ecto_bp2.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "tmp_ecto_bp2.word.pdf")

# tmp  ---------------------------------------------------------------
ClassDimPlot(srt_ecto, "CSSclusters")
srt_ecto <- RunDEtest(srt_ecto,
  group_by = "CSSclusters", test.use = "wilcox",
  BPPARAM = BiocParallel::MulticoreParam(workers = 30), force = TRUE
)

de_filter <- filter(srt_ecto@tools$DEtest_CSSclusters$AllMarkers_wilcox, p_val_adj < 0.05)
res <- SCP::RunEnrichment(geneID = de_filter$gene, geneID_groups = de_filter$group1, db_version = "3.13")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "bar", combine = TRUE, ncol = 3)
p <- panel_fix(p, height = 1.5, width = 2, save = "tmp_ecto_bp.bar.pdf")
p <- EnrichmentPlot(res$enrichment, db = "GO_BP", plot_type = "wordcloud", combine = TRUE, ncol = 6)
p <- panel_fix(p, height = 2.5, width = 2.5, save = "tmp_ecto_bp.word.pdf")

srt_ecto <- RunKNNMap(
  srt_query = srt_ecto, srt_ref = srt_esc, ref_umap = "CSSUMAP2D", ref_group = "Annotation"
)
ProjectionPlot(
  srt_query = srt_ecto, query_group = "CSSclusters",
  srt_ref = srt_esc, ref_group = "Annotation"
)
