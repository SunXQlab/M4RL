#############
#  library  #
#############

library(spatstat)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggsci)

setwd("/ST_CCC_results/")
rm(list = ls())
gc()

# Seurat-obj --------------------------------------------------------------

# choose your directory
my_working_dir = paste0(getwd())

#InstallData("ssHippo")
seurat.obj <- readRDS(paste0("/prior_data/UKF262_seurat_object_CT5.rds"))

## color
celltype <- unique(seurat.obj$cell_type)
scales::show_col(pal_igv(palette = "default", alpha = 0.8)(20))
mycolors_nejm <- pal_igv(palette = "default", alpha = 0.8)(20)

mycolor_ct <- mycolors_nejm[1:length(celltype)]
names(mycolor_ct) <- celltype
scales::show_col(mycolor_ct)
mycolor_ct

## plot
seurat.obj@meta.data$sdimx <- seurat.obj@reductions$spatial@cell.embeddings[,1]
seurat.obj@meta.data$sdimy <- seurat.obj@reductions$spatial@cell.embeddings[,2]
plot_df <- seurat.obj@meta.data

p <- ggplot(plot_df, aes(x = sdimx, y = sdimy, color = cell_type)) +
  geom_point(size = 1.0) +
  scale_color_manual(values = mycolor_ct) +
  coord_fixed() +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# save cell type distribution
pdf("./Spatial_Plot.pdf", width = 6, height = 4)
print(p)
dev.off()

# Inputs --------------------------------------------------------------------

## filtering
table(seurat.obj$cell_type)
seurat.obj_filter = subset(seurat.obj, subset = cell_type != "Unassigned")

# annotation
table(seurat.obj_filter$cell_type)
df_anno = data.frame(Barcode=colnames(seurat.obj_filter),
                     Cluster=seurat.obj_filter$cell_type)
head(df_anno)
unique(celltype)

# expression
df_count = as.data.frame(as.matrix(seurat.obj_filter@assays$RNA@counts))
df_count[1:10,1:10]
dim(df_count)

df_norm = as.data.frame(as.matrix(seurat.obj_filter@assays$RNA@data))
df_norm[1:10,1:10]
dim(df_norm)

if (identical(df_count, df_norm)) {
  message("Need Normalization!")
  seurat.obj_filter <- SCTransform(seurat.obj_filter, verbose = TRUE)
  message("Run SCTransform")
  df_norm = as.data.frame(as.matrix(seurat.obj_filter@assays$SCT@data))
  df_norm[1:10,1:10]
}
# location
df_loca <- as.data.frame(Embeddings(seurat.obj_filter[["spatial"]]))
head(df_loca)

# signals
Idents(seurat.obj_filter) <- "cell_type"

st_markers <- FindMarkers(
  seurat.obj_filter,
  assay = "SCT",
  slot = "data",
  ident.1 = "Malignant",
  ident.2 = "Macrophages",
  logfc.threshold = 0.6,
  min.pct = 0.05
)
st_markers$ident.1 <- "Malignant"
st_markers$gene <- rownames(st_markers)
df_markers <- st_markers[st_markers$p_val_adj<=0.05,]
ICGs_list <- split(df_markers$gene,df_markers$ident.1)
str(ICGs_list)

Databases <- readRDS('/prior_data/Databases.rds')
ligs_in_db <- Databases$LigRec.DB$source %>% unique() 
ligs_in_db <- intersect(ligs_in_db, rownames(seurat.obj_filter))
recs_in_db <- Databases$LigRec.DB$target %>% unique() 
recs_in_db <- intersect(recs_in_db, rownames(seurat.obj_filter))

expr.ct <- 0.001
pct.ct <- 0.001
data <- as.matrix(df_norm)
clusters <- celltype %>% as.character() %>% unique()

abundant.cutoff = expr.ct
all_mean <- rowMeans(data)
hist(log10(all_mean), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))
abline(v=log10(abundant.cutoff), col="red", lwd=2, lty=2)

meanExpr_of_LR <- lapply(clusters, function(cluster){
  
  cluster.ids <- df_anno$Barcode[celltype == cluster]
  source_mean <- rowMeans(data[,cluster.ids])
  names(source_mean) <- rownames(data)
  source_mean
  
}) %>% do.call('cbind',.) %>% as.data.frame()
colnames(meanExpr_of_LR) <- clusters

pct_of_LR <- lapply(clusters, function(cluster){
  
  cluster.ids <- df_anno$Barcode[celltype == cluster]
  dat <- data[,cluster.ids]
  pct <- rowSums(dat>0)/ncol(dat)
  names(pct) <- rownames(data)
  pct
  
}) %>% do.call('cbind',.) %>% as.data.frame()
colnames(pct_of_LR) <- clusters

Recs_expr_list <- lapply(clusters, function(cluster){
  
  recs <- rownames(data)[meanExpr_of_LR[,cluster] >= expr.ct & pct_of_LR[,cluster] >= pct.ct]
  intersect(recs, recs_in_db)
  
})
names(Recs_expr_list) <- clusters
str(Recs_expr_list)

Ligs_expr_list <- lapply(clusters, function(cluster){
  
  ligs <- rownames(data)[meanExpr_of_LR[,cluster] >= expr.ct & pct_of_LR[,cluster] >= pct.ct]
  intersect(ligs, ligs_in_db)
  
})
names(Ligs_expr_list) <- clusters
str(Ligs_expr_list)

rownames(df_count) <- toupper(rownames(df_count))
rownames(df_norm) <- toupper(rownames(df_norm))
ICGs_list <- lapply(ICGs_list, toupper)
Ligs_expr_list <- lapply(Ligs_expr_list, toupper)
Recs_expr_list <- lapply(Recs_expr_list, toupper)

# save
save(ICGs_list,Ligs_expr_list,Recs_expr_list,
     file = paste0(my_working_dir, "/seurat_output.rda"))

## Reference
## L. Yan, J. Cheng, Q. Nie, X. Sun, Dissecting multilayer cell-cell communications 
## with signaling feedback loops from spatial transcriptomics data. 
## Genome Research, gr.279857.124 [Online ahead of print] (2025). https://doi.org/10.1101/gr.279857.124
