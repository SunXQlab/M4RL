#############
#  library  #
#############

library(Matrix)
library(dplyr)
library(Seurat)

setwd("/ST_CCC_results/")
sapply("/prior/codes_from_stMLnet_Rpackage.R", source)


###############
## get MLnet ##
###############

## load ####
Databases <- Databases

## data

GCMat <- df_norm
BarCluTable <- df_anno
clusters <- BarCluTable$Cluster %>% as.character() %>% unique()

#load("./st_decompostion_output.rda")
Ligs_up_list <- Ligs_expr_list

str(Ligs_up_list)
str(Recs_expr_list)
str(ICGs_list)

## parameters

wd <- paste0("./runscMLnet/")
dir.create(wd,recursive = T)

## database ####

quan.cutoff = 0.98

RecTF.DB <- Databases$RecTF.DB %>% 
  .[.$score > quantile(.$score, quan.cutoff),] %>%
  dplyr::distinct(source, target)

LigRec.DB <- Databases$LigRec.DB %>%
  dplyr::distinct(source, target) %>%
  dplyr::filter(target %in% RecTF.DB$source)

TFTG.DB <- Databases$TFTG.DB %>%
  dplyr::distinct(source, target) %>%
  dplyr::filter(source %in% RecTF.DB$target)

## get Macrophages-Malignant multi-layer ####

LigClu <- "Macrophages"
RecClu <- "Malignant"
message(paste(LigClu,RecClu,sep = "-"))

MLnet <- mainfunc_with_pval(LigClu, RecClu, wd,RecTF.method = "Fisher",TFTG.method = "Fisher")
MLnet_list <- list()
MLnet_list[[1]] <- MLnet
Output <- matrix(ncol = 1, nrow = 10) %>% as.data.frame()
Output[,1] <- c(Ligs_up_list[[LigClu]] %>% length(),
                Recs_expr_list[[RecClu]] %>% length(),
                ICGs_list[[RecClu]] %>% length(),
                nrow(MLnet$LigRec),nrow(MLnet$RecTF),nrow(MLnet$TFTar),
                ifelse(nrow(MLnet$LigRec)==0,0,MLnet$LigRec$source %>% unique() %>% length()),
                ifelse(nrow(MLnet$LigRec)==0,0,MLnet$LigRec$target %>% unique() %>% length()),
                ifelse(nrow(MLnet$TFTar)==0,0,MLnet$TFTar$source %>% unique() %>% length()),
                ifelse(nrow(MLnet$TFTar)==0,0,MLnet$TFTar$target %>% unique() %>% length()))



names(MLnet_list) <- paste(LigClu,RecClu,sep = "_")
colnames(Output) <- paste(LigClu,RecClu,sep = "_")
rownames(Output) <- c("Lig_bk","Rec_bk","ICG_bk",
                      "LRpair","RecTFpair","TFTGpair",
                      "Ligand", "Receptor", "TF", "TG")

write.csv(Output, file = paste0(wd,"/TME_",RecClu,".csv"))

## Reference
## L. Yan, J. Cheng, Q. Nie, X. Sun, Dissecting multilayer cell-cell communications 
## with signaling feedback loops from spatial transcriptomics data. 
## Genome Research, gr.279857.124 [Online ahead of print] (2025). https://doi.org/10.1101/gr.279857.124