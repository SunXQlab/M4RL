#############
#  library  #
#############

library(dplyr)
library(ggsci)
library(ggplot2)
library(igraph)
library(ggraph)

setwd("/ST_CCC_results/")

sapply("/prior/codes_from_stMLnet_Rpackage.R", source)
###########
## color ##
###########

## load ####

res_path <- "/ST_CCC_results/"
load(paste0(res_path,'seurat_output.rda'))

# celltype

celltype <- unique(df_anno$Cluster)

scales::show_col(pal_igv(palette = "default", alpha = 0.8)(15))
mycolors_nejm <- pal_igv(palette = "default", alpha = 0.8)(15)

mycolor_ct <- mycolors_nejm[1:length(celltype)]
names(mycolor_ct) <- celltype
scales::show_col(mycolor_ct)

# nodekey

scales::show_col(pal_locuszoom(palette = "default", alpha = 0.8)(7))
mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.8)(7)

nodekey <- c("Ligand","Receptor","TF","Target")
mycolor_key <- mycolors_locus[1:4]
names(mycolor_key) <- nodekey
scales::show_col(mycolor_key)

# nodetype

scales::show_col(pal_locuszoom(palette = "default", alpha = 0.8)(7))
mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.8)(7)

nodetype <- c("cell","Sender","Receiver")
mycolor_nt <- mycolors_locus[1:3]
names(mycolor_nt) <- nodetype
scales::show_col(mycolor_nt)

#############
## workdir ##
#############

plotdir = './visualize_CCC/'
dir.create(plotdir,recursive = T)

## selecting ligands receptors TFs TGs of interest
key_lig <- unique(c("FN1", "EGF", "L1CAM", "IGF1", "FGF2",
                    "TGFB1","ANGPTL4","PLAU","PTN","VEGFA","GNAI2"))

key_rec <- unique(c("EGFR", "ERBB3", "ERBB4", "IGF1R", "ITGAV", "SDC2", "FGFR2",
                    "IL17RC", "TGFBR3","SDC4","PLAUR"))

key_tf <- unique(c("CREB1", "NFKB1", "RELA", "EP300","NME2", "AR","ERF1",
                   "CREBBP","APBB2","ARHGAP35","E2F1"))

key_tg <- unique(c("CCK", "CCL4", "GFAP", "CD74", "FAM107A", "PTMA","S100A10",
            "SNAP25", "RTN1","APP","THY1","YWHAN","EIF4A2","FBLN1",
            "TSPAN7","APOD","BEX1","HLA−DPA1","NEFL","NPTXR","TMEM59L",
            "TUBB2A","PCSK1N"))

mlnet <- MLnet
mlnet$LigRec <- mlnet$LigRec[mlnet$LigRec$source %in% key_lig,]
mlnet$LigRec <- mlnet$LigRec[mlnet$LigRec$target %in% key_rec,]
mlnet$RecTF <- mlnet$RecTF[mlnet$RecTF$source %in% mlnet$LigRec$target,]
mlnet$RecTF <- mlnet$RecTF[mlnet$RecTF$target %in% key_tf,]
mlnet$TFTar <- mlnet$TFTar[mlnet$TFTar$source %in% mlnet$RecTF$target,]

mlnet_check <- mlnet
mlnet_check$TFTar <- mlnet_check$TFTar[mlnet_check$TFTar$target %in% key_tg,]
mlnet_check$RecTF <- mlnet_check$RecTF[mlnet_check$RecTF$target %in% mlnet_check$TFTar$source,]
mlnet_check$LigRec <- mlnet_check$LigRec[mlnet_check$LigRec$target %in% mlnet_check$RecTF$source,]

colodb = pal_locuszoom(palette = "default", alpha = 0.5)(4)
names(colodb) <- c("Ligand", "Receptor","TF", "Target")
scales::show_col(colodb)

downstream <- 'Target'
gtitle <- paste0('TAM-TC_IGF1_IGF1R_select')
drawMLnetworkPlot_V5(mlnet=mlnet_check,downstream=downstream,
                     colodb=colodb,gtitle=gtitle,wd=plotdir,
                     p_height = 6,p_width = 10)

## Reference
## L. Yan, J. Cheng, Q. Nie, X. Sun, Dissecting multilayer cell-cell communications 
## with signaling feedback loops from spatial transcriptomics data. 
## Genome Research, gr.279857.124 [Online ahead of print] (2025). https://doi.org/10.1101/gr.279857.124