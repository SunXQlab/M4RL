###################
#  main function  #
###################

mainfunc <- function(LigClu, Rebclu, workdir, RecTF.method, TFTG.method)
{
  
  ## LigRec
  
  source_abundant <- Ligs_expr_list[[LigClu]]
  cat("source_background:",length(source_abundant),"\n")
  
  target_abundant <- Recs_expr_list[[Rebclu]]
  cat("target_background:",length(target_abundant),"\n")

  tryCatch({
    LigRecTab <- getLigRec(LigRec.DB, source_abundant, target_abundant)
  }, error = function(e){
    cat(conditionMessage(e),"\n")
  })
  tag1 = exists("LigRecTab")
  if(!tag1) LigRecTab = data.frame()

  ## TFTG
  
  target_gene <- rownames(GCMat)[rowMeans(GCMat)>0]
  cat("target_background:",length(target_gene),"\n")
  target_icg <- ICGs_list[[Rebclu]]
  
  if(LigClu %in% names(target_icg)){
    target_icg <- target_icg[[LigClu]]
  }
  cat("target_icg:",length(target_icg),"\n")
  

  tryCatch({
    TFTGTab <- getTFTG(TFTG.DB, target_icg, target_gene, TFTG.method)
  },error = function(e){
    cat(conditionMessage(e),"\n")
  })
  tag2 = exists("TFTGTab")
  if(!tag2) TFTGTab = data.frame()

  ## RecTF
  
  if(tag1 & tag2){
    Rec.list <- getNodeList(LigRecTab, "target")
    TF.list <- getNodeList(TFTGTab, "source")
    target.tfs <- TF.list

    tryCatch({
      RecTFTab <- getRecTF(RecTF.DB, Rec.list, TF.list, RecTF.method)
    },error=function(e){
      cat(conditionMessage(e),"\n")
    })
  }

  tag3 = exists("RecTFTab")
  if(!tag3) RecTFTab = data.frame()
  
  ## merge
  
  if(tag3 & tag1){
    Receptors_in_Tab <- getNodeList(RecTFTab, "source")
    LigRecTab_new <- LigRecTab[LigRecTab[,2] %in% Receptors_in_Tab,]
  }else{
    LigRecTab_new <- data.frame()
  }

  cat("LR pairs:",nrow(LigRecTab_new),"\n")

  if(tag3 & tag2){
    TFs_in_Tab <- getNodeList(RecTFTab, "target")
    TFTGTab_new <- TFTGTab[TFTGTab[,1] %in% TFs_in_Tab,]
  }else{
    TFTGTab_new = data.frame()
  }

  cat("TFTG pairs:",nrow(TFTGTab_new),"\n")
  
  ## result
  
  result <- list("LigRec" = LigRecTab_new,
                 "RecTF" = RecTFTab,
                 "TFTar" = TFTGTab_new)
  
  foldername = paste(LigClu,RecClu,sep = "_")
  workdir = paste(workdir,foldername,sep = "/")
  if(!dir.exists(workdir)) dir.create(workdir,recursive = TRUE)
  saveRDS(result, file = paste0(workdir,"/scMLnet.rds"))
  
  return(result)
}

# main function with calculate pval
mainfunc_with_pval <- function(LigClu, Rebclu, workdir, RecTF.method, TFTG.method)
{
  
  ## LigRec
  
  source_abundant <- Ligs_expr_list[[LigClu]]
  cat("source_background:",length(source_abundant),"\n")
  target_abundant <- Recs_expr_list[[Rebclu]]
  cat("target_background:",length(target_abundant),"\n")
  
  tryCatch({
    LigRecTab <- getLigRec(LigRec.DB, source_abundant, target_abundant)
  }, error = function(e){
    cat(conditionMessage(e),"\n")
  })
  tag1 = exists("LigRecTab")
  if(!tag1) LigRecTab = data.frame()
  
  ## TFTG
  
  target_gene <- rownames(GCMat)[rowMeans(GCMat)>0]
  cat("target_background:",length(target_gene),"\n")
  target_icg <- ICGs_list[[Rebclu]]
  #target_icg <- target_gene
  if(LigClu %in% names(target_icg)){
    target_icg <- target_icg[[LigClu]]
  }
  cat("target_icg:",length(target_icg),"\n")
  
  tryCatch({
    TFTGTab <- getTFTG(TFTG.DB, target_icg, target_gene, TFTG.method)
  },error = function(e){
    cat(conditionMessage(e),"\n")
  })
  tag2 = exists("TFTGTab")
  if(!tag2) TFTGTab = data.frame()
  
  ## RecTF
  
  if(tag1 & tag2){
    Rec.list <- getNodeList(LigRecTab, "target")
    TF.list <- getNodeList(TFTGTab, "source")
    target.tfs <- TF.list
    tryCatch({
      RecTFTab <- getRecTF(RecTF.DB, Rec.list, TF.list, RecTF.method)
    },error=function(e){
      cat(conditionMessage(e),"\n")
    })
  }
  
  tag3 = exists("RecTFTab")
  if(!tag3) RecTFTab = data.frame()
  
  ## merge
  
  if(tag3 & tag1){
    Receptors_in_Tab <- getNodeList(RecTFTab, "source")
    LigRecTab_new <- LigRecTab[LigRecTab[,2] %in% Receptors_in_Tab,]
  }else{
    LigRecTab_new <- data.frame()
  }
  cat("LR pairs:",nrow(LigRecTab_new),"\n")
  
  if(tag3 & tag2){
    TFs_in_Tab <- getNodeList(RecTFTab, "target")
    TFTGTab_new <- TFTGTab[TFTGTab[,1] %in% TFs_in_Tab,]
  }else{
    TFTGTab_new = data.frame()
  }
  cat("TFTG pairs:",nrow(TFTGTab_new),"\n")
  
  ## result
  
  result <- list("LigRec" = LigRecTab_new,
                 "RecTF" = RecTFTab,
                 "TFTar" = TFTGTab_new)
  
  foldername = paste(LigClu,RecClu,sep = "_")
  workdir = paste(workdir,foldername,sep = "/")
  if(!dir.exists(workdir)) dir.create(workdir,recursive = TRUE)
  saveRDS(result, file = paste0(workdir,"/scMLnet.rds"))
  
  # calculate P_value
  if (tag2 & RecTF.method == 'Fisher'){
    if(tag1 & tag2){
      Rec.list <- getNodeList(LigRecTab, "target")
      TF.list <- getNodeList(TFTGTab, "source")
      target.tfs <- TF.list
      tryCatch({
        RecTFpval <- getRecpval(RecTF.DB, Rec.list, TF.list, RecTF.method)
      },error=function(e){
        cat(conditionMessage(e),"\n")
      })
    }
    
    if (length(LigRecTab)==0){
      RecTFpval = data.frame()
    }else{
      RecTFpval <- as.data.frame(RecTFpval) %>% rownames_to_column(.,var = "target")
      LigRecTab$pval <- NA
      for (rec in LigRecTab$target){
        pos1 <- which(LigRecTab$target %in% rec)
        pos2 <- which(RecTFpval$target %in% rec)
        LigRecTab$pval[pos1] <- RecTFpval$RecTFpval[pos2]
      }
    }
    saveRDS(LigRecTab, file = paste0(workdir,"/cellpair_LRI_pval.rds")) 
  }else{
    
    LigRecTab = data.frame()
    saveRDS(LigRecTab, file = paste0(workdir,"/cellpair_LRI_pval.rds")) 
  }
  
  return(result)
}

# get nodes in prior Database
getNodeList <- function(Database, Nodetype)
{
  
  NodeList <- Database[,Nodetype] %>% unlist() %>% unique()
  
  return(NodeList)
  
}

getRecpval <- function(RecTF.DB, Rec.list, TF.list, method = c('Fisher','Search'))
{
  
  if(method=='Search'){
    RecTFTable <- getRecTFSearch(RecTF.DB, Rec.list, TF.list)
  }else if(method=='Fisher'){
    RecTFpval <- getRecFisherpval(RecTF.DB, Rec.list, TF.list)
  }
  
  return(RecTFpval)
}

getRecFisherpval <- function(RecTF.DB, Rec.list, TF.list)
{
  
  if (!is.data.frame(RecTF.DB))
    stop("RecTF.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(RecTF.DB))
    stop("RecTF.DB must contain a column named 'source'")
  if (!"target" %in% colnames(RecTF.DB))
    stop("RecTF.DB must contain a column named 'target'")
  
  # make sure Rec.list in RecTF.DB
  Rec.list <- Rec.list[Rec.list %in% RecTF.DB$source]
  Rec.list <- as.vector(Rec.list)
  
  # make sure TF.list in RecTF.DB
  TF.list <- TF.list[TF.list %in% RecTF.DB$target]
  TF.list <- as.vector(TF.list)
  
  # get TF activated by Receptors
  TFofRec <- lapply(Rec.list, function(x){
    RecTF.DB %>% dplyr::filter(source == x)  %>% dplyr::select(target) %>% unlist() %>% unique()
  })
  names(TFofRec) <- Rec.list
  
  # get all TF
  TFofALL <- RecTF.DB %>% dplyr::select(target) %>% unlist() %>% unique()
  
  # perform fisher test
  Recs <- lapply(TFofRec, function(x){
    fisher_test(subset1 = x, subset2 = TF.list, backgrond = TFofALL)
  })
  
  pval <- unlist(Recs)
  return(pval)
  # return(Recs)
}

getLigRec <- function(LigRec.DB, source_up, target_up)
{
  
  require(dplyr)
  
  if (!is.data.frame(LigRec.DB)) 
    stop("LigRec.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(LigRec.DB)) 
    stop("LigRec.DB must contain a column named 'source'")
  if (!"target" %in% colnames(LigRec.DB)) 
    stop("LigRec.DB must contain a column named 'target'")
  
  # get ligand and receptor list
  LigGene <- LigRec.DB %>% dplyr::select(source) %>% unlist() %>% unique()
  RecGene <- LigRec.DB %>% dplyr::select(target) %>% unlist() %>% unique()
  TotLigRec <- paste(LigRec.DB$source, LigRec.DB$target, sep = "_") %>% unique()

  # get high expressed ligand and receptor
  LigHighGene <- intersect(LigGene,source_up)
  RecHighGene <- intersect(RecGene,target_up)

  # get activated LR pairs
  LRList <- paste(rep(LigHighGene,each = length(RecHighGene)),RecHighGene,sep = "_")
  LRList <- intersect(LRList,TotLigRec)
  
  # check result
  if(length(LRList)==0)
    stop("Error: No significant LigRec pairs")
  
  # get result
  LRTable <- LRList %>% strsplit(.,split = "_") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(LRTable) <- c("source","target")

  cat(paste0("get ",length(LRList)," activated LR pairs\n"))
  return(LRTable)
  
}

getTFTG <- function(TFTG.DB, target.degs, target.genes, method = c('Fisher','Search'))
{
  
  if(method=='Search'){
    TFTGTable <- getTFTGSearch(TFTG.DB, target.degs)
  }else if(method=='Fisher'){
    TFTGTable <- getTFTGFisher(TFTG.DB, target.degs, target.genes)
  }
  
  return(TFTGTable)
}

getTFTGSearch <- function(TFTG.DB, target.degs)
{
  
  require(dplyr)
  
  if (!is.data.frame(TFTG.DB)) 
    stop("TFTG.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(TFTG.DB)) 
    stop("TFTG.DB must contain a column named 'source'")
  if (!"target" %in% colnames(TFTG.DB)) 
    stop("TFTG.DB must contain a column named 'target'")
  
  # get TF and Target list
  TFGene <- TFTG.DB %>% dplyr::select(source) %>% unlist() %>% unique()
  TargetGene <- TFTG.DB %>% dplyr::select(target) %>% unlist() %>% unique()
  TotTFTG <- paste(TFTG.DB$source, TFTG.DB$target, sep = "_") %>% unique()
  
  # get activated target
  TargetHighGene <- intersect(TargetGene,target.degs)
  
  # get potential activated TFTG pairs
  TFTGList <- paste(TFGene, TargetHighGene,sep = "_")
  TFTGList <- intersect(TFTGList,TotTFTG)
  
  # check result
  if(length(TFTGList)==0)
    stop("Error: No significant TFTG pairs")
  
  # get result
  TFTGTable <- TFTGList %>% strsplit(.,split = "_") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(TFTGTable) <- c("source","target")
  
  cat(paste0("get ",length(TFTGList)," activated TFTG pairs\n"))
  return(TFTGTable)
  
}

getTFTGFisher <- function(TFTG.DB, target.degs, target.genes)
{
  
  if (!is.data.frame(TFTG.DB))
    stop("TFTG.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(TFTG.DB))
    stop("TFTG.DB must contain a column named 'source'")
  if (!"target" %in% colnames(TFTG.DB))
    stop("TFTG.DB must contain a column named 'target'")
  
  # get TF list
  TF.list <- TFTG.DB %>% dplyr::select(source) %>% unlist() %>% unique()
  
  # get Target list
  TG.list <- lapply(TF.list, function(x){
    TFTG.DB %>% filter(source == x)  %>% dplyr::select(target) %>% unlist() %>% unique() %>% intersect(.,target.genes)
  })
  names(TG.list) <- TF.list
  
  # get target differently expressed genes
  DEGs <- target.degs
  
  # perform fisher test
  TFs <- lapply(TG.list, function(x){fisher_test(subset1 = x, subset2 = DEGs, backgrond = target.genes)})
  TFs <- unlist(TFs)
  TFs <- names(TFs)[TFs <= 0.05]
  TFs <- TFs[TFs %in% target.genes]
  
  # get activated LR pairs
  TFTGList <- TG.list[TFs]
  TFTGList <- lapply(TFTGList, function(x){intersect(x, DEGs)})
  TFTGList <- paste(rep(TFs, times = lengths(TFTGList)), unlist(TFTGList), sep = "_")
  
  # check result
  if(length(TFTGList)==0)
    stop("Error: No significant TFTG pairs")
  
  # get result
  TFTGTable <- TFTGList %>% strsplit(.,split = "_") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(TFTGTable) <- c("source","target")
  
  cat(paste0("get ",length(TFTGList)," activated TFTG pairs\n"))
  return(TFTGTable)
  
}

# perform fisher test to get activate TFs or pathways
fisher_test <- function(subset1,subset2,backgrond)
{
  a=length(intersect(subset1,subset2))
  b=length(subset1)-a
  c=length(subset2)-a
  d=length(backgrond)-a-b-c
  matrix=matrix(c(a,c,b,d),nrow=2)
  fisher.test(matrix,alternative="greater")$p.value
}

getRecTF <- function(RecTF.DB, Rec.list, TF.list, method = c('Fisher','Search'))
{
  
  if(method=='Search'){
    RecTFTable <- getRecTFSearch(RecTF.DB, Rec.list, TF.list)
  }else if(method=='Fisher'){
    RecTFTable <- getRecTFFisher(RecTF.DB, Rec.list, TF.list)
  }
  
  return(RecTFTable)
}
getRecTFFisher <- function(RecTF.DB, Rec.list, TF.list)
{
  
  if (!is.data.frame(RecTF.DB))
    stop("RecTF.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(RecTF.DB))
    stop("RecTF.DB must contain a column named 'source'")
  if (!"target" %in% colnames(RecTF.DB))
    stop("RecTF.DB must contain a column named 'target'")
  
  # make sure Rec.list in RecTF.DB
  Rec.list <- Rec.list[Rec.list %in% RecTF.DB$source]
  Rec.list <- as.vector(Rec.list)
  
  # make sure TF.list in RecTF.DB
  TF.list <- TF.list[TF.list %in% RecTF.DB$target]
  TF.list <- as.vector(TF.list)
  
  # get TF activated by Receptors
  TFofRec <- lapply(Rec.list, function(x){
    RecTF.DB %>% dplyr::filter(source == x)  %>% dplyr::select(target) %>% unlist() %>% unique()
  })
  names(TFofRec) <- Rec.list
  
  # get all TF
  TFofALL <- RecTF.DB %>% dplyr::select(target) %>% unlist() %>% unique()
  
  # perform fisher test
  Recs <- lapply(TFofRec, function(x){
    fisher_test(subset1 = x, subset2 = TF.list, backgrond = TFofALL)
  })
  
  Recs <- unlist(Recs)
  Recs <- names(Recs)[Recs <= 0.05]
  # Recs <- Recs[Recs %in% target_gene]
  
  # get activated RecTF pairs
  RecTFList <- TFofRec[Recs]
  RecTFList <- lapply(RecTFList, function(x){intersect(x, TF.list)})
  RecTFList <- paste(rep(Recs, times = lengths(RecTFList)), unlist(RecTFList), sep = "_")
  
  # check result
  if(length(RecTFList)==0)
    stop("Error: No significant RecTF pairs")
  
  # get result
  RecTFTable <- RecTFList %>% strsplit(.,split = "_") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(RecTFTable) <- c("source","target")
  
  cat(paste0("get ",length(RecTFList)," activated RecTF pairs\n"))
  return(RecTFTable)
  # return(Recs)
}

getRecTFSearch <- function(RecTF.DB, Rec.list, TF.list)
{
  
  if (!is.data.frame(RecTF.DB))
    stop("RecTF.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(RecTF.DB))
    stop("RecTF.DB must contain a column named 'source'")
  if (!"target" %in% colnames(RecTF.DB))
    stop("RecTF.DB must contain a column named 'target'")
  
  # make sure Rec.list in RecTF.DB
  Rec.list <- Rec.list[Rec.list %in% RecTF.DB$source]
  Rec.list <- as.vector(Rec.list)
  
  # make sure TF.list in RecTF.DB
  TF.list <- TF.list[TF.list %in% RecTF.DB$target]
  TF.list <- as.vector(TF.list)
  
  # RecTF pairs in RecTF.DB
  TotRecTF <- paste(RecTF.DB$source, RecTF.DB$target, sep = "_") %>% unique()
  
  # get activated LR pairs
  RecTFList <- paste(rep(Rec.list,each = length(TF.list)),TF.list,sep = "_")
  RecTFList <- intersect(RecTFList,TotRecTF)
  
  # check result
  if(length(RecTFList)==0)
    stop("Error: No significant RecTF pairs")
  
  # get result
  RecTFTable <- RecTFList %>% strsplit(.,split = "_") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(RecTFTable) <- c("source","target")
  
  cat(paste0("get ",length(RecTFList)," activated RecTF pairs\n"))
  return(RecTFTable)
}
################
## getLRscore ##
################

## 对表达谱进行填补
run_Imputation <- function(exprMat, use.seed = TRUE, seed = 2021)
{
  
  expr.Impute <- CreateSeuratObject(exprMat,verbose=F)
  if(use.seed) set.seed(seed)
  message('Using imputation method ALRA wrapped in Seurat')
  expr.Impute <- RunALRA(expr.Impute)
  exprMat.Impute <- expr.Impute@assays$alra@data
  
  return(exprMat.Impute)
  
}

## 获取指定sender和Receiver的LRscore和TG表达量

calculate_LRTG_score_V2 <- function(exprMat, distMat, annoMat, group = NULL,
                                    LRpairs, TGs, Receiver, Sender = NULL, 
                                    far.ct = 0.75, close.ct = 0.25, 
                                    downsample = FALSE)
{
  
  # 获取Sender和Receiver
  receBars = annoMat %>% dplyr::filter(Cluster == Receiver) %>% 
    dplyr::select(Barcode) %>% unlist() %>% as.character()
  if(is.character(Sender)){
    sendBars = annoMat %>% dplyr::filter(Cluster == Sender) %>% 
      dplyr::select(Barcode) %>% unlist() %>% as.character()
  }else{
    sendBars = annoMat %>% dplyr::filter(Cluster != Receiver) %>% 
      dplyr::select(Barcode) %>% unlist() %>% as.character()
  }
  
  # 获取Lig和Rec列表
  Receptors = lapply(LRpairs, function(lr){stringr::str_split(lr,"_",simplify = T)[,2]})
  Ligands = lapply(LRpairs, function(lr){stringr::str_split(lr,"_",simplify = T)[,1]})
  
  # get exprMat of Ligand
  LigMats = lapply(TGs, function(tg){
    # print(tg)
    ligands = Ligands[[tg]]
    if(length(ligands)==1){
      lig_count = exprMat[ligands, sendBars]
      lig_count = matrix(lig_count,nrow = 1)
    }else{
      lig_count = exprMat[ligands, sendBars] %>% as.matrix()
    }
    rownames(lig_count) = LRpairs[[tg]]
    colnames(lig_count) = sendBars
    lig_count
  })
  names(LigMats) = TGs
  
  # get exprMat of Receptor
  RecMats = lapply(TGs, function(tg){
    receptors = Receptors[[tg]]
    if(length(receptors)==1){
      rec_count = exprMat[receptors, receBars]
      rec_count = matrix(rec_count,nrow = 1)
    }else{
      rec_count = exprMat[receptors, receBars] %>% as.matrix()
    }
    rownames(rec_count) = LRpairs[[tg]]
    colnames(rec_count) = receBars
    rec_count
  })
  names(RecMats) = TGs
  
  # 获取距离矩阵
  distMat = distMat[sendBars, receBars]
  
  # 获取细胞配对
  if(!is.null(group)){
    cpMat <- get_cell_pairs(group, distMat, far.ct, close.ct)
  }else{
    cpMat <- NULL
  }
  distMat = 1/distMat
  
  # 获取LR打分矩阵
  t1 <- Sys.time(); message(paste0('Start at: ',as.character(t1)))
  LRs_score = lapply(TGs, function(tg){
    # 对每个TG计算LRscore
    # print(tg)
    LigMat = LigMats[[tg]]
    RecMat = RecMats[[tg]]
    lr = LRpairs[[tg]] 
    
    if(is.null(cpMat)){
      
      LR_score = RecMat*(LigMat%*%distMat)
      LR_score = t(LR_score)
      colnames(LR_score) = lr
      rownames(LR_score) = receBars
      
    }else{
      
      LR_score = lapply(unique(cpMat$Receiver), function(j){
        # j = unique(cpMat$Receiver)[1]
        is <- cpMat$Sender[cpMat$Receiver == j] %>% unique()
        if(length(is)==1){
          RecMat[,j]*(LigMat[,is]*distMat[is,j])
        }else{
          RecMat[,j]*(LigMat[,is]%*%distMat[is,j])
        }
      }) %>% do.call('cbind',.) %>% t()
      colnames(LR_score) = lr
      rownames(LR_score) = unique(cpMat$Receiver)
      
    }
    LR_score
    
  })
  names(LRs_score) = TGs
  t2 <- Sys.time(); message(paste0('End at: ',as.character(t2)))
  t2-t1 # 17 sec
  
  # 获取Target表达谱
  if(is.null(cpMat)){
    
    TGs_expr = lapply(TGs, function(tg){ exprMat[tg, receBars] })
    
  }else{
    
    TGs_expr = lapply(TGs, function(tg){ exprMat[tg, unique(cpMat$Receiver)] })
    
  }
  names(TGs_expr) = TGs
  
  # downsample
  if(length(receBars)>500){
    if(isTRUE(downsample)){
      set.seed(2021)
      
      if(is.null(cpMat)){
        keep_cell = sample(receBars, size = 500, replace = F)
      }else{
        keep_cell = sample(unique(cpMat$Receiver), size = 500, replace = F)
      }
      
      LRs_score = lapply(LRs_score, function(LR_score){ LR_score = LR_score[keep_cell,] })
      TGs_expr = lapply(TGs_expr, function(TG_count){ TG_count = TG_count[keep_cell] })
    }
  }
  
  # LR打分矩阵~靶基因表达量
  LRTG_score = list(LRs_score = LRs_score, TGs_expr = TGs_expr)
  
  return(LRTG_score)
  
}

calculate_LRTG_score_V3 <- function(exprMat, distMat, annoMat, group = NULL,
                                    LRpairs, TGs, Receiver, Sender = NULL, 
                                    far.ct = 0.75, close.ct = 0.25, 
                                    downsample = FALSE)
{
  
  # get Sender and Receiver
  receBars = annoMat %>% dplyr::filter(Cluster == Receiver) %>% 
    dplyr::select(Barcode) %>% unlist() %>% as.character()
  if(is.character(Sender)){
    sendBars = annoMat %>% dplyr::filter(Cluster == Sender) %>% 
      dplyr::select(Barcode) %>% unlist() %>% as.character()
  }else{
    sendBars = annoMat %>% dplyr::filter(Cluster != Receiver) %>% 
      dplyr::select(Barcode) %>% unlist() %>% as.character()
  }
  
  # get Ligand and Receptor
  Receptors = lapply(LRpairs, function(lr){stringr::str_split(lr,"_",simplify = T)[,2]})
  Ligands = lapply(LRpairs, function(lr){stringr::str_split(lr,"_",simplify = T)[,1]})
  
  # get exprMat of Ligand
  LigMats = lapply(TGs, function(tg){
    # print(tg)
    ligands = Ligands[[tg]]
    if(length(ligands)==1){
      lig_count = exprMat[ligands, sendBars]
      lig_count = matrix(lig_count,nrow = 1)
    }else{
      lig_count = exprMat[ligands, sendBars] %>% as.matrix()
    }
    rownames(lig_count) = LRpairs[[tg]]
    colnames(lig_count) = sendBars
    lig_count
  })
  names(LigMats) = TGs
  
  # get exprMat of Receptor
  RecMats = lapply(TGs, function(tg){
    receptors = Receptors[[tg]]
    if(length(receptors)==1){
      rec_count = exprMat[receptors, receBars]
      rec_count = matrix(rec_count,nrow = 1)
    }else{
      rec_count = exprMat[receptors, receBars] %>% as.matrix()
    }
    rownames(rec_count) = LRpairs[[tg]]
    colnames(rec_count) = receBars
    rec_count
  })
  names(RecMats) = TGs
  
  # get and transform distMat
  l.scale = 100
  distMat = distMat[sendBars, receBars]
  
  # get cell pair
  if(!is.null(group)){
    cpMat <- get_cell_pairs(group, distMat, far.ct, close.ct)
  }else{
    cpMat <- NULL
  }
  distMat = exp(-distMat^2/(2*l.scale^2))
  
  # get LRscore
  t1 <- Sys.time(); message(paste0('Start at: ',as.character(t1)))
  LRs_score = lapply(TGs, function(tg){
    # print(tg)
    LigMat = LigMats[[tg]]
    RecMat = RecMats[[tg]]
    lr = LRpairs[[tg]] 
    
    if(is.null(cpMat)){
      
      LR_score = RecMat*(LigMat%*%distMat)
      LR_score = t(LR_score)
      colnames(LR_score) = lr
      rownames(LR_score) = receBars
      
    }else{
      
      LR_score = lapply(unique(cpMat$Receiver), function(j){
        # j = 'GATATTTCCTACATGG.1'
        is <- cpMat$Sender[cpMat$Receiver == j] %>% unique()
        RecMat[,j]*(LigMat[,is]%*%distMat[is,j])
      }) %>% do.call('rbind',.) 
      colnames(LR_score) = lr
      rownames(LR_score) = unique(cpMat$Receiver)
      
    }
    LR_score
    
  })
  names(LRs_score) = TGs
  t2 <- Sys.time(); message(paste0('End at: ',as.character(t2)))
  t2-t1 # 17 sec
  
  # get expression of each TG
  if(is.null(cpMat)){
    
    TGs_expr = lapply(TGs, function(tg){ exprMat[tg, receBars] })
    
  }else{
    
    TGs_expr = lapply(TGs, function(tg){ exprMat[tg, unique(cpMat$Receiver)] })
    
  }
  names(TGs_expr) = TGs
  
  # downsample
  if(length(receBars)>500){
    if(isTRUE(downsample)){
      set.seed(2021)
      
      if(is.null(cpMat)){
        keep_cell = sample(receBars, size = 500, replace = F)
      }else{
        keep_cell = sample(unique(cpMat$Receiver), size = 500, replace = F)
      }
      
      LRs_score = lapply(LRs_score, function(LR_score){ LR_score = LR_score[keep_cell,] })
      TGs_expr = lapply(TGs_expr, function(TG_count){ TG_count = TG_count[keep_cell] })
    }
  }
  
  LRTG_score = list(LRs_score = LRs_score, TGs_expr = TGs_expr)
  
  return(LRTG_score)
  
}

calculate_LRTG_score_V4 <- function(exprMat, distMat, annoMat, group = NULL,
                                    LRpairs, TGs, Receiver, Sender = NULL, 
                                    far.ct = 0.75, close.ct = 0.25, 
                                    downsample = FALSE)
{
  
  # get Sender and Receiver
  receBars = annoMat %>% dplyr::filter(Cluster == Receiver) %>% 
    dplyr::select(Barcode) %>% unlist() %>% as.character()
  if(is.character(Sender)){
    sendBars = annoMat %>% dplyr::filter(Cluster == Sender) %>% 
      dplyr::select(Barcode) %>% unlist() %>% as.character()
  }else{
    sendBars = annoMat %>% dplyr::filter(Cluster != Receiver) %>% 
      dplyr::select(Barcode) %>% unlist() %>% as.character()
  }
  
  # get Ligand and Receptor
  Receptors = lapply(LRpairs, function(lr){stringr::str_split(lr,"_",simplify = T)[,2]})
  Ligands = lapply(LRpairs, function(lr){stringr::str_split(lr,"_",simplify = T)[,1]})
  
  # get exprMat of Ligand
  LigMats = lapply(TGs, function(tg){
    # print(tg)
    ligands = Ligands[[tg]]
    if(length(ligands)==1){
      lig_count = exprMat[ligands, sendBars]
      lig_count = matrix(lig_count,nrow = 1)
    }else{
      lig_count = exprMat[ligands, sendBars] %>% as.matrix()
    }
    rownames(lig_count) = LRpairs[[tg]]
    colnames(lig_count) = sendBars
    lig_count
  })
  names(LigMats) = TGs
  
  # get exprMat of Receptor
  RecMats = lapply(TGs, function(tg){
    receptors = Receptors[[tg]]
    if(length(receptors)==1){
      rec_count = exprMat[receptors, receBars]
      rec_count = matrix(rec_count,nrow = 1)
    }else{
      rec_count = exprMat[receptors, receBars] %>% as.matrix()
    }
    rownames(rec_count) = LRpairs[[tg]]
    colnames(rec_count) = receBars
    rec_count
  })
  names(RecMats) = TGs
  
  # get and transform distMat
  distMat = distMat[sendBars, receBars]
  
  # get cell pair
  if(!is.null(group)){
    cpMat <- get_cell_pairs(group, distMat, far.ct, close.ct)
  }else{
    cpMat <- NULL
  }
  distMat[sendBars, receBars] = 1
  
  # get LRscore
  t1 <- Sys.time(); message(paste0('Start at: ',as.character(t1)))
  LRs_score = lapply(TGs, function(tg){
    print(tg)
    LigMat = LigMats[[tg]]
    RecMat = RecMats[[tg]]
    lr = LRpairs[[tg]] 
    
    if(is.null(cpMat)){
      
      LR_score = RecMat*(LigMat%*%distMat)/length(sendBars)
      LR_score = t(LR_score)
      colnames(LR_score) = lr
      rownames(LR_score) = receBars
      
    }else{
      
      LR_score = lapply(unique(cpMat$Receiver), function(j){
        # j = 'GATCCCTTTATACTGC.1'
        is <- cpMat$Sender[cpMat$Receiver == j] %>% unique()
        RecMat[,j]*(LigMat[,is]%*%distMat[is,])/length(is)
      }) %>% do.call('rbind',.) 
      colnames(LR_score) = lr
      rownames(LR_score) = unique(cpMat$Receiver)
      
    }
    LR_score
    
  })
  names(LRs_score) = TGs
  t2 <- Sys.time(); message(paste0('End at: ',as.character(t2)))
  t2-t1 # 17 sec
  
  # get expression of each TG
  if(is.null(cpMat)){
    
    TGs_expr = lapply(TGs, function(tg){ exprMat[tg, receBars] })
    
  }else{
    
    TGs_expr = lapply(TGs, function(tg){ exprMat[tg, unique(cpMat$Receiver)] })
    
  }
  names(TGs_expr) = TGs
  
  # downsample
  if(length(receBars)>500){
    if(isTRUE(downsample)){
      set.seed(2021)
      
      if(is.null(cpMat)){
        keep_cell = sample(receBars, size = 500, replace = F)
      }else{
        keep_cell = sample(unique(cpMat$Receiver), size = 500, replace = F)
      }
      
      LRs_score = lapply(LRs_score, function(LR_score){ LR_score = LR_score[keep_cell,] })
      TGs_expr = lapply(TGs_expr, function(TG_count){ TG_count = TG_count[keep_cell] })
    }
  }
  
  LRTG_score = list(LRs_score = LRs_score, TGs_expr = TGs_expr)
  
  return(LRTG_score)
  
}


## 获取微环境细胞间的LRscore和TG表达量

calculate_LRTG_allscore_V2 <- function(exprMat, distMat, annoMat, group = NULL,
                                       mulNetList, Receiver, Sender = NULL, 
                                       far.ct = 0.75, close.ct = 0.25,
                                       downsample = FALSE){
  
  if(is.null(Sender)){
    
    mulNetList = mulNetList[grep(paste0("_",Receiver),names(mulNetList))]
    mulNet_tab = lapply(mulNetList, function(mlnet){
      
      ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
      rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
      tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
      
      res = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
        merge(., tftg, by = 'TF') %>% 
        dplyr::select(Ligand, Receptor, TF, Target) %>% 
        arrange(Ligand, Receptor)
      
    })
    mulNet_tab = do.call("rbind", mulNet_tab)
    
    LRpairs = by(mulNet_tab, as.character(mulNet_tab$Target), function(x){paste(x$Ligand, x$Receptor, sep = "_")})
    LRpairs = lapply(LRpairs, function(lrtg){lrtg[!duplicated(lrtg)]})
    TGs = names(LRpairs)
    
    # cat("calculate the regulatory score of LR pairs from microenvironment")
    cat(paste0("calculate the regulatory score of LR pairs from microenvironment to ",Receiver))
    LRTG_allscore = calculate_LRTG_score_V2(exprMat, distMat, annoMat, group,
                                            LRpairs, TGs, Receiver, Sender, 
                                            far.ct, close.ct, downsample)
    
  }else if(length(Sender)==1){
    
    cellpair = paste(Sender,Receiver,sep = "_")
    mulNet = mulNetList[[cellpair]]
    
    TGs = mulNet %>% .[['TFTar']] %>% dplyr::select(target) %>% unlist() %>% as.character() %>% unique()
    LRpairs = get_LRTG_link(mulNet, TGs)
    
    cat(paste0("calculate the regulatory score of LR pairs from ",Sender,' to ',Receiver))
    LRTG_allscore = calculate_LRTG_score_V2(exprMat, distMat, annoMat, group,
                                            LRpairs, TGs, Receiver, Sender, 
                                            far.ct, close.ct, downsample)
    
  }else{
    
    cellpair = paste(Sender,Receiver,sep = "_")
    cellpair <- intersect(names(mulNetList),cellpair)
    if(length(cellpair)==0){
      LRTG_allscore = NA
      return(LRTG_allscore)
    }
    
    mulNetlist = mulNetList[cellpair]
    mulNet_tab = lapply(mulNetlist, function(mlnet){
      
      ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
      rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
      tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
      
      res = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
        merge(., tftg, by = 'TF') %>% 
        dplyr::select(Ligand, Receptor, TF, Target) %>% 
        arrange(Ligand, Receptor)
      
    })
    
    LRpairs_TGs_list <- lapply(mulNet_tab, function(ml_tab){
      
      lrpairs = by(ml_tab, as.character(ml_tab$Target), function(x){
        paste(x$Ligand, x$Receptor, sep = "_")
      })
      lrpairs = lapply(lrpairs, function(lrtg){lrtg[!duplicated(lrtg)]})
      tgs = names(lrpairs)
      
      list(LRpairs = lrpairs, TGs = tgs)
      
    })
    names(LRpairs_TGs_list) <- names(mulNet_tab)
    
    LRTG_allscore <- list()
    for (cp in cellpair) {
      
      receiver <- gsub('.*_','',cp)
      sender <- gsub('_.*','',cp)
      
      LRpairs <- LRpairs_TGs_list[[cp]]$LRpairs
      TGs <- LRpairs_TGs_list[[cp]]$TGs
      
      cat(paste0("calculate the regulatory score of LR pairs from ",sender,' to ',receiver))
      LRTG_allscore[[cp]] = calculate_LRTG_score_V2(exprMat, distMat, annoMat, group,
                                                    LRpairs, TGs, receiver, sender, 
                                                    far.ct, close.ct, downsample)
      
    }
  }
  
  return(LRTG_allscore)
}

calculate_LRTG_allscore_V3 <- function(exprMat, distMat, annoMat, group = NULL,
                                       mulNetList, Receiver, Sender = NULL, 
                                       far.ct = 0.75, close.ct = 0.25,
                                       downsample = FALSE){
  
  if(is.null(Sender)){
    
    mulNetList = mulNetList[grep(paste0("_",Receiver),names(mulNetList))]
    mulNet_tab = lapply(mulNetList, function(mlnet){
      
      ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
      rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
      tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
      
      res = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
        merge(., tftg, by = 'TF') %>% 
        dplyr::select(Ligand, Receptor, TF, Target) %>% 
        arrange(Ligand, Receptor)
      
    })
    mulNet_tab = do.call("rbind", mulNet_tab)
    
    LRpairs = by(mulNet_tab, as.character(mulNet_tab$Target), function(x){paste(x$Ligand, x$Receptor, sep = "_")})
    LRpairs = lapply(LRpairs, function(lrtg){lrtg[!duplicated(lrtg)]})
    TGs = names(LRpairs)
    
    # cat("calculate the regulatory score of LR pairs from microenvironment")
    cat(paste0("calculate the regulatory score of LR pairs from microenvironment to ",Receiver))
    LRTG_allscore = calculate_LRTG_score_V3(exprMat, distMat, annoMat, group,
                                            LRpairs, TGs, Receiver, Sender, 
                                            far.ct, close.ct, downsample)
    
  }else if(length(Sender)==1){
    
    cellpair = paste(Sender,Receiver,sep = "_")
    mulNet = mulNetList[[cellpair]]
    
    TGs = mulNet %>% .[['TFTar']] %>% dplyr::select(target) %>% unlist() %>% as.character() %>% unique()
    LRpairs = get_LRTG_link(mulNet, TGs)
    
    cat(paste0("calculate the regulatory score of LR pairs from ",Sender,' to ',Receiver))
    LRTG_allscore = calculate_LRTG_score_V3(exprMat, distMat, annoMat, group,
                                            LRpairs, TGs, Receiver, Sender, 
                                            far.ct, close.ct, downsample)
    
  }else{
    
    cellpair = paste(Sender,Receiver,sep = "_")
    cellpair <- intersect(names(mulNetList),cellpair)
    if(length(cellpair)==0){
      LRTG_allscore = NA
      return(LRTG_allscore)
    }
    
    mulNetlist = mulNetList[cellpair]
    mulNet_tab = lapply(mulNetlist, function(mlnet){
      
      ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
      rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
      tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
      
      res = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
        merge(., tftg, by = 'TF') %>% 
        dplyr::select(Ligand, Receptor, TF, Target) %>% 
        arrange(Ligand, Receptor)
      
    })
    
    LRpairs_TGs_list <- lapply(mulNet_tab, function(ml_tab){
      
      lrpairs = by(ml_tab, as.character(ml_tab$Target), function(x){
        paste(x$Ligand, x$Receptor, sep = "_")
      })
      lrpairs = lapply(lrpairs, function(lrtg){lrtg[!duplicated(lrtg)]})
      tgs = names(lrpairs)
      
      list(LRpairs = lrpairs, TGs = tgs)
      
    })
    names(LRpairs_TGs_list) <- names(mulNet_tab)
    
    LRTG_allscore <- list()
    for (cp in cellpair) {
      
      receiver <- gsub('.*_','',cp)
      sender <- gsub('_.*','',cp)
      
      LRpairs <- LRpairs_TGs_list[[cp]]$LRpairs
      TGs <- LRpairs_TGs_list[[cp]]$TGs
      
      cat(paste0("calculate the regulatory score of LR pairs from ",sender,' to ',receiver))
      LRTG_allscore[[cp]] = calculate_LRTG_score_V3(exprMat, distMat, annoMat, group,
                                                    LRpairs, TGs, receiver, sender, 
                                                    far.ct, close.ct, downsample)
      
    }
  }
  
  return(LRTG_allscore)
}

calculate_LRTG_allscore_V4 <- function(exprMat, distMat, annoMat, group = NULL,
                                       mulNetList, Receiver, Sender = NULL, 
                                       far.ct = 0.75, close.ct = 0.25,
                                       downsample = FALSE){
  
  if(is.null(Sender)){
    
    mulNetList = mulNetList[grep(paste0("_",Receiver),names(mulNetList))]
    mulNet_tab = lapply(mulNetList, function(mlnet){
      
      ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
      rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
      tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
      
      res = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
        merge(., tftg, by = 'TF') %>% 
        dplyr::select(Ligand, Receptor, TF, Target) %>% 
        arrange(Ligand, Receptor)
      
    })
    mulNet_tab = do.call("rbind", mulNet_tab)
    
    LRpairs = by(mulNet_tab, as.character(mulNet_tab$Target), function(x){paste(x$Ligand, x$Receptor, sep = "_")})
    LRpairs = lapply(LRpairs, function(lrtg){lrtg[!duplicated(lrtg)]})
    TGs = names(LRpairs)
    
    # cat("calculate the regulatory score of LR pairs from microenvironment")
    cat(paste0("calculate the regulatory score of LR pairs from microenvironment to ",Receiver))
    LRTG_allscore = calculate_LRTG_score_V4(exprMat, distMat, annoMat, group,
                                            LRpairs, TGs, Receiver, Sender, 
                                            far.ct, close.ct, downsample)
    
  }else if(length(Sender)==1){
    
    cellpair = paste(Sender,Receiver,sep = "_")
    mulNet = mulNetList[[cellpair]]
    
    TGs = mulNet %>% .[['TFTar']] %>% dplyr::select(target) %>% unlist() %>% as.character() %>% unique()
    LRpairs = get_LRTG_link(mulNet, TGs)
    
    cat(paste0("calculate the regulatory score of LR pairs from ",Sender,' to ',Receiver))
    LRTG_allscore = calculate_LRTG_score_V4(exprMat, distMat, annoMat, group,
                                            LRpairs, TGs, Receiver, Sender, 
                                            far.ct, close.ct, downsample)
    
  }else{
    
    cellpair = paste(Sender,Receiver,sep = "_")
    cellpair <- intersect(names(mulNetList),cellpair)
    if(length(cellpair)==0){
      LRTG_allscore = NA
      return(LRTG_allscore)
    }
    
    mulNetlist = mulNetList[cellpair]
    mulNet_tab = lapply(mulNetlist, function(mlnet){
      
      ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
      rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
      tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
      
      res = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
        merge(., tftg, by = 'TF') %>% 
        dplyr::select(Ligand, Receptor, TF, Target) %>% 
        arrange(Ligand, Receptor)
      
    })
    
    LRpairs_TGs_list <- lapply(mulNet_tab, function(ml_tab){
      
      lrpairs = by(ml_tab, as.character(ml_tab$Target), function(x){
        paste(x$Ligand, x$Receptor, sep = "_")
      })
      lrpairs = lapply(lrpairs, function(lrtg){lrtg[!duplicated(lrtg)]})
      tgs = names(lrpairs)
      
      list(LRpairs = lrpairs, TGs = tgs)
      
    })
    names(LRpairs_TGs_list) <- names(mulNet_tab)
    
    LRTG_allscore <- list()
    for (cp in cellpair) {
      
      receiver <- gsub('.*_','',cp)
      sender <- gsub('_.*','',cp)
      
      LRpairs <- LRpairs_TGs_list[[cp]]$LRpairs
      TGs <- LRpairs_TGs_list[[cp]]$TGs
      
      cat(paste0("calculate the regulatory score of LR pairs from ",sender,' to ',receiver))
      LRTG_allscore[[cp]] = calculate_LRTG_score_V4(exprMat, distMat, annoMat, group,
                                                    LRpairs, TGs, receiver, sender, 
                                                    far.ct, close.ct, downsample)
      
    }
  }
  
  return(LRTG_allscore)
}


## 获取指定TG上游LRpair
get_LRTG_link <- function(mulNet, TGs)
{
  
  ligrec = data.frame(Ligand = mulNet$LigRec$source, Receptor = mulNet$LigRec$target) 
  rectf = data.frame(Receptor = mulNet$RecTF$source, TF = mulNet$RecTF$target)
  tftg = data.frame(TF = mulNet$TFTar$source, Target = mulNet$TFTar$target)
  
  mulNet_tab = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
    merge(., tftg, by = 'TF') %>% 
    dplyr::select(Ligand, Receptor, TF, Target) %>% 
    arrange(Ligand, Receptor) %>%
    filter(Target %in% TGs)
  
  LRTG_link = by(mulNet_tab, as.character(mulNet_tab$Target), function(x){paste(x$Ligand, x$Receptor, sep = "_")})
  LRTG_link = lapply(LRTG_link, function(lrtg){lrtg[!duplicated(lrtg)]})
  LRTG_link = LRTG_link[TGs]
  
  return(LRTG_link)
}

## get_cell_pairs：获取与特定细胞类型的相近/相离的细胞配对
get_cell_pairs <- function(group=NULL, distMat, far.ct = 0.75, close.ct = 0.25)
{
  
  distMat_long <- reshape2::melt(distMat)
  colnames(distMat_long) <- c('Sender','Receiver','Distance')
  distMat_long$Sender <- as.character(distMat_long$Sender)
  distMat_long$Receiver <- as.character(distMat_long$Receiver)
  
  if(is.null(group)){
    respon_cellpair <- distMat_long[,1:2]
    group = 'all'
  }
  if(group == 'close') 
    respon_cellpair <- distMat_long[distMat_long$Distance <= quantile(distMat_long$Distance,close.ct),1:2]
  if(group == 'far') 
    respon_cellpair <- distMat_long[distMat_long$Distance >= quantile(distMat_long$Distance,far.ct),1:2]
  
  return(respon_cellpair)
}

## 求互信息
getMI <- function(LRTG_score)
{
  
  require(infotheo)
  
  # 获取df_keys
  df_keys = data.frame(
    LRpair = lapply(LRTG_score$LRs_score, function(df){colnames(df)}) %>% unlist(),
    TG = lapply(1:length(LRTG_score$LRs_score), function(i){
      rep(names(LRTG_score$LRs_score)[i],ncol(LRTG_score$LRs_score[[i]]))
    }) %>% unlist()
  )
  
  ## 计算相关性
  t1 <- Sys.time()
  df_mi <- lapply(1:nrow(df_keys), function(i){
    
    LRexpr = LRTG_score$LRs_score[[df_keys[i,'TG']]][,df_keys[i,'LRpair']]
    TGexpr = LRTG_score$TGs_expr[[df_keys[i,'TG']]]
    
    dat <- data.frame(LRexpr,TGexpr)
    dat <- discretize(dat)
    mi <- mutinformation(dat)
    mi[1,2]
    
  }) %>% unlist()
  t2 <- Sys.time()
  t2-t1
  df_mi <- do.call('cbind',list(df_keys,df_mi)) %>% as.data.frame()
  colnames(df_mi) <- c(colnames(df_keys),"MI")
  
  return(df_mi)
  
}

## 求相关性
getPCC <- function(LRTG_score)
{
  
  # 获取df_keys
  df_keys = data.frame(
    LRpair = lapply(LRTG_score$LRs_score, function(df){colnames(df)}) %>% unlist(),
    TG = lapply(1:length(LRTG_score$LRs_score), function(i){
      rep(names(LRTG_score$LRs_score)[i],ncol(LRTG_score$LRs_score[[i]]))
    }) %>% unlist()
  )
  
  ## 计算相关性
  t1 <- Sys.time()
  df_cor <- lapply(1:nrow(df_keys), function(i){
    
    LRexpr = LRTG_score$LRs_score[[df_keys[i,'TG']]][,df_keys[i,'LRpair']]
    TGexpr = LRTG_score$TGs_expr[[df_keys[i,'TG']]]
    
    cor <- cor.test(LRexpr,TGexpr,method = 'pearson')
    pval <- signif(cor$p.value,digits = 3)
    R <- signif(cor$estimate,digits = 2)
    c(R,pval)
    
  }) %>% do.call('rbind',.)
  t2 <- Sys.time()
  t2-t1
  colnames(df_cor) <- c("R","pval")
  df_cor <- cbind(df_keys,df_cor) %>% as.data.frame()
  
  return(df_cor)
  
}

####################
## train RF model ##
####################

# mian
get_pim_auto = function(trainx, trainy, ncores = 1, auto_para = TRUE, 
                        n.trees = 500, node.feature = 'sqrt', n.trys = 10, 
                        tree.method = 'variance', node.size = 5,
                        nPrem = 10, verbose = TRUE)
{
  
  # data
  LRpairs = colnames(trainx)
  data = as.data.frame(cbind(trainx, trainy))
  colnames(data) = c(paste0('LR',seq(ncol(trainx))), "Target")
  
  # get Lig, Rec
  LRTab = data.frame(LRpair = LRpairs)
  LRTab$Ligand = strsplit(LRpairs,"_") %>% do.call('rbind',.) %>% .[,1]
  LRTab$Receptor = strsplit(LRpairs,"_") %>% do.call('rbind',.) %>% .[,2]
  
  # remove sample which TG expr = 0
  cat("\nRemove zero Target")
  zeroTarget = which(trainy==0)
  if(length(zeroTarget)>1) data = data[-zeroTarget,]
  
  # remove LRpair which LRscore = 0
  cat("\nRemove zero LR")
  zeroLR = which(colSums(data[,-ncol(data)]==0)==nrow(data))
  zeroLR = c(zeroLR,which(is.na(colSums(data[,-ncol(data)]==0))))
  if(length(zeroLR)>1){
    data = data[,-zeroLR]
    LRTab = LRTab[-zeroLR,]
  }
  
  # normalize
  cat("\nNormalize data")
  data = scale(data)
  data = as.data.frame(data)
  
  # get optimal parameters of RF model
  cat("\nParameter Tuning")
  if(verbose){
    
    t1 <- Sys.time()
    cat(paste0("\nStart at: ",as.character(t1)))
    
  }
  if(auto_para == TRUE){
    
    # parameters of automatic tuning parameter
    fitControl <- trainControl(
      method = "cv",
      number = 5,
      search = 'random',
      allowParallel = TRUE)
    
    # model
    set.seed(2021)
    rfFit <- caret::train(Target ~ ., data = data, 
                          method = 'ranger', 
                          trControl = fitControl,
                          tuneLength = 20)
    
    # optimal parameter
    sel_mtry = as.numeric(rfFit$bestTune$mtry)
    sel_splitrule = as.character(rfFit$bestTune$splitrule)
    sel_min.node.size = as.numeric(rfFit$bestTune$min.node.size)
    sel_num.trees = n.trees
    
  }else{
    
    # sel_mtry = floor(sqrt(ncol(trainx)))
    sel_mtry = n.trys
    sel_splitrule = tree.method
    sel_min.node.size = node.size
    sel_num.trees = n.trees
    rfFit <- list()
    
  }
  parameters = list(sel_mtry, sel_splitrule, sel_min.node.size, sel_num.trees)
  names(parameters) = c('mtry','splitrule','min.node.size','n.trees')
  cat(paste0("\nThe final parameters used for the model: mtry = ",sel_mtry,
             ', splitrule = ',sel_splitrule, ', min.node.size = ',sel_min.node.size, 
             ' ,num.trees = ',sel_num.trees))
  if(verbose){
    
    t2 <- Sys.time()
    message(paste0("\nEnd at: ",as.character(t2)))
    message(paste0('\nAbout ',signif(t2-t1,digits = 4),' ',units(t2-t1)))
    
  } 
  
  # train RF model
  cat("\nTrain final model")
  finalFit <- ranger(formula = Target ~ ., data = data, 
                     num.trees = sel_num.trees, seed = 2021, 
                     splitrule = sel_splitrule, mtry = sel_mtry, 
                     min.node.size = sel_min.node.size,
                     importance = 'permutation', keep.inbag = TRUE, 
                     oob.error = TRUE, num.threads = ncores)
  
  # get permutation importance based on LRpair
  cat("\nobtain variable importance")
  df_IM = finalFit$variable.importance
  df_IM = data.frame(IM = df_IM)
  df_IM$LRpair = LRTab$LRpair
  
  # get permutation importance based on Lig/Rec
  cat("\nobtain permutation importance")
  if(verbose){
    
    t1 <- Sys.time()
    message(paste0("\nStart at: ",as.character(t1)))
    
  }
  df_pIM <- lapply(seq(nPrem), function(j){
    
    forest = finalFit
    data_X=data[,-ncol(data)]
    data_y=data$Target
    rf_oob_pim(forest, data_X, data_y, LRTab)
    
  })
  df_pIM <- Reduce("+", df_pIM)/nPrem
  if(verbose){
    
    t2 <- Sys.time()
    message(paste0("\nEnd at: ",as.character(t2)))
    message(paste0('\nAbout ',signif(t2-t1,digits = 4),' ',units(t2-t1)))
    
  }
  
  result <- list(model = finalFit,
                 df_IM = df_IM,
                 df_pIM = df_pIM,
                 rfFit = rfFit)
  
  return(result)
}

# get pIM from shuffle and OOB data
rf_oob_pim = function(forest, data_X, data_y, LRTab)
{
  
  # get pred value based on OOB data
  oob_pred = rf_oob_pred(forest, data_X)
  oob_pred = as.data.frame(oob_pred)
  
  # MSE before shuffle (500 x 1)
  oob_mse_bef = lapply(seq(ncol(oob_pred)), function(k){
    
    oob_pred_k = oob_pred[,k]
    mean((oob_pred_k-data_y)^2,na.rm = T)
    
  })
  oob_mse_bef = unlist(oob_mse_bef)
  
  # get Lig/Rec
  Vars = unique(c(LRTab$Ligand,LRTab$Receptor))
  
  # MSE after shuffle (500 x len(Vars))
  oob_mse_aft = lapply(seq(length(Vars)), function(i){
    
    # shuffle one Lig/Rec
    newX = shuffle_LigRec(data_X,i,LRTab)
    
    # pred value based on OOB data after shuffle one Lig/Rec
    oob_pred_i = rf_oob_pred(forest, newX)
    oob_pred_i = as.data.frame(oob_pred_i)
    
    # MSE after shuffle one Lig/Rec
    oob_mse_aft_i = lapply(seq(ncol(oob_pred_i)), function(k){
      
      oob_pred_i_k = oob_pred_i[,k]
      mean((oob_pred_i_k-data_y)^2,na.rm = T)
      
    })
    oob_mse_aft_i = unlist(oob_mse_aft_i)
    
  })
  oob_mse_aft = do.call("cbind",oob_mse_aft)
  
  # mean MSE of each Lig/Rec
  df_pIM = apply(oob_mse_aft, 2, function(oob_mse_aft_i){ mean(oob_mse_aft_i-oob_mse_bef) })
  df_pIM = data.frame(pIM = df_pIM)
  rownames(df_pIM) = paste0('shuffle_',Vars)
  
  return(df_pIM)
  
}

# OOB data
# get pred value of OOB data from each tree(n.trees = 500 => 500 trees)
# for each tree, trainx was split into OOB data (used for evaluation) and train data (used for training RF)
rf_oob_pred <- function(forest, data_X) 
{
  
  preds <- predict(forest, data_X, predict.all=TRUE)
  oob <- forest$inbag.counts
  oob <- do.call("cbind",oob)
  oob <- oob==0
  oob[which(!oob)] = NA
  preds.oob = oob*preds$predictions 
  return (preds.oob)
  
}

# shuffle
shuffle_LigRec <- function(data_X,i,LRTab)
{
  
  Vars = unique(c(LRTab$Ligand,LRTab$Receptor))
  var <- Vars[i]
  var_ids <- union(which(LRTab$Ligand %in% var),which(LRTab$Receptor %in% var))
  for(id in var_ids) {
    data_X[,id] <- sample(x = unlist(data_X[,id]), size = nrow(data_X), replace = FALSE)
  }
  return(data_X)
  
}

###################
## calculate_auc ##
###################

## get different metrices and ROC/PRC preformance object
get_evaluate_metrics <- function(pred,label)
{
  
  find_optimal_cutoff <- function(TPR, FPR, threshold){
    
    y = TPR - FPR
    Youden_index = which.max(y)
    # optimal_threshold = threshold[Youden_index]
    return(Youden_index)
    
  }
  
  if(length(which(label==TRUE))!=0 & length(which(label==FALSE))!=0){
    
    pred <- prediction(pred, label)
    
    perf_ROC <- performance(pred, measure = "tpr", x.measure = "fpr")
    ind_ROC <- find_optimal_cutoff(perf_ROC@y.values[[1]],perf_ROC@x.values[[1]],perf_ROC@alpha.values[[1]])
    cutoff_ROC <- perf_ROC@alpha.values[[1]][ind_ROC]
    
    perf_PRC <- performance(pred, measure = "prec", x.measure = "rec")
    ind_PRC <- find_optimal_cutoff(perf_PRC@y.values[[1]],perf_PRC@x.values[[1]],perf_PRC@alpha.values[[1]])
    cutoff_PRC <- perf_PRC@alpha.values[[1]][ind_PRC]
    
    ACC <- performance(pred, measure = "acc")@y.values[[1]] %>% .[ind_ROC] %>% signif(.,4)
    ERR <- performance(pred, measure = "err")@y.values[[1]] %>% .[ind_ROC]  %>% signif(.,4)
    PPV <- performance(pred, measure = "ppv")@y.values[[1]] %>% .[ind_ROC]  %>% signif(.,4)
    MCC <- performance(pred, measure = "mat")@y.values[[1]] %>% .[ind_ROC]  %>% signif(.,4)
    
    AUC <- performance(pred, measure = "auc")@y.values[[1]] %>% signif(.,4)
    AUCPR <- performance(pred, measure = "aucpr")@y.values[[1]] %>% signif(.,4)
    
    res = list(perf_ROC = perf_ROC,
               perf_PRC = perf_PRC,
               perf_metrics = c(ROC_AUC=AUC,PRC_AUC=AUCPR,
                                ACC=ACC,ERR=ERR,PPV=PPV,MCC=MCC))
    
  }else{
    
    res <- list(perf_ROC = NA, perf_PRC = NA, 
                perf_metrics = rep(0,6))
    names(res$perf_metrics) <- c('ROC_AUC','PRC_AUC','ACC','ERR','PPV','MCC')
    
  }
  
  
  return(res)
  
}

## get points from ROC/PRC preformance object to draw ROC/PRC curve
get_curve_input <- function(res_metrics, curve_type = 'ROC')
{
  
  if(curve_type == 'ROC'){
    
    model <- res_metrics$perf_ROC
    AUC <- res_metrics$perf_metrics['ROC_AUC'] %>% signif(.,4)
    
  }else{
    
    model <- res_metrics$perf_PRC
    AUC <- res_metrics$perf_metrics['PRC_AUC'] %>% signif(.,4)
    
  }
  
  if(class(model)[1] != 'performance'){
    v_x <- 0
    v_y <- 0
  }else{
    v_x <- model@x.values[[1]]
    v_y <- model@y.values[[1]]
  }
  
  res <- data.frame(x = v_x, y = v_y, AUC = AUC, row.names = NULL)
  res[is.na(res)] <- 0
  if(curve_type == 'ROC') colnames(res) <- c('FPR','TPR','AUC')
  if(curve_type == 'PRC') colnames(res) <- c('Recall','Precision','AUC')
  
  return(res)
}

## organize the result from 'get_evaluate_metrics' function
clean_res_metrcs <- function(res_metrics, sheetID, method){
  
  df_metrics <- res_metrics$perf_metrics
  df_metrics <- matrix(c(df_metrics,sheetID)) %>% t() %>% as.data.frame()
  rownames(df_metrics) <- NULL
  colnames(df_metrics) <- c(names(res_metrics$perf_metrics),'sheet')
  df_metrics$method <- method
  
  df_ROC <- get_curve_input(res_metrics = res_metrics, curve_type = 'ROC')
  df_ROC$method <- method
  df_ROC$sheet <- sheetID 
  
  df_PRC <- get_curve_input(res_metrics = res_metrics, curve_type = 'PRC')
  df_PRC$method <- method
  df_PRC$sheet <- sheetID
  
  res <- list(df_metrics=df_metrics,df_ROC=df_ROC,df_PRC=df_PRC)
}

##########################
## search feedback loop ##
##########################

check_feedback_loop <- function(celltypes, workdir, savepath, cts_of_interest){
  
  cp_of_inter <- data.frame(matrix(ncol = 4,dimnames = list(c(),c('ct1','ct2','keys_of_ct1','keys_of_ct2'))))
  for (i in 1:length(celltypes)) {
    
    ct1 <- celltypes[i]
    for (j in (i+1):length(celltypes)) {
      
      ct2 <- celltypes[j]
      cp1 <- paste0(ct1,'_',ct2)
      cp2 <- paste0(ct2,'_',ct1)
      if(cp1 %in% list.files(paste0(workdir,'/runscMLnet/')) & cp2 %in% list.files(paste0(workdir,'/runscMLnet/'))){
        
        cat('check in ',cp1,' and ',cp2,'\n')
        mlnet1 <- readRDS(paste0(workdir,cp1,'/scMLnet.rds'))
        mlnet2 <- readRDS(paste0(workdir,cp2,'/scMLnet.rds'))
        
        ct2_tgs <- unique(mlnet1$TFTar$target)
        ct1_tgs <- unique(mlnet2$TFTar$target)
        
        ct1_ligs <- unique(mlnet1$LigRec$source)
        ct2_ligs <- unique(mlnet2$LigRec$source)
        
        ct1_keys <- intersect(ct1_ligs,ct1_tgs)
        ct2_keys <- intersect(ct2_ligs,ct2_tgs)
        
        if(length(ct2_keys)>0|length(ct1_keys)>0){
          cp_of_inter <- rbind(cp_of_inter,c(ct1,ct2,length(ct1_keys),length(ct2_keys)))
        }
      }
    }
    
    cp_of_inter <- na.omit(cp_of_inter)
    cp_of_inter <- cp_of_inter[cp_of_inter$keys_of_ct1 != 0 & cp_of_inter$keys_of_ct2 != 0,]
    rownames(cp_of_inter) <- 1:nrow(cp_of_inter)
    
    ## get genes list ####
    
    key_of_inter <- list()
    for (k in 1:nrow(cp_of_inter)) {
      
      ct1 <- cp_of_inter$ct1[k]
      ct2 <- cp_of_inter$ct2[k]
      
      mlnet1 <- readRDS(paste0(workdir,paste(ct1,ct2,sep = '_'),'/scMLnet.rds'))
      mlnet2 <- readRDS(paste0(workdir,paste(ct2,ct1,sep = '_'),'/scMLnet.rds'))
      
      ct2_tgs <- unique(mlnet1$TFTar$target)
      ct1_tgs <- unique(mlnet2$TFTar$target)
      
      ct1_ligs <- unique(mlnet1$LigRec$source)
      ct2_ligs <- unique(mlnet2$LigRec$source)
      
      ct1_keys <- intersect(ct1_ligs,ct1_tgs)
      ct2_keys <- intersect(ct2_ligs,ct2_tgs)
      
      key_of_inter[[paste(ct1,ct2,sep = '_')]] <- list(
        ct1_keys = ct1_keys,
        ct2_keys = ct2_keys
      )
      
    }
  }
  
  ## transform to table ####
  
  fbloop <- lapply(1:length(key_of_inter), function(i){
    
    ls <- key_of_inter[[i]]
    cp <- names(key_of_inter)[i]
    ct1 <- strsplit(cp,'_')[[1]][1]
    ct2 <- strsplit(cp,'_')[[1]][2]
    ct1_keys <- paste(ls$ct1_keys,collapse = ' ')
    ct2_keys <- paste(ls$ct2_keys,collapse = ' ')
    rbind(c(ct1,ct2,ct2_keys),c(ct2,ct1,ct1_keys))
    
  }) %>% do.call('rbind',.) %>% as.data.frame()
  colnames(fbloop) <- c('Sender','Receiver','Siganl')
  
  if (length(cts_of_interest) == 0){
    write.csv(fbloop,paste0(savepath,'feddback_loop_detail.csv'))
  }else{
    cts_of_interest <- cts_of_interest
    fbloop <- fbloop[fbloop$Sender %in% cts_of_interest & fbloop$Receiver %in% cts_of_interest,]
    write.csv(fbloop,paste0(savepath,'feddback_loop_detail.csv'))
  }
  
  return(fbloop)
}

##################
## visualizeCCI ##
##################

# NetworkPlot
DrawLeg <- function(alltype,allcolos,xstart,ystart,cirr,jiange)
{
  ThisX <- xstart
  ThisY <- ystart
  for(i in seq(1,length(alltype)))
  {
    ThisType <- alltype[i]
    draw.circle(ThisX,ThisY,cirr,col = allcolos[ThisType])
    text(ThisX+cirr+0.1,ThisY,ThisType,adj = 0)
    
    ThisY <- ThisY - (2*cirr + jiange)
  }
}

DrawCellComm <- function(CellTab,colodb,gtitle = 'TME')
{
  
  aaa <- CellTab
  g <- graph.data.frame(aaa,directed=TRUE)
  alltype <- c(aaa$cell_from,aaa$cell_to)
  alltype <- unique(alltype)
  # allcolos <- rainbow(length(alltype))
  ChoColorNum <- alltype# 1:length(alltype) # sample(1:length(colodb),length(alltype), replace = FALSE)
  
  allcolos <- colodb[ChoColorNum]
  names(allcolos) <- alltype
  
  edge.start <- ends(g,es=E(g),names=FALSE)
  
  layout <- in_circle()  #igraph packages
  coords <- layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  
  loop.angle <- ifelse(coords_scale[V(g),1]>0,-atan(coords_scale[V(g),2]/coords_scale[V(g),1]),pi-atan(coords_scale[V(g),2]/coords_scale[V(g),1]))
  
  vertex.label.color <- 'black'
  V(g)$size <- 20
  V(g)$color <- allcolos[V(g)]
  V(g)$label.color <- vertex.label.color
  V(g)$label.cex <- 0
  
  label <- FALSE
  if(label){
    E(g)$label<-E(g)$n
  }
  
  edge.max.width = 10
  if(max(E(g)$n)==min(E(g)$n)){
    E(g)$width <- 2
  }else{
    E(g)$width <- 1+edge.max.width/(max(E(g)$n)-min(E(g)$n))*(E(g)$n-min(E(g)$n))
  }
  
  edge.label.color <- "black"
  E(g)$arrow.width<- 1
  E(g)$arrow.size <- 0.1
  E(g)$label.color <- edge.label.color
  E(g)$label.cex <- 1
  E(g)$color <- V(g)$color[edge.start[,1]]
  
  if(sum(edge.start[,2]==edge.start[,1])!=0)
  {
    E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  
  #draw 
  shape = 'circle'
  margin = 0.2
  edgecurved = 0.2
  vertexlabelcex = 1
  plot(
    g,
    edge.curved=edgecurved,
    vertex.label = "",
    vertex.shape=shape,
    layout=coords_scale,
    margin=margin,
    vertex.shape="fcircle",
    vertex.label.cex = vertexlabelcex,
    axes = FALSE
  )
  
  #draw legend
  xstart <- 1.5
  ystart <- 0.5
  jiange <- 0.02
  cirr <- 0.04
  DrawLeg(alltype,allcolos,xstart,ystart,cirr,jiange)
  
  #title
  title(main = list(gtitle, cex=2))
  
}

# MLnetPlot

prepareMLnetworkPlotData_V3 <- function(
  mlnet = NULL, # 画图依据
  lrtg_im = NULL,  #过滤依据
  Key,Type, do.check = FALSE
)
{
  
  # Key = 'IGF1'
  # Type = 'Ligand'
  # mlnet = MLnet
  # lrtg_im = LRTG_im
  
  if(Type == 'Ligand'){
    
    mlnet$LigRec <- mlnet$LigRec[mlnet$LigRec$source %in% Key,]
    mlnet$RecTF <- mlnet$RecTF[mlnet$RecTF$source %in% mlnet$LigRec$target,]
    mlnet$TFTar <- mlnet$TFTar[mlnet$TFTar$source %in% mlnet$RecTF$target,]
    lrtg_im = lrtg_im %>% dplyr::select(.,Ligand,Receptor,Target,Score = im_norm) %>% dplyr::filter(Ligand %in% Key)
    
  }else{
    
    mlnet$LigRec <- mlnet$LigRec[mlnet$LigRec$target %in% Key,]
    mlnet$RecTF <- mlnet$RecTF[mlnet$RecTF$source %in% mlnet$LigRec$target,]
    mlnet$TFTar <- mlnet$TFTar[mlnet$TFTar$source %in% mlnet$RecTF$target,]
    lrtg_im = lrtg_im %>% dplyr::select(Ligand,Receptor,Target,Score = im_norm) %>% dplyr::filter(Receptor %in% Key)
    
  }
  
  if(do.check & nrow(lrtg_im)>10){
    
    message('using the score of LRTG to narrow down for better visualization')
    tg_im_check <- lrtg_im %>% group_by(Target) %>% 
      summarise(sum_Score=sum(Score)) %>% arrange(desc(sum_Score)) %>% 
      dplyr::select(Target) %>% unlist() %>% head(10)
    
    mlnet_check <- mlnet
    mlnet_check$TFTar <- mlnet_check$TFTar[mlnet_check$TFTar$target %in% tg_im_check,]
    mlnet_check$RecTF <- mlnet_check$RecTF[mlnet_check$RecTF$target %in% mlnet_check$TFTar$source,]
    mlnet_check$LigRec <- mlnet_check$LigRec[mlnet_check$LigRec$target %in% mlnet_check$RecTF$source,]
    mlnet <- mlnet_check
    
  }
  
  return(mlnet)
}

drawMLnetworkPlot_V4 <- function(mlnet,downstream = c('TF','Target'),colodb,
                                 gtitle = 'TME',wd = './',p_height=4,p_width=6)
{
  
  # mlnet <- MLnet_key
  Ligs <- unique(mlnet$LigRec$source)
  Recs <- unique(mlnet$LigRec$target)
  TFs <- unique(mlnet$TFTar$source)
  TGs <- unique(mlnet$TFTar$target)
  TGs <- TGs[!TGs %in% c(Ligs,Recs,TFs)]
  df_nodes <- data.frame(name = c(Ligs, Recs, TFs, TGs),
                         key = c(rep('Ligand', length(Ligs)),
                                 rep('Receptor', length(Recs)),
                                 rep('TF', length(TFs)),
                                 rep('Target', length(TGs))),
                         stringsAsFactors = FALSE)
  
  df_nodes <- df_nodes[!duplicated(df_nodes$name), ]
  
  # 边表也检查 NA
  df_edges <- do.call("rbind", list(mlnet$LigRec, mlnet$RecTF, mlnet$TFTar))
  df_edges <- df_edges[!is.na(df_edges$source) & !is.na(df_edges$target), ]
  
  # 构建 igraph 对象
  subnet <- graph_from_data_frame(d = df_edges, vertices = df_nodes)
  df_nodes <- df_nodes[match(names(V(subnet)),df_nodes$node),]
  root_index <- grep('Ligand',df_nodes$key)
  set.seed(4)
  coords <- layout_(subnet,layout = as_tree(root = root_index)) %>% as.data.frame()
  coords <- cbind(coords,df_nodes$key)
  colnames(coords) <- c('dim_x','dim_y','type')
  
  dist_TG <- 1;len_TG <- table(coords$type)[['Target']]
  dist_TF <- 1.5;len_TF <- table(coords$type)[['TF']]
  dist_Rec <- 2.5;len_Rec <- table(coords$type)[['Receptor']]
  dist_Lig <- 2.5;len_Lig <- table(coords$type)[['Ligand']]
  if(len_Lig==1){
    coords$dim_x[coords$type == 'Ligand'] = 0
  }else{
    dim_x_1 = seq(to = -dist_Lig/2,by = dist_Lig,length.out = ceiling(len_Lig/2))
    dim_x_2 = seq(from = dist_Lig/2,by = dist_Lig,length.out = len_Lig-ceiling(len_Lig/2))
    coords$dim_x[coords$type == 'Ligand'] = c(dim_x_1,dim_x_2) %>% scale(.,scale = F) %>% as.vector()
  }
  if(len_Rec==1){
    coords$dim_x[coords$type == 'Receptor'] = 0
  }else{
    dim_x_1 = seq(to = -dist_Rec/2,by = dist_Rec,length.out = ceiling(len_Rec/2))
    dim_x_2 = seq(from = dist_Rec/2,by = dist_Rec,length.out = len_Rec-ceiling(len_Rec/2))
    coords$dim_x[coords$type == 'Receptor'] = c(dim_x_1,dim_x_2) %>% scale(.,scale = F) %>% as.vector()
  }
  if(len_TG<len_TF & len_TF>10){
    
    dist_TG <- 1
    dist_TF <- 0.5
    
  }else if(len_TG>len_TF & len_TG>10){
    
    dist_TG <- 0.5
    dist_TF <- 1
    
  }
  if(len_TF==1){
    coords$dim_x[coords$type == 'TF'] = 0
  }else{
    dim_x_1 = seq(to = -dist_TF/2,by = dist_TF,length.out = ceiling(len_TF/2))
    dim_x_2 = seq(from = dist_TF/2,by = dist_TF,length.out = len_TF-ceiling(len_TF/2))
    coords$dim_x[coords$type == 'TF'] = c(dim_x_1,dim_x_2) %>% scale(.,scale = F) %>% as.vector()
  }
  if(len_TG==1){
    coords$dim_x[coords$type == 'Target'] = 0
  }else{
    dim_x_1 = seq(to = -dist_TG/2,by = dist_TG,length.out = ceiling(len_TG/2))
    dim_x_2 = seq(from = dist_TG/2,by = dist_TG,length.out = len_TG-ceiling(len_TG/2))
    coords$dim_x[coords$type == 'Target'] = c(dim_x_1,dim_x_2) %>% scale(.,scale = F) %>% as.vector()
  }
  coords$dim_y <- lapply(coords$type,switch,'Ligand'=0.9,'Receptor'=0.6,'TF'=0.3,'Target'=0) %>% unlist()
  # plot(subnet, layout = as.matrix(coords[,1:2]), vertex.color=df_nodes$color)
  
  layout <- create_layout(subnet, layout = 'tree')
  layout[,1:2] <- coords[,1:2]
  # head(layout)
  
  temp <- function(key){
    
    max(coords$dim_x[coords$type==key])+0.5
    
  }
  df_anno <- data.frame(x = lapply(unique(layout$key), temp) %>% unlist(),
                        y = c(0.92,0.62,0.32,0.02),
                        lab = unique(layout$key),
                        key = unique(layout$key))
  if(downstream != 'Target'){
    
    layout <- layout[layout$key != 'Target',]
    df_anno <- df_anno[df_anno$key != 'Target',]
    
  }
  
  # windowsFonts(HEL=windowsFont("Helvetica CE 55 Roman"),
  #              RMN=windowsFont("Times New Roman"),
  #              ARL=windowsFont("Arial"))
  
  pt <- ggraph(layout) +
    geom_edge_link(aes(),color="grey",
                   arrow = arrow(length = unit(1.5, 'mm')), 
                   start_cap = circle(3, 'mm'),
                   end_cap = circle(3, 'mm')) + 
    geom_node_point(aes(fill = key,color = key),shape=21,size = 8) +
    geom_node_text(aes(label=name),size=2) + # ,fontface='bold',family='Arial'
    xlim(c(min(layout$x)-2,max(layout$x)+2)) +
    geom_text(data = df_anno, aes(x,y,label=lab,color=key),
              vjust = 1, hjust = 0, size = 4,fontface='bold') + # ,family='ARL'
    scale_fill_manual(values = colodb) + scale_color_manual(values = colodb) +
    guides(fill='none',color='none') +
    theme_graph()
  print(pt)
  
  ## save
  
  pdf(paste0(wd,'MLnetPlot3-',gtitle,'.pdf'),height = p_height,width = p_width)
  print(pt)
  dev.off()
  
  png(paste0(wd,'MLnetPlot3-',gtitle,'.png'),height = p_height,width = p_width, units = 'in', res = 300)
  print(pt)
  dev.off()
  
}

drawMLnetworkPlot_V5 <- function(mlnet, downstream = c('TF','Target'), colodb,
                                 gtitle = 'TME', wd = './', p_height = 4, p_width = 6) {
  
  Ligs <- unique(mlnet$LigRec$source)
  Recs <- unique(mlnet$LigRec$target)
  TFs  <- unique(mlnet$TFTar$source)
  TGs_all <- unique(mlnet$TFTar$target)
  
  # 保留 Target 中不属于 Ligand 或 Receptor 的，允许与 TF 重叠
  TGs <- setdiff(TGs_all, union(Ligs, Recs))
  if (length(TGs) == 0) {
    message("⚠️ 无独立 Target 节点，图中将不包含 Target 层")
    downstream <- "TF"
  }
  
  # 构建 node 表
  df_nodes <- data.frame(
    name = c(Ligs, Recs, TFs, TGs),
    key = c(rep('Ligand', length(Ligs)),
            rep('Receptor', length(Recs)),
            rep('TF', length(TFs)),
            rep('Target', length(TGs))),
    stringsAsFactors = FALSE
  )
  df_nodes <- df_nodes[!duplicated(df_nodes$name), ]
  df_nodes$color <- colodb[df_nodes$key]
  
  # 构建 edge 表
  df_edges <- do.call("rbind", list(mlnet$LigRec, mlnet$RecTF, mlnet$TFTar))
  df_edges <- df_edges[!is.na(df_edges$source) & !is.na(df_edges$target), ]
  colnames(df_edges)[1:2] <- c("from", "to")  # Ensure consistent naming
  
  # 创建 igraph 对象
  library(igraph)
  subnet <- graph_from_data_frame(d = df_edges, vertices = df_nodes)
  
  # 匹配节点顺序
  df_nodes <- df_nodes[match(names(V(subnet)), df_nodes$name), ]
  root_index <- grep('Ligand', df_nodes$key)
  
  # 布局
  set.seed(4)
  coords <- layout_(subnet, layout = as_tree(root = root_index)) %>% as.data.frame()
  coords <- cbind(coords, df_nodes$key)
  colnames(coords) <- c('dim_x', 'dim_y', 'type')
  
  ## 布局位置调整函数
  adjust_x <- function(type_label, dist, len) {
    if (len == 1) {
      return(0)
    } else {
      dim_x_1 <- seq(to = -dist / 2, by = dist, length.out = ceiling(len / 2))
      dim_x_2 <- seq(from = dist / 2, by = dist, length.out = len - ceiling(len / 2))
      return(c(dim_x_1, dim_x_2) %>% scale(., scale = FALSE) %>% as.vector())
    }
  }
  
  # 各层数量
  len_Lig <- sum(coords$type == 'Ligand')
  len_Rec <- sum(coords$type == 'Receptor')
  len_TF  <- sum(coords$type == 'TF')
  len_TG  <- sum(coords$type == 'Target')
  
  coords$dim_x[coords$type == 'Ligand']   <- adjust_x('Ligand', 1, len_Lig)
  coords$dim_x[coords$type == 'Receptor'] <- adjust_x('Receptor', 1, len_Rec)
  
  # 根据 Target 和 TF 相对数量决定纵向间距
  dist_TG <- if (len_TG > len_TF & len_TG > 10) 1 else 1
  dist_TF <- if (len_TF > len_TG & len_TF > 10) 1 else 1
  
  coords$dim_x[coords$type == 'TF']     <- adjust_x('TF', dist_TF, len_TF)
  coords$dim_x[coords$type == 'Target'] <- adjust_x('Target', dist_TG, len_TG)
  
  coords$dim_y <- sapply(coords$type, switch,
                         'Ligand' = 0.9,
                         'Receptor' = 0.6,
                         'TF' = 0.3,
                         'Target' = 0)
  
  # 生成 layout 用于 ggraph
  library(ggraph)
  layout <- create_layout(subnet, layout = 'tree')
  layout[, 1:2] <- coords[, 1:2]
  layout$key <- coords$type
  
  # 构建注释信息
  df_anno <- data.frame(
    x = sapply(unique(layout$key), function(k) max(coords$dim_x[coords$type == k]) + 0.5),
    y = c(0.92, 0.62, 0.32, 0.02)[1:length(unique(layout$key))],
    lab = unique(layout$key),
    key = unique(layout$key)
  )
  
  # 移除 Target 层（如果用户不想显示）
  if (downstream != 'Target') {
    layout <- layout[layout$key != 'Target', ]
    df_anno <- df_anno[df_anno$key != 'Target', ]
  }
  
  ## 绘图
  pt <- ggraph(layout) +
    geom_edge_link(aes(), color = "grey",
                   arrow = arrow(length = unit(1.5, 'mm')),
                   start_cap = circle(3, 'mm'),
                   end_cap = circle(3, 'mm')) +
    geom_node_point(aes(fill = key, color = key), shape = 21, size = 8) +
    geom_node_text(aes(label = name), size = 2) +
    xlim(c(min(layout$x) - 2, max(layout$x) + 2)) +
    geom_text(data = df_anno, aes(x, y, label = lab, color = key),
              vjust = 1, hjust = 0, size = 4, fontface = 'bold') +
    scale_fill_manual(values = colodb) +
    scale_color_manual(values = colodb) +
    guides(fill = 'none', color = 'none') +
    theme_graph()
  
  print(pt)
  
  ## 保存图像
  pdf(paste0(wd, 'MLnetPlot3-', gtitle, '.pdf'), height = p_height, width = p_width)
  print(pt)
  dev.off()
  
  png(paste0(wd, 'MLnetPlot3-', gtitle, '.png'), height = p_height, width = p_width, units = 'in', res = 300)
  print(pt)
  dev.off()
}


# EdgeBundlingPlot

prepareEdgeBundlingPlotData_V2 <- function(df_LRTGscore, do.check = FALSE){
  
  # edge
  df_LRTGscore$from = paste(df_LRTGscore$source,df_LRTGscore$source_group,sep = "_")
  df_LRTGscore$to = paste(df_LRTGscore$target,df_LRTGscore$target_group,sep = "_")
  df_edge <- df_LRTGscore[,c('from','to')]
  
  # hierarchy
  d1 <- data.frame(from="cell", to=c('Sender','Receiver'))
  d2 <- data.frame(from=c(rep('Sender',length(unique(df_LRTGscore$source_group))),'Receiver'), 
                   to=c(unique(df_LRTGscore$source_group),unique(df_LRTGscore$target_group)))
  d3 <- data.frame(from=c(df_LRTGscore$source_group,df_LRTGscore$target_group), 
                   to=c(df_edge$from,df_edge$to))
  d3 <- d3[!duplicated(d3),]
  df_hierarchy <- do.call('rbind',list(d1,d2,d3))
  
  # node
  df_node  <-  data.frame(
    name = unique(c(as.character(df_hierarchy$from), as.character(df_hierarchy$to)))
  ) 
  df_node$label <- gsub("_.*","",df_node$name)
  df_node$group <- df_hierarchy$from[match(df_node$name, df_hierarchy$to)]
  df_node$group <- factor(df_node$group,levels=unique(df_node$group))
  df_node$score <- df_LRTGscore$count[match(df_node$name,df_LRTGscore$from)]
  df_node$score <- replace(df_node$score,is.na(df_node$score),1e-5)
  
  # check size
  if(do.check & nrow(df_edge)>50){
    
    message('try to narrow down the number of LR for better visualization')
    df_node$rank <- Inf
    df_node$rank[df_node$group %in% d2$to[d2$from == 'Sender']] <- df_node$score[df_node$group %in% d2$to[d2$from == 'Sender']] %>% rank()
    
    df_node_check <- df_node
    df_node_check <- df_node_check[df_node_check$rank >= max(df_node_check$rank[is.finite(df_node_check$rank)])-49,]
    
    df_edge_check <- df_edge[df_edge$from %in% df_node_check$name,]
    
    df_node_check <- df_node_check[df_node_check$name %in% 
                                     c(d1$from,d1$to,
                                       df_node_check$group %>% as.character(),
                                       df_edge_check$from,df_edge_check$to),]
    
    df_hierarchy_check <- df_hierarchy[df_hierarchy$to %in% df_node_check$name,]
    
    df_node <- df_node_check
    rownames(df_node) <- seq(nrow(df_node))
    df_edge <- df_edge_check
    rownames(df_edge) <- seq(nrow(df_edge))
    df_hierarchy <- df_hierarchy_check
    rownames(df_hierarchy) <- seq(nrow(df_hierarchy))
    
  }
  
  # attr
  
  index_myleaves <- grep('^cell$|Sender|Receiver',df_node$group,invert = T)[-1]
  index_midleaves <- floor(length(index_myleaves)/2) # floor(nrow(df_node)/2)
  index_myleaves <- c(index_midleaves:max(index_myleaves),min(index_myleaves):(index_midleaves-1))
  nleaves <- length(index_myleaves)
  
  df_node$id <- 0
  df_node$id[index_myleaves] <- seq(1:length(index_myleaves))
  
  df_node$angle <- 90 - 360 * df_node$id / nleaves
  df_node$hjust <- lapply(df_node$angle,function(x){
    if(x < -90 & x != 90) hjust <- 0.5
    if(x >= -90 & x != 90) hjust <- 2
    if(x == 90) hjust <- 0
    hjust
  }) %>% unlist()
  df_node$angle <- ifelse(df_node$angle < -90, df_node$angle+180, df_node$angle)
  
  # merge
  df_input <- list(df_node = df_node, df_edge = df_edge, df_hierarchy = df_hierarchy)
  return(df_input)
  
}

drawEdgeBundlingPlot <-function(df_input,colodb,gtitle = 'TME',wd = './', p_height=7.5, p_width=7){
  
  # graph
  
  df_node <- df_input$df_node # node attr
  df_hierarchy <- df_input$df_hierarchy # edge
  mygraph <- graph_from_data_frame(df_hierarchy, vertices = df_node)
  
  node_color <- colodb[names(colodb) %in% unique(df_node$label)]
  
  df_edge <- df_input$df_edge # 
  from  <-  match(df_edge$from, df_node$name)
  to  <-  match(df_edge$to, df_node$name)
  
  ## plot
  
  pt <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
    #alpha线的透明度，width线的宽度，tension是线的“密集”程度
    geom_conn_bundle(data = get_con(from = from, to = to), 
                     alpha=1, aes(colour=group, width=score),tension=0.5) +
    # 设置边颜色
    scale_edge_color_manual(values = node_color,guide='none')+
    # 设置边粗细
    scale_edge_width(range = c(0,1)) + 
    # 设置节点标签，字体大小，文本注释信息
    geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label=label,
                       angle = angle, hjust=hjust*0.5, colour=group), size=3, alpha=1, show.legend=F) +
    # 设置节点的大小，颜色和透明度
    geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, 
                        colour=group, size=score, alpha=1)) +
    # 设置节点颜色
    scale_colour_manual(values = node_color,name="celltype",
                        breaks=names(node_color)[!names(node_color) %in% c('cell','Sender','Receiver')],) +
    # 设置节点大小的范围
    scale_size_continuous(range = c(1,5)) + theme_void() +
    # 去除alpha图例
    scale_alpha(guide = 'none') +
    theme(
      legend.position = 'bottom',
      # legend.direction = 'vertical',
      legend.box = 'vertical'
    ) 
  print(pt)
  
  ## save
  
  pdf(paste0(wd,'EdgeBundlingPlot-',gtitle,'.pdf'),height = p_height, width = p_width)
  print(pt)
  dev.off()
  
  png(paste0(wd,'EdgeBundlingPlot-',gtitle,'.png'),
      height = p_height, width = p_width, units = 'in', res = 300)
  print(pt)
  dev.off()
  
  
}

# AlluviumPlot

prepareAlluviumPlotData_V2 <- function(
  mlnet = NULL, # list of mlnet consisting of LigRec, RecTF, TFTar
  lrtg_im = NULL, # replace of mlnet.list when using score
  color.by = 'Nodekey', # color: Nodekey/Node/cellpair/Ligand/...
  do.check = FALSE
){
  
  if(!is.null(mlnet)){
    
    # mlnet = MLnet
    df_ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
    df_rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
    df_tftar = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
    
    df_mlnet = df_ligrec %>% merge(., df_rectf, by = 'Receptor') %>% 
      merge(., df_tftar, by = 'TF') %>% 
      dplyr::select(Ligand, Receptor, TF, Target) %>% 
      arrange(Ligand, Receptor)
    df_mlnet$Score <- 1
    df_mlnet$ID <- seq(nrow(df_mlnet))
    df_plot <- df_mlnet
    
  }else if(!is.null(lrtg_im)){
    
    # lrtg_im = LRTG_im_merge
    df_lrtg_im = lrtg_im
    colnames(df_lrtg_im)[grep('im_norm',colnames(df_lrtg_im))] <- 'Score'
    df_lrtg_im$ID <- seq(nrow(df_lrtg_im))
    df_plot <- df_lrtg_im
    
  }
  
  if(do.check & nrow(df_plot)>30){
    
    if(sd(df_plot$Score)!=0){
      
      message('using the score of LRTG to narrow down for better visualization')
      df_plot_check <- df_plot
      df_plot_check$rank <- df_plot_check$Score %>% rank()
      df_plot_check <- df_plot_check[df_plot_check$rank >= max(df_plot_check$rank)-29,]
      df_plot <- df_plot_check
      
    }else{
      
      message('using the number of Target to narrow down for better visualization')
      df_plot_check <- df_plot
      df_plot_check <- df_plot_check %>% group_by(Ligand,Receptor,TF) %>% 
        summarise(count=n()) %>% arrange(desc(count)) %>% ungroup() %>%
        inner_join(df_plot_check,., by = c('Ligand','Receptor','TF')) 
      df_plot_check$rank <- df_plot_check$count %>% rank()
      df_plot_check <- df_plot_check[df_plot_check$rank >= max(df_plot_check$rank)-29,]
      df_plot <- df_plot_check
      
    }
    
  }
  
  df_plot <- df_plot[,-grep('LRpair|IM|rank',colnames(df_plot))]
  df_plot_long <- reshape2::melt(df_plot, id = c('Score','ID','Sender'))
  colnames(df_plot_long) <- c("Score",'ID','Sender',"Nodekey","Node")
  
  # df_plot_long$Score <- (df_plot_long$Score-min(df_plot_long$Score))/(max(df_plot_long$Score)-min(df_plot_long$Score))
  
  if(color.by == 'Nodekey'){
    df_plot_long$color <- df_plot_long$Nodekey
  }else if(color.by == 'Node'){
    df_plot_long$color <- df_plot_long$Node
  }else{
    df_plot_long$color <- df_plot[[color.by]][match(df_plot_long$ID,df_plot$ID)]
  }
  
  return(df_plot_long)
  
}

drawAlluviumPlot <- function(
  dat, colodb, 
  gtitle = 'groundtrue', # vector
  wd = './visualize_CCI/AlluviumPlot/',
  p_height=9, p_width=18
){
  
  # dat <- df_plot_long
  # dat <- df_MLnet_long_check
  
  pt_color <- colodb
  if(!any(unique(dat$color) %in% names(pt_color))){
    pt_color <- rainbow(length(unique(dat$color)))
    names(pt_color) <- unique(dat$color)
  }else{
    pt_color <- pt_color[names(pt_color) %in% unique(dat$color)]
  }
  
  pt <- ggplot(data = dat,
               aes(x = Nodekey,y = Score,
                   stratum = Node, alluvium = ID,label = Node)) +
    geom_flow(aes(fill = color)) + geom_stratum() + 
    theme_minimal() + geom_text(stat = "stratum", size = 3,check_overlap = F) +
    ggtitle("Multilayer Network", gtitle) +
    scale_fill_manual(values = pt_color,  name="color",breaks=names(pt_color)) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 11),
      plot.title = element_text(size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      panel.grid = element_blank(),
      legend.position = 'none'
    )
  pt
  
  ## save
  
  pdf(paste0(wd,'AlluviumPlot-',gtitle,'.pdf'),height = p_height,width = p_width)
  print(pt)
  dev.off()
  
  png(paste0(wd,'AlluviumPlot-',gtitle,'.png'),height = p_height,width = p_width, units = 'in', res = 300)
  print(pt)
  dev.off()
  
  return(pt)
}

# Heatmap plot pval 
# Heatmap plot pval 
DrawHeatmapPlot <- function(InputDir, Sender, Receiver, outputdir,p_height = 9,p_width = 8.4){
  
  inputdir <- InputDir
  sender <- Sender
  receiver <- Receiver
  
  wd <- paste0(InputDir,"/runscMLnet/")
  files <- list.files(wd)[grep('TME_',list.files(wd),invert = TRUE)]
  
  LRI_allpval <- lapply(files, function(f) {
    cat("Processing file:", f, "\n")
    
    LRI_pval <- readRDS(paste0(wd, f, "/cellpair_LRI_pval.rds"))
    print(str(LRI_pval))  # Print the structure of LRI_pval
    
    if (length(LRI_pval) > 0) {
      LRI_pval$Sender <- strsplit(f, "_")[[1]][1]
      LRI_pval$Receiver <- strsplit(f, "_")[[1]][2]
      return(LRI_pval)
    } else {
      # If LRI_pval is empty, return a placeholder or handle it as needed
      cat("Warning: Empty LRI_pval for file", f, "\n")
      return(NULL)  # or return an empty data frame, depending on your needs
    }
  }) %>% do.call('rbind', .)
  
  if (length(sender) != 0 ){
    
    LRI_pval_Inter <- LRI_allpval[LRI_allpval$Sender == sender,]
    df_plot <- data.frame(cellpair = paste0(LRI_pval_Inter$Sender,"_",LRI_pval_Inter$Receiver),
                          LRpair = paste0(LRI_pval_Inter$source,"_",LRI_pval_Inter$target),
                          pval = LRI_pval_Inter$pval)
    
    df_plot <- df_plot[order(df_plot$pval,decreasing=F),] 
    
    # load LR signaling score
    ct <- sender
    workpath <- paste0(inputdir,'/runModel/')
    files = list.files(workpath)
    files = files[grep(paste0("LRTG_allscore_",ct),files)]  
    
    df_LRTGscore = lapply(files, function(file){
      
      print(file)
      LRS_score = readRDS(paste0(workpath,file))[[1]]
      LRS_score_merge = do.call('cbind',LRS_score)
      if (length(unique(colnames(LRS_score_merge))) == 1){
        LRpair <- unique(colnames(LRS_score_merge))
        LRS_score_merge = LRS_score_merge[,1] 
        
        # file <- gsub('-','_',file)
        df_LigRec <- data.frame(
          source = LRpair %>% gsub('_.*','',.),
          target = LRpair %>% gsub('.*_','',.),
          LRpair = LRpair,
          count = mean(LRS_score_merge),
          source_group = strsplit(file,'[_\\.]')[[1]][3],
          target_group = strsplit(file,'[_\\.]')[[1]][4])
        
      }else{
        LRS_score_merge = do.call('cbind',LRS_score) %>% .[,!duplicated(colnames(.))]
        
        # file <- gsub('-','_',file)
        df_LigRec <- data.frame(
          source = colnames(LRS_score_merge) %>% gsub('_.*','',.),
          target = colnames(LRS_score_merge) %>% gsub('.*_','',.),
          LRpair = colnames(LRS_score_merge),
          count = colMeans(LRS_score_merge),
          source_group = strsplit(file,'[_\\.]')[[1]][3],
          target_group = strsplit(file,'[_\\.]')[[1]][4])
      }
      
    }) %>% do.call('rbind',.)
    
    df_LRTGscore$pval <- NA
    for (lr in df_LRTGscore$LRpair){
      pos1 <- which(df_plot$LRpair %in% lr)
      pos2 <- which(df_LRTGscore$LRpair %in% lr)
      df_LRTGscore$pval[pos2] <- df_plot$pval[pos1]
    }
    
    df <- data.frame(cellpair = paste0(df_LRTGscore$source_group,"_",df_LRTGscore$target_group),
                     LRpair = df_LRTGscore$LRpair,
                     pval = df_LRTGscore$pval,
                     count = df_LRTGscore$count)
    
  }
  
  if(length(receiver) != 0 ){
    
    LRI_pval_Inter <- LRI_allpval[LRI_allpval$Receiver == receiver,]
    df_plot <- data.frame(cellpair = paste0(LRI_pval_Inter$Sender,"_",LRI_pval_Inter$Receiver),
                          LRpair = paste0(LRI_pval_Inter$source,"_",LRI_pval_Inter$target),
                          pval = LRI_pval_Inter$pval)
    
    df_plot <- df_plot[order(df_plot$pval,decreasing=F),] 
    
    # load LR signaling score
    ct <- receiver
    workpath <- paste0(inputdir,'/runModel/')
    files = list.files(workpath)
    files = files[grep(paste0(ct,".rds"),files)]  
    
    df_LRTGscore = lapply(files, function(file){
      
      print(file)
      LRS_score = readRDS(paste0(workpath,file))[[1]]
      LRS_score_merge = do.call('cbind',LRS_score)
      if (length(unique(colnames(LRS_score_merge))) == 1){
        LRpair <- unique(colnames(LRS_score_merge))
        LRS_score_merge = LRS_score_merge[,1] 
        
        # file <- gsub('-','_',file)
        df_LigRec <- data.frame(
          source = LRpair %>% gsub('_.*','',.),
          target = LRpair %>% gsub('.*_','',.),
          LRpair = LRpair,
          count = mean(LRS_score_merge),
          source_group = strsplit(file,'[_\\.]')[[1]][3],
          target_group = strsplit(file,'[_\\.]')[[1]][4])
        
      }else{
        LRS_score_merge = do.call('cbind',LRS_score) %>% .[,!duplicated(colnames(.))]
        
        # file <- gsub('-','_',file)
        df_LigRec <- data.frame(
          source = colnames(LRS_score_merge) %>% gsub('_.*','',.),
          target = colnames(LRS_score_merge) %>% gsub('.*_','',.),
          LRpair = colnames(LRS_score_merge),
          count = colMeans(LRS_score_merge),
          source_group = strsplit(file,'[_\\.]')[[1]][3],
          target_group = strsplit(file,'[_\\.]')[[1]][4])
      }
      
    }) %>% do.call('rbind',.)
    
    df_LRTGscore$pval <- NA
    for (lr in df_LRTGscore$LRpair){
      pos1 <- which(df_plot$LRpair %in% lr)
      pos2 <- which(df_LRTGscore$LRpair %in% lr)
      df_LRTGscore$pval[pos2] <- df_plot$pval[pos1]
    }
    
    df <- data.frame(cellpair = paste0(df_LRTGscore$source_group,"_",df_LRTGscore$target_group),
                     LRpair = df_LRTGscore$LRpair,
                     pval = df_LRTGscore$pval,
                     count = df_LRTGscore$count)
    
  }
  
  p1 <- ggplot(df, aes(x = cellpair, y = LRpair, color = pval, size = count)) +
    geom_point(pch = 16) +
    scale_color_gradient(low = "red", high = "yellow")+
    #theme_linedraw() + 
    #theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, size = 10,hjust= NULL, vjust = NULL),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_discrete(position = "bottom")
  print(p1)
  
  ## save
  
  pdf(paste0(plotdir,ct,'_bubble_pval.pdf'),height = p_height,width = p_width)
  print(p1)
  dev.off()
  
}

# Enrichment
Perf_ORA <- function(ls_tg,DB=c('GO','KEGG'))
{
  
  # ls_tg <- LRTG_pim$C5AR2$Target
  geneList <- ls_tg
  g2s <- toTable(org.Hs.egSYMBOL)
  geneList <- g2s$gene_id[match(geneList,g2s$symbol)]
  geneList <- na.omit(geneList)
  
  tryCatch({
    if(DB == 'GO'){
      
      res_enrich <- enrichGO(geneList, 'org.Hs.eg.db', ont = 'ALL', keyType = 'ENTREZID', 
                             minGSSize = 1, pvalueCutoff = 0.99)@result
      res_enrich$GeneRatio <- res_enrich$Count/length(enrich_gobp@gene)
      res_enrich <- res_enrich[,c("ID","Description","GeneRatio","p.adjust",'ONTOLOGY')]
      
    }else{
      
      res_enrich <- enrichKEGG(geneList, minGSSize = 1, pvalueCutoff = 0.99)@result
      res_enrich$GeneRatio <- res_enrich$Count/length(enrich_kegg@gene)
      res_enrich <- res_enrich[,c("ID","Description","GeneRatio","p.adjust")]
      res_enrich$ONTOLOGY <- 'KEGG'
      
    }
  }, error = function(e){
    res_enrich <- as.data.frame(matrix(
      ncol = 5,dimnames = list(c(),c("ID","Description","GeneRatio","p.adjust",'ONTOLOGY'))))
  })
  
  if(!exists('res_enrich')){
    res_enrich <- as.data.frame(matrix(
      ncol = 5,dimnames = list(c(),c("ID","Description","GeneRatio","p.adjust",'ONTOLOGY'))))
  }
  
  return(res_enrich)
}

Perf_GSEA <- function(df_tg,DB=c('KEGG','BIOCARTA'))
{
  
  # df_tg <- LRTG_pim[[4]]
  geneList <- df_tg$pIM
  names(geneList) = df_tg$Target
  geneList <- geneList[order(geneList, decreasing = T)]
  
  cp.gmt <- read.gmt("./vaild_scRNAseq/visualize_CCI/input/c2.cp.v7.4.symbols.gmt")
  cp.gmt$database <- as.character(cp.gmt$term) %>% lapply(., function(term){strsplit(term,"_")[[1]][1]}) %>% unlist() 
  cp.gmt.sel <- cp.gmt[cp.gmt$database == DB,]
  
  tryCatch(expr = {
    res_gsea <- GSEA(geneList, TERM2GENE=cp.gmt.sel, minGSSize = 1, pvalueCutoff = 0.99, verbose=FALSE, seed = 1)
  },error = function(e){
    res_gsea <- NA
  })
  if(!exists('res_gsea')){
    res_gsea <- NA
  }
  
  return(res_gsea)
}

Perf_Enrich <- function(ls_im,Type=c('ORA','GSEA'),DB=c('GO','KEGG'))
{
  
  # ls_im <- LRTG_im_spl
  Recs <- names(ls_im)
  if(Type == 'ORA'){
    
    res_enrich <- lapply(1:length(Recs), function(i){
      
      # rec <- 'ACVR2A'
      print(paste0(i,"/",length(Recs)))
      tgs <- ls_im[[i]]$Target
      
      res_ORA <- Perf_ORA(tgs,DB=DB)
      if(nrow(res_ORA)==0){
        res_ORA <- as.data.frame(matrix(
          ncol = 6,dimnames = list(c(),
                                   c("ID","Description","GeneRatio","p.adjust","ONTOLOGY","Regulator"))))
        res_ORA
      }else{
        res_ORA$Regulator <- Recs[i]
        res_ORA
      }
      
    }) %>% do.call('rbind',.)
    
  }else{
    
    res_enrich <- lapply(1:length(Recs), function(i){
      
      # i=4
      print(paste0(i,"/",length(Recs)))
      df_tg <- ls_im[[i]]
      
      res_GSEA <- Perf_GSEA(df_tg,DB=DB)
      
    })
    names(res_enrich) <- Recs
    
  }
  
  return(res_enrich)
}

####################
## other function ##
####################

## RCTD function: decompose a doublet into two cells
decompose_doublet_fast <- function(barcode){
  
  gene_list = myRCTD@internal_vars$gene_list_bulk
  cell_type_info = myRCTD@cell_type_info$renorm
  type1 = results$results_df[barcode, "first_type"]
  type2 = results$results_df[barcode, "second_type"]
  weights = results$weights_doublet[barcode, ]
  # counts <- as.matrix(myRCTD@spatialRNA@counts)
  # bead = counts[gene_list, barcode]
  bead = myRCTD@spatialRNA@counts[gene_list, barcode]
  N_genes = length(gene_list)
  epsilon = 1e-10
  
  de_result = lapply(1:N_genes, function(ind){
    gene = gene_list[ind]
    if(bead[ind] == 0){
      c(expect_1 = 0,
        expect_2 = 0,
        variance = 0) 
    }else{
      denom = weights[1] * cell_type_info[[1]][gene,type1] + 
        weights[2] * cell_type_info[[1]][gene,type2] + 
        epsilon
      posterior_1 = (weights[1] * cell_type_info[[1]][gene,type1] + epsilon / 2) / denom
      c(expect_1 = posterior_1 * bead[gene],
        expect_2 = bead[gene] - posterior_1 * bead[gene],
        variance = posterior_1 * bead[gene] * (1 - posterior_1)) 
    }
  })
  
  l1 = do.call('c', de_result)
  l2 = split(l1,stringr::str_extract(names(l1),"(expect_2)|(expect_1)|(variance)"))
  
  return(l2)
}

## RCTD function: remap celltypes
remap_celltypes <- function(cell_dict_file, cell_ident) {
  cell_type_dict <- read.csv(file=cell_dict_file, header=TRUE, sep=",")
  cell_type_dict$Name <- factor(cell_type_dict$Name)
  rownames(cell_type_dict) = cell_type_dict[,'Cluster']
  true_type_names = lapply(cell_ident, function(x) cell_type_dict[as.character(x),"Name"])
  true_type_names = unlist(true_type_names)
}

## RCTD function: create downsampled data for deconvolution
create_downsampled_data <- function(reference, refdir, cell_type_import = NULL,n_samples = 1000, each_cell_type = T,save.file = T) {
  if(!each_cell_type)
    index_keep = sample(which(reference@meta.data$liger_ident_coarse != ""),n_samples,replace=FALSE)
  else {
    if(is.null(cell_type_import))
      cell_type_import = levels(reference@meta.data$liger_ident_coarse)
    cell_type_import = cell_type_import[unlist(lapply(cell_type_import, function(x) nchar(x) > 0))]
    index_keep = c(); i = 1
    repeat{
      new_index = which(reference@meta.data$liger_ident_coarse == cell_type_import[i])
      new_samples = min(n_samples, length(new_index))
      index_keep = c(index_keep, sample(new_index,new_samples,replace=FALSE))
      if((i = i + 1) > length(cell_type_import))
        break
    }
  }
  reference@assays$RNA@counts = reference@assays$RNA@counts[,index_keep]
  reference@meta.data = reference@meta.data[index_keep,]
  reference@meta.data$liger_ident_coarse = droplevels(reference@meta.data$liger_ident_coarse)
  reference@assays$RNA@data <-matrix(0,nrow=2,ncol=2)
  reference@assays$RNA@scale.data <- matrix(0,nrow=2,ncol=2)
  if(save.file)
    saveRDS(reference, paste(refdir,"/scRefSubsampled", n_samples, ".RDS",sep=""))
  return(reference)
}

## load simulation data
prepare_input_data <- function(sheetID)
{
  
  exprMat <- openxlsx::read.xlsx("./data/ExpressionMatrix_100_slides.xlsx", 
                                 sheet = sheetID, colNames = F, startRow = 2)
  colnames(exprMat) <- paste0('Cell_',1:ncol(exprMat))
  rownames(exprMat) <- c(paste0('Lig',1:5),paste0('Rec',1:2),paste0('TF',1:3),paste0('TG',1:4))
  
  locaMat <- openxlsx::read.xlsx("./data/Lable_Coordinates_100_slides.xlsx", 
                                 sheet = sheetID, colNames = F, )
  locaMat <- locaMat[,-1]
  rownames(locaMat) <- colnames(exprMat)
  colnames(locaMat) <- c('dim_x','dim_y')
  
  annoMat <- openxlsx::read.xlsx("./data/ExpressionMatrix_100_slides.xlsx", 
                                 sheet = sheetID, colNames = F, rows = 1)
  annoMat <- rbind(colnames(exprMat),annoMat) %>% t() %>% as.data.frame()
  colnames(annoMat) <- c('Barcode','Cluster')
  rownames(annoMat) <- NULL
  annoMat$Cluster <- paste0("CT_",annoMat$Cluster)
  
  locaMat <- locaMat[!duplicated(locaMat),]
  annoMat <- annoMat[annoMat$Barcode %in% rownames(locaMat),]
  exprMat <- exprMat[colnames(exprMat) %in% rownames(locaMat)]
  
  res = list(exprMat=exprMat, locaMat=locaMat, annoMat=annoMat)
  return(res)
  
}
                                                      
#=======================================added by yll 2025-01-08================================
DrawBubblePlot <- function(respath,outputdir,ct,p_width,p_height){
  
  inputdir0 <- paste0(res_path,'getPIM/')
  files <- list.files(inputdir0)[grep('_im_',list.files(inputdir0))]
  CPs <- gsub('LRTG_im_clean_|LRTG_pim_clean_|\\.rds','',files)
  
  
  inputdir <- paste0(res_path,'runModel/')
  files = list.files(inputdir)
  files <- files[lapply(CPs, function(cp){grep(cp,files)}) %>% unlist() %>% unique()]
  files = files[grep(paste0('_',ct[1],'.rds',"|",'_',ct[2],'.rds'),files)]
  
  df_LRTGscore = lapply(files, function(file){
    
    print(file)
    LRS_score = readRDS(paste0(inputdir,file))[[1]]
    LRS_score_merge = do.call('cbind',LRS_score) %>% .[,!duplicated(colnames(.))]
    
    # file <- gsub('-','_',file)
    df_LigRec <- data.frame(
      source = colnames(LRS_score_merge) %>% gsub('_.*','',.),
      target = colnames(LRS_score_merge) %>% gsub('.*_','',.),
      LRpair = colnames(LRS_score_merge),
      count = colMeans(LRS_score_merge),
      source_group = strsplit(file,'[_\\.]')[[1]][3],
      target_group = strsplit(file,'[_\\.]')[[1]][4]
    )
    
  }) %>% do.call('rbind',.)
  df_LRTGscore$cellpair <- paste0(df_LRTGscore$source_group,"->",df_LRTGscore$target_group)
  df_LRTGscore$LRpair <- gsub("_","-",df_LRTGscore$LRpair)
  
  
  bubble_plot_LRscore(df_LRTGscore,outputdir,p_width,p_height,save_name=paste0(ct[1],"_",ct[2]))
  
}
                                                      
bubble_plot_LRscore <- function(res,wd_path,save_name){
  
  df_LRTGscore = res
  n.colors = 10
  color.heatmap = c("Spectral")  # ,"viridis"
  direction = -1
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
    })
  }else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }
  
  g <- ggplot(df_LRTGscore, aes(x = cellpair, y = LRpair,size=count,color = count)) +
    geom_point() +  # Color points by 'count'
    # scale_color_viridis(option="B") +  # Apply viridis color scale (optional)
    scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                           na.value = "white", limits = c(quantile(df_LRTGscore$count,0, na.rm = T), 
                                                          quantile(df_LRTGscore$count, 1, na.rm = T)), 
                           breaks = c(quantile(df_LRTGscore$count, 0, na.rm = T), 
                                      quantile(df_LRTGscore$count,1, na.rm = T)), labels = c("min", "max")) + 
    guides(color = guide_colourbar(barwidth = 0.5, title = "LR signaling score"))+
    theme_linedraw() + 
    theme(panel.grid.major = element_blank())+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+  # Rotate x-axis labels
    geom_vline(xintercept = seq(1.5, length(unique(df_LRTGscore$cellpair)) - 
                                  0.5, 1), lwd = 0.1, colour = "grey90")+
    geom_hline(yintercept = seq(1.5, length(unique(df_LRTGscore$LRpair)) - 
                                  0.5, 1), lwd = 0.1, colour = "grey90")
  
  plotdir = paste0(wd_path,'/visualize/')
  dir.create(plotdir,recursive = T)
  
  pdf(paste0(plotdir,"bubble_LRsignalingscore_",save_name,".pdf"),
      width = 6,height = 5.5)
  print(g)
  dev.off()
  
}

  DrawHeatmapNetworkPlot <- function(InputDir,OutputDir, Metric, mycolor_ct,p_width,p_height){
  suppressMessages(library(CellChat))
  suppressMessages(library(ComplexHeatmap))
  
  inputdir <- InputDir
  key <- Metric
  
  plotdir <- paste0(OutputDir,"/visualize/NetworkPlot/")
  dir.create(plotdir, recursive = T, showWarnings = F)
  
  files <- list.files(inputdir)[grep("LRTG_im_clean_", list.files(inputdir))]
  
  LRTG_detail <- lapply(files, function(f){
    
    LRTG_im <- readRDS(paste0(inputdir,f))
    LRTG_im <- na.omit(LRTG_im)
    c(length(unique(LRTG_im$LRpair)),length(unique(LRTG_im$Target)),
      sum(LRTG_im$IM),sum(LRTG_im$im_norm),
      sum(LRTG_im$IM)/nrow(LRTG_im),sum(LRTG_im$im_norm)/nrow(LRTG_im))
    
  }) %>% do.call('rbind',.) %>% as.data.frame()
  df_cellpair <- gsub('LRTG_im_clean_|.rds','',files) %>% strsplit(.,"_") %>% do.call('rbind',.) %>% as.data.frame()
  LRTG_detail <- cbind(df_cellpair,LRTG_detail)
  LRTG_detail <- na.omit(LRTG_detail)
  colnames(LRTG_detail) <- c('cell_from','cell_to','n_LRs','n_TGs','IM','IM_norm','mean_IM','mean_IM_norm')
  
  
  celltype = unique(c(LRTG_detail$cell_from,LRTG_detail$cell_to))
  LRTG_detail[[key]] <- as.numeric(LRTG_detail[[key]])
  
  tmeTab0 <- LRTG_detail[,c('cell_from','cell_to',key)] 
  
  mat2 <- matrix(0,nrow = length(celltype),ncol = length(celltype))
  rownames(mat2) <- celltype
  colnames(mat2) <- celltype
  for (i in 1:nrow(tmeTab0)){
    c1 <- tmeTab0$cell_from[i]
    c2 <- tmeTab0$cell_to[i]
    val <- tmeTab0$n_LRs[i]
    mat2[c1,c2] <- val
    
  }
  
  pt <- Heatmap_plot_Network(net = mat2,key=key,color.use=mycolor_ct[rownames(mat2)], color.heatmap = "Reds")
  
  pdf(paste0(plotdir,"NetwrokHeatmap_",key,".pdf"),height = p_height,width = p_width)
  print(pt)
  dev.off()
  
}
                                                      
Heatmap_plot_Network <- function(net = tmeTab,
                               key=NULL,
                               comparison = c(1, 2),
                               color.use = NULL,
                               color.heatmap = c("#2166ac", "#b2182b"),
                               title.name = NULL,
                               width = NULL,
                               height = NULL,
                               font.size = 8, 
                               font.size.title = 10,
                               cluster.rows = FALSE,
                               cluster.cols = FALSE,
                               sources.use = NULL,
                               targets.use = NULL,
                               remove.isolate = FALSE,
                               row.show = NULL,
                               col.show = NULL){
  
  if (key=='n_LRs'){
    legend.name <- "Number"
    title.name = paste0(key, " number")
  }else if(key=='IM'){
    legend.name <- "Score"
    title.name = paste0(key, " Score")
  }else if(key=='n_TGs'){
    legend.name <- "Number"
    title.name = paste0(key, " number")
  }else if(key=='IM_norm'){
    legend.name <- "Score"
    title.name = paste0(key, " Score")
  }else if(key=='mean_IM'){
    legend.name <- "Score"
    title.name = paste0(key, " Score")
  }else if(key=='mean_IM_norm'){
    legend.name <- "Score"
    title.name = paste0(key, " Score")
  }
  
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    if (length(idx) > 0) {
      net <- net[-idx, ]
      net <- net[, -idx]
    }
  }
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show]
    color.use <- color.use[col.show]
  }
  if (min(mat) < 0) {
    color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), 
                                   c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                      "\\1", min(mat, na.rm = T))) + 1), 0, round(max(mat, 
                                                                                                                      na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1",                                                                                                                                                 max(mat, na.rm = T))) + 1))
  }else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 1) {
      color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                                 name = color.heatmap))))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                      "\\1", min(mat, na.rm = T))) + 1), round(max(mat, 
                                                                                                                   na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1", 
                                                                                                                                                  max(mat, na.rm = T))) + 1))
  }
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), 
                                              border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                      show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), 
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                          show_annotation_name = FALSE)
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  }
  else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                name = legend.name, bottom_annotation = col_annotation, 
                left_annotation = row_annotation, top_annotation = ha2, 
                right_annotation = ha1, cluster_rows = cluster.rows, 
                cluster_columns = cluster.rows, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                column_names_gp = gpar(fontsize = font.size), column_title = title.name, 
                column_title_gp = gpar(fontsize = font.size.title), 
                column_names_rot = 90, row_title = "Sources (Sender)",
                row_title_gp = gpar(fontsize = font.size.title), row_title_rot = 90, 
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                            fontface = "plain"), title_position = "leftcenter-rot", 
                                            border = NA, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8), 
                                            grid_width = unit(2, "mm")))
  
  return(ht1)
}

