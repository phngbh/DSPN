library(tidyverse)
library(caret, lib.loc = "./Rlib/")


extract_corr = function(filepath,
                        phenotype = NULL,
                        omics = NULL,
                        type = c("cv","tetsing"),
                        cor = c("pearson","spearman")){
  ###Extract individual omics result and save as a df
  
  if(file.exists(filepath)){
    cat("Extracting ",omics,"\n")
    result = readRDS(filepath)
    
    if (type == "cv"){
      if( "corr" %in% names(result[[1]])){
        corr = list()
        #rmse = list()
        
        for(k in 1:length(result)){
          
          if (cor == "pearson"){
            corr[[k]] = result[[k]]$corr$pearson
          } else if (cor == "spearman"){
            corr[[k]] = result[[k]]$corr$spearman
          }
          
          #rmse[[k]] = sqrt(mean(result[[k]]$diff[,1]^2))
          
        }
        
        df = data.frame(phenotype = rep(phenotype, each = length(result)),
                        omics = rep(omics, each = length(result)),
                        corr_acc = unlist(corr), 
                        #rmse = unlist(rmse),
                        type = rep("regression", each = length(result))) %>%
          rownames_to_column("cv") 
        
      } else if ("acc" %in% names(result[[1]])){
        acc = list()
        #rmse = list()
        
        for(k in 1:length(result)){
          
          acc[[k]] = result[[k]]$acc
          #rmse[[k]] <- NULL
          
        }
        
        df = data.frame(phenotype = rep(phenotype, each = length(result)),
                        omics = rep(omics, each = length(result)),
                        corr_acc = unlist(acc), 
                        #rmse = unlist(rmse),
                        type = rep("classification",each = length(result))) %>%
          rownames_to_column("cv") 
      }
      
    } else if (type =="testing"){
      if ("corr" %in% names(result)){
        
        if (cor == "pearson"){
          corr = result$corr$pearson
        } else if (cor =="spearman"){
          corr = result$corr$spearman
        }
        
        #rmse = sqrt(mean(result$diff[,1]^2))
        df = data.frame(phenotype = phenotype,
                        omics = omics, 
                        corr_acc = corr, 
                        #rmse = rmse,
                        type = "regression") %>% mutate(cv = "None") 
      } else if ("acc" %in% names(result)){
        acc = result$acc
        #rmse = NULL
        df = data.frame(phenotype = phenotype, 
                        omics = omics, 
                        corr_acc = acc, 
                        #rmse = rmse,
                        type = "classification") %>% mutate(cv = "None")
      }
    }
    
    return(df)
  } else {
    df <- NULL
    return(df)
  }
  
}

make_corr_df = function(dir,
                        omics = c("Genomics","Transcriptomics","Proteomics","Metabolomics","Methylomics","Olink"),
                        phenotype = NULL,
                        type = c("cv","testing"),
                        cor = c("pearson","spearman")){
  ###Extract all results and save as a df
  
  comb_list = list()
  for (i in 1:length(omics)){
    tmp = combn(omics,i, simplify = F)
    comb_list = c(comb_list, tmp)
  } 
  comb_list = lapply(comb_list, function(x) paste(x, collapse = ''))
  
  res_list = list()
  if (type == "cv"){
    for(p in phenotype){
      cat("Extracting ",p,"\n")
      res_list[[p]] = lapply(comb_list, function(x) extract_corr(paste0(dir,x,",",p,",cv.rds"),
                                                                 phenotype = p,
                                                                 omics = x,
                                                                 type = type,
                                                                 cor = cor))
      res_list[[p]] = do.call("rbind", res_list[[p]])
    }
  } else if (type == "testing"){
    for(p in phenotype){
      cat("Extracting ",p,"\n")
      res_list[[p]] = lapply(comb_list, function(x) extract_corr(paste0(dir,x,",",p,",mod.rds"),
                                                                 phenotype = p,
                                                                 omics = x,
                                                                 type = type,
                                                                 cor = cor))
      res_list[[p]] = do.call("rbind", res_list[[p]])
    }
  }
  
  res_df = do.call("rbind",res_list)
  return(res_df)
}

extract_pred = function(filepath,
                        phenotype = NULL,
                        omics = NULL){
  ###Extract individual omics result and save as a df
  
  if(file.exists(filepath)){
    cat("Extracting ",omics,"\n")
    result = readRDS(filepath)
    
    if("corr" %in% names(result)){
      pred = result$pred
      obs = pred - result$diff[,1]
      
      df = data.frame(sample = rownames(result$diff),
                      omics = rep(omics, each = length(pred)),
                      phenotype = rep(phenotype, each = length(pred)),
                      prediction = pred,
                      observation = obs,
                      type = rep("regression", each = length(pred))) 
      
    } else if ("acc" %in% names(result)){
      df = data.frame(sample = rownames(result$obs),
                      omics = rep(omics, each = length(result$pred)),
                      phenotype = rep(phenotype, each = length(result$pred)),
                      prediction = as.numeric(as.character(result$pred)),
                      observation = as.numeric(as.character(result$obs[,1])),
                      type = rep("classification", each = length(result$pred))) 
    }
    
    return(df)
  } else{
    df <- NULL
    return(df)
  }
  
}

make_pred_df = function(dir,
                        omics = c("Genomics","Transcriptomics","Proteomics","Metabolomics","Methylomics","Olink"),
                        phenotype = NULL){
  ###Extract all results and save as a df
  
  comb_list = list()
  for (i in 1:length(omics)){
    tmp = combn(omics,i, simplify = F)
    comb_list = c(comb_list, tmp)
  } 
  comb_list = lapply(comb_list, function(x) paste(x, collapse = ''))
  
  res_list = list()
  for(p in phenotype){
    cat("Extracting ",p,"\n")
    res_list[[p]] = lapply(comb_list, function(x) extract_pred(paste0(dir,x,",",p,",mod.rds"),
                                                               phenotype = p,
                                                               omics = x))
    res_list[[p]] = do.call("rbind", res_list[[p]])
  }
  
  res_df = do.call("rbind",res_list)
  return(res_df)
}

extract_features = function(filepath,
                            phenotype = NULL,
                            omics = NULL){
  ### Extract feature importance 
  
  cat("Extracting ",omics,"\n")
  result = readRDS(filepath)
  
  fit = result$fit
  
  if("corr" %in% names(result)){
    imp = varImp(fit)$importance %>% rownames_to_column("Feature") %>% 
      filter(Overall != 0) %>% mutate(Phenotype = phenotype, Omics = omics)
  } else if ("acc" %in% names(result)){
    imp = varImp(fit)$importance
    rn = rownames(imp)
    imp_med = apply(imp,1, function(x) mean(x, na.rm = T))
    names(imp_med) = rn
    imp_med = imp_med[imp_med > 0]
    imp = data.frame(Feature = names(imp_med), Overall = imp_med, 
                     Phenotype = rep(phenotype, each = length(imp_med)),
                     Omics = rep(omics, each = length(imp_med)))
  }
  
  
  return(imp)
  
}

extract_features_rf = function(filepath,
                            phenotype = NULL,
                            omics = NULL){
  ### Extract feature importance 
  
  cat("Extracting ",omics,"\n")
  result = readRDS(filepath)
  
  fit = result$fit
  imp = varImp(fit)$importance %>% rownames_to_column("Feature") %>% 
    filter(Overall != 0) %>% mutate(Phenotype = phenotype, Omics = omics)
  
  
  return(imp)
  
}

make_feature_df = function(dir,
                           omics = c("Genomics","Transcriptomics","Proteomics","Metabolomics","Methylomics","Olink"),
                           phenotype = NULL){
  ###Extract all results and save as a df
  
  comb_list = list()
  for (i in 1:length(omics)){
    tmp = combn(omics,i, simplify = F)
    comb_list = c(comb_list, tmp)
  } 
  comb_list = lapply(comb_list, function(x) paste(x, collapse = ''))
  
  res_list = list()
  for(p in phenotype){
    cat("Extracting ",p,"\n")
    res_list[[p]] = lapply(comb_list, function(x) extract_features(paste0(dir,x,",",p,",mod.rds"),
                                                                   phenotype = p,
                                                                   omics = x))
    res_list[[p]] = do.call("rbind", res_list[[p]])
  }
  
  res_df = do.call("rbind",res_list)
  return(res_df)
}

calc_mnsi = function(mat,age){
  mat = na.omit(mat)
  age = age[rownames(mat)]
  mnsi = vector(mode = "numeric", length = nrow(mat))
  for (i in 1:nrow(mat)){
    
    app_r = ifelse(any(mat[i,paste0("un12_",c(1,2,5:7))] == 1), 1, 0)
    app_l = ifelse(any(mat[i,paste0("un13_",c(1,2,5:7))] == 1), 1, 0)
    app = app_r + app_l
    
    ank_r = ifelse(mat[i,"un19"] == 0,0,ifelse(mat[i,"un19"] == 1,0.5,1))
    ank_l = ifelse(mat[i,"un20"] == 0,0,ifelse(mat[i,"un20"] == 1,0.5,1))
    ank = ank_r + ank_l
    
    vib_tmp_r = ifelse(mat[i,"un17"] < 5.75 - 0.026*age[i], 1, 0)
    if (vib_tmp_r == 0) {vib_r <- 0}
    if (vib_tmp_r == 1 & mat[i,"un17"] > 0) {vib_r <- 0.5}
    if (mat[i,"un17"] == 0) {vib_r <- 1}
    vib_tmp_l = ifelse(mat[i,"un18"] < 5.75 - 0.026*age[i], 1, 0)
    if (vib_tmp_l == 0) vib_l <- 0
    if (vib_tmp_l == 1 & mat[i,"un18"] > 0) vib_l <- 0.5
    if (mat[i,"un18"] == 0) vib_l <- 1
    vib = vib_r + vib_l
    
    mon_r = ifelse(mat[i,"un14"] == 0,1,ifelse(mat[i,"un14"] > 0 & mat[i,"un14"] <= 7,0.5,0))
    mon_l = ifelse(mat[i,"un15"] == 0,1,ifelse(mat[i,"un15"] > 0 & mat[i,"un15"] <= 7,0.5,0))
    mon = mon_r + mon_l
    
    mnsi[i] = app + ank + vib + mon
    rm(app_r,app_l,app,ank_r,ank_l,ank,vib_tmp_r,vib_tmp_l,vib_r,vib_l,vib,mon_r,mon_l,mon)
  }
  
  names(mnsi) = rownames(mat)
  return(mnsi)  
  
}

mnsi_corr = function(pred = NULL,
                     info = NULL
                     ){
  
  omics_comb = levels(as.factor(pred$omics))
  corr = list()
  for (o in omics_comb){
    pred_o = filter(pred, omics == o) %>% 
      select(sample,phenotype,prediction) %>% 
      pivot_wider(names_from = "phenotype", values_from = "prediction") %>%
      column_to_rownames("sample")
    info = info[rownames(pred_o),]
    age = info$utalter
    names(age) = rownames(info)
    mnsi_p = calc_mnsi(pred_o,age)
    mnsi_o = info$mnsi
    names(mnsi_o) = rownames(info)
    mnsi_o = mnsi_o[names(mnsi_p)]
    
    corr[[o]] = cor.test(mnsi_p,mnsi_o, method = "pearson", exact = F)$estimate
  }
  
  return(corr)
  
}

