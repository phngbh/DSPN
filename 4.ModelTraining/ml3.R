library(tidyverse)
library(snpStats, lib.loc = "./Rlib/")
library(caret, lib.loc = "./Rlib/")
library(caretEnsemble, lib.loc = "./Rlib/")
library(haven)
library(limma)
library(ChAMP)
library(glmnet)
library(pROC, lib.loc = "./Rlib/")
library(glmnetUtils, lib.loc = "./Rlib/")
library(DMwR, lib.loc = "Rlib")
library(MLmetrics, lib.loc = "Rlib")
library(randomForest)
library(kernlab, lib.loc = "./Rlib/")
library(nnls, lib.loc = "Rlib")
#library(SuperLearner, lib.loc = "Rlib")
library(caTools)

myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

fil_col = function(
  ###Filter uniform variables
  mat = NULL, 
  type = "remove"
){
  if(type=="remove"){
    keep = apply(mat, 2, function(x) nlevels(as.factor(x)) >= 2)
    mat = mat[,keep, drop=F]
  }
  
  if(type=="convert"){
    conv = apply(mat, 2, function(x) nlevels(as.factor(x)) < 2)
    tmp = mat[,conv] %>% apply(.,2,as.numeric)
    cn = colnames(tmp)
    for( i in cn){
      mat[,i] = as.numeric(mat[,i])
    }
  }
  
  return(mat)
}

bed_to_df = function(
  ###Load in plink bed file and make df
  bed = NULL #path to bed file
){
  
  dat = read.plink(bed = bed)
  dn = dimnames(dat$genotypes)
  df = as.data.frame(dat$genotypes)
  df.char = apply(df,2,as.character)
  df.char.fil = fil_col(mat = df.char) %>% as.data.frame %>% model.matrix(~ 0 + .,.)
  row.names(df.char.fil) <- dn[[1]]
  return(df.char.fil)
  
}

bed_to_dimnames = function(
  ###Extract dimnames from bed file
  bed = NULL #path to bed file
){
  dat = read.plink(bed = bed)
  dn = dimnames(dat$genotypes)
  names(dn) = c("samples","features")
  return(dn)
}

mat_to_dimnames = function(
  ###Extract dimnames from a list of (omics) matrices
  list = NULL #name of omics list
){
  omics = list()
  for (o in names(list)){
    omics[[o]] = dimnames(list[[o]])
    names(omics[[o]]) = c("samples","features")
  }
  return(omics)
}

feature_sel = function(
  ###Select most relevant and least redundant features
  feature_mat = NULL,
  phenotype = NULL,
  o = NULL, #String, name of omics comb
  p = NULL # name of the outcome
  
){

  feature_mat = as.data.frame(feature_mat)

  if (p %in% c("un14","un15","un17","un18","utmnsi","prog")){
    corr = sapply(feature_mat, function(x) cor.test(x,phenotype, method = "pearson",exact = F)$estimate) %>%
      abs() %>% sort(decreasing = T)
    names(corr) = gsub(".cor","",names(corr)) %>%
      gsub("Covar.","",.) %>%
      gsub(paste0(o,"."),"",.)
    corr_top = head(corr, n=1000)
    feature_mat_fil = feature_mat[,names(corr_top)]

    ft.cor = cor(feature_mat_fil)
    hc = findCorrelation(ft.cor, cutoff = 0.9)
    if (length(hc) > 0){
      feature_mat_fil = feature_mat_fil[,-hc]
      feature_mat_final = cbind(feature_mat_fil,feature_mat[,names(corr)[1001:(1000+length(hc))],drop = F])
    } else {
      feature_mat_final = feature_mat_fil
    }

  } else {
    
    if(nlevels(as.factor(phenotype)) > 2){
      phenotype[phenotype == 2] <- 1
    }
    
    phenotype = as.character(phenotype)
    
    if(o == "Transcriptomics"){
      mod = model.matrix(~ 0 + phenotype, data.frame(phenotype = phenotype))
      fit = lmFit(t(feature_mat),mod)
      contrast = makeContrasts(phenotype1 - phenotype0, levels = colnames(coef(fit)))
      tmp <- contrasts.fit(fit, contrast)
      tmp <- eBayes(tmp)
      topde <- topTable(tmp, sort.by = "P", n = 2000)
      feature_mat_fil = feature_mat[,rownames(topde)[1:1000]]
      
      ft.cor = cor(feature_mat_fil)
      hc = findCorrelation(ft.cor, cutoff = 0.9)
      if (length(hc) > 0){
        feature_mat_fil = feature_mat_fil[,-hc]
        feature_mat_final = cbind(feature_mat_fil,feature_mat[,rownames(topde)[1001:(1000+length(hc))],drop = F])
      } else {
        feature_mat_final = feature_mat_fil
      }

    } 
    
    if(o == "Methylomics"){
      
      try = myTryCatch(champ.DMP(beta = t(feature_mat),pheno=phenotype, adjPVal = 0.3))
      
      if( !is.null(try$error)){
        feature_mat_final <- NULL
        
      } else {
        myDMP <- champ.DMP(beta = t(feature_mat),pheno=phenotype, adjPVal = 0.3)
        topdm = myDMP[[1]] %>% rownames_to_column("probe") %>% arrange(.,adj.P.Val) %>% column_to_rownames("probe") %>% rownames()
        
        if (length(topdm) <= 1000){
          feature_mat_fil = feature_mat[,topdm]
          ft.cor = cor(feature_mat_fil)
          hc = findCorrelation(ft.cor, cutoff = 0.9)
          if (length(hc) > 0){
            feature_mat_final = feature_mat_fil[,-hc]
          } else {
            feature_mat_final = feature_mat_fil
          }
          
        } else{
          feature_mat_fil = feature_mat[,topdm[1:1000]]
          
          ft.cor = cor(feature_mat_fil)
          hc = findCorrelation(ft.cor, cutoff = 0.9)
          if (length(hc) > 0){
            feature_mat_fil = feature_mat_fil[,-hc]
            if (length(topdm) >= 1000 + length(hc)){
              feature_mat_final = cbind(feature_mat_fil,feature_mat[,topdm[1001:(1000+length(hc))],drop = F])
            } else {
              feature_mat_final = cbind(feature_mat_fil,feature_mat[,topdm[1001:length(topdm)],drop = F])
            }
            
          } else {
            feature_mat_final = feature_mat_fil
          }
        }
      }
  
    }

  }
  
  return(feature_mat_final)
  
}

dimnames_to_mat = function(
  ###Retrieve feature matrix from list of dimnames
  feature_list = NULL,
  sample_list = NULL,
  genomics = NULL,
  phenotype = NULL,
  omics = NULL,
  covar = NULL,
  lab = NULL,
  features = NULL, #Character vector of features to select (only provide if feature_sel=F)
  feature_sel = T, #Whether to run feature_sel()
  p = NULL #name of the outcome
){
  
  mat_list = vector(mode = "list", length = length(feature_list))
  names(mat_list) <- names(feature_list)
  
  
  for(i in names(feature_list)) {
    if (i == "Covar"){
      mat = covar[sample_list[["Clinical"]], feature_list[["Covar"]]]
    
    } else if (i == "Lab") {
      mat = lab[sample_list[["Clinical"]], feature_list[["Lab"]]]
      
    } else if(i== "Genomics") {
      mat = genomics[sample_list[["Genomics"]],]
      
    } else if (i %in% c("Transcriptomics","Methylomics")) {
      if (feature_sel){
        ftr_mat = omics[[i]][sample_list[[i]],feature_list[[i]]]
        mat = feature_sel(feature_mat = ftr_mat, phenotype = phenotype[,1], o = i, p = p)
        rm(ftr_mat)
      } else{
        mat = omics[[i]][sample_list[[i]],intersect(feature_list[[i]],features)]
      }
      
    } else {
      mat = omics[[i]][sample_list[[i]],feature_list[[i]]]
    }
    mat_list[[i]] <- mat
  }
  
  names(mat_list) <- NULL
  mat = do.call("cbind", mat_list)
  rownames(mat) = sample_list[["Clinical"]]
  return(mat)
}

match_samples_features = function(
  ###Return a list of samples with matched IDs and features for each omics combination
  omics = NULL, #a list of dimnames of multiple omics
  genomics = genomics, #a list of dimnames of genomics data 
  covar = NULL, #dataframe of dummy covariates
  lab = NULL, #df of laboratory measurements
  omics_id = NULL, #dataframe of matched IDs for omics and clinical data
  omics_type = NULL #a character string or character vector of omics names 
){
  
  omics_full = c(list(Genomics = genomics),omics)
  clin_id = intersect(rownames(covar), rownames(lab))
  covar = covar[clin_id,]
  lab = lab[clin_id,]
  
  omics_id = omics_id[,c(omics_type,"Clinical"), drop=FALSE]
  for (i in omics_type){
    ind = match(as.character(as.data.frame(omics_id)[,i]),omics_full[[i]][["samples"]])
    omics_id[,i] = ind
  }
  omics_id = na.omit(omics_id) %>% filter(Clinical %in% as.numeric(rownames(covar)))
  covar = covar[as.character(omics_id$Clinical),]

  omics_sample = lapply(omics_type, function(x)
    return(omics_full[[x]][["samples"]][as.data.frame(omics_id)[,x]]))
  names(omics_sample) = omics_type
  omics_sample = c(omics_sample, list(Clinical = as.character(omics_id$Clinical)))
  
  omics_feature = lapply(omics_type, function(x)
    return(omics_full[[x]][["features"]]))
  names(omics_feature) = omics_type
  
  data_list = list(Features = c(omics_feature,list(Covar = colnames(covar), Lab = colnames(lab))), 
                   Samples = omics_sample)
  
  return(data_list)
}

make_comb_list = function(
  ###Return a list of combinations of omics, each includes matched samples and features
  omics_list = c("Genomics","Transcriptomics","Proteomics",
                 "Metabolomics","Methylomics","Olink"),
  omics = NULL, #a list of multiple omics
  covar = NULL, #dataframe of covariates
  lab = NULL, #df of laboratory measurements
  omics_id = NULL #dataframe of matched IDs for omics and clinical data
){

  omics = mat_to_dimnames(omics)
  
  comb_list = list()
  for (i in 1:length(omics_list)){
    tmp = combn(omics_list,i, simplify = F)
    comb_list = c(comb_list, tmp)
  }
  
  genomics = bed_to_dimnames("./small_bed/uh_adipo.bed")
  genomics[["features"]] <- NULL
  
  data_list = lapply(comb_list, function(x) match_samples_features(omics = omics,
                                                                    genomics = genomics, 
                                                                    covar = covar, 
                                                                    lab = lab,
                                                                    omics_id = omics_id, 
                                                                    omics_type = x))
  name_comb = lapply(comb_list, function(x) paste(x, collapse = '')) %>% unlist()
  names(data_list) = name_comb
  
  return(data_list)
}

make_train_test = function(
  ###Split 80% data for training and 20% for testing
  ###and return a list of the two subsets
  omics_list = NULL, #list of matched samples and features for each omics comb 
  seed = 993
){
  
  data_list = list(train = list(), test = list())
  
  for (i in names(omics_list)){
    for(j in names(omics_list[[i]][["Samples"]])){
      
      samples = omics_list[[i]][["Samples"]][[j]]
      train_size = floor(0.8*length(samples))
      
      set.seed(seed)
      train_ind = sample(length(samples), size = train_size)
      
      train_set = samples[train_ind]
      test_set = samples[-train_ind]
      
      data_list[["train"]][[i]][[j]] = train_set
      data_list[["test"]][[i]][[j]] = test_set
      
    }
  }
  
  return(data_list)
}

extract_train_test = function(
  ###Return a train or test set based on predefined train/test IDs
  omics_list = NULL, # list of #list of matched samples and features for each omics comb
  split_list = NULL, # list of train and test IDs
  train = T # whether train or test set
){
  
  if(train){
    for(c in names(omics_list)){
      for(o in names(omics_list[[c]][["Samples"]])){
        omics_list[[c]][["Samples"]][[o]] <- split_list[["train"]][[c]][[o]]
      }
    }
  } else{
    for(c in names(omics_list)){
      for(o in names(omics_list[[c]][["Samples"]])){
        omics_list[[c]][["Samples"]][[o]] <- split_list[["test"]][[c]][[o]]
      }
    }
  }
  
  return(omics_list)
}

make_cv_list <- function(
  ###Create a list of cross-validation subsets for each omics combination
  omics_list = NULL, #list of training data 
  kfold = 5, #how many folds
  times = 2, #number of k-fold cross validations 
  seed = c(993,999) #setting a fixed seed
){
  
  cv_list = list()

  for (c in names(omics_list)){
    
    cv_list[[c]] = list()
      
    samples = omics_list[[c]][["Samples"]][["Clinical"]]
    
    fold_ids <- cut(1:length(samples), breaks = kfold, labels = F)
    
    size <- floor(((kfold-1)/kfold)*length(samples))
    print(paste0("Fold size is ",as.character(size)))
    print(paste0("Test size is ", as.character(length(samples)-size)))
    
    cv_list[[c]] = list(train = list(), test = list())
    
    for (j in 1:length(seed)){
      
      set.seed(seed[j])
      samples <- samples[sample(length(samples))]
      
      test_id <- lapply(1:kfold, function(x) return(samples[fold_ids == x]))
      train_id <- lapply(1:kfold, function(x) return(samples[fold_ids != x]))
      
      cv_list[[c]][["train"]] = c(cv_list[[c]][["train"]], train_id)
      cv_list[[c]][["test"]] = c(cv_list[[c]][["test"]], test_id)
        
    }
    
  }
  return(cv_list)
}

use_rf <- function(
  ### Train random forest
  x_train=x_train,
  y_train=y_train,
  x_test=x_test,
  y_test=y_test,
  #hyperparam=hyperparam,
  seed = 993,
  p = p, #name of the outcome
  resampling = NULL
){
  
  
  if (p %in% c("un14","un15","un17","un18","utmnsi", "prog")){
    
    control <- trainControl(method="none")
    
    set.seed(seed)
    tune = tuneRF(x_train, y_train[,1], stepFactor = 1.5, improve = 1e-5, tree = 500,trace = F, plot = F) %>% as.data.frame()
    mtry = tune$mtry[which(tune$OOBError == min(tune$OOBError))] %>% as.numeric() %>% max()
    tunegrid <- expand.grid(.mtry=mtry)
    
    set.seed(seed)
    fit <- caret::train(x = x_train,
                        y = as.double(y_train),
                        method="rf", 
                        metric="RMSE",
                        #tuneLength = 20,
                        tuneGrid=tunegrid, 
                        trControl=control,
                        importance = TRUE)
    
    pred <- predict(fit, x_test)
    pred <- as.data.frame(pred)[,1]
    rmse <- sqrt(sum((pred - y_test)^2))
    pearson <- cor.test(pred,y_test, method = "pearson", exact = F)$estimate
    spearman <- cor.test(pred,y_test, method = "spearman", exact = F)$estimate
    perf = list(rmse = rmse, pearson = pearson, spearman = spearman)
    
    return(list(pred=pred, obs = y_test, perf = perf, fit=fit))
  } 
  
  if (!p %in% c("un14","un15","un17","un18","utmnsi","prog")){
    
    if (nlevels(as.factor(y_train)) == 2 & nlevels(as.factor(y_test)) == 2) {
      
      y_train = ifelse(y_train == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
      y_test = ifelse(y_test == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
      
      if (p == "utmnsi2"){
        perc.over = 300
      } else if (p == "utmnsi3"){
        perc.over = 800
      } else if (p == "inc3"){
        perc.over = 400
      } else {
        perc.over = 200
      }
      
      smotest <- list(name = "SMOTE with more upsampling",
                      func = function (x, y) {
                        library(DMwR)
                        dat <- if (is.data.frame(x)) x else as.data.frame(x)
                        dat$.y <- y
                        dat <- SMOTE(.y ~ ., data = dat, perc.over = perc.over, perc.under = 100)
                        list(x = dat[, !grepl(".y", colnames(dat), fixed = TRUE)], 
                             y = dat$.y)
                      },
                      first = TRUE)
      
      if (resampling == "up"){
        control <- trainControl(method="repeatedcv", number=5, repeats=2, summaryFunction = twoClassSummary, classProbs = TRUE, sampling = "up")
      } else if (resampling == "smote"){
        control <- trainControl(method="repeatedcv", number=5, repeats=2, summaryFunction = twoClassSummary, classProbs = TRUE, sampling = smotest)
      }
      
      # weights <- ifelse(y_train[,1] == "One",
      #                   (1/table(as.factor(y_train[,1]))[1]) * 0.5,
      #                   (1/table(as.factor(y_train[,1]))[2]) * 0.5)
      
      set.seed(seed)
      tune = tuneRF(x_train, y_train, stepFactor = 1.5, improve = 1e-5, tree = 500,trace = F, plot = F) %>% as.data.frame()
      mtry = tune$mtry[which(tune$OOBError == min(tune$OOBError))] %>% as.numeric() %>% max()
      tunegrid <- expand.grid(.mtry=mtry)
      
      set.seed(seed)
      fit <- caret::train(x = x_train,
                          y = y_train,
                          method="rf", 
                          metric="ROC",
                          #tuneLength = 20,
                          tuneGrid=tunegrid, 
                          trControl=control,
                          importance = TRUE)
      
      pred <- predict(fit, x_test, s = 'lambda.min')
      pred.prob <- predict(fit, x_test, s = 'lambda.min', type = "prob")
      roc <- roc(response = y_test, predictor = pred.prob$One, levels = c("Zero","One"))
      auc = auc(roc)
      opt.thres = coords(roc, x = "best", ret = "threshold", best.method = "closest.topleft", transpose = T)
      obs = y_test 
      positive =  ifelse(pred.prob$One > opt.thres,1,0) 
      sens = sum(positive)/length(positive) #caret::sensitivity(pred,obs)
      negative = ifelse(pred.prob$Zero < opt.thres,1,0) 
      spec = sum(negative)/length(negative)  #caret::specificity(pred,obs)
      perf = list(roc = roc, opt.thres = opt.thres, auc = auc, sens = sens, spec = spec)
      
      return(list(pred=pred, obs=y_test, perf = perf, fit=fit))
      
    } else {
      
      message("Labels of train or test set has only one level")
      res = NULL
      return(res)
      
    }
  }

}

use_glm <- function(
  ###Train a linear model for each omics combination and phenotype 
  x_train=NULL,
  y_train=NULL,
  x_test=NULL,
  y_test=NULL,
  seed = 993,
  p = NULL, # Name of the outcome
  resampling = NULL
){
  
  
  if(p %in% c("un14","un15","un17","un18","utmnsi","prog")){
  
    control <- trainControl(method="repeatedcv", number=5, repeats=2)
    set.seed(seed)
    fit <- caret::train(x = x_train, 
                        y = as.double(y_train), 
                        method = "glmnet",
                        metric = "RMSE",
                        tuneLength = 20, 
                        trControl = control,
                        na.action = na.omit)
    
    pred <- predict(fit, x_test, s = 'lambda.min')
    rmse <- sqrt(sum((pred - y_test)^2))
    pearson <- cor.test(pred,y_test, method = "pearson", exact = F)$estimate
    spearman <- cor.test(pred,y_test, method = "spearman", exact = F)$estimate
    perf = list(rmse = rmse, pearson = pearson, spearman = spearman)
    
    return(list(pred=pred, obs = y_test, perf = perf, fit=fit))
    
  }
  
  if(!p %in% c("un14","un15","un17","un18","utmnsi","prog")){
    
    if (nlevels(as.factor(y_train)) == 2 & nlevels(as.factor(y_test)) == 2) {
      
      y_train = ifelse(y_train == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
      y_test = ifelse(y_test == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
      
      if (p == "utmnsi2"){
        perc.over = 300
      } else if (p == "utmnsi3"){
        perc.over = 800
      } else if (p == "inc3"){
        perc.over = 400
      } else {
        perc.over = 200
      }
      
      smotest <- list(name = "SMOTE with more upsampling",
                      func = function (x, y) {
                        library(DMwR)
                        dat <- if (is.data.frame(x)) x else as.data.frame(x)
                        dat$.y <- y
                        dat <- SMOTE(.y ~ ., data = dat, perc.over = perc.over, perc.under = 100)
                        list(x = dat[, !grepl(".y", colnames(dat), fixed = TRUE)], 
                             y = dat$.y)
                      },
                      first = TRUE)
      
      if (resampling == "up"){
        control <- trainControl(method="repeatedcv", number=5, repeats=2, summaryFunction = twoClassSummary, classProbs = TRUE, sampling = "up")
      } else if (resampling == "smote"){
        control <- trainControl(method="repeatedcv", number=5, repeats=2, summaryFunction = twoClassSummary, classProbs = TRUE, sampling = smotest)
      }
      
      # weights <- ifelse(y_train[,1] == "One",
      #                   (1/table(as.factor(y_train[,1]))[1]) * 0.5,
      #                   (1/table(as.factor(y_train[,1]))[2]) * 0.5)
      set.seed(seed)
      fit <- caret::train(x = x_train,
                          y = y_train,
                          method = "glmnet",
                          metric = "ROC",
                          #weights = weights,
                          tuneLength = 20,
                          trControl = control,
                          na.action = na.omit)
      
      pred <- predict(fit, x_test, s = 'lambda.min')
      pred.prob <- predict(fit, x_test, s = 'lambda.min', type = "prob")
      roc <- roc(response = y_test, predictor = pred.prob$One, levels = c("Zero","One"))
      auc = auc(roc)
      opt.thres = coords(roc, x = "best", ret = "threshold", best.method = "closest.topleft", transpose = T)
      obs = y_test 
      positive =  ifelse(pred.prob$One > opt.thres,1,0) 
      sens = sum(positive)/length(positive) #caret::sensitivity(pred,obs)
      negative = ifelse(pred.prob$Zero < opt.thres,1,0) 
      spec = sum(negative)/length(negative)  #caret::specificity(pred,obs)
      perf = list(roc = roc, opt.thres = opt.thres, auc = auc, sens = sens, spec = spec)
      
      return(list(pred=pred, obs=y_test, perf = perf, fit=fit))
      
    } else {
      
      message("Labels of train or test set has only one level")
      res = NULL
      return(res)
      
    }
    
  }
  
}

use_svm <- function(
  ###Train a linear model for each omics combination and phenotype 
  x_train=NULL,
  y_train=NULL,
  x_test=NULL,
  y_test=NULL,
  seed = 993,
  p = NULL, # Name of the outcome
  resampling = NULL
){
  
  
  if(p %in% c("un14","un15","un17","un18","utmnsi","prog")){
    
    control <- trainControl(method="repeatedcv", number=5, repeats=2)
    set.seed(seed)
    fit <- caret::train(x = x_train, 
                        y = as.double(y_train), 
                        method = "svmLinear",
                        metric = "RMSE",
                        tuneLength = 20, 
                        trControl = control,
                        na.action = na.omit)
    
    pred <- predict(fit, x_test, s = 'lambda.min')
    rmse <- sqrt(sum((pred - y_test)^2))
    pearson <- cor.test(pred,y_test, method = "pearson", exact = F)$estimate
    spearman <- cor.test(pred,y_test, method = "spearman", exact = F)$estimate
    perf = list(rmse = rmse, pearson = pearson, spearman = spearman)
    
    return(list(pred=pred, obs = y_test, perf = perf, fit=fit))
    
  }
  
  if(!p %in% c("un14","un15","un17","un18","utmnsi","prog")){
    
    if (nlevels(as.factor(y_train)) == 2 & nlevels(as.factor(y_test)) == 2) {
      
      y_train = ifelse(y_train == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
      y_test = ifelse(y_test == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
      
      if (p == "utmnsi2"){
        perc.over = 300
      } else if (p == "utmnsi3"){
        perc.over = 800
      } else if (p == "inc3"){
        perc.over = 400
      } else {
        perc.over = 200
      }
      
      smotest <- list(name = "SMOTE with more upsampling",
                      func = function (x, y) {
                        library(DMwR)
                        dat <- if (is.data.frame(x)) x else as.data.frame(x)
                        dat$.y <- y
                        dat <- SMOTE(.y ~ ., data = dat, perc.over = perc.over, perc.under = 100)
                        list(x = dat[, !grepl(".y", colnames(dat), fixed = TRUE)], 
                             y = dat$.y)
                      },
                      first = TRUE)
      
      if (resampling == "up"){
        control <- trainControl(method="repeatedcv", number=5, repeats=2, summaryFunction = twoClassSummary, classProbs = TRUE, sampling = "up")
      } else if (resampling == "smote"){
        control <- trainControl(method="repeatedcv", number=5, repeats=2, summaryFunction = twoClassSummary, classProbs = TRUE, sampling = smotest)
      }
      
      # weights <- ifelse(y_train[,1] == "One",
      #                   (1/table(as.factor(y_train[,1]))[1]) * 0.5,
      #                   (1/table(as.factor(y_train[,1]))[2]) * 0.5)
      
      #df = data.frame(x_train, outcome = y_train)
      
      set.seed(seed)
      fit <- caret::train(x = x_train,
                          y = y_train,
                          method = "svmRadial",
                          metric = "ROC",
                          #weights = weights,
                          tuneLength = 10,
                          trControl = control,
                          na.action = na.omit)
      
      pred <- predict(fit, x_test, s = 'lambda.min')
      
      if (any(is.na(pred))) {
        
        return(list(pred=pred, obs=y_test, perf = NULL, fit=fit))
        
      } else {
        
        pred.prob <- predict(fit, x_test, s = 'lambda.min', type = "prob")
        roc <- roc(response = y_test, predictor = pred.prob$One, levels = c("Zero","One"))
        auc = auc(roc)
        opt.thres = coords(roc, x = "best", ret = "threshold", best.method = "closest.topleft", transpose = T)
        obs = y_test
        positive =  ifelse(pred.prob$One > opt.thres,1,0)
        sens = sum(positive)/length(positive) #caret::sensitivity(pred,obs)
        negative = ifelse(pred.prob$Zero < opt.thres,1,0)
        spec = sum(negative)/length(negative)  #caret::specificity(pred,obs)
        perf = list(roc = roc, opt.thres = opt.thres, auc = auc, sens = sens, spec = spec)
        
        return(list(pred=pred, obs=y_test, perf = perf, fit=fit))
      }
      
    } else {
      
      message("Labels of train or test set has only one level")
      res = NULL
      return(res)
      
    }
    
  }
  
}

make_fit = function(
  ###Train and return CV score if cv==T, linear models if cv==F
  train_list = NULL, #list of training data
  test_list = NULL, #list of testing data
  phenotype = NULL, #df of responses
  cv_list = NULL, #list of CV subsets
  omics = NULL,
  covar = NULL,
  lab = NULL,
  method = c("elnet","rf","svm","xgboost"),
  cv = TRUE,
  n = 1, #index of omics combination
  overlapped = F, #use overlapping dataset of all omics or not?
  useProteomics = T, #if overlapped == T, use proteomic features or not?
  useClinical = c("yes","noLab","onlyClin", "onlyCovar"), #use of clinical data as features?
  resampling = NULL,
  featureList = NULL, #Preprovide a list of features
  out.dir = NULL 
){
  
  # ncores = detectCores() - 1
  # cl = makeCluster(ncores)
  # registerDoParallel(cl)
  # clusterCall(cl, function(x) .libPaths(x), .libPaths())
  
  # tmp = apply(neuropathy, 2, as.character) %>% 
  #   as.data.frame() %>% 
  #   transform(., un14 = as.numeric(un14),
  #             un15 = as.numeric(un15),
  #             un17 = as.numeric(un17),
  #             un18 = as.numeric(un18),
  #             utmnsi = as.numeric(utmnsi)) %>%
  #   dplyr::select(-c(un12_3,un12_4,un12_8,un12_9,un13_3,un13_4,un13_8,un13_9))
  # tmp = apply(neuropathy, 2, as.character) %>% as.data.frame()
  # rownames(tmp) <- rownames(neuropathy)
  # neuropathy <- tmp
  
  o = names(train_list)[n]
  
  for(p in colnames(phenotype))  {
    print(paste0("Processing omics ",o," and phenotype ",p))
    
    #Process train set
    pheno = phenotype[,p,drop=F] %>% as.matrix()
    rownames(pheno) = rownames(phenotype)
    pheno = na.omit(pheno)
    clin_id = train_list[[o]][["Samples"]][["Clinical"]]
    pheno = pheno[intersect(rownames(pheno), as.character(clin_id)),,drop = F]
    ind = match(rownames(pheno), as.character(clin_id))
    sample_list = train_list[[o]][["Samples"]]
    sample_list = lapply(names(sample_list), function(x) return(sample_list[[x]][ind]))
    names(sample_list) = names(train_list[[o]][["Samples"]])
    feature_list = train_list[[o]][["Features"]]
    
    if (useClinical == "onlyClin"){
      feature_list = feature_list[c("Covar","Lab")]
    } else if (useClinical == "noLab"){
      feature_list["Lab"] <- NULL
    } else if (useClinical == "onlyCovar"){
      feature_list = feature_list["Covar"]
    }
    
    if (cv){
      cvglm_list = vector(mode = "list", length = length(cv_list[[o]]$train)) 
      
      print("...Cross-validating")
      for(k in 1:length(cvglm_list)) {
        
        print(paste("...Fold",k))
        
        if(length(grep("Genomics",o)) > 0){
          
          if (overlapped) {
            if (useProteomics) {
              genomics_k = bed_to_df(paste0("./Genomics/small_bed_neuro/GenomicsTranscriptomicsProteomicsMetabolomicsMethylomicsOlink","_",p,"_",k,".bed"))
            } else {
              genomics_k = bed_to_df(paste0("./Genomics/small_bed_neuro/GenomicsTranscriptomicsMetabolomicsMethylomicsOlink","_",p,"_",k,".bed"))
            }
          } else {
            genomics_k = bed_to_df(paste0("./Genomics/small_bed_neuro/",o,"_",p,"_",k,".bed"))
          }
          
        } else {
          genomics_k <- NULL
        }
        
        train_ind = cv_list[[o]][["train"]][[k]]
        test_ind = cv_list[[o]][["test"]][[k]]
        
        ind = match(train_ind, sample_list$Clinical) %>% na.omit()
        sample_list_k = lapply(names(sample_list), function(x) return(sample_list[[x]][ind]))
        names(sample_list_k) = names(sample_list)
        pheno_k = pheno[sample_list_k$Clinical,,drop=F]
        
        
        feature_mat_k = dimnames_to_mat(feature_list, sample_list_k, 
                                        genomics = genomics_k,
                                        omics = omics,
                                        covar = covar,
                                        lab = lab,
                                        phenotype = pheno_k,
                                        p = p) %>% na.omit()
        
        if (!is.null(featureList)){
          feature_mat_k = feature_mat_k[,intersect(colnames(feature_mat_k), featureList[[p]][[o]])]
        }
        
        if (length(grep("Methylomics",o)) > 0){
          if (length(intersect(feature_list[["Methylomics"]], colnames(feature_mat_k))) == 0){
            print("There is no methylomics features!")
            cvglm_list[[k]] <- NULL
            next 
          }
        }
        
        pheno_k = pheno_k[rownames(feature_mat_k),,drop=F]
        #pheno_k = ifelse(pheno_k == "1","One","Zero") 
        
        if (any(rownames(feature_mat_k) != rownames(pheno_k))) stop("Rownames of x and y do not match!")
        
        ind = match(test_ind, sample_list$Clinical) %>% na.omit()
        sample_list_k_test = lapply(names(sample_list), function(x) return(sample_list[[x]][ind]))
        names(sample_list_k_test) = names(sample_list)
        pheno_k_test = pheno[sample_list_k_test$Clinical,,drop=F]
        
        feature_mat_k_test = dimnames_to_mat(feature_list, sample_list_k_test, 
                                             genomics = genomics_k,
                                             omics = omics,
                                             covar = covar,
                                             lab = lab,
                                             phenotype = pheno_k_test,
                                             feature_sel = F,
                                             features = colnames(feature_mat_k),
                                             p = p) %>% na.omit()
        pheno_k_test = pheno_k_test[rownames(feature_mat_k_test),,drop=F]
        feature_mat_k_test = feature_mat_k_test[,match(colnames(feature_mat_k), colnames(feature_mat_k_test))]
        
        if (any(rownames(feature_mat_k_test) != rownames(pheno_k_test))) stop("Rownames of x and y do not match!")
        if (any(colnames(feature_mat_k_test) != colnames(feature_mat_k))) stop("Colnames of train and test do not match!")
        
        #pheno_k_test = ifelse(pheno_k_test == "1","One","Zero") 
        
        if (method == "elnet") {
          cvglm_list[[k]] = use_glm( x_train=feature_mat_k,
                                     y_train=pheno_k,
                                     x_test=feature_mat_k_test,
                                     y_test=pheno_k_test,
                                     seed = 993,
                                     p = p,
                                     resampling = resampling)
        }

        if (method == "rf") {
          cvglm_list[[k]] = use_rf( x_train=feature_mat_k,
                                   y_train=pheno_k,
                                   x_test=feature_mat_k_test,
                                   y_test=pheno_k_test,
                                   seed = 993,
                                    p = p,
                                   resampling = resampling)
        }
        
        if (method == "svm") {
          cvglm_list[[k]] = use_svm( x_train=feature_mat_k,
                                     y_train=pheno_k,
                                     x_test=feature_mat_k_test,
                                     y_test=pheno_k_test,
                                     seed = 993,
                                     p = p,
                                     resampling = resampling)
        }

        rm(feature_mat_k)
        rm(feature_mat_k_test)
      }
      
      saveRDS(cvglm_list, file = paste0(out.dir,o,",",p,",cv.rds"))
      rm(cvglm_list)
     
    } else {
      
      if(length(grep("Genomics",o)) > 0){
        
        if (overlapped) {
          if (useProteomics) {
            genomics = bed_to_df(paste0("./Genomics/small_bed_neuro/GenomicsTranscriptomicsProteomicsMetabolomicsMethylomicsOlink","_",p,".bed"))
          } else {
            genomics = bed_to_df(paste0("./Genomics/small_bed_neuro/GenomicsTranscriptomicsMetabolomicsMethylomicsOlink","_",p,".bed"))
          }
        } else {
          genomics = bed_to_df(paste0("./Genomics/small_bed_neuro/",o,"_",p,".bed"))
        }
      
      } else {
        genomics <- NULL
      }
      
      feature_mat = dimnames_to_mat(feature_list, sample_list, 
                                    genomics = genomics,
                                    omics = omics,
                                    covar = covar,
                                    lab = lab,
                                    phenotype = pheno,
                                    p = p) %>% na.omit()
      
      if (length(grep("Methylomics",o)) > 0){
        if(length(intersect(feature_list[["Methylomics"]], colnames(feature_mat))) == 0){
          print("There is no methylomics features!")
          next 
        }
      }
      
      pheno = pheno[rownames(feature_mat),,drop=F]
      #pheno = ifelse(pheno == "1","One","Zero") 
      
      #Process test set
      test_pheno = phenotype[,p,drop=F] %>% as.matrix()
      rownames(test_pheno) = rownames(phenotype)
      test_pheno = na.omit(test_pheno)
      test_clin = test_list[[o]][["Samples"]][["Clinical"]]
      test_pheno = test_pheno[intersect(rownames(test_pheno), as.character(test_clin)),,drop = F]
      ind = match(rownames(test_pheno), as.character(test_clin))
      sample_list_test = test_list[[o]][["Samples"]]
      sample_list_test = lapply(names(sample_list_test), function(x) return(sample_list_test[[x]][ind]))
      names(sample_list_test) = names(test_list[[o]][["Samples"]])
      feature_test = test_list[[o]][["Features"]]
      test_omics = dimnames_to_mat(feature_test, sample_list_test, 
                                   genomics = genomics,
                                   omics = omics,
                                   covar = covar,
                                   lab = lab,
                                   feature_sel = F,
                                   features = colnames(feature_mat),
                                   p = p) %>% na.omit()
      test_omics = test_omics[,match(colnames(feature_mat), colnames(test_omics))]
      test_pheno = test_pheno[rownames(test_omics),,drop=F]
      #test_pheno = ifelse(test_pheno == "1","One","Zero") 
      
      if (any(rownames(feature_mat) != rownames(pheno))) stop("Rownames of x and y do not match!")
      if (any(rownames(test_omics) != rownames(test_pheno))) stop("Rownames of x and y do not match!")
      if (any(colnames(test_omics) != colnames(feature_mat))) stop("Features of train and test do not match!")
      
      print("...Fitting")
      if (method == "elnet"){
        pred = use_glm(x_train=feature_mat,
                       y_train=pheno,
                       x_test = test_omics,
                       y_test = test_pheno,
                       seed = 993,
                       p = p)

        saveRDS(pred, file = paste0(out.dir,o,",",p,",mod.rds"))

      } else if (method == "rf"){
        pred = use_rf(x_train=feature_mat,
                       y_train=pheno,
                       x_test = test_omics,
                       y_test = test_pheno,
                       seed = 993,
                      p = p)

        saveRDS(pred, file = paste0(out.dir,o,",",p,",mod.rds"))
      
      } else if (method == "svm"){
        pred = use_svm(x_train=feature_mat,
                      y_train=pheno,
                      x_test = test_omics,
                      y_test = test_pheno,
                      seed = 993,
                      p = p)
        
        saveRDS(pred, file = paste0(out.dir,o,",",p,",mod.rds"))
      }

      rm(pred)
      rm(test_omics)
    }
    
    #rm(genomics)
    
  }
  
  #stopCluster(cl)
}




########## FOR ENSEMBLE TRAINING #################

standardise = function(
  ###Scale datasets into standard distribution
  df = NULL #dataframe or matrix to be transformed
){
  rm = apply(as.data.frame(df), 2, function(x) sd(x) == 0)
  df = df[,!rm]
  df = apply(as.data.frame(df), 2, function(x) (x - mean(x))/sd(x)) %>% as.matrix()
  return(df)
}

bed_to_df_copyNumber = function(
  ###Load in plink bed file and make df
  bed = NULL #path to bed file
){
  
  dat = read.plink(bed = bed)
  df = as.data.frame(dat$genotypes)
  rn = rownames(df)
  df.char = apply(df,2,as.numeric)
  rownames(df.char) = rn
  for (i in 1:nrow(df.char)){
    for (j in 1:ncol(df.char)){
      if (df.char[i,j] == 1){
        df.char[i,j] = 2
      } else if (df.char[i,j] == 2){
        df.char[i,j] = 1
      } else if (df.char[i,j] == 3){
        df.char[i,j] = 0
      }
    }
  }
  return(df.char)
}

make_overlap_set = function(
  ###Create a list of datasets of the same corresponding samples and different features
  data_list = NULL,
  id_table = NULL,
  y = NULL #phenotype df or matrix
){
  list_full = na.omit(id_table) %>% filter(as.character(Clinical) %in% rownames(y)) %>% lapply(.,as.character)
  id_fil = list()
  for (i in names(list_full)){
    id_fil[[i]] = which(list_full[[i]] %in% rownames(data_list[[i]]))
  }
  #id_fil = lapply(list_full, function(x, data_list = data_list) which(x %in% rownames(data_list[[x]])))
  id_final = Reduce(intersect, id_fil)
  for (i in names(data_list)){
    data_list[[i]] = data_list[[i]][list_full[[i]][id_final],]
    rownames(data_list[[i]]) = rownames(data_list[["Clinical"]][list_full[["Clinical"]][id_final],])
  }
  #data_list = lapply(data_list, function(x, id_final = id_final) return(x[id_final,]))
  return(data_list)
}

make_cv_list_ens <- function(
  ###Create a list of cross-validation subsets, same for all feature groups
  outcome = NULL, #named vector of outcome (factor) 
  #kfold = 5, #how many folds
  times = 20, #number of repeated splits
  partition = 0.8,
  kfold_inner = 5, #how many inner folds
  times_inner = 1, #number of k-fold cross validations for inner training
  seed = 993 #setting fixed seeds 
){
  
  set.seed(seed)
  outer_train = createDataPartition(y = outcome, times = times, p = partition)
  outer_test = lapply(outer_train, function(x) setdiff(1:length(outcome), x))
  # outer_validate_test = lapply(outer_train, function(x) names(outcome[setdiff(1:length(outcome), x)]))
  # outer_validate_ind = lapply(names(outer_validate_test), function(x) {set.seed(seed); createDataPartition(y = outcome[outer_validate_test[[x]]], times = 1, p = 0.5)[[1]]})
  # outer_validate = lapply(1:length(outer_validate_ind), function(x) outer_validate_test[[x]][outer_validate_ind[[x]]])
  # outer_test = lapply(1:length(outer_validate), function(x) setdiff(outer_validate_test[[x]], outer_validate[[x]]))
  
  outer = list(train = lapply(outer_train, function(x) names(outcome)[x]), test = lapply(outer_test, function(x) names(outcome)[x]))
  #outer = list(train = lapply(outer_train, function(x) names(outcome)[x]), validate = outer_validate, test = outer_test)
  
  inner_train = lapply(outer$train, function(x) {set.seed(seed); createMultiFolds(y = outcome[x], k = kfold_inner, times = times_inner)})
  outcome_outer = lapply(outer$train, function(x) names(outcome[x]))
  inner_test = list()
  for (i in 1:length(inner_train)){
    inner_test[[i]] = list()
    for ( j in 1: length(inner_train[[i]])){
      inner_train[[i]][[j]] = outcome_outer[[i]][inner_train[[i]][[j]]]
      inner_test[[i]][[j]] = setdiff(outcome_outer[[i]], inner_train[[i]][[j]])
    }
  }
  inner = list(train = inner_train, test = inner_test)
  return(list(outer = outer, inner = inner))
}

LogLoss <- function(pred, true, eps = 1e-15, weights = NULL) {

  pred = pmin(pmax(pred, eps), 1 - eps) # Bound the results
  
  if (is.null(weights)) {
    return(-(sum(
      true * log(pred) + (1 - true) * log(1 - pred)
    )) / length(true))
  } else{
    return(-weighted.mean(true * log(pred) + (1 - true) * log(1 - pred), weights))
  }
}

caretLogLoss <- function(data, lev = NULL, model = NULL) {
  cls <- levels(data$obs) #find class names
  loss <- LogLoss(
    pred = data[, cls[2]],
    true = as.numeric(data$obs) - 1,
    weights = data$weights
  )
  names(loss) <- c('myLogLoss')
  loss
}

bestSE = function(x, metric, num, maximize){
  #Customize selection function for parametter tunning (mean + 1SE) 
  
  perf <- x[, metric] + (x[, paste(metric, "SD", sep = "")])/sqrt(num)
  if (maximize){
    bestIter = which.max(perf)
  } else {
    bestIter = which.min(perf)
  }
  bestIter
}

fit_ensemble = function(
  data_list = NULL,
  y = NULL, #Named vector of outcome
  cv_list = NULL,
  n_cv = NULL, #Numer of inner folds for cross-validation
  seed = 993,
  n = NULL, #Iteration number in paralelizing 
  p_metric = c("AUROC","wLogLoss")
){
  
  df = do.call(cbind, data_list) %>% as.data.frame()
  i = n
  
  cat("####Iteration ",i," ####\n\n")
  
  x_train_i = df[cv_list$outer$train[[i]],]
  x_test_i = df[cv_list$outer$test[[i]],]
  
  y_train_i = y[cv_list$outer$train[[i]]]
  y_test_i = y[cv_list$outer$test[[i]]]
  
  # weights = ifelse(y_train == "One", table(y_train)[[2]]/table(y_train)[[1]], 1)
  # weights_test = ifelse(y_test == "One", table(y_test)[[2]]/table(y_test)[[1]], 1)
  
  if (nlevels(as.factor(y_train_i)) < 2 | nlevels(as.factor(y_test_i)) < 2){
    stop("Train or test set has only one class label!")
  }
  
  if (any(rownames(x_train_i) != names(y_train_i)) | any(rownames(x_test_i) != names(y_test_i))){
    stop("Samples in train and test sets do not match!")
  }
  
  if (p_metric == "wLogLoss"){
    sampling = NULL
    sumFunc = caretLogLoss
    metric = "myLogLoss"
    maximize = F
    weights_train_i = ifelse(y_train_i == "One", table(y_train_i)[[2]]/table(y_train_i)[[1]], 1)
    weights_test_i = ifelse(y_test_i == "One", table(y_test_i)[[2]]/table(y_test_i)[[1]], 1)
  } else if (p_metric == "AUROC"){
    sampling = "smote"
    sumFunc = twoClassSummary
    metric = "ROC"
    maximize = T
    weights_test_i = ifelse(y_test_i == "One", table(y_test_i)[[2]]/table(y_test_i)[[1]], 1)
  }
  
  my_control <- trainControl(
    method="repeatedcv",
    number=5,
    repeats = 4,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=sumFunc,
    sampling = sampling,
    allowParallel = T
  )
  
  cat("Fit single models for each modality\n")
  
  single_mods = list()
  for (m in names(data_list)){
    
    cat("...Data type: ",m,"\n")
    x_train_m = data_list[[m]][cv_list$outer$train[[i]],]
    
    if (p_metric == "AUROC"){
      set.seed(seed)
      fit_m <- caret::train(x = x_train_m,
                            y = y_train_i,
                            method="glmnet", 
                            metric=metric,
                            tuneLength = 20,
                            #weights = weights_n, #Got error when including weights, investigate later
                            maximize = maximize,
                            trControl=my_control,
                            importance = TRUE)
      single_mods[[m]] = fit_m
      
    } else if (p_metric == "wLogLoss"){
      set.seed(seed)
      fit_m <- caret::train(x = x_train_m,
                            y = y_train_i,
                            method="glmnet", 
                            metric=metric,
                            tuneLength = 20,
                            weights = weights_train_i, #Got error when including weights, investigate later
                            maximize = maximize,
                            trControl=my_control,
                            importance = TRUE)
      single_mods[[m]] = fit_m
    }
  }
  
  cat("Make meta model\n")
  
  cat("...Do cross-validation and make predictions on inner test sets\n")
  
  pred = list()
  for (m in names(data_list)){
    
    pred[[m]] = c()
    cat("......Data type: ",m,"\n")
    
    for (n in 1:n_cv){
      
      x_train_n = data_list[[m]][cv_list$inner$train[[i]][[n]],]
      x_test_n = data_list[[m]][cv_list$inner$test[[i]][[n]],]
      y_train_n = y[cv_list$inner$train[[i]][[n]]]
      
      if (p_metric == "AUROC"){
        set.seed(seed)
        fit_n <- caret::train(x = x_train_n,
                              y = y_train_n,
                              method="glmnet", 
                              metric=metric,
                              tuneLength = 20,
                              #weights = weights_n, #Got error when including weights, investigate later
                              maximize = maximize,
                              trControl=my_control,
                              importance = TRUE)
        pred_n = predict(fit_n, x_test_n, s = "lambda.min", type = "prob")$One
        
      } else if (p_metric == "wLogLoss"){
        weights_train_n = ifelse(y_train_n == "One", table(y_train_n)[[2]]/table(y_train_n)[[1]], 1)
        set.seed(seed)
        fit_n <- caret::train(x = x_train_n,
                              y = y_train_n,
                              method="glmnet", 
                              metric=metric,
                              tuneLength = 20,
                              weights = weights_train_n, #Got error when including weights, investigate later
                              maximize = maximize,
                              trControl=my_control,
                              importance = TRUE)
        pred_n = predict(fit_n, x_test_n, s = "lambda.min", type = "prob")$One
      }
      
      names(pred_n) = rownames(x_test_n)
      pred[[m]] = c(pred[[m]], pred_n)
      
    }
  }
  
  cat("...Fit a meta model\n")
  pred = as.data.frame(pred)
  y_inner_tests = y[rownames(pred)]
  weights_inner_tests = ifelse(y_inner_tests == "One", table(y_inner_tests)[[2]]/table(y_inner_tests)[[1]], 1)
  if (p_metric == "AUROC"){
    set.seed(seed)
    fit_meta <- caret::train(x = pred,
                          y = y_inner_tests,
                          method="glmnet", 
                          metric=metric,
                          tuneLength = 20,
                          #weights = weights_n, #Got error when including weights, investigate later
                          maximize = maximize,
                          trControl=my_control,
                          importance = TRUE)
    
  } else if (p_metric == "wLogLoss"){
    set.seed(seed)
    fit_meta <- caret::train(x = pred,
                          y = y_inner_tests,
                          method="glmnet", 
                          metric=metric,
                          tuneLength = 20,
                          weights = weights_inner_tests, #Got error when including weights, investigate later
                          maximize = maximize,
                          trControl=my_control,
                          importance = TRUE)
  }
  
  cat("Evaluate performannce on the outer test set\n")
  single_preds = predict(single_mods, x_test_i, type = "prob") %>% as.data.frame()
  keep_col = grep("One", colnames(single_preds))
  single_preds = single_preds[,keep_col, drop = F]
  colnames(single_preds) = names(single_mods)
  final_pred = predict(fit_meta, single_preds, type = "prob")$One
  roc <- roc(response = y_test_i, predictor = final_pred, levels = c("Zero","One"))
  auroc = auc(roc)[[1]]
  auprc = MLmetrics::PRAUC(final_pred, ifelse(y_test_i == "One",1,0))
  ll = LogLoss(final_pred, ifelse(y_test_i == "One",1,0), weights = weights_test_i)
  
  perf = data.frame(AUROC = auroc, AUPRC = auprc, WeightedLogLoss = ll)
  
  return(list(perf = perf, single_models = single_mods, meta_model = fit_meta, pred = final_pred, test_set = names(y_test_i)))
}

fit_concatenation = function(
  data_list = NULL,
  y = NULL, #Named vector of outcome
  cv_list = NULL,
  n_cv = NULL, #Numer of inner folds for cross-validation
  seed = 993,
  n = NULL, #Iteration number in paralelizing 
  p_metric = c("AUROC","wLogLoss")
){
  
  df = do.call(cbind, data_list) %>% as.data.frame()
  i = n
  
  cat("####Iteration ",i," ####\n\n")
  
  x_train_i = df[cv_list$outer$train[[i]],]
  x_test_i = df[cv_list$outer$test[[i]],]
  
  y_train_i = y[cv_list$outer$train[[i]]]
  y_test_i = y[cv_list$outer$test[[i]]]
  
  if (nlevels(as.factor(y_train_i)) < 2 | nlevels(as.factor(y_test_i)) < 2){
    stop("Train or test set has only one class label!")
  }
  
  if (any(rownames(x_train_i) != names(y_train_i)) | any(rownames(x_test_i) != names(y_test_i))){
    stop("Samples in train and test sets do not match!")
  }
  
  if (p_metric == "wLogLoss"){
    sampling = NULL
    sumFunc = caretLogLoss
    metric = "myLogLoss"
    maximize = F
    weights_train_i = ifelse(y_train_i == "One", table(y_train_i)[[2]]/table(y_train_i)[[1]], 1)
    weights_test_i = ifelse(y_test_i == "One", table(y_test_i)[[2]]/table(y_test_i)[[1]], 1)
  } else if (p_metric == "AUROC"){
    sampling = "smote"
    sumFunc = twoClassSummary
    metric = "ROC"
    maximize = T
    weights_test_i = ifelse(y_test_i == "One", table(y_test_i)[[2]]/table(y_test_i)[[1]], 1)
  }
  
  my_control <- trainControl(
    method="repeatedcv",
    number=5,
    repeats = 4,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=sumFunc,
    sampling = sampling,
    allowParallel = T
  )
  
  cat("Fit concatened dataset\n")
  if (p_metric == "AUROC"){
    set.seed(seed)
    fit_cat <- caret::train(x = x_train_i,
                          y = y_train_i,
                          method="glmnet", 
                          metric=metric,
                          tuneLength = 20,
                          #weights = weights_n, #Got error when including weights, investigate later
                          maximize = maximize,
                          trControl=my_control,
                          importance = TRUE)
    
  } else if (p_metric == "wLogLoss"){
    set.seed(seed)
    fit_cat <- caret::train(x = x_train_i,
                          y = y_train_i,
                          method="glmnet", 
                          metric=metric,
                          tuneLength = 20,
                          weights = weights_train_i, #Got error when including weights, investigate later
                          maximize = maximize,
                          trControl=my_control,
                          importance = TRUE)
  }
  
  cat("Evaluate performannce on the outer test set\n")
  pred_cat = predict(fit_cat, x_test_i, s = "lambda.min", type = "prob")$One
  roc <- roc(response = y_test_i, predictor = pred_cat, levels = c("Zero","One"))
  auroc = auc(roc)[[1]]
  auprc = MLmetrics::PRAUC(pred_cat, ifelse(y_test_i == "One",1,0))
  ll = LogLoss(pred_cat, ifelse(y_test_i == "One",1,0), weights = weights_test_i)
  
  perf = data.frame(AUROC = auroc, AUPRC = auprc, WeightedLogLoss = ll)
  
  return(list(perf = perf, model = fit_cat, pred = pred_cat, test_set = names(y_test_i)))
  
}

fit_forwardSelect = function(
  #Function to do CV and identify best model through foward feature selection
  data_list = NULL,
  y = NULL, #Named vector of outcome
  cv_list = NULL,
  n_cv = NULL, #Number of inner CVs
  n_datasets = NULL, #Number of datasets to evaluate
  p_metric = c("wLogLoss","AUROC"),
  seed = 993,
  n = NULL #Fold number in parallelizing 
){
  
  df = do.call(cbind, data_list) %>% as.data.frame()
  i = n
  
  cat("Fit train set",i,"\n")
  
  y_train_i = y[cv_list$outer$train[[i]]]
  y_test_i = y[cv_list$outer$test[[i]]]
  
  if (nlevels(as.factor(y_train_i)) < 2 | nlevels(as.factor(y_test_i)) < 2){
    stop("Train or test set has only one class label!")
  }
  
  if (p_metric == "wLogLoss"){
    sampling = NULL
    sumFunc = caretLogLoss
    metric = "myLogLoss"
    maximize = F
    weights_test_i = ifelse(y_test_i == "One", table(y_test_i)[[2]]/table(y_test_i)[[1]], 1)
  } else if (p_metric == "AUROC"){
    sampling = "smote"
    sumFunc = twoClassSummary
    metric = "ROC"
    maximize = T
    weights_test_i = ifelse(y_test_i == "One", table(y_test_i)[[2]]/table(y_test_i)[[1]], 1)
  }
  
  my_control <- trainControl(
    method="repeatedcv",
    number=5,
    repeats = 4,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=sumFunc,
    sampling = sampling,
    allowParallel = T
  )
  
  comb_list = lapply(1:n_datasets, function(x) combn(names(data_list),x, simplify = F))
  
  perf_validate = list()
  perf_test = list()
  pred = list()
  var_list = list()
  best_d = NULL
  
  for (j in 1:length(comb_list)){
    
    perf_j = list(roc = list(), pr = list(), ll = list())
    
    if (!is.null(best_d)){
      ind = lapply(comb_list[[j]], function(x) all(best_d %in% x)) %>% unlist()
      comb_list_fil = comb_list[[j]][ind]
    } else {
      comb_list_fil = comb_list[[j]]
    }
    
    for (d in 1:length(comb_list_fil)){
      
      comb = paste0(comb_list_fil[[d]], collapse = "")
      cat("Fit ",comb,"\n")
      
      roc_d = list()
      pr_d = list()
      ll_d = list()
      for (n in 1:n_cv){
        x_train_n = do.call(cbind, data_list[comb_list_fil[[d]]])[cv_list$inner$train[[i]][[n]],] %>% as.data.frame()
        y_train_n = y[rownames(x_train_n) ]
        x_test_n = do.call(cbind, data_list[comb_list_fil[[d]]])[cv_list$inner$test[[i]][[n]],] %>% as.data.frame()
        y_test_n = y[rownames(x_test_n) ]
        
        if (p_metric == "AUROC"){
          set.seed(seed)
          weights_test_n = ifelse(y_test_n == "One", table(y_test_n)[[2]]/table(y_test_n)[[1]], 1)
          fit_n <- caret::train(x = x_train_n,
                                y = y_train_n,
                                method="glmnet", 
                                metric=metric,
                                tuneLength = 20,
                                #weights = weights_n, #Got error when including weights, investigate later
                                maximize = maximize,
                                trControl=my_control,
                                importance = TRUE)
          pred_n = predict(fit_n, x_test_n, s = "lambda.min", type = "prob")$One
          roc <- roc(response = y_test_n, predictor = pred_n, levels = c("Zero","One"))
          roc_n = auc(roc)[[1]]
          pr_n = MLmetrics::PRAUC(pred_n, ifelse(y_test_n == "One",1,0))
          ll_n = LogLoss(pred_n, ifelse(y_test_n == "One",1,0), weights = weights_test_n)
          #perf_n = caTools::colAUC(pred_n, ifelse(y_test_n == "One",1,0))
    
        } else if (p_metric == "wLogLoss"){
          weights_n = ifelse(y_train_n == "One", table(y_train_n)[[2]]/table(y_train_n)[[1]], 1)
          weights_test_n = ifelse(y_test_n == "One", table(y_test_n)[[2]]/table(y_test_n)[[1]], 1)
          set.seed(seed)
          fit_n <- caret::train(x = x_train_n,
                                y = y_train_n,
                                method="glmnet", 
                                metric=metric,
                                tuneLength = 20,
                                weights = weights_n, #Got error when including weights, investigate later
                                maximize = maximize,
                                trControl=my_control,
                                importance = TRUE)
          pred_n = predict(fit_n, x_test_n, s = "lambda.min", type = "prob")$One
          roc <- roc(response = y_test_n, predictor = pred_n, levels = c("Zero","One"))
          roc_n = auc(roc)[[1]]
          pr_n = MLmetrics::PRAUC(pred_n, ifelse(y_test_n == "One",1,0))
          ll_n = LogLoss(pred_n, ifelse(y_test_n == "One",1,0), weights = weights_test_n)
        }
        roc_d[[n]] = roc_n
        pr_d[[n]] = pr_n
        ll_d[[n]] = ll_n
      }
      
      # perf_d = unlist(perf_d) %>% median()
      # perf_j[[d]] = perf_d
      perf_j$roc[[d]] = unlist(roc_d) %>% median()
      perf_j$pr[[d]] = unlist(pr_d) %>% median()
      perf_j$ll[[d]] = unlist(ll_d) %>% median()
      
    }
    
    if (p_metric == "AUROC"){
      best_ind = which.max(unlist(perf_j$roc))
    } else if (p_metric == "wLogLoss") {
      best_ind = which.min(unlist(perf_j$ll))
    }
    
    best_d = comb_list_fil[[best_ind]]
    cat("The best ",j, " base model is ",best_d,"\n")
    
    perf_validate[[j]] = data.frame(Complexity = paste0(j," Dataset(s)"), Model = paste0(best_d, collapse = ""), Value = c(perf_j$roc[[best_ind]],perf_j$pr[[best_ind]],perf_j$ll[[best_ind]]), Type = c("AUROC","AUCPR","Weighted LogLoss"))
    
    x_train_i = do.call(cbind, data_list[comb_list_fil[[best_ind]]])[cv_list$outer$train[[i]],] %>% as.data.frame()
    x_test_i = do.call(cbind, data_list[comb_list_fil[[best_ind]]])[cv_list$outer$test[[i]],] %>% as.data.frame()
    
    if (any(rownames(x_train_i) != names(y_train_i)) | any(rownames(x_test_i) != names(y_test_i))){
      stop("Samples in train and test sets do not match!")
    }
    
    if (p_metric == "AUROC"){
      set.seed(seed)
      fit_best <- caret::train(x = x_train_i,
                               y = y_train_i,
                               method="glmnet", 
                               metric=metric,
                               tuneLength = 20,
                               #weights = weights_train_i, #Got error when including weights, investigate later
                               maximize = maximize,
                               trControl=trainControl(method="repeatedcv",
                                                      number = 5,
                                                      repeats = 4,
                                                      savePredictions="final",
                                                      classProbs=TRUE,
                                                      summaryFunction=sumFunc,
                                                      sampling = sampling
                               ),
                               importance = TRUE)
    } else if (p_metric == "wLogLoss"){
      weights_train_i = ifelse(y_train_i == "One", table(y_train_i)[[2]]/table(y_train_i)[[1]], 1)
      set.seed(seed)
      fit_best <- caret::train(x = x_train_i,
                               y = y_train_i,
                               method="glmnet", 
                               metric=metric,
                               tuneLength = 20,
                               weights = weights_train_i,
                               maximize = maximize,
                               trControl=trainControl(method="repeatedcv",
                                                      number = 5,
                                                      repeats = 4,
                                                      savePredictions="final",
                                                      classProbs=TRUE,
                                                      summaryFunction=sumFunc,
                                                      sampling = sampling
                               ),
                               importance = TRUE)
    }
    
    var_imp = varImp(fit_best)$importance
    var = var_imp$Overall
    names(var) = rownames(var_imp)
    var = var[order(var, decreasing = T)]
    #var = var[var != 0]
    var_list[[j]] = var
    
    pred_test = predict(fit_best, x_test_i, s = "lambda.min", type = "prob")$One
    pred[[j]] = pred_test
    #roc_test = caTools::colAUC(pred_test, ifelse(y_test_i == "One",1,0))[,1]
    roc_test = auc(roc(response = y_test_i, predictor = pred_test, levels = c("Zero","One")))[[1]]
    pr_test = MLmetrics::PRAUC(pred_test, ifelse(y_test_i == "One",1,0))
    ll_test = LogLoss(pred_test, ifelse(y_test_i == "One",1,0), weights = weights_test_i)
    
    perf_test[[j]] = data.frame(Complexity = paste0(j," Dataset(s)"), Model = paste0(best_d, collapse = ""), Value = c(roc_test,pr_test,ll_test), Type = c("AUROC","AUCPR","Weighted LogLoss"))
    
  }
  
  perf_validate = do.call(rbind, perf_validate)
  perf_test = do.call(rbind, perf_test)
  return(list(perf_validate = perf_validate, perf_test = perf_test, pred = pred, var = var_list, train_samples = cv_list$outer$train[[i]], test_samples = cv_list$outer$test[[i]]))
}

fit_forwardSelectFromClinical = function(
  #Function to do CV and identify best model through foward feature selection
  data_list = NULL,
  y = NULL, #Named vector of outcome
  cv_list = NULL,
  n_cv = NULL, #Number of inner CVs
  n_datasets = NULL, #Number of datasets to evaluate
  p_metric = c("wLogLoss","AUROC"),
  seed = 993,
  n = NULL #Fold number in parallelizing
){
  
  df = do.call(cbind, data_list) %>% as.data.frame()
  i = n
  
  cat("Fit train set",i,"\n")
  
  y_train_i = y[cv_list$outer$train[[i]]]
  y_test_i = y[cv_list$outer$test[[i]]]
  
  if (nlevels(as.factor(y_train_i)) < 2 | nlevels(as.factor(y_test_i)) < 2){
    stop("Train or test set has only one class label!")
  }
  
  if (p_metric == "wLogLoss"){
    sampling = NULL
    sumFunc = caretLogLoss
    metric = "myLogLoss"
    maximize = F
    weights_test_i = ifelse(y_test_i == "One", table(y_test_i)[[2]]/table(y_test_i)[[1]], 1)
  } else if (p_metric == "AUROC"){
    sampling = "smote"
    sumFunc = twoClassSummary
    metric = "ROC"
    maximize = T
    weights_test_i = ifelse(y_test_i == "One", table(y_test_i)[[2]]/table(y_test_i)[[1]], 1)
  }
  
  my_control <- trainControl(
    method="repeatedcv",
    number=5,
    repeats = 4,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=sumFunc,
    sampling = sampling,
    allowParallel = T
  )
  
  cat("Cross-validation with clinical as baseline\n")
  roc_clin = list()
  pr_clin = list()
  ll_clin = list()
  for (n in 1:n_cv){
    x_train_clin = data_list$Clinical[cv_list$inner$train[[i]][[n]],] %>% as.data.frame()
    y_train_clin = y[rownames(x_train_clin) ]
    x_test_clin = data_list$Clinical[cv_list$inner$test[[i]][[n]],] %>% as.data.frame()
    y_test_clin = y[rownames(x_test_clin) ]
    
    if (p_metric == "AUROC"){
      set.seed(seed)
      weights_test_clin = ifelse(y_test_clin == "One", table(y_test_clin)[[2]]/table(y_test_clin)[[1]], 1)
      fit_clin <- caret::train(x = x_train_clin,
                            y = y_train_clin,
                            method="glmnet", 
                            metric=metric,
                            tuneLength = 20,
                            #weights = weights_n, #Got error when including weights, investigate later
                            maximize = maximize,
                            trControl=my_control,
                            importance = TRUE)
      pred_clin = predict(fit_clin, x_test_clin, s = "lambda.min", type = "prob")$One
      roc <- roc(response = y_test_clin, predictor = pred_clin, levels = c("Zero","One"))
      roc_n = auc(roc)[[1]]
      pr_n = MLmetrics::PRAUC(pred_clin, ifelse(y_test_clin == "One",1,0))
      ll_n = LogLoss(pred_clin, ifelse(y_test_clin == "One",1,0), weights = weights_test_clin)
      #perf_n = caTools::colAUC(pred_n, ifelse(y_test_n == "One",1,0))
      
    } else if (p_metric == "wLogLoss"){
      weights_clin = ifelse(y_train_clin == "One", table(y_train_clin)[[2]]/table(y_train_clin)[[1]], 1)
      weights_test_clin = ifelse(y_test_clin == "One", table(y_test_clin)[[2]]/table(y_test_clin)[[1]], 1)
      set.seed(seed)
      fit_clin <- caret::train(x = x_train_clin,
                            y = y_train_clin,
                            method="glmnet", 
                            metric=metric,
                            tuneLength = 20,
                            weights = weights_clin, #Got error when including weights, investigate later
                            maximize = maximize,
                            trControl=my_control,
                            importance = TRUE)
      pred_clin = predict(fit_clin, x_test_clin, s = "lambda.min", type = "prob")$One
      roc <- roc(response = y_test_clin, predictor = pred_clin, levels = c("Zero","One"))
      roc_n = auc(roc)[[1]]
      pr_n = MLmetrics::PRAUC(pred_clin, ifelse(y_test_clin == "One",1,0))
      ll_n = LogLoss(pred_clin, ifelse(y_test_clin == "One",1,0), weights = weights_test_clin)
    }
    roc_clin[[n]] = roc_n
    pr_clin[[n]] = pr_n
    ll_clin[[n]] = ll_n
  }
  
  perf_clin_validate = data.frame(Complexity = "1 Dataset(s)", Model = "Clinical", Value = c(median(unlist(roc_clin)), median(unlist(pr_clin)), median(unlist(ll_clin))), Type = c("AUROC","AUPRC","Weighted LogLoss"))
  
  cat("Testing with clinical as baseline\n")
  x_train_i_clin = data_list$Clinical[cv_list$outer$train[[i]],] %>% as.data.frame()
  x_test_i_clin = data_list$Clinical[cv_list$outer$test[[i]],] %>% as.data.frame()
  
  if (any(rownames(x_train_i_clin) != names(y_train_i)) | any(rownames(x_test_i_clin) != names(y_test_i))){
    stop("Samples in train and test sets do not match!")
  }
  
  if (p_metric == "AUROC"){
    set.seed(seed)
    fit_clin <- caret::train(x = x_train_i_clin,
                             y = y_train_i,
                             method="glmnet", 
                             metric=metric,
                             tuneLength = 20,
                             #weights = weights_train_i, #Got error when including weights, investigate later
                             maximize = maximize,
                             trControl=trainControl(method="repeatedcv",
                                                    number = 5,
                                                    repeats = 4,
                                                    savePredictions="final",
                                                    classProbs=TRUE,
                                                    summaryFunction=sumFunc,
                                                    sampling = sampling
                             ),
                             importance = TRUE)
  } else if (p_metric == "wLogLoss"){
    weights_train_i = ifelse(y_train_i == "One", table(y_train_i)[[2]]/table(y_train_i)[[1]], 1)
    set.seed(seed)
    fit_clin <- caret::train(x = x_train_i_clin,
                             y = y_train_i,
                             method="glmnet", 
                             metric=metric,
                             tuneLength = 20,
                             weights = weights_train_i,
                             maximize = maximize,
                             trControl=trainControl(method="repeatedcv",
                                                    number = 5,
                                                    repeats = 4,
                                                    savePredictions="final",
                                                    classProbs=TRUE,
                                                    summaryFunction=sumFunc,
                                                    sampling = sampling
                             ),
                             importance = TRUE)
  }
  
  var_imp_clin = varImp(fit_clin)$importance
  var_clin = var_imp_clin$Overall
  names(var_clin) = rownames(var_imp_clin)
  var_clin = var_clin[order(var_clin, decreasing = T)]
  
  pred_test_clin = predict(fit_clin, x_test_i_clin, s = "lambda.min", type = "prob")$One
  #roc_test = caTools::colAUC(pred_test, ifelse(y_test_i == "One",1,0))[,1]
  roc_test_clin = auc(roc(response = y_test_i, predictor = pred_test_clin, levels = c("Zero","One")))[[1]]
  pr_test_clin = MLmetrics::PRAUC(pred_test_clin, ifelse(y_test_i == "One",1,0))
  ll_test_clin = LogLoss(pred_test_clin, ifelse(y_test_i == "One",1,0), weights = weights_test_i)
  
  perf_clin_test = data.frame(Complexity = "1 Dataset(s)", Model = "Clinical", Value = c(roc_test_clin,pr_test_clin,ll_test_clin), Type = c("AUROC","AUPRC","Weighted LogLoss"))
  
  comb_list = lapply(1:(n_datasets - 1), function(x) combn(setdiff(names(data_list),"Clinical"),x, simplify = F))
  
  perf_validate = list()
  perf_test = list()
  pred = list()
  var_list = list()
  best_d = NULL
  
  for (j in 1:length(comb_list)){
    
    perf_j = list(roc = list(), pr = list(), ll = list())
    
    if (!is.null(best_d)){
      ind = lapply(comb_list[[j]], function(x) all(best_d %in% x)) %>% unlist()
      comb_list_fil = comb_list[[j]][ind]
    } else {
      comb_list_fil = comb_list[[j]]
    }
    
    for (d in 1:length(comb_list_fil)){
      
      comb = paste0(c("Clinical",comb_list_fil[[d]]), collapse = "")
      cat("Fit ",comb,"\n")
      
      roc_d = list()
      pr_d = list()
      ll_d = list()
      for (n in 1:n_cv){
        x_train_n = do.call(cbind, data_list[c("Clinical",comb_list_fil[[d]])])[cv_list$inner$train[[i]][[n]],] %>% as.data.frame()
        y_train_n = y[rownames(x_train_n) ]
        x_test_n = do.call(cbind, data_list[c("Clinical",comb_list_fil[[d]])])[cv_list$inner$test[[i]][[n]],] %>% as.data.frame()
        y_test_n = y[rownames(x_test_n) ]
        
        if (p_metric == "AUROC"){
          set.seed(seed)
          weights_test_n = ifelse(y_test_n == "One", table(y_test_n)[[2]]/table(y_test_n)[[1]], 1)
          fit_n <- caret::train(x = x_train_n,
                                y = y_train_n,
                                method="glmnet", 
                                metric=metric,
                                tuneLength = 20,
                                #weights = weights_n, #Got error when including weights, investigate later
                                maximize = maximize,
                                trControl=my_control,
                                importance = TRUE)
          pred_n = predict(fit_n, x_test_n, s = "lambda.min", type = "prob")$One
          roc <- roc(response = y_test_n, predictor = pred_n, levels = c("Zero","One"))
          roc_n = auc(roc)[[1]]
          pr_n = MLmetrics::PRAUC(pred_n, ifelse(y_test_n == "One",1,0))
          ll_n = LogLoss(pred_n, ifelse(y_test_n == "One",1,0), weights = weights_test_n)
          #perf_n = caTools::colAUC(pred_n, ifelse(y_test_n == "One",1,0))
          
        } else if (p_metric == "wLogLoss"){
          weights_n = ifelse(y_train_n == "One", table(y_train_n)[[2]]/table(y_train_n)[[1]], 1)
          weights_test_n = ifelse(y_test_n == "One", table(y_test_n)[[2]]/table(y_test_n)[[1]], 1)
          set.seed(seed)
          fit_n <- caret::train(x = x_train_n,
                                y = y_train_n,
                                method="glmnet", 
                                metric=metric,
                                tuneLength = 20,
                                weights = weights_n, #Got error when including weights, investigate later
                                maximize = maximize,
                                trControl=my_control,
                                importance = TRUE)
          pred_n = predict(fit_n, x_test_n, s = "lambda.min", type = "prob")$One
          roc <- roc(response = y_test_n, predictor = pred_n, levels = c("Zero","One"))
          roc_n = auc(roc)[[1]]
          pr_n = MLmetrics::PRAUC(pred_n, ifelse(y_test_n == "One",1,0))
          ll_n = LogLoss(pred_n, ifelse(y_test_n == "One",1,0), weights = weights_test_n)
        }
        roc_d[[n]] = roc_n
        pr_d[[n]] = pr_n
        ll_d[[n]] = ll_n
      }
      
      perf_j$roc[[d]] = unlist(roc_d) %>% median()
      perf_j$pr[[d]] = unlist(pr_d) %>% median()
      perf_j$ll[[d]] = unlist(ll_d) %>% median()
    }
    
    if (p_metric == "AUROC"){
      best_ind = which.max(unlist(perf_j$roc))
    } else if (p_metric == "wLogLoss") {
      best_ind = which.min(unlist(perf_j$ll))
    }
    
    best_d = comb_list_fil[[best_ind]]
    cat("The best ",j+1, " base model is ",paste0(c("Clinical",best_d), collapse = ""),"\n")
    
    perf_validate[[j]] = data.frame(Complexity = paste0(j+1," Dataset(s)"), Model = paste0(c("Clinical",best_d), collapse = ""), Value = c(perf_j$roc[[best_ind]],perf_j$pr[[best_ind]],perf_j$ll[[best_ind]]), Type = c("AUROC","AUPRC","Weighted LogLoss"))
    
    x_train_i = do.call(cbind, data_list[comb_list_fil[[best_ind]]])[cv_list$outer$train[[i]],] %>% as.data.frame()
    x_test_i = do.call(cbind, data_list[comb_list_fil[[best_ind]]])[cv_list$outer$test[[i]],] %>% as.data.frame()
    
    if (any(rownames(x_train_i) != names(y_train_i)) | any(rownames(x_test_i) != names(y_test_i))){
      stop("Samples in train and test sets do not match!")
    }
    
    if (p_metric == "AUROC"){
      set.seed(seed)
      fit_best <- caret::train(x = x_train_i,
                               y = y_train_i,
                               method="glmnet", 
                               metric=metric,
                               tuneLength = 20,
                               #weights = weights_train_i, #Got error when including weights, investigate later
                               maximize = maximize,
                               trControl=trainControl(method="repeatedcv",
                                                      number = 5,
                                                      repeats = 4,
                                                      savePredictions="final",
                                                      classProbs=TRUE,
                                                      summaryFunction=sumFunc,
                                                      sampling = sampling
                               ),
                               importance = TRUE)
    } else if (p_metric == "wLogLoss"){
      weights_train_i = ifelse(y_train_i == "One", table(y_train_i)[[2]]/table(y_train_i)[[1]], 1)
      set.seed(seed)
      fit_best <- caret::train(x = x_train_i,
                               y = y_train_i,
                               method="glmnet", 
                               metric=metric,
                               tuneLength = 20,
                               weights = weights_train_i,
                               maximize = maximize,
                               trControl=trainControl(method="repeatedcv",
                                                      number = 5,
                                                      repeats = 4,
                                                      savePredictions="final",
                                                      classProbs=TRUE,
                                                      summaryFunction=sumFunc,
                                                      sampling = sampling
                               ),
                               importance = TRUE)
    }
    
    var_imp = varImp(fit_best)$importance
    var = var_imp$Overall
    names(var) = rownames(var_imp)
    var = var[order(var, decreasing = T)]
    #var = var[var != 0]
    var_list[[j]] = var
    
    pred_test = predict(fit_best, x_test_i, s = "lambda.min", type = "prob")$One
    pred[[j]] = pred_test
    #roc_test = caTools::colAUC(pred_test, ifelse(y_test_i == "One",1,0))[,1]
    roc_test = auc(roc(response = y_test_i, predictor = pred_test, levels = c("Zero","One")))[[1]]
    pr_test = MLmetrics::PRAUC(pred_test, ifelse(y_test_i == "One",1,0))
    ll_test = LogLoss(pred_test, ifelse(y_test_i == "One",1,0), weights = weights_test_i)
    
    perf_test[[j]] = data.frame(Complexity = paste0(j+1," Dataset(s)"), Model = paste0(c("Clinical",best_d), collapse = ""), Value = c(roc_test,pr_test,ll_test), Type = c("AUROC","AUPRC","Weighted LogLoss"))
    
  }
  
  pred = c(list(pred_test_clin), pred)
  perf_validate = do.call(rbind, perf_validate) %>% rbind(perf_clin_validate,.)
  perf_test = do.call(rbind, perf_test) %>% rbind(perf_clin_test,.)
  return(list(perf_validate = perf_validate, perf_test = perf_test, pred = pred, var = c(list(var_clin),var_list), train_samples = cv_list$outer$train[[i]], test_samples = cv_list$outer$test[[i]])) 
}

fit_forwardSelectEnsemble = function(
  #Function to do CV and identify best model through foward feature selection
  data_list = NULL,
  y = NULL, #Named vector of outcome
  cv_list = NULL,
  n_datasets = NULL, #Number of datasets to evaluate
  metric = c("wLogLoss","AUROC"),
  seed = 993,
  n = NULL #Fold number in parallelizing 
){
  
  df = do.call(cbind, data_list) %>% as.data.frame()
  i = n
  
  cat("Folds ",i,"\n")
  
  x_train_i = df[cv_list$outer$train[[i]],]
  x_test_i = df[cv_list$outer$test[[i]],]
  
  y_train = y[cv_list$outer$train[[i]]]
  y_test = y[cv_list$outer$test[[i]]]
  
  if (nlevels(as.factor(y_train)) < 2 | nlevels(as.factor(y_test)) < 2){
    stop("Train or test set has only one class label!")
  }
  
  if (any(rownames(x_train_i) != names(y_train)) | any(rownames(x_test_i) != names(y_test))){
    stop("Samples in train and test sets do not match!")
  }
  
  if (metric == "wLogLoss"){
    sampling = NULL
    sumFunc = caretLogLoss
    metric = "myLogLoss"
    maximize = F
    weights = ifelse(y_train == "One", table(y_train)[[2]]/table(y_train)[[1]], 1)
    weights_test = ifelse(y_test == "One", table(y_test)[[2]]/table(y_test)[[1]], 1)
  } else if (metric == "AUROC"){
    sampling = "smote"
    sumFunc = twoClassSummary
    metric = "ROC"
    maximize = T
    weights = rep(1, length(y_train))
    weights_test = ifelse(y_test == "One", table(y_test)[[2]]/table(y_test)[[1]], 1)
  }
  
  my_control <- trainControl(
    method="repeatedcv",
    #number=5,
    #repeats = 4,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=sumFunc,
    sampling = sampling,
    index = lapply(cv_list$inner_train[[i]],function(x) match(x, names(y_train))),
    allowParallel = T
  )
  
  cat("Fit base models\n")
  base = list()
  for (d in names(data_list)){
    
    cat("...Data type: ",d,"\n")
    
    x_train_d = data_list[[d]][cv_list$outer$train[[i]],]
    
    set.seed(seed)
    fit_d <- caret::train(x = x_train_d,
                          y = y_train,
                          method="glmnet", 
                          metric=metric,
                          tuneLength = 20,
                          weights = weights, 
                          maximize = maximize,
                          trControl=my_control,
                          importance = TRUE)
    base[[d]] = fit_d
  }
  
  cat("Evaluating base models performances\n")
  pred_base_test <- predict(base, x_test_i, type = "prob") %>% as.data.frame()
  keep_col = grep("One", colnames(pred_base_test))
  pred_base_test = pred_base_test[,keep_col, drop = F]
  colnames(pred_base_test) = names(base)
  roc_base = caTools::colAUC(pred_base_test, ifelse(y_test == "One",1,0))
  ll_base = lapply(pred_base_test, function(x) LogLoss(x, ifelse(y_test == "One",1,0), weights = weights_test)) %>% as.data.frame()
  pr_base = lapply(pred_base_test, function(x) MLmetrics::PRAUC(x, ifelse(y_test == "One",1,0))) %>% as.data.frame()
  best_ind_roc_base = apply(roc_base,1,which.max)
  best_ind_ll_base = apply(ll_base,1,which.min)
  best_ind_pr_base = apply(pr_base,1,which.max)
  best_base = colnames(pred_base_test)[best_ind_roc_base]
  cat("...The best base model is ",best_base,"\n")
  
  return(list(perf = data.frame(Complexity = "1 Dataset", 
                                Model = c(colnames(roc_base)[best_ind_roc_base],colnames(pr_base)[best_ind_pr_base], colnames(ll_base)[best_ind_ll_base]), 
                                Value = c(roc_base[1,best_ind_roc_base],pr_base[1,best_ind_pr_base],ll_base[1,best_ind_ll_base]), 
                                Type = c("AUROC","AUCPR","Weighted LogLoss")), mod = base[[best_base]]))
  
  cat("Prepare for ensemble learning\n")
  pred_base_train = lapply(base, function(x) arrange(x$pred, Resample, rowIndex)$One) %>% as.data.frame()
  y_train_truth = y_train[arrange(base[[1]]$pred, Resample, rowIndex)$rowIndex]
  y_train_truth = ifelse(y_train_truth == "One", 1, 0)
  comb_list = lapply(2:n_datasets, function(x) combn(names(data_list),x, simplify = F))

  perf_all = list()
  fit_all = list()
  best_c = NULL

  for (j in 1:length(comb_list)){

    ind_j = lapply(comb_list[[j]], function(x) best_base %in% x) %>% unlist()
    comb_list_j = comb_list[[j]][ind_j]

    perf = list(roc = list(), pr = list(), ll = list())
    fit = list()

    if (!is.null(best_c)){
      ind = lapply(comb_list_j, function(x) all(best_c %in% x)) %>% unlist()
      comb_list_fil = comb_list_j[ind]
    } else {
      comb_list_fil = comb_list_j
    }

    for (c in 1:length(comb_list_fil)){

      comb = paste0(comb_list_fil[[c]], collapse = "")
      cat("Fit ensemble of ",comb,"\n")

      pred_base_c = pred_base_train[,comb_list_fil[[c]]]
      y_train_best_ind = abs(pred_base_c - y_train_truth) %>% apply(.,1,which.min)
      y_train_best = colnames(pred_base_c)[y_train_best_ind]
      # cat("...Best datasets include:\n")
      # print(table(y_train_best))

      set.seed(seed)
      fit_meta = caret::train(x = x_train_i[arrange(base[[1]]$pred, Resample, rowIndex)$rowIndex,],
                              y = y_train_best %>% as.factor(),
                              method="glmnet",
                              metric="logLoss",
                              tuneLength = 20,
                              trControl=trainControl(
                                method="boot",
                                number=10,
                                savePredictions="final",
                                classProbs=TRUE,
                                summaryFunction=multiClassSummary
                              ),
                              importance = TRUE)

      pred_weight = predict(fit_meta, x_test_i, s = "lambda.min", type = "prob")
      pred_base_fil = pred_base_test[,colnames(pred_weight)]
      pred_ensemble = apply(pred_base_fil*pred_weight,1,sum)
      perf_roc = caTools::colAUC(pred_ensemble, y_test)[,1]
      perf_pr = MLmetrics::PRAUC(pred_ensemble, ifelse(y_test == "One",1,0))
      perf_ll = LogLoss(pred_ensemble, ifelse(y_test == "One",1,0), weights = weights_test)

      fit[[comb]] = fit_meta
      perf$roc[[comb]] = perf_roc
      perf$pr[[comb]] = perf_pr
      perf$ll[[comb]] = perf_ll
    }

    best_ind_roc = unlist(perf$roc) %>% which.max()
    best_ind_pr = unlist(perf$pr) %>% which.max()
    best_ind_ll = unlist(perf$ll) %>% which.min()
    perf_all[[j]] = data.frame(Complexity = paste0(j+1," Datasets"), Model = c(names(fit)[best_ind_roc],names(fit)[best_ind_pr], names(fit)[best_ind_ll]),
                               Value = c(perf$roc[[best_ind_roc]], perf$pr[[best_ind_pr]], perf$ll[[best_ind_ll]]), Type = c("AUROC","AUCPR","Weighted LogLoss"))
    fit_all[[j]] = fit[[best_ind_ll]]

    best_c = comb_list_fil[[best_ind_ll]]
    cat("The best ensemble of",best_base,"and", j, "base model is",best_c,"\n")
  }

  perf_all = do.call(rbind, perf_all) %>% rbind(data.frame(Complexity = "1 Dataset", Model = c(colnames(roc_base)[best_ind_roc_base],colnames(pr_base)[best_ind_pr_base], colnames(ll_base)[best_ind_ll_base]),
                                                           Value = c(roc_base[1,best_ind_roc_base],pr_base[1,best_ind_pr_base],ll_base[1,best_ind_ll_base]), Type = c("AUROC","AUCPR","Weighted LogLoss")),.)
  return(list(perf = perf_all, mod = fit_all))
  
}

#########################################################################
######################## DEPRECATED #####################################
#########################################################################

# use_glm <- function(
#   ###Train a linear model for each omics combination and phenotype 
#   x_train=NULL,
#   y_train=NULL,
#   x_test=NULL,
#   y_test=NULL,
#   cv = TRUE,
#   seed = 993
# ){
#   
#   control <- trainControl(method="repeatedcv", number=5, repeats=1)
#   
#   if(cv){
#     set.seed(seed)
#     fit <- caret::train(x = x_train, 
#                         y = as.double(y_train), 
#                         method = "glmnet",
#                         metric = "RMSE",
#                         tuneLength = 10, 
#                         trControl = control,
#                         na.action = na.omit)
#     
#     pred <- predict(fit, x_test, s = 'lambda.min')
#     diff <- pred - y_test
#     pearson <- cor.test(pred,y_test, method = "pearson", exact = F)$estimate
#     spearman <- cor.test(pred,y_test, method = "spearman", exact = F)$estimate
#     corr <- list(pearson = pearson, spearman = spearman)
#     
#   }else{
#     set.seed(seed)
#     fit <- caret::train(x = x_train, 
#                         y = as.double(y_train), 
#                         method = "glmnet",
#                         metric = "RMSE",
#                         tuneLength = 10, 
#                         trControl = control,
#                         na.action = na.omit)
#     
#     pred <- predict(fit, x_test, s = 'lambda.min')
#     diff <- pred - y_test
#     pearson <- cor.test(pred,y_test, method = "pearson", exact = F)$estimate
#     spearman <- cor.test(pred,y_test, method = "spearman", exact = F)$estimate
#     corr <- list(pearson = pearson, spearman = spearman)
#   }
#   
#   return(list(pred=pred, diff=diff, corr = corr, fit=fit))
# }
# 
# 
# make_fit = function(
#   ###Train and return CV score if cv==T, linear models if cv==F
#   omics_list = NULL, #list of training data
#   test_list = NULL, #list of testing data
#   phenotype = NULL, #df of phenotypes
#   cv_list = NULL, #list of CV subsets
#   omics = NULL,
#   covar = NULL,
#   cv = TRUE,
#   n = 1 #index of omics combination  
# ){
#   
#   # ncores = detectCores() - 1
#   # cl = makeCluster(ncores)
#   # registerDoParallel(cl)
#   # clusterCall(cl, function(x) .libPaths(x), .libPaths())
#   
#   o = names(omics_list)[n]
#   n_phenotype = colnames(dplyr::select(phenotype,-zz_nr))
#   
#   for(p in n_phenotype)  {
#     
#     print(paste0("Processing omics ",o," for phenotype ",p))
#     
#     #Process train set
#     pheno = phenotype[,p,drop=F] %>% as.matrix()
#     rownames(pheno) = phenotype$zz_nr
#     pheno = na.omit(pheno)
#     clin_id = omics_list[[o]][["Samples"]][["Clinical"]]
#     pheno = pheno[intersect(rownames(pheno), as.character(clin_id)),,drop = F]
#     ind = match(rownames(pheno), as.character(clin_id))
#     sample_list = omics_list[[o]][["Samples"]]
#     sample_list = lapply(names(sample_list), function(x) return(sample_list[[x]][ind]))
#     names(sample_list) = names(omics_list[[o]][["Samples"]])
#     feature_list = omics_list[[o]][["Features"]]
#     
#     if (cv){
#       cvglm_list = vector(mode = "list", length = 10) 
#       
#       print("...Cross-validating")
#       for(k in 1:10) {
#         
#         print(paste("...Fold",k))
#         
#         if(length(grep("Genomics",o)) > 0){
#           genomics_k = bed_to_df(paste0("./small_bed/",p,"_",k,".bed"))
#         } else {
#           genomics_k <- NULL
#         }
#         
#         train_ind = cv_list[[o]][["train"]][[k]]
#         test_ind = cv_list[[o]][["test"]][[k]]
#         
#         ind = match(train_ind, sample_list$Clinical) %>% na.omit()
#         sample_list_k = lapply(names(sample_list), function(x) return(sample_list[[x]][ind]))
#         names(sample_list_k) = names(sample_list)
#         pheno_k = pheno[sample_list_k$Clinical,,drop=F]
#         
#         feature_mat_k = dimnames_to_mat(feature_list, sample_list_k, 
#                                         genomics = genomics_k,
#                                         omics = omics,
#                                         covar = covar,
#                                         phenotype = pheno_k)
#         
#         if (any(rownames(feature_mat_k) != rownames(pheno_k))) stop("Rownames of x and y do not match!")
#         
#         ind = match(test_ind, sample_list$Clinical) %>% na.omit()
#         sample_list_k_test = lapply(names(sample_list), function(x) return(sample_list[[x]][ind]))
#         names(sample_list_k_test) = names(sample_list)
#         pheno_k_test = pheno[sample_list_k_test$Clinical,,drop=F]
#         
#         feature_mat_k_test = dimnames_to_mat(feature_list, sample_list_k_test, 
#                                              genomics = genomics_k,
#                                              omics = omics,
#                                              covar = covar,
#                                              phenotype = pheno_k_test,
#                                              feature_sel = F,
#                                              features = setdiff(colnames(feature_mat_k), colnames(covar)))
#         
#         if (any(rownames(feature_mat_k_test) != rownames(pheno_k_test))) stop("Rownames of x and y do not match!")
#         
#         cvglm_list[[k]] = use_glm( x_train=feature_mat_k,
#                                    y_train=pheno_k,
#                                    x_test=feature_mat_k_test,
#                                    y_test=pheno_k_test,
#                                    cv = TRUE,
#                                    seed = 993)
#         rm(feature_mat_k)
#         rm(feature_mat_k_test)
#       }
#       
#       saveRDS(cvglm_list, file = paste0("./Result1/",o,",",p,",cv.rds"))
#       rm(cvglm_list)
#       
#     } else {
#       
#       if(length(grep("Genomics",o)) > 0){
#         genomics = bed_to_df(paste0("./small_bed/",p,".bed"))
#       } else {
#         genomics <- NULL
#       }
#       
#       feature_mat = dimnames_to_mat(feature_list, sample_list, 
#                                     genomics = genomics,
#                                     omics = omics,
#                                     covar = covar,
#                                     phenotype = pheno)
#       
#       #Process test set
#       test_pheno = phenotype[,p,drop=F] %>% as.matrix()
#       rownames(test_pheno) = phenotype$zz_nr
#       test_pheno = na.omit(test_pheno)
#       test_clin = test_list[[o]][["Samples"]][["Clinical"]]
#       test_pheno = test_pheno[intersect(rownames(test_pheno), as.character(test_clin)),,drop = F]
#       ind = match(rownames(test_pheno), as.character(test_clin))
#       sample_list_test = test_list[[o]][["Samples"]]
#       sample_list_test = lapply(names(sample_list_test), function(x) return(sample_list_test[[x]][ind]))
#       names(sample_list_test) = names(test_list[[o]][["Samples"]])
#       feature_test = test_list[[o]][["Features"]]
#       test_omics = dimnames_to_mat(feature_test, sample_list_test, 
#                                    genomics = genomics,
#                                    omics = omics,
#                                    covar = covar,
#                                    phenotype = test_pheno,
#                                    feature_sel = F,
#                                    features = setdiff(colnames(feature_mat), colnames(covar)))
#       test_omics = test_omics[,match(colnames(feature_mat), colnames(test_omics))]
#       
#       if (any(rownames(feature_mat) != rownames(pheno))) stop("Rownames of x and y do not match!")
#       if (any(rownames(test_omics) != rownames(test_pheno))) stop("Rownames of x and y do not match!")
#       if (any(colnames(test_omics) != colnames(feature_mat))) stop("Features of train and test do not match!")
#       
#       print("...Fitting")
#       pred = use_glm(x_train=feature_mat,
#                      y_train=pheno,
#                      x_test = test_omics,
#                      y_test = test_pheno,
#                      cv = F,
#                      seed = 993)
#       saveRDS(pred, file = paste0("./Result1/",o,",",p,",mod.rds"))
#       rm(pred)
#       rm(test_omics)
#     }
#     
#     #rm(genomics)
#     
#   }
#   
#   #stopCluster(cl)
# }
# 
# make_fit_neuro = function(
#   ###Train and return CV score if cv==T, linear models if cv==F
#   omics_list = NULL, #list of training data
#   test_list = NULL, #list of testing data
#   phenotype = NULL, #df of phenotypes
#   neuropathy = NULL, #df of neuropathy vars
#   cv_list = NULL, #list of CV subsets
#   covar = NULL,
#   cv = TRUE,
#   n = 1 #index of omics combination  
# ){
#   
#   # ncores = detectCores() - 1
#   # cl = makeCluster(ncores)
#   # registerDoParallel(cl)
#   # clusterCall(cl, function(x) .libPaths(x), .libPaths())
#   
#   o = names(omics_list)[n]
#   tmp = apply(neuropathy, 2, as.character) %>% 
#     as.data.frame() %>% 
#     transform(., un14 = as.numeric(un14),
#               un15 = as.numeric(un15),
#               un17 = as.numeric(un17),
#               un18 = as.numeric(un18)) %>%
#     dplyr::select(-c(un12_3,un12_4,un12_8,un12_9,un13_3,un13_4,un13_8,un13_9))
#   rownames(tmp) <- rownames(neuropathy)
#   neuropathy <- tmp
#   ne = colnames(neuropathy)
#   phenotype = column_to_rownames(phenotype,"zz_nr")
#   
#   print(paste0("Processing omics ",o))
#   
#   for(n in ne)  {
#     
#     print(paste0("Processing data for neuropathy variable ",n))
#     
#     #Process train set
#     
#     neuro = neuropathy[,n,drop=F]  %>% as.matrix()
#     rownames(neuro) = rownames(neuropathy)
#     neuro = na.omit(neuro)
#     clin_id = omics_list[[o]][["Samples"]][["Clinical"]]
#     neuro = neuro[intersect(rownames(neuro), as.character(clin_id)),,drop = F]
#     ind = match(rownames(neuro), as.character(clin_id))
#     clin_id = clin_id[ind]
#     train_phenotype = phenotype[clin_id,]
#     train_covar = covar[clin_id,]
#     train_phenotype = cbind(train_phenotype, train_covar) %>% na.omit()
#     neuro = neuro[rownames(train_phenotype),,drop=F]
#     
#     if (cv){
#       cvglm_list = vector(mode = "list", length = 10) 
#       
#       print("...Cross-validating")
#       for(k in 1:10) {
#         
#         train_ind = cv_list[[o]][["train"]][[k]]
#         test_ind = cv_list[[o]][["test"]][[k]]
#         
#         pheno_train = train_phenotype[intersect(train_ind, rownames(train_phenotype)),,drop=F]
#         neuro_train = neuro[intersect(train_ind, rownames(neuro)),,drop=F] 
#         pheno_test = train_phenotype[intersect(test_ind, rownames(train_phenotype)),,drop=F]
#         neuro_test = neuro[intersect(test_ind, rownames(neuro)),,drop=F] 
#         
#         if (class(neuro[,1]) == "numeric"){
#           
#           cvglm_list[[k]] = use_glm_neuro( x_train=pheno_train,
#                                            y_train=neuro_train,
#                                            x_test=pheno_test,
#                                            y_test=neuro_test,
#                                            seed = 993,
#                                            type = "regression")
#         }
#         
#         if (class(neuro[,1]) == "character"){
#           
#           cvglm_list[[k]] = use_glm_neuro( x_train=pheno_train,
#                                            y_train=neuro_train,
#                                            x_test=pheno_test,
#                                            y_test=neuro_test,
#                                            seed = 993,
#                                            type = "classification")
#         } 
#         
#       }
#       
#       saveRDS(cvglm_list, file = paste0("./Result2/",o,",",n,",cv.rds"))
#       rm(cvglm_list)
#       
#     } else {
#       
#       #Process test set
#       test_neuro = neuropathy[,n,drop=F] %>% as.matrix()
#       rownames(test_neuro) = rownames(neuropathy)
#       test_neuro = na.omit(test_neuro)
#       test_clin = test_list[[o]][["Samples"]][["Clinical"]]
#       test_neuro = test_neuro[intersect(rownames(test_neuro), as.character(test_clin)),,drop = F]
#       ind = match(rownames(test_neuro), as.character(test_clin))
#       test_clin = test_clin[ind]
#       test_pheno = phenotype[test_clin,]
#       test_covar = covar[test_clin,]
#       test_pheno = cbind(test_pheno, test_covar)
#       
#       print("...Fitting")
#       
#       if(class(test_neuro[,1]) == "numeric"){
#         
#         pred = use_glm_neuro(x_train=train_phenotype,
#                              y_train=neuro,
#                              x_test = test_pheno,
#                              y_test = test_neuro,
#                              seed = 993,
#                              type = "regression")
#       }
#       
#       if(class(test_neuro[,1]) == "character"){
#         
#         pred = use_glm_neuro(x_train=train_phenotype,
#                              y_train=neuro,
#                              x_test = test_pheno,
#                              y_test = test_neuro,
#                              seed = 993,
#                              type = "classification")
#       }
#       
#       saveRDS(pred, file = paste0("./Result2/",o,",",n,",mod.rds"))
#       rm(pred)
#       
#     }
#     
#   }
#   
#   #stopCluster(cl)
# }
# 
# fit_ensemble = function(
#   data_list = NULL,
#   y = NULL, #Named vector of outcome
#   cv_list = NULL,
#   seed = 993,
#   n = NULL #Fold number in parallelizing 
# ){
#   
#   cat("Fit concatened dataset\n")
#   set.seed(seed)
#   fit_cat <- caret::train(x = df[cv_list$outer$train[[i]],],
#                           y = y_train,
#                           method="glmnet", 
#                           metric="ROC",
#                           tuneLength = 20,
#                           weights = weights, 
#                           maximize = T,
#                           trControl=my_control,
#                           importance = TRUE)
#   
#   pred_cat = predict(fit_cat, x_test_i, s = "lambda.min", type = "prob")
#   
#   cat("Fit individual datasets\n")
#   ef[[i]] = list()
#   
#   for (d in names(data_list)){
#     
#     cat("...Data type: ",d,"\n")
#     
#     x_train = data_list[[d]][cv_list$outer$train[[i]],]
#     #x_test = data_list[[d]][cv_list$outer$test[[i]],]
#     
#     set.seed(seed)
#     fit <- caret::train(x = x_train,
#                         y = y_train,
#                         method="glmnet", 
#                         metric="ROC",
#                         tuneLength = 20,
#                         #weights = weights, 
#                         maximize = T,
#                         trControl=my_control,
#                         importance = TRUE)
#     ef[[i]][[d]] = fit
#     
#     
#   }
#   
#   cat("Fit meta model and computing performance\n")
#   preds = list()
#   vars = list()
#   for (a in names(ef[[i]])){
#     preds[[a]] = arrange(ef[[i]][[a]]$pred, Resample, rowIndex)$One
#     # vars[[a]] = varImp(ef[[i]][[a]])$importance #%>% mutate(Importance = (.$One + .$Zero)/2) %>% dplyr::select(-One,-Zero) 
#     # rownames(vars[[a]]) = rownames(varImp(ef[[i]][[a]])$importance)
#   }
#   preds = as.data.frame(preds)
#   y_train_truth = y_train[arrange(ef[[i]][[1]]$pred, Resample, rowIndex)$rowIndex] 
#   y_train_truth = ifelse(y_train_truth == "One", 1, 0)
#   y_train_best = abs(preds - y_train_truth) %>% apply(.,1,which.min)
#   y_train_best = colnames(preds)[y_train_best]
#   cat("...Best datasets include:\n")
#   print(table(y_train_best))
#   set.seed(seed)
#   fit_meta = caret::train(#x = preds,
#     #y = y_train[arrange(ef[[i]][[1]]$pred, Resample, rowIndex)$rowIndex],
#     x = x_train_i[arrange(ef[[i]][[1]]$pred, Resample, rowIndex)$rowIndex,],
#     y = y_train_best %>% as.factor(),
#     method="glmnet", 
#     metric="AUC",
#     tuneLength = 20, 
#     trControl=trainControl(
#       method="boot",
#       number=10,
#       savePredictions="final",
#       classProbs=TRUE,
#       summaryFunction=multiClassSummary
#     ),
#     importance = TRUE)
#   
#   pred_mods <- predict(ef[[i]], x_test_i, type = "prob") %>% as.data.frame()
#   keep_col = grep("One", colnames(pred_mods))
#   pred_mods = pred_mods[,keep_col, drop = F]
#   colnames(pred_mods) = names(ef[[i]])
#   # pred_ensemble = predict(fit_meta, pred_mods, s = "lambda.min", type = "prob")
#   # pred_mods$ensemble = pred_ensemble$One
#   pred_weight = predict(fit_meta, x_test_i, s = "lambda.min", type = "prob")
#   pred_mods_fil = pred_mods[,colnames(pred_weight)]
#   # pred_best = predict(fit_meta, x_test_i, s = "lambda.min", type = "raw") %>% as.character()
#   # pred_mods$ensemble = lapply(1:nrow(pred_mods), function (x) pred_mods[x,pred_best[x]]) %>% unlist()
#   pred_mods$concatenation = pred_cat$One
#   pred_mods$ensemble = apply(pred_mods_fil*pred_weight,1,sum)
#   # y_test = ifelse(y_test == "One",1,0)
#   perf_roc = caTools::colAUC(pred_mods, y_test)
#   perf_pr = lapply(pred_mods, function(x) MLmetrics::PRAUC(x, ifelse(y_test == "One",1,0))) %>% as.data.frame()
#   perf_ll = lapply(pred_mods, function(x) LogLoss(x, ifelse(y_test == "One",1,0), weights = weights_test)) %>% as.data.frame()
#   # cf = coef(fit_meta$finalModel, fit_meta$bestTune$mtry)
#   # weight = cf[-1]
#   # names(weight) = names(vars)
#   # vars_df = lapply(names(vars), function(x) vars[[x]]*abs(weight[x])) %>% do.call(rbind, .)
#   # var_imp = varImp(fit_meta)$importance #%>% mutate(Importance = (.$One + .$Zero)/2) %>% dplyr::select(-One,-Zero)
#   # rownames(var_imp) = names(vars)
#   # vars_df = lapply(names(vars), function(x) vars[[x]]*abs(var_imp[x,1])) %>% do.call(rbind, .)
#   # vars_df$Importance = vars_df$Importance/sum(vars_df$Importance)*100
#   # vars_df$Overall = vars_df$Overall/sum(vars_df$Overall)*100
#   
#   # cat("...Meta model and computing performance\n")
#   # class(ef[[i]]) <- c("caretList")
#   # x_test_i = df[cv_list$outer$test[[i]],]
#   # pred_mods <- predict(ef[[i]], x_test_i) %>% as.data.frame()
#   # colnames(pred_mods) = names(ef[[i]])
#   # model_ensemble <- caretStack(
#   #   ef[[i]],
#   #   method = "glmnet",
#   #   metric="ROC",
#   #   trControl=trainControl(
#   #     method = "boot",
#   #     number=10,
#   #     savePredictions="final",
#   #     summaryFunction=twoClassSummary,
#   #     classProbs=TRUE
#   #   ))
#   # pred_ensemble = predict(model_ensemble$ens_model, pred_mods, type = "prob")
#   # pred_mods$ensemble = pred_ensemble$One 
#   # perf = caTools::colAUC(pred_mods, y_test)
#   # var_imp = varImp(model_ensemble)
#   
#   ef[[i]] = list(mod_list = ef[[i]] ,
#                  mod_meta = fit_meta,
#                  pred_mods = pred_mods,
#                  perf = list(AUROC = perf_roc, PR = perf_pr, logloss = perf_ll),
#                  #vars = vars_df,
#                  weight = pred_weight
#                  #coef = cf
#   )
#   
#   
#   
#   cat("Done\n")
#   return(ef)
# }