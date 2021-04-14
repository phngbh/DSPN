#if(!"snpStats" %in% rownames(installed.packages())) BiocManager::install("snpStats", lib = "./Rlib/")
library(tidyverse)
library(snpStats, lib.loc = "./Rlib/")
library(caret, lib.loc = "./Rlib/")
library(haven)
library(limma)
library(ChAMP)
library(glmnet)
#library(glmnetUtils, lib.loc = "./Rlib/")
library(DMwR, lib.loc = "Rlib")
library(MLmetrics, lib.loc = "Rlib")
library(randomForest)

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
  o = NULL #String, name of omics comb
  
){

  feature_mat = as.data.frame(feature_mat)

  if (class(phenotype) == "numeric"){
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

  } else if (class(phenotype) == "character") {
    
    if(nlevels(as.factor(phenotype)) > 2){
      phenotype[phenotype == "2"] <- "1"
    }
    
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
      
      try = myTryCatch(champ.DMP(beta = t(feature_mat),pheno=phenotype, adjPVal = 0.5))
      
      if( !is.null(try$error)){
        feature_mat_final <- NULL
        
      } else {
        myDMP <- champ.DMP(beta = t(feature_mat),pheno=phenotype, adjPVal = 0.5)
        topdm = rownames(myDMP[[1]])
        
        if (length(topdm) < 1000){
          feature_mat_final = feature_mat[,rownames(topdm)]
        } else{
          feature_mat_fil = feature_mat[,rownames(topdm)[1:1000]]
          
          ft.cor = cor(feature_mat_fil)
          hc = findCorrelation(ft.cor, cutoff = 0.9)
          if (length(hc) > 0){
            feature_mat_fil = feature_mat_fil[,-hc]
            if (nrow(topdm) >= 1000 + length(hc)){
              feature_mat_final = cbind(feature_mat_fil,feature_mat[,rownames(topdm)[1001:(1000+length(hc))],drop = F])
            } else {
              feature_mat_final = cbind(feature_mat_fil,feature_mat[,rownames(topdm)[1001:nrow(topdm)],drop = F])
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
  feature_sel = T #Whether to run feature_sel()
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
        mat = feature_sel(feature_mat = ftr_mat, phenotype = phenotype[,1], o = i)
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
  seed = 993
){
  
  control <- trainControl(method="none")
  
  if (class(y_train[,1]) == "numeric"){
    
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
    diff <- pred - y_test
    pearson <- cor.test(pred,y_test, method = "pearson", exact = F)$estimate
    spearman <- cor.test(pred,y_test, method = "spearman", exact = F)$estimate
    corr = list(pearson = pearson, spearman = spearman)
    
    return(list(pred=pred, diff=diff, corr = corr, fit=fit))
  } 
  
  if (class(y_train[,1]) == "character"){
    
    set.seed(seed)
    tune = tuneRF(x_train, as.factor(y_train[,1]), stepFactor = 1.5, improve = 1e-5, tree = 500,trace = F, plot = F) %>% as.data.frame()
    mtry = tune$mtry[which(tune$OOBError == min(tune$OOBError))] %>% as.numeric() %>% max()
    tunegrid <- expand.grid(.mtry=mtry)
    
    set.seed(seed)
    fit <- caret::train(x = x_train,
                        y = y_train,
                        method="rf", 
                        metric="Accuracy",
                        #tuneLength = 20,
                        tuneGrid=tunegrid, 
                        trControl=control,
                        importance = TRUE)
    
    pred <- predict(fit, x_test)
    pred <- as.data.frame(pred)[,1]
    acc = mean(pred == y_test, na.rm = T)
    
    return(list(pred=pred, obs=y_test, acc = acc, fit=fit))
  }

}

use_glm <- function(
  ###Train a linear model for each omics combination and phenotype 
  x_train=NULL,
  y_train=NULL,
  x_test=NULL,
  y_test=NULL,
  seed = 993
){
  
  
  if(class(y_train[,1]) == "numeric"){
    print("...Numeric response")
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
    diff <- pred - y_test
    pearson <- cor.test(pred,y_test, method = "pearson", exact = F)$estimate
    spearman <- cor.test(pred,y_test, method = "spearman", exact = F)$estimate
    corr = list(pearson = pearson, spearman = spearman)
    
    return(list(pred=pred, diff=diff, corr = corr, fit=fit))
    
  }
  
  if(class(y_train[,1]) == "character"){
    control <- trainControl(method="repeatedcv", number=5, repeats=2, summaryFunction = twoClassSummary, classProbs = TRUE, sampling = "smote")
    # weights <- ifelse(y_train[,1] == "One",
    #                   (1/table(as.factor(y_train[,1]))[1]) * 0.5,
    #                   (1/table(as.factor(y_train[,1]))[2]) * 0.5)
    set.seed(seed)
    fit <- caret::train(x = x_train,
                        y = y_train,
                        method = "glmnet",
                        metric = "Sens",
                        #weights = weights,
                        tuneLength = 20,
                        trControl = control,
                        na.action = na.omit)
    
    # set.seed(seed)
    # fit <- glmnet::cv.glmnet(x = as.matrix(x_train),
    #                          y = y_train[,1] %>% factor(levels = c("One","Zero")),
    #                          family = "binomial",
    #                          type.measure = "auc",
    #                          relax = T)
    
    pred <- predict(fit, x_test, s = 'lambda.min') %>% factor(levels = c("One","Zero"))
    obs = y_test %>% factor(levels = c("One","Zero"))
    sens = caret::sensitivity(pred,obs)
    spec = caret::specificity(pred,obs)
    perf = list(sens = sens, spec = spec)
    
    return(list(pred=pred, obs=y_test, perf = perf, fit=fit))
    
  }
  
}


make_fit = function(
  ###Train and return CV score if cv==T, linear models if cv==F
  train_list = NULL, #list of training data
  test_list = NULL, #list of testing data
  neuropathy = NULL, #df of responses
  cv_list = NULL, #list of CV subsets
  omics = NULL,
  covar = NULL,
  lab = NULL,
  method = c("elnet","rf","xgboost"),
  cv = TRUE,
  n = 1 #index of omics combination  
){
  
  # ncores = detectCores() - 1
  # cl = makeCluster(ncores)
  # registerDoParallel(cl)
  # clusterCall(cl, function(x) .libPaths(x), .libPaths())
  
  o = names(train_list)[n]
  tmp = apply(neuropathy, 2, as.character) %>% 
    as.data.frame() %>% 
    transform(., un14 = as.numeric(un14),
              un15 = as.numeric(un15),
              un17 = as.numeric(un17),
              un18 = as.numeric(un18),
              utmnsi = as.numeric(utmnsi)) %>%
    dplyr::select(-c(un12_3,un12_4,un12_8,un12_9,un13_3,un13_4,un13_8,un13_9))
  rownames(tmp) <- rownames(neuropathy)
  neuropathy <- tmp
  
  #for(p in colnames(neuropathy))  {
   p = "un14" 
    print(paste0("Processing omics ",o," for neuropathy variable ",p))
    
    #Process train set
    pheno = neuropathy[,p,drop=F] %>% as.matrix()
    rownames(pheno) = rownames(neuropathy)
    pheno = na.omit(pheno)
    clin_id = train_list[[o]][["Samples"]][["Clinical"]]
    pheno = pheno[intersect(rownames(pheno), as.character(clin_id)),,drop = F]
    ind = match(rownames(pheno), as.character(clin_id))
    sample_list = train_list[[o]][["Samples"]]
    sample_list = lapply(names(sample_list), function(x) return(sample_list[[x]][ind]))
    names(sample_list) = names(train_list[[o]][["Samples"]])
    feature_list = train_list[[o]][["Features"]]
    
    if (cv){
      cvglm_list = vector(mode = "list", length = 10) 
      
      print("...Cross-validating")
      for(k in 1:10) {
        
        print(paste("...Fold",k))
        
        if(length(grep("Genomics",o)) > 0){
          genomics_k = bed_to_df(paste0("./Genomics/small_bed_neuro/",o,"_",p,"_",k,".bed"))
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
                                        phenotype = pheno_k) %>% na.omit()
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
                                             features = colnames(feature_mat_k)) %>% na.omit()
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
                                     seed = 993)
        }

        if (method == "rf") {
          cvglm_list[[k]] = use_rf( x_train=feature_mat_k,
                                     y_train=pheno_k,
                                     x_test=feature_mat_k_test,
                                     y_test=pheno_k_test,
                                     seed = 993)
        }

        rm(feature_mat_k)
        rm(feature_mat_k_test)
       }

      if (method == "elnet"){
        saveRDS(cvglm_list, file = paste0("./test/",o,",",p,",cv.rds"))
      } else if (method == "rf"){
        saveRDS(cvglm_list, file = paste0("./test/",o,",",p,",cv.rds"))
      }

      rm(cvglm_list)
     
    } else {
      
      if(length(grep("Genomics",o)) > 0){
        genomics = bed_to_df(paste0("./Genomics/small_bed_neuro/",o,"_",p,".bed"))
      } else {
        genomics <- NULL
      }
      
      feature_mat = dimnames_to_mat(feature_list, sample_list, 
                                    genomics = genomics,
                                    omics = omics,
                                    covar = covar,
                                    lab = lab,
                                    phenotype = pheno) %>% na.omit()
      pheno = pheno[rownames(feature_mat),,drop=F]
      #pheno = ifelse(pheno == "1","One","Zero") 
      
      #Process test set
      test_pheno = neuropathy[,p,drop=F] %>% as.matrix()
      rownames(test_pheno) = rownames(neuropathy)
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
                                   features = colnames(feature_mat)) %>% na.omit()
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
                       seed = 993)

        saveRDS(pred, file = paste0("./test/",o,",",p,",mod.rds"))

      } else if (method == "rf"){
        pred = use_rf(x_train=feature_mat,
                       y_train=pheno,
                       x_test = test_omics,
                       y_test = test_pheno,
                       seed = 993)

        saveRDS(pred, file = paste0("./test/",o,",",p,",mod.rds"))
      }

      rm(pred)
      rm(test_omics)
    }
    
    #rm(genomics)
    
  #}
  
  #stopCluster(cl)
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
