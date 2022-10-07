#!/usr/bin/env Rscript 

#Get sample IDs and phenotype files
clin <- read_sas("/data/gesamt_k05720g_v2007.sas7bdat", NULL) %>% 
  filter(!utdm100y15 %in% c(6,7,8)) %>% droplevels()
neuro = clin[,c(grep("un12",colnames(clin)),grep("un13",colnames(clin)))] %>%
  cbind(dplyr::select(clin, "un14","un15","un17","un18","un19","un20","utmnsi","u3tmnsi"))
rownames(neuro) = clin$zz_nr
mnsi = neuro %>% 
  rownames_to_column("id") %>%
  dplyr::select(id, utmnsi, u3tmnsi) %>% 
  filter(!is.na(utmnsi)) %>% 
  mutate(utmnsi2 = ifelse(utmnsi <= 2, 0, 1), 
         utmnsi3 = ifelse(utmnsi <= 3, 0, 1),
         u3tmnsi2 = ifelse(is.na(u3tmnsi), NA, ifelse(u3tmnsi <= 2, 0, 1)),
         u3tmnsi3 = ifelse(is.na(u3tmnsi), NA, ifelse(u3tmnsi <= 3, 0, 1))) %>%
  mutate(inc2 = ifelse(is.na(u3tmnsi2) | utmnsi2 == 1, NA, ifelse(utmnsi2 == 0 & u3tmnsi2 == 0, 0, 1)),
         inc3 = ifelse(is.na(u3tmnsi3) | utmnsi3 == 1, NA, ifelse(utmnsi3 == 0 & u3tmnsi3 == 0, 0, 1)),
         prog = u3tmnsi - utmnsi) %>%
  remove_rownames() %>%
  column_to_rownames("id") %>%
  dplyr::select(-c(utmnsi,u3tmnsi, u3tmnsi2, u3tmnsi3))
mnsi = filter(mnsi, !is.na(inc3))
dat = readRDS("/1.Processing/Clinical/clinical_processed.rds")
inc = mnsi$inc3
names(inc) = rownames(mnsi)
inc = na.omit(inc)
inc = ifelse(inc == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
sample_list = readRDS("/2.Partitioning/sample_partition.rds")
resamples = sample_list$featureSel$Clinical

#Do elastic net
sampling = NULL
sumFunc = caretLogLoss
metric = "myLogLoss"
maximize = F
  
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

perf_list = list()
var = list()
for (i in 1:length(resamples)){
  cat("Fold",i,"\n")
  test_id = setdiff(1:nrow(dat), resamples[[i]])
  x_train = dat[resamples[[i]],]
  y_train = inc[resamples[[i]]]
  x_test = dat[test_id,] 
  y_test = inc[test_id]
  
  weights = ifelse(y_train == "One", table(y_train)[[2]]/table(y_train)[[1]], 1)
  
  if (any(rownames(x_train) != names(y_train)) | any(rownames(x_test) != names(y_test))){
      message("Samples in train and test sets do not match!")
      next
    }
  
  set.seed(993)
  fit <- caret::train(x = x_train,
                        y = y_train,
                        method="glmnet", 
                        metric=metric,
                        tuneLength = 20,
                        weights = weights,
                        maximize = maximize,
                        trControl=my_control,
                        importance = TRUE)
  pred = predict(fit, x_test, s = "lambda.min", type = "prob")$One
  roc <- roc(response = y_test, predictor = pred, levels = c("Zero","One"))
  auc = auc(roc)
  perf_list[[i]] = auc
  var_imp = varImp(fit)$importance
  var = var_imp$Overall
  names(var) = rownames(var_imp)
  var = var[order(var, decreasing = T)]
  #var = var[var != 0]
  var[[i]] = var
}

#Select best features
var_fil = lapply(var, function(x) names(x[x>0]))
var_rra = aggregateRanks(var_fil)
var_rra$adjP = var_rra$Score*100
var_rra$adjP = p.adjust(var_rra$adjP, "fdr")
var_sel = var_rra$adjP[var_rra$adjP < 0.05]
#var_sel = var_sel[-8]

saveRDS(dat[,names(var_sel)],"clinical_selected.rds")