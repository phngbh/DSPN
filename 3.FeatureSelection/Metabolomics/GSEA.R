#!/usr/bin/env Rscript 

#Get sample IDs and phenotype files
clin <- read_sas("/data/gesamt_k05720g_v2007.sas7bdat", NULL)
clin = filter(clin, !utdm100y15 %in% c(6,7,8)) %>% droplevels()
sampleIDs = clin %>% dplyr::select(lg_dnaAxiom_s4f4, 
                           un_expr_f4ogtt,
                           un_protSoma_f4,
                           un_metabMetabolon_f4,
                           un_meth450k_f4,
                           zz_nr)
sampleIDs$Clinical = clin$zz_nr
colnames(sampleIDs) = c("Genomics", "Transcriptomics", "Proteomics",
                        "Metabolomics", "Methylomics", "Olink", "Clinical")
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
sample_list = readRDS("/2.Partitioning/sample_partition.rds")
resamples = sample_list$featureSel$Metabolomics
samples = unlist(resamples) %>% unique()
info = mnsi[as.character(sampleIDs$Clinical)[match(samples, as.character(sampleIDs$Metabolomics))],inc3,drop = F]
rownames(info) = samples

#Do gene set enrichment analysis
dat = readRDS("/1.Processing/Metabolomics/metabolomics_processed.rds")
pathway_list = readRDS("pathway_list.rds")
lowest_level_pathways = readRDS("lowest_level_pathways.rds")
meta_anno_fil = readRDS("meta_anno_fil.rds")
keep_meta = meta_anno_fil$Mnumber[!is.na(meta_anno_fil$ChEBI)]
modmatrix = model.matrix(~ 0 + ., data=info)
set.seed(993)
gsea = list(edge = list(), probe = list(), auc = list(), pathway = list())

for (i in 1:length(resamples)){
  cat("Iter ",i,"\n")
  dat_tmp = dat[keep_meta,resamples[[i]]]

  cat("...DE analysis\n")
  mod_tmp = modmatrix[resamples[[i]],]
  if( any(table(mod_tmp[,1]) < 2) | any(table(mod_tmp[,2]) < 2)){
    cat("......One of the factor levels has less than 2 observations => Stop!")
    next
  }
  fit = lmFit(dat_tmp, mod_tmp)
  contrast = makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contrast)
  tmp <- eBayes(tmp)
  topde <- topTable(tmp, sort.by = "P", n = Inf) %>% mutate(ID = rownames(.)) %>%
    mutate(Name = meta_anno_fil$Biochemical_F4[match(.$ID,meta_anno_fil$Mnumber)],
           ChEBI = meta_anno_fil$ChEBI[match(.$ID, meta_anno_fil$Mnumber)])
  cat("...GSEA\n")
  ranklist = topde$t
  names(ranklist) = topde$ChEBI
  ranklist = sort(ranklist)
  set.seed(993)
  fgseaRes_tmp <- fgsea(pathways = pathway_list,
                  stats    = ranklist,
                  minSize  = 1,
                  maxSize  = 200,
                  eps = 0) %>% arrange(pval) %>% filter(padj < 0.1)
  if (nrow(fgseaRes_tmp) == 0){
    message("No significant pathway found")
    next
  }
  edge_tmp = fgseaRes_tmp$leadingEdge %>% unlist() %>% unique()
  gsea$edge[[i]] = edge_tmp
  gsea$probe[[i]] = meta_anno_fil$Mnumber[match(edge_tmp, meta_anno_fil$ChEBI)]
  gsea$pathway[[i]] = fgseaRes_tmp
  
  cat("...Elastic net\n")
  
  if(is.null(gsea$ilmn[[i]])){
    message("No significant pathway found")
    next
  }
  
  probelist_tmp = gsea$probe[[i]]
  test_ind = setdiff(1:ncol(dat), resamples[[i]])
  x_train = dat[probelist_tmp,resamples[[i]]] %>% t()
  y_train = info[resamples[[i]],"inc3"]
  y_train = ifelse(y_train == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
  x_test = dat[probelist_tmp,test_ind] %>% t()
  y_test = info[test_ind,"inc3"]
  y_test = ifelse(y_test == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
  my_control <- trainControl(
    method="repeatedcv",
    number=5,
    repeats = 2,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=twoClassSummary,
    sampling = "smote",
    allowParallel = T
  )
  set.seed(993)
  fit <- caret::train(x = x_train,
                        y = y_train,
                        method="glmnet", 
                        metric="ROC",
                        tuneLength = 20,
                        maximize = T,
                        trControl=my_control,
                        importance = TRUE)
  pred = predict(fit, x_test, s = "lambda.min", type = "prob")$One
  roc <- roc(response = y_test, predictor = pred, levels = c("Zero","One"))
  auc = auc(roc)
  gsea$auc[[i]] = auc

}
saveRDS(gsea,"gsea_list.rds")

#Select the top significant pathways
keep = vector("logical",length = length(gsea$auc))
for (i in 1:length(keep)){
  if (is.null(gsea$auc[[i]])){
    next
  } else if (gsea$auc[[i]] <= 0.5){
    next
  } else {
    keep[i] = TRUE
  }
}

pwlist = gsea$pathway[keep]
pathways = list(pw = list(), pval = list(), NES = list())
for (i in 1: length(pwlist)){
  if(is.null(pwlist[[i]])){
    next
  }
  pathways$pw[[i]] = pwlist[[i]]$pathway
  pathways$pval[[i]] = pwlist[[i]]$pval
  pathways$NES[[i]] = pwlist[[i]]$NES
}
pwlist_agg = aggregateRanks(pathways$pw)
pwlist_agg$adjP = pwlist_agg$Score*length(pathways$pw)
pwlist_agg$adjP = p.adjust(pwlist_agg$adjP, method = "fdr")
toppw = rownames(filter(pwlist_agg, adjP < 0.05))

toppw_pval = list()
for (p in toppw){
  tmplist = list()
  for (i in 1:length(pwlist)){
    ind = which(pathways$pw[[i]] == p)
    if (length(ind) > 0){
      tmplist[[i]] = pathways$pval[[i]][ind]
    } else {
      tmplist[[i]] = NULL
    }
  }
  tmpvec = unlist(tmplist)
  toppw_pval[[p]] = mean(tmpvec)
}

#Final gene set enrichment analysis on the selected pathways
fit = lmFit(dat, modmatrix)
contrast = makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contrast)
tmp <- eBayes(tmp)
topde <- topTable(tmp, sort.by = "P", n = Inf)
topde <- topTable(tmp, sort.by = "P", n = Inf) %>% mutate(ID = rownames(.)) %>%
    mutate(Name = meta_anno_fil$Biochemical_F4[match(.$ID,meta_anno_fil$Mnumber)],
           ChEBI = meta_anno_fil$ChEBI[match(.$ID, meta_anno_fil$Mnumber)])
ranklist = topde$t
topde$ChEBI = meta_anno_fil$ChEBI[match(topde$ID, meta_anno_fil$Mnumber)]
names(ranklist) = topde$ChEBI
ranklist = sort(ranklist)
set.seed(993)
fgseaRes <- fgsea(pathways = geneset_reactome,
                stats    = ranklist,
                minSize  = 5,
                maxSize  = 200) %>% arrange(pval) #%>% filter(padj < 0.3)
fgseaRes[, leadingEdge := mapIdsList(
                                     x=org.Hs.eg.db, 
                                     keys=leadingEdge,
                                     keytype="ENTREZID", 
                                     column="SYMBOL")]
kable(fgseaRes)
saveRDS(fgseaRes, "gsea_final.rds")

#Extract selected features for training
edge = fgseaRes$leadingEdge %>% unlist() %>% unique()
probe = meta_anno_fil$Mnumber[meta_anno_fil$ChEBI %in% edge]
dat_selected = dat[probe,,drop=F]

saveRDS(dat_selected,"metabolomics_selected.rds")
