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
resamples = sample_list$featureSel$Methylomics
samples = unlist(resamples) %>% unique()
info = mnsi[as.character(sampleIDs$Clinical)[match(samples, as.character(sampleIDs$Methylomics))],inc3,drop = F]
rownames(info) = samples

#Do gene set enrichment analysis
dat = readRDS("/1.Processing/Methylomics/methylomics_processed.rds")
lowest_level_pathways = readRDS("lowest_level_pathways.rds")
modmatrix = model.matrix(~ 0 + ., data=info)
fit = lmFit(dat, modmatrix)
contrast = makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contrast)
tmp <- eBayes(tmp)
topde <- topTable(tmp, sort.by = "P", n = Inf)
data("dualmap450kEID")
geneset_reactome = reactomePathways(names(mapEIDto450k.lv))
geneset_reactome = geneset_reactome[intersect(names(geneset_reactome), lowest_level_pathways)]
for (p in names(geneset_reactome)){
  probes = mapEIDto450k.lv[intersect(geneset_reactome[[p]], names(mapEIDto450k.lv))] %>% unlist() %>% unique()
  na_genes = setdiff(geneset_reactome[[p]], names(mapEIDto450k.lv))
  if(length(probes) > 0){
    geneset_reactome[[p]] = c(probes, na_genes)
  } else {
    geneset_reactome[[p]] = NULL
  }
}
ranklist = topde$t
names(ranklist) = rownames(topde)
ranklist = sort(ranklist)
set.seed(993)
fgseaRes <- fgsea(pathways = geneset_reactome,
                stats    = ranklist,
                minSize  = 15,
                maxSize  = 200) %>% arrange(pval) %>% filter(padj < 0.05)
fgseaRes[, leadingEdge := mapIdsList(
                                     x=org.Hs.eg.db, 
                                     keys=leadingEdge,
                                     keytype="ENTREZID", 
                                     column="SYMBOL")]
kable(fgseaRes)
saveRDS(fgseaRes, "gsea_final.rds")

#Extract selected features for training
edge_mnsi3 = fgseaRes$leadingEdge %>% unlist() %>% unique()
dat_selected = dat[edge_mnsi3,]
saveRDS(dat_selected,"methylomics_selected.rds")
