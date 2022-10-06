#!/usr/bin/env Rscript

gene_snp = read.csv("gene_annot.txt", sep = "", check.names = F, header = F)
geneset_reactome = reactomePathways(levels(as.factor(as.character(gene_snp$V1))))
lowest_level_pathways = readRDS("lowest_level_pathways.rds")
geneset_reactome = geneset_reactome[intersect(names(geneset_reactome), lowest_level_pathways)]
for (p in names(geneset_reactome)){
  gen = intersect(levels(as.factor(as.character(gene_snp$V1))), geneset_reactome[[p]])
  df = filter(gene_snp, V1 %in% as.integer(gen))
  na_genes = setdiff(geneset_reactome[[p]], levels(as.factor(as.character(gene_snp$V1))))
  if(length(gen) > 0){
    geneset_reactome[[p]] = c(df$V2 %>% unique(), na_genes)
  } else {
    geneset_reactome[[p]] = NULL
  }
}

gsea = list()
for (i in 1:100){
  cat("Iter",i,"\n")
  gwas_tmp = read.csv(paste0("./gwas/samples_",i,".inc3.assoc.logistic"), sep = "") %>% 
    filter(SNP %in% gene_snp$V2)
  gwas_tmp$GENE = gene_snp$V1[match(gwas_tmp$SNP,gene_snp$V2)]
  gwas_tmp_unq = gwas_tmp %>% distinct(STAT, GENE, .keep_all = T)
  ranklist = gwas_tmp_unq$STAT
  names(ranklist) = gwas_tmp_unq$SNP
  ranklist = sort(ranklist)
  set.seed(993)
  fgseaRes_tmp <- fgsea(pathways = geneset_reactome,
                        stats    = ranklist,
                        eps = 0,
                        minSize  = 15,
                        maxSize  = 500) %>% arrange(pval) %>% filter(padj < 5e-8)
  if(nrow(fgseaRes_tmp) == 0){
    message("No significant pathway found\n")
    next
  }
  fgseaRes_tmp$leadingEdge_gen = lapply(fgseaRes_tmp$leadingEdge, 
                                          function(x)
                                            return(gwas_tmp_unq$GENE[match(x,gwas_tmp_unq$SNP)] %>%
                                                     unique()))
  gsea[[i]] = fgseaRes_tmp
}

pathways = list()
for (i in 1: length(gsea)){
  if(is.null(gsea[[i]])){
    next
  }
  pathways[[i]] = gsea[[i]]$pathway
}
pwlist_agg = aggregateRanks(pathways)
pwlist_agg$adjP = pwlist_agg$Score*length(pathways) 
pwlist_agg$adjP = p.adjust(pwlist_agg$adjP, method = "fdr")
toppw = rownames(pwlist_agg[pwlist_agg$adjP < 0.05,])

gwas = read.csv(paste0("./gwas/gwas_all.inc3.assoc.logistic"), sep = "") %>% 
    filter(SNP %in% gene_snp$V2)
gwas$GENE = gene_snp$V1[match(gwas$SNP,gene_snp$V2)]
gwas_unq = gwas_mnsi3 %>% distinct(STAT, GENE, .keep_all = T)
ranklist = gwas_unq$STAT
names(ranklist) = gwas_unq$SNP
ranklist = sort(ranklist)
set.seed(993)
fgseaRes <- fgsea(pathways = geneset_reactome[toppw],
                    stats    = ranklist,
                    eps = 0,
                    minSize  = 15,
                    maxSize  = 500) %>% arrange(pval) 
fgseaRes$leadingEdge_gen = lapply(fgseaRes$leadingEdge, 
                                      function(x)
                                        return(gwas_mnsi3_unq$GENE[match(x,gwas_mnsi3_unq$SNP)] %>%
                                                 unique()))
kable(fgseaRes)
fgseaRes$leadingEdge_s = lapply(fgseaRes$leadingEdge_gen, 
                                      function(x) sort(x) %>% paste0(collapse = " ")) %>% unlist()
fgseaRes$similarity = sapply(fgseaRes$leadingEdge_s,
                                   function(x) agrep(x,fgseaRes$leadingEdge_s, max=0.2))
fgseaRes_fil = distinct(fgseaRes, similarity, .keep_all = T)

edge = fgseaRes$leadingEdge %>% unlist() %>% unique()
write.table(data.frame(V1 = edge), file = "edge_snp.txt", col.names = F, row.names = F, quote = F)