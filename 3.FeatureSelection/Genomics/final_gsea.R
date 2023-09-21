#!/usr/bin/env Rscript

library(dplyr, lib.loc = "/lustre/groups/cbm01/workspace/phong.nguyen/DSPN/Rlib")
library(RobustRankAggreg, lib.loc = "/lustre/groups/cbm01/workspace/phong.nguyen/DSPN/Rlib")
library(fgsea, lib.loc = "/lustre/groups/cbm01/workspace/phong.nguyen/DSPN/Rlib")
library(reactome.db, lib.loc = "/lustre/groups/cbm01/workspace/phong.nguyen/DSPN/Rlib")

args = commandArgs(trailingOnly=TRUE)

pwlist <- list()
for (i in 1:100){
  res <- read.table(paste0("Genomics/magma_geneset/",as.character(args[[1]]),"_",i,".gsa.out"), header = T, sep = "", skip = 4) 
  res$adjP <- p.adjust(res$P, method = "fdr")
  res <- arrange(res, adjP) %>% filter(adjP < 0.1)
  pwlist[[i]] <- as.character(res$FULL_NAME) 
}

pwlist <- pwlist[lapply(pwlist,length)>0]
set.seed(993)
pwlist_agg = aggregateRanks(pwlist)
pwlist_agg$adjP = pwlist_agg$Score*length(pwlist)
pwlist_agg$adjP = p.adjust(pwlist_agg$adjP, method = "fdr")
toppw = filter(pwlist_agg, adjP < 0.05)$Name %>% as.character()

write.table(data.frame(V1 = toppw), paste0("toppw_",as.character(args[[1]]),".txt"), col.names = F, row.names = F, quote = F)

gene_df <- read.table(paste0("Genomics/magma_gene/",as.character(args[[1]]),".genes.out"), header = T, sep = "") %>%
  arrange(ZSTAT) 
ranklist <- gene_df$ZSTAT
names(ranklist) <- gene_df$GENE
geneset_reactome = reactomePathways(names(ranklist))
names(geneset_reactome) <- gsub(" ","_", names(geneset_reactome))
geneset_reactome = geneset_reactome[intersect(names(geneset_reactome),toppw)]
set.seed(993)
fgseaRes <- fgsea(pathways = geneset_reactome,
                        stats    = ranklist,
                        #minSize  = 15,
                        maxSize  = 200,
                  nperm = 5000) %>% arrange(pval)
edge_genes <- fgseaRes$leadingEdge %>% unlist() %>% unique()

write.table(data.frame(V1 = edge_genes), paste0("Genomics/leadingEdge_",as.character(args[[1]]),".txt"), col.names = F, row.names = F, quote = F)

gene_snps <- read.table("Genomics/gene_annot.txt") %>% filter(V1 %in% edge_genes)
write.table(data.frame(V1 = unique(gene_snps$V2)), paste0("Genomics/leadingEdge_snps_",as.character(args[[1]]),".txt"), col.names = F, row.names = F, quote = F)

