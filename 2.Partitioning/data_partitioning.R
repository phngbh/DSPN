#!/usr/bin/env Rscript

#Load processed data
gen = read.csv("/1.Processing/Genomics/genomics_processed.fam", header = F, sep = "\t")
tra = readRDS("/1.Processing/Transcriptomics/transcriptomics_processed.rds")
pro = readRDS("/1.Processing/Proteomics/proteomics_processed.rds")
meta = readRDS("/1.Processing/Metabolomics/metabolomics_processed.rds")
meth = readRDS("/1.Processing/Methylomics/methylomics_processed.rds")
oli = readRDS("/1.Processing/Clinical/olink_processed.rds")
cli = readRDS("/1.Processing/Clinical/clinical_processed.rds")

#Extract sample IDs
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

#Extract and compute outcome labels
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

##Example data partitioning for incident DSPN prediction (using MNSI threshold of 3)

#Make training sample list
overlapped_ids = na.omit(sampleIDs) %>% dplyr::select(Clinical) %>% as.character()
inc3_samples = filter(mnsi, !is.na(inc3))
train_samples = intersect(overlapped_ids, rownames(inc3_samples))

#Make feature selection sets
#Genomics
gen_ids_tot = as.character(gen$V1)
gen_ids_inc3 = sampleIDs$Genomics[match(rownames(inc3_samples), as.character(sampleIDs$Clinical))] %>% na.omit() %>% as.character()
gen_ids_train = sampleIDs$Genomics[match(train_samples, as.character(sampleIDs$Clinical))] %>% as.character()
gen_ids_ftrsel = intersect(gen_ids_inc3, setdiff(gen_ids_tot, gen_ids_train))
inc3_gen = mnsi$inc3[as.character(sampleIDs$Clinical)[match(gen_ids_ftrsel, as.character(sampleIDs$Genomics))]]
rownames(inc3_gen) = gen_ids_ftrsel
gen_ids_resamples = createDataPartition(y = inc3_gen, times = 100, p = 0.8)

#Transcriptomics
tra_ids_tot = rownames(tra)
tra_ids_inc3 = sampleIDs$Transcriptomics[match(rownames(inc3_samples), as.character(sampleIDs$Clinical))] %>% na.omit() %>% as.character()
tra_ids_train = sampleIDs$Transcriptomics[match(train_samples, as.character(sampleIDs$Clinical))] %>% as.character()
tra_ids_ftrsel = intersect(tra_ids_inc3, setdiff(tra_ids_tot, tra_ids_train))
inc3_tra = mnsi$inc3[as.character(sampleIDs$Clinical)[match(tra_ids_ftrsel, as.character(sampleIDs$Transcriptomics))]]
rownames(inc3_tra) = tra_ids_ftrsel
tra_ids_resamples = createDataPartition(y = inc3_tra, times = 100, p = 0.8)

#Proteomics
pro_ids_tot = rownames(pro)
pro_ids_inc3 = sampleIDs$Proteomics[match(rownames(inc3_samples), as.character(sampleIDs$Clinical))] %>% na.omit() %>% as.character()
pro_ids_train = sampleIDs$Proteomics[match(train_samples, as.character(sampleIDs$Clinical))] %>% as.character()
pro_ids_ftrsel = intersect(pro_ids_inc3, setdiff(pro_ids_tot, pro_ids_train))
inc3_pro = mnsi$inc3[as.character(sampleIDs$Clinical)[match(pro_ids_ftrsel, as.character(sampleIDs$Proteomics))]]
rownames(inc3_pro) = pro_ids_ftrsel
pro_ids_resamples = createDataPartition(y = inc3_pro, times = 100, p = 0.8)

#Metabolomics
meta_ids_tot = rownames(meta)
meta_ids_inc3 = sampleIDs$Metabolomics[match(rownames(inc3_samples), as.character(sampleIDs$Clinical))] %>% na.omit() %>% as.character()
meta_ids_train = sampleIDs$Metabolomics[match(train_samples, as.character(sampleIDs$Clinical))] %>% as.character()
meta_ids_ftrsel = intersect(meta_ids_inc3, setdiff(meta_ids_tot, meta_ids_train))
inc3_meta = mnsi$inc3[as.character(sampleIDs$Clinical)[match(meta_ids_ftrsel, as.character(sampleIDs$Metabolomics))]]
rownames(inc3_meta) = meta_ids_ftrsel
meta_ids_resamples = createDataPartition(y = inc3_meta, times = 100, p = 0.8)

#Methylomics
meth_ids_tot = rownames(meth)
meth_ids_inc3 = sampleIDs$Methylomics[match(rownames(inc3_samples), as.character(sampleIDs$Clinical))] %>% na.omit() %>% as.character()
meth_ids_train = sampleIDs$Methylomics[match(train_samples, as.character(sampleIDs$Clinical))] %>% as.character()
meth_ids_ftrsel = intersect(meth_ids_inc3, setdiff(meth_ids_tot, meth_ids_train))
inc3_meth = mnsi$inc3[as.character(sampleIDs$Clinical)[match(meth_ids_ftrsel, as.character(sampleIDs$Methylomics))]]
rownames(inc3_meth) = meth_ids_ftrsel
meth_ids_resamples = createDataPartition(y = inc3_meth, times = 100, p = 0.8)

#Clinical data
cli_ids_tot = rownames(cli)
cli_ids_ftrsel = intersect(inc3_samples, setdiff(cli_ids_tot, train_samples))
inc3_cli = inc3_samples$inc3[cli_ids_ftrsel]
cli_ids_resamples = createDataPartition(y = inc3_cli, times = 100, p = 0.8)

sample_partition = list(train = train_samples, featureSel = list(Genomics = gen_ids_ftrsel, Transcriptomics = tra_ids_ftrsel, Proteomics = pro_ids_ftrsel, Metabolomics = meta_ids_ftrsel, Methylomics = meth_ids_ftrsel))
saveRDS(sample_partition, file = "sample_partition.rds")