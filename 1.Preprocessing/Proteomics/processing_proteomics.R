#!/usr/bin/env Rscript

dat = read.csv("/data/somalogicKF4_somamer_data_060617.csv", check.names = F, sep = ";", row.names = 1)
info_s = read.csv("data/somalogicKF4_sample_info_060617.csv", check.names = F, sep = ";")
info_p = read.csv("somamer_info_edited.csv", check.names = F, sep = ";")

drop = c("2795-23_3", "3590-8_3", "5071-3_3", "5073-30_2", "5118-74_2")
info_p = dplyr::filter(info_p, ColCheck == "PASS" & !SeqId %in% drop)
info_s = dplyr::filter(info_s, RowCheck == "PASS")
info_p$Target[info_p$Target == "14.03.2003"] = "14-3-3"

clin <- read_sas("data/gesamt_k05720g_v2007.sas7bdat", NULL)
info_s$mnsi = clin$utmnsi[match(as.numeric(info_s$SampID), as.numeric(clin$un_protSoma_f4))]

dat = dat[as.character(info_s$SampID), info_p$SeqId] %>% t() %>% log2()
rownames(info_p) = rownames(dat)
saveRDS(dat, file = "proteomics_processed.rds")