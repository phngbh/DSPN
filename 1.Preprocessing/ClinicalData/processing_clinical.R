#!/usr/bin/env Rscript

clin <- read_sas("/data/gesamt_k05720g_v2007.sas7bdat", NULL) %>% 
  filter(!utdm100y15 %in% c(6,7,8)) %>% droplevels()
covar = clin %>% dplyr::select(utalter, ucsex, utgroe, uttumf, uthyact, ul_chola, utmdiu, utmata, utmglyk, utmepilep_vitd,
                        ul_hbava, utcigreg, utalkcat, utphys, utmi, utschl, utmoadi, utmcablo, utmhypot, utmadins,
                        utgfr_ckd_crcc, un05_1, utmansaid, utbmi, utdm100y15,utmeddia,utmbbl, utmace, 
                        utmhypol, utmstati, utmfibra, utmahypol, utlipi, utmantco, utmantpl, utmnitro,zz_nr) %>%
  apply(.,2, as.factor) %>%
  as.data.frame() %>%
  transform(., utalter = as.numeric(utalter), utgroe = as.numeric(utgroe), utbmi = as.numeric(utbmi),
            uttumf = as.numeric(uttumf), zz_nr = as.numeric(zz_nr), ul_chola = as.numeric(ul_chola),
            ul_hbava = as.numeric(ul_hbava), utgfr_ckd_crcc = as.numeric(utgfr_ckd_crcc)) %>%
  column_to_rownames("zz_nr") %>%
  model.matrix(~ 0 + .,.)
phenotype = clin %>% dplyr::select(utglukfast_a, utglukfast_n, utglukrand_a, utglukrand_n, 
                            utgluk2a, utgluk2n, utgfr_ckd_cr, utgfr_ckd_cc, 
                            utgfr_ckd_crcc, ul_chola, ul_choln, ul_hdla, ul_hdln, ul_ldla, 
                            ul_ldln, ul_tria, ul_trin, ul_hbava, ul_hsrea, ul_hsren, ul_kreaa, 
                            ul_krean, ul_vitd, ul_hsrea, ul_hsren, ul_kreaa, ul_krean, uh_ins, 
                            uh_ins_2h, uh_adipo, uh_adipo_2h, uh_il1ra, uh_il1ra_2h, uh_tnfa, 
                            uh_il6, uh_sicam1, uh_omentin, uh_crp, uh_il18, uh_cystatinc, 
                            uh_rbp, uh_leptin_edta, uh_sfrp5, uh_wnt5a, uh_vitd_f, uh_alb, 
                            uh_pth1_84, uh_dbp, uh_crp, uh_il18, uh_sod3, uh_mpo, uh_sod3, zz_nr)
keep_pheno = apply(phenotype, 2, function(x) sum(!is.na(x)) > nrow(phenotype)/10*9)
phenotype = phenotype[,keep_pheno] %>% column_to_rownames("zz_nr") %>% na.omit()
dat = cbind(covar[intersect(rownames(covar), rownames(phenotype)),], phenotype[intersect(rownames(covar), rownames(phenotype)),])
saveRDS(dat, file = "clinical_processed.rds")

olink_inf = clin[, grep("uh_o_", colnames(clin))] %>% as.data.frame()
rownames(olink_inf) = clin$zz_nr
olink_inf = na.omit(olink_inf)
saveRDS(olink_inf, file = "olink_processed.rds")

