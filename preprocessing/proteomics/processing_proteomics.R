library(dplyr)
library(tidyverse)
library(haven)
library(ComplexHeatmap)

### Load data ######

dat = read.csv("somalogicKF4_somamer_data_060617.csv", check.names = F, sep = ";", row.names = 1)
info_s = read.csv("somalogicKF4_sample_info_060617.csv", check.names = F, sep = ";")
info_p = read.csv("somalogicKF4_somamer_info_060617_archive.csv", check.names = F, sep = ";")

drop = c("2795-23_3", "3590-8_3", "5071-3_3", "5073-30_2", "5118-74_2")
info_p = filter(info_p, ColCheck == "PASS" & !SeqId %in% drop)
info_s = filter(info_s, RowCheck == "PASS")
info_p$Target[info_p$Target == "14.03.2003"] = "14-3-3"

clin <- read_sas("gesamt_k05720g_v2007.sas7bdat", NULL)
#info_s$mnsi = clin$utmnsi[match(as.numeric(info_s$SampID), as.numeric(clin$un_protSoma_f4))]

dat = dat[as.character(info_s$SampID), info_p$SeqId] %>% t() %>% log2()

### Visualize ####

meansd = data.frame(mean = apply(dat,1,mean), sd = apply(dat,1,sd), 
                    iqr = apply(dat,1,IQR), med = apply(dat, 1, median))
ggplot(meansd, aes(x=mean,y=sd))+
  geom_point(col = "black",alpha=0.5)+
  stat_smooth(method = "loess",se=F,col="red")

pca = prcomp(t(dat), scale. = T)
pca_df = data.frame(pca$x, info_s)
pca_var = summary(pca)$importance
varp = round(pca_var[2,]*100,2)
ggplot(pca_df,aes(x=PC1,y=PC2,col = PlateId)) + 
  geom_point(size = 3, alpha = 0.5)+
  xlab(paste("PC1 (",varp[1],"%)")) +
  ylab(paste("PC2 (",varp[2],"%)")) +
  theme(panel.background = element_rect(fill = "white"))

ggplot(pca_df,aes(x=PC1,y=PC3,col = PlateId)) + 
  geom_point(size = 3, alpha = 0.5)+
  xlab(paste("PC1 (",varp[1],"%)")) +
  ylab(paste("PC3 (",varp[3],"%)")) +
  theme(panel.background = element_rect(fill = "white"))

dat.z = apply(dat,2,function(x) (x-meansd$mean)/meansd$sd)
# column_an = HeatmapAnnotation(MNSI = info_s$mnsi)
# Heatmap(dat.z, cluster_rows = T, cluster_columns = T, show_column_names = F,
#         show_row_names = F, name = "Z-score",row_title = "Proteins", top_annotation = column_an,
#         column_title = "Samples", column_title_rot = 0)

saveRDS(dat, file = "proteins_processed.rds")
