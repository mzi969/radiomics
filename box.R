rm(list=ls())
library(ggpubr)
setwd("D:/Desktop/1")
getwd() 
box<-read.csv("primary_box.csv")
group <- box$group
radsocre <-box$radscore
tiff(file = "primary_box.tiff", res = 600, width = 4800, height = 3600, compression = "lzw")
ggboxplot(box, x = "group", y = "radsocre",
          color = "group",add="jitter", palette = "jama"
                  )+
  stat_compare_means()+theme(text=element_text(size=16,  family="serif"))
dev.off()
###########2##########
box<-read.csv("validation_box.csv")
group <- box$group
radsocre <-box$radscore
tiff(file = "validation_box.tiff", res = 600, width = 4800, height = 3600, compression = "lzw")
ggboxplot(box, x = "group", y = "radsocre",
          color = "group",add="jitter", palette = "jama"
                  )+
  stat_compare_means()+theme(text=element_text(size=16,  family="serif"))
dev.off()
