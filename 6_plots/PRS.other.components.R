library(ggplot2)
library(metafor)
library(patchwork)

##### load NC or SCZ data
load("sda_quantile/Heatmap_data/LM_detect.final.RData")
PRS_assoc = PRS_assoc[c("V32","V78","V87","V102","V104","V115")]
PRS_assoc = lapply(PRS_assoc,function(x) x$PRS6.PGC3$model)
stopifnot(all(sapply(PRS_assoc,function(x) identical(row.names(PRS_assoc$V32),row.names(x)))))
df = lapply(PRS_assoc,function(x) x[,c(1,6)])
df = Reduce(merge,df)
df$PRS6.PGC3 = scale(df$PRS6.PGC3)
df.melt = reshape2::melt(df[,-1])
df.melt$PRS = rep(df$PRS6.PGC3,6)
PRS_p = PRS_p[unique(levels(df.melt$variable)),"PRS6.PGC3"]
PRS_p = signif(PRS_p,3)
df.melt$variable = as.character(df.melt$variable)
df.melt$variable = gsub("V","C",df.melt$variable)
df.melt$comp = paste0(df.melt$variable,"\n","p=",rep(PRS_p,each = nrow(df)))
df.melt$comp = gsub("0.0",".0",df.melt$comp)
df.melt$comp = factor(df.melt$comp, levels = unique(df.melt$comp))

my_PRS <- ggplot(df.melt, aes(x = PRS, y = value)) +
  geom_point(size = 4, 
             color = "black",alpha = 0.6) +
  stat_smooth(method = "lm",
              se = T,
              linewidth = 1,
              color = "black",
              fullrange = T) +
  ylim(-2,2) +
  xlim(-3,3) +
  xlab("Whole genome PRS") +
  ylab("Individual score") +
  facet_wrap(.~comp) +
  theme_classic(base_size = 30) 

my_PRS

ggsave(filename = "sda_quantile/Paper/PRS.other.components.svg",plot = my_PRS, width=7.5, height=7,dpi = 600)


##### load metanalaysis of NC and SCZ 
load("sda_quantile/Paper/other.components_meta.RData")
# pdf(file = paste0("sda_quantile/Paper/other.components_meta.PRS1.whole.pdf"),width=20, height=10)
# par(mfrow = c(3,3))
purrr::iwalk(meta,function(x,y){
  pdf(file = paste0("sda_quantile/Paper/",y,"_meta.PRS1.whole.pdf"),width=10, height=3.9)
  forest(x,
         header = c("Diagnosis"),annotate = T,slab = c("NC","SCZ"), 
         cex = 2, cex.lab = 1.8,cex.axis = 2.2, level = 99.5,
        xlim=c(-2.5,3.5))
  dev.off()
})





