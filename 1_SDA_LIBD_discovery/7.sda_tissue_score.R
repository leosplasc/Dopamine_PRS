###############################################################################################################
######################## SDA tissue specificity ##################################
###############################################################################################################

my_scale = function(x){return(x/max(abs(x)))}

load(paste0("final_SDA.RData"))
tissue.score = as.data.frame(apply(final_SDA$B,2,my_scale))

save(tissue.score, file = paste0("tissue.score.RData"))



