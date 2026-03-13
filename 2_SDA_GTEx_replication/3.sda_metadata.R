###############################################################################################################
######################## SDA metadata  ####################################
###############################################################################################################

load("sda_preprocessing.RData")
rm(list = setdiff(ls(),"sample_list"))
source("3a.sda_metadata_functions.R")

data = metadata_GTEx(data_dir = "",
              tissues = c("caudate","dlpfc","hippo"),common_metadata = c("SAMPID","Sex","Age","PMI"),
              tissue_specific_metadata = c("RIN","TotalAssignedGene","rRNA_rate"),
              data_name = "final_SDA",output_dir = "",
              sample_list = sample_list)

save(data,file = paste0("data.SDA.RData"))




