#!/usr/bin/Rscript

###########################################################################################
###### model run script - take model code and config and run an instance of the model #####
###########################################################################################


source("core_model.R")
run_file<-commandArgs(trailingOnly=TRUE)[1]
experiment_name<-commandArgs(trailingOnly=TRUE)[2]
run_name<-unlist(strsplit(unlist(strsplit(run_file,"/"))[3],"\\."))[1]

print(run_file)
print(experiment_name)

#overwrite default config values with run file
source(run_file)

#run the model
model<-model_run()

# write output
write.csv(model,file=paste("experiments/",experiment_name,"/",run_name,".csv",sep=""))

