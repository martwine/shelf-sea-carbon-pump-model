#!/usr/bin/Rscript

###########################################################################################
###### model run script - take model code and config and run an instance of the model #####
###########################################################################################


source("core_model.R")
run_name<-commandArgs(trailingOnly=TRUE)

#overwrite default config values with run file
config_file=paste(run_name,".R",sep="")


source(paste("configs/",config_file,sep=""))


#run the model
model<-model_run()

print(cumsum(model$airseaFlux))
print(warnings())

#create a run directory
dir.create(paste("results/",run_name,sep=""))

# write output
write.csv(model,file=paste(out_dir,"/",run_name,"/",run_name,".csv",sep=""))

#move config
file.copy(paste("configs/",config_file,sep=""),paste(out_dir,"/",run_name,"/",config_file,sep=""),overwrite=TRUE)
