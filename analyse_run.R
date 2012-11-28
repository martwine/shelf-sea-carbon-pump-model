##!/usr/bin/Rscript


############################################################
# analyse a single run                                     #
############################################################



#run_name<-commandArgs(trailingOnly=TRUE)
config_file=paste(run_name,".R",sep="")
source(paste("configs/",config_file,sep=""))
run_data<-read.csv(file=paste(out_dir,"/",run_name,"/",run_name,".csv",sep=""))


#check internal consistency

if(sum(run_data$airseaFlux)-sum(diff(run_data$total_C))<1e-6){print("carbon checks out")}else{print(paste("Something might be wrong - the discrepancy between as-flux is:",sum(run_data$airseaFlux)-sum(diff(run_data$total_C))))}

plot(run_data$temperature,run_data$airseaFlux)
