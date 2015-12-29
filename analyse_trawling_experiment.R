#!/usr/bin/Rscript
#############################################################
####### analyse data from lhc experiment ####################
#############################################################

experiment_name=commandArgs(trailingOnly=TRUE)[1]
experiment_nruns=as.numeric(commandArgs(trailingOnly=TRUE)[2])
experiment_path=paste("experiments/",experiment_name,"/",sep="")

print(experiment_nruns)

trim_param<-function(param){
	param<-gsub("<-","",param)
	param<-gsub("  "," ",param)
	param<-gsub("^\\s+|\\s+$","",param)
	
	unlist(strsplit(param," "))
	
	#c<-list()
	#c[[psplit[1]]]<-psplit[2]
	#c
}

get_params<-function(run_config_file){
	params<-readLines(run_config_file)
	x<-sapply(params,trim_param)
	y<- t(as.data.frame(as.numeric(x[2,])))
	colnames(y)<-x[1,]
	rownames(y)<-NULL
	y	
	#d<-as.list(unlist(lapply(params,trim_param)))
	#print(d$run_number)
	#d$run_number<-NULL
	#as.data.frame(as.matrix((unlist(d))))	
}

get_paramnames<-function(run_config_file){
	params<-readLines(run_config_file)
	sapply(params,trim_param)[1,]
}


parse_data<-function(trawlclim_data_file,notrawl_data_file,noclim_data_file,baseline_data_file){
	print(trawlclim_data_file)
	trawlclim_data<-read.csv(trawlclim_data_file)
	print(notrawl_data_file)
	notrawl_data<-read.csv(notrawl_data_file)
	print(noclim_data_file)
	noclim_data<-read.csv(noclim_data_file)
	print(baseline_data_file)
	baseline_data<-read.csv(baseline_data_file)
	print("finished")

	#phystotalc=phys_data$total_C,
	#redtotalc=red_data$total_C,
	#nonredtotalc=nonred_data$totalC,
	
	trawlclimtotalc<-tail(trawlclim_data$total_C,365)
	notrawltotalc<-tail(notrawl_data$total_C,365)
	noclimtotalc<-tail(noclim_data$total_C,365)
	baselinetotalc<-tail(baseline_data$total_C,365)
	
	delta_total_carbon_trawlclim=tail(trawlclimtotalc,n=1)-trawlclimtotalc[1]
	delta_total_carbon_notrawl=tail(notrawltotalc,n=1)-notrawltotalc[1]
	delta_total_carbon_noclim=tail(noclimtotalc,n=1)-noclimtotalc[1]
	delta_total_carbon_baseline=tail(baselinetotalc,n=1)-baselinetotalc[1]


	mean_storage_trawlclim=mean(trawlclimtotalc)
	mean_storage_notrawl=mean(notrawltotalc)
	mean_storage_noclim=mean(noclimtotalc)
	mean_storage_baseline=mean(baselinetotalc)
	

	trawlclim_delta_cumulativeasflux<-tail(cumsum(tail(baseline_data$airseaFlux,365)),1)-tail(cumsum(tail(trawlclim_data$airseaFlux,356)),1)
	trawlonly_delta_cumulativeasflux<-tail(cumsum(tail(baseline_data$airseaFlux,365)),1)-tail(cumsum(tail(noclim_data$airseaFlux,365)),1)
	climonly_delta_cumulativeasflux<-tail(cumsum(tail(baseline_data$airseaFlux,365)),1)-tail(cumsum(tail(notrawl_data$airseaFlux,365)),1)
	
	t(as.data.frame(c(delta_total_carbon_trawlclim,
	                  delta_total_carbon_notrawl,
	                  delta_total_carbon_noclim,
	                  delta_total_carbon_baseline,
	                  mean_storage_trawlclim,
	                  mean_storage_notrawl,
	                  mean_storage_noclim,
	                  mean_storage_baseline,
	                  trawlclim_delta_cumulativeasflux,
	                  climonly_delta_cumulativeasflux,
	                  trawlonly_delta_cumulativeasflux),
	                row.names=c("delta_total_carbon_trawlclim","delta_total_carbon_notrawl","delta_total_carbon_noclim","delta_total_carbon_baseline","mean_storage_trawlclim","mean_storage_notrawl","mean_storage_noclim","mean_storage_baseline","trawlclim_delta_cumulativeasflux","climonly_delta_cumulativeasflux","trawlonly_delta_cumulativeasflux"),colnames=NULL))
	
}

#data frame containing summary data for each run

trawlclim_number<-seq(to=experiment_nruns)
print(trawlclim_number)
notrawl_number<-seq(from=experiment_nruns+1,to=2*experiment_nruns,by=1)
noclim_number<-seq(from=2*experiment_nruns+1,to=3*experiment_nruns,by=1)
baseline_number<-seq(from=3*experiment_nruns+1,to=4*experiment_nruns,by=1)

experiment_data<-as.data.frame(trawlclim_number)

experiment_data$baseline_config_file<-paste(experiment_path,experiment_name,"_",baseline_number,".R",sep="")
experiment_data$noclim_config_file<-paste(experiment_path,experiment_name,"_",noclim_number,".R",sep="")
experiment_data$notrawl_config_file<-paste(experiment_path,experiment_name,"_",notrawl_number,".R",sep="")
experiment_data$trawlclim_config_file<-paste(experiment_path,experiment_name,"_",trawlclim_number,".R",sep="")

experiment_data$baseline_data_file<-paste(experiment_path,experiment_name,"_",baseline_number,".csv",sep="")
experiment_data$noclim_data_file<-paste(experiment_path,experiment_name,"_",noclim_number,".csv",sep="")
experiment_data$notrawl_data_file<-paste(experiment_path,experiment_name,"_",notrawl_number,".csv",sep="")
experiment_data$trawlclim_data_file<-paste(experiment_path,experiment_name,"_",trawlclim_number,".csv",sep="")

param_vals<-t(sapply(experiment_data$trawlclim_config_file,get_params))
colnames(param_vals)<-get_paramnames(experiment_data$trawlclim_config_file[1])

print("still working")
carbon_uptakevals<-t(mapply(parse_data,experiment_data$trawlclim_data_file,experiment_data$notrawl_data_file,experiment_data$noclim_data_file,experiment_data$baseline_data_file))
colnames(carbon_uptakevals)<-c("delta_total_carbon_trawlclim","delta_total_carbon_notrawl","delta_total_carbon_noclim","delta_total_carbon_baseline","mean_storage_trawlclim","mean_storage_notrawl","mean_storage_noclim","mean_storage_baseline","trawlclim_delta_cumulativeasflux","climonly_delta_cumulativeasflux","trawlonly_delta_cumulativeasflux")

experiment_data<-cbind(experiment_data, param_vals,carbon_uptakevals)

save(experiment_data,file=paste(experiment_path,"analysis",experiment_name,"_summary.Rdata",sep=""))


