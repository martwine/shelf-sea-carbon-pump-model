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


parse_data<-function(trawl_data_file,notrawl_data_file){
	print(trawl_data_file)
	trawl_data<-read.csv(trawl_data_file)
	print(notrawl_data_file)
	notrawl_data<-read.csv(notrawl_data_file)
	print("finished")

	#phystotalc=phys_data$total_C,
	#redtotalc=red_data$total_C,
	#nonredtotalc=nonred_data$totalC,
	
	trawltotalc<-tail(trawl_data$total_C,365)
	notrawltotalc<-tail(notrawl_data$total_C,365)
		
		delta_total_carbon_trawl=tail(trawltotalc,n=1)-trawltotalc[1]
		delta_total_carbon_notrawl=tail(notrawltotalc,n=1)-notrawltotalc[1]


		mean_storage_trawl=mean(trawltotalc)
		mean_storage_notrawl=mean(notrawltotalc)

	delta_cumulativeasflux<-tail(cumsum(notrawl_data$airseaFlux),1)-tail(cumsum(trawl_data$airseaFlux),1)
print(tail(cumsum(notrawl_data$airseaFlux),1))
print("tha was the number")
	t(as.data.frame(c(delta_total_carbon_trawl,delta_total_carbon_notrawl,mean_storage_trawl,mean_storage_notrawl,delta_cumulativeasflux),row.names=c("delta_total_carbon_trawl","delta_total_carbon_notrawl","mean_storage_trawl","mean_storage_notrawl","delta_cumulativeasflux"),colnames=NULL))
	
}

#data frame containing summary data for each run

trawl_number<-seq(to=experiment_nruns)
print(trawl_number)
notrawl_number<-seq(from=experiment_nruns+1,to=2*experiment_nruns,by=1)

experiment_data<-as.data.frame(trawl_number)

experiment_data$notrawl_config_file<-paste(experiment_path,experiment_name,"_",notrawl_number,".R",sep="")
experiment_data$trawl_config_file<-paste(experiment_path,experiment_name,"_",trawl_number,".R",sep="")

experiment_data$notrawl_data_file<-paste(experiment_path,experiment_name,"_",notrawl_number,".csv",sep="")
experiment_data$trawl_data_file<-paste(experiment_path,experiment_name,"_",trawl_number,".csv",sep="")

param_vals<-t(sapply(experiment_data$trawl_config_file,get_params))
colnames(param_vals)<-get_paramnames(experiment_data$trawl_config_file[1])

print("still working")
carbon_uptakevals<-t(mapply(parse_data,experiment_data$notrawl_data_file,experiment_data$trawl_data_file))
colnames(carbon_uptakevals)<-c("delta_total_carbon_trawl","delta_total_carbon_notrawl","mean_storage_trawl","mean_storage_notrawl","delta_cumulativeasflux")

experiment_data<-cbind(experiment_data, param_vals,carbon_uptakevals)

save(experiment_data,file=paste(experiment_path,"analysis",experiment_name,"_summary.Rdata",sep=""))


