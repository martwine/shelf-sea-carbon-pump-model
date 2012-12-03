#############################################################
####### analyse data from lhc experiment ####################
#############################################################

experiment_name="physical_pump_lhc"
experiment_nruns=10000
experiment_path=paste("experiments/",experiment_name,"/",sep="")



trim_param<-function(param){
	param<-gsub(" ","",param)
	param<-gsub("<","",param)
	psplit<-unlist(strsplit(param,"-"))
	c<-list()
	c[[psplit[1]]]<-psplit[2]
	c
}

get_params<-function(run_config_file){
	params<-readLines(run_config_file)
	d<-as.list(unlist(lapply(params,trim_param)))
	d$run_number<-NULL
	d	
}

parse_data<-function(run_data_file){
	run_data<-read.csv(run_data_file)
	list(
		asflux=tail(cumsum(run_data$airseaFlux),n=1),
		delta_total_carbon=tail(run_data$total_C,n=1)-run_data$total_C[2]
	)
	
}

#data frame containing summary data for each run

run_number<-seq(to=experiment_nruns)
experiment_data<-as.data.frame(run_number)
experiment_data$run_config_file<-paste(experiment_path,experiment_name,"_",experiment_data$run_number,".R",sep="")
experiment_data$run_data_file<-paste(experiment_path,experiment_name,"_",experiment_data$run_number,".csv",sep="")
experiment_data<-cbind(experiment_data, t(sapply(experiment_data$run_config_file,get_params)),t(sapply(experiment_data$run_data_file,parse_data)))

cbind(experiment_data,get_params(experiment_data$run_config_file),parse_data(experiment_data$run_data_file))





