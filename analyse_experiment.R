#############################################################
####### analyse data from lhc experiment ####################
#############################################################

experiment_name="Physics_lhc"
experiment_nruns=10
experiment_path=paste("experiments/",experiment_name,sep="")



trim_param<-function(param){
	param<-gsub(" ","",param)
	param<-gsub("<","",param)
	psplit<-unlist(strsplit(param,"-"))
	c<-list()
	c[[psplit[1]]]=psplit[2]
	c
}

get_params<-function(run_config_file){
	params<-readLines(run_config_file)
	unlist(lapply(params,trim_param))	
}

parse_data<-function(run_data_file){
	run_data<-read.csv(run_data_file)
	asflux<-tail(cumsum(run_data$airseaFlux),n=1)
	delta_total_carbon<-tail(run_data$total_C,n=1)-run_data$total_C[2]
}

#data frame containing summary data for each run


