#! /usr/bin/RScript
library(lhs)
#############################################################
#### set-up multi-run experiments from a single config file #
#############################################################

max_runs<-10000
n_vals<-20 #numver of values for each variable

### READ IN CONFIG FILE AND MUNGE THE PARAMETER VALUES

rm(list=ls())

experiment_name<-commandArgs(trailingOnly=TRUE)

config_file=paste(experiment_name,".R",sep="")

source(paste("configs/",config_file,sep=""))

vars<-ls()
vars<-vars[-which(vars=="config_file")][-which(vars="experiment_name")]

#get a list of key (parameter) - value pairs
varlist<-list()
for(variable in vars){
	varlist[[variable]]=get(variable)
}

#count the length of the elements of varlist
varlength<-lapply(varlist,length)

#divide into fixed and varying values
fixedvars<-names(varlength[which(varlength==1)])
fixedvarlist<-lapply(fixedvars,get)
names(fixedvarlist)<-fixedvars

changingvars<-names(varlength[which(varlength>1)])
changingvarlist<-lapply(changingvars,get)
names(changingvarlist)<-changingvars

cvars<-new.env()
counter=1
for(variable in changingvarlist){
	varname<-changingvars[counter]
	if(variable[3]==1) #logarthmic distribution case
		{cvars[[varname]]<-exp(seq(from=log(variable[1]), to=log(variable[2]), length.out=n_vals))}
	else {cvars[[varname]]<-seq(from=variable[1],to=variable[2],length.out=n_vals)}
	counter<-counter+1
}

#construct the matrix of parameter values for each run (by lhc sampling if too may combinations)
cvarlist<-as.list(cvars)
if(nrow(expand.grid(cvarlist))<max_runs){
		run_matrix<-as.data.frame(expand.grid(cvarlist))
	} else {
		lh<-round((n_vals-1)*improvedLHS(max_runs,length(changingvars)))+1
		run_matrix=data.frame(run_number=seq(from=1,to=max_runs))
		counter=1

		for(vari in names(cvarlist)){
			run_matrix[[vari]]<-cvarlist[[vari]][lh[,counter]]
			counter<-counter+1
		}
		
}

#add in the non-varying parameters
for(vari in names(fixedvarlist))
	{
	run_matrix[[vari]]<-fixedvarlist[[vari]]
}

#make directory
dir.create(paste("experiments/",out_dir,sep=""))

#write config files
for(i in 1:nrow(run_matrix)) {
    row <- run_matrix[i,]
    # do stuff with row
	textout<-c()
	for(item in names(row)){
		slug<-paste(item,"<-",row[item],"\n")
		textout<-c(textout,slug)
		print(textout)
	}
	cat(	textout,
		file=paste("experiments/",out_dir,"/",row$RUN_NAME,"_",i,".R"sep=""),
		append=FALSE
		)	
}










