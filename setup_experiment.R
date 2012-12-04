#! /usr/bin/Rscript
library(lhs)
#############################################################
#### set-up multi-run experiments from a single config file #
#############################################################

#usage Rscript seup_experiment.R <name of experiment> <force_lhc?>



### READ IN CONFIG FILE AND MUNGE THE PARAMETER VALUES

rm(list=ls())

experiment_name<-commandArgs(trailingOnly=TRUE)[1]

config_file=paste(experiment_name,".R",sep="")

source(paste("configs/",config_file,sep=""))

vars<-ls()
vars<-vars[-which(vars=="config_file")]
vars<-vars[-which(vars=="experiment_name")]
vars<-vars[-which(vars=="RUN_NAME")]
vars<-vars[-which(vars=="out_dir")]

#reload to get out_dir back

experiment_name<-commandArgs(trailingOnly=TRUE)[1]
max_runs<-as.integer(commandArgs(trailingOnly=TRUE)[2])
chunk_size<-as.integer(commandArgs(trailingOnly=TRUE)[3])
print(chunk_size)
source(paste("configs/",config_file,sep=""))



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

#how many runs and variable values?
# if 3 values for variable in config file, this specifies min, max and dsitrbution (log vs linear)
# if any other number of variables, each value is used in the lhc (needs to match n_vals)
n_vals<-10 #number of values for each variable where specific values not set


cvars<-new.env()
counter=1
for(variable in changingvarlist){
	varname<-changingvars[counter]
	if(length(variable)==3) # min, max and dist_type specified
		{if(variable[3]==1) #logarthmic distribution case
			{cvars[[varname]]<-exp(seq(from=log(variable[1]), to=log(variable[2]), length.out=n_vals))}
		else {cvars[[varname]]<-seq(from=variable[1],to=variable[2],length.out=n_vals)}
	} else #take each values specified
		{cvars[[varname]]<-variable
	}
	counter<-counter+1
}

#construct the matrix of parameter values for each run (by lhc sampling if too may combinations)
cvarlist<-as.list(cvars)

lh<-round((n_vals-1)*randomLHS(max_runs,length(changingvars)))+1
run_matrix=data.frame(run_number=seq(from=1,to=max_runs))
counter=1

for(vari in names(cvarlist)){
	run_matrix[[vari]]<-cvarlist[[vari]][lh[,counter]]
	counter<-counter+1
}
		


#add in the non-varying parameters
for(vari in names(fixedvarlist))
	{
	run_matrix[[vari]]<-fixedvarlist[[vari]]
}
print(run_matrix)
#make directory
dir.create(paste("experiments/",out_dir,sep=""))

#write config files
for(i in 1:nrow(run_matrix)) {
    row <- run_matrix[i,]
    # do stuff with row
	textout<-c()
	for(item in names(row)){
		if(is.character(row[item])){
			slug<-paste(item,'<-\"',row[item],'\"\n')
		} else {
			slug<-paste(item,"<-",row[item],"\n")
		}		
		textout<-c(textout,slug)
		print(textout)
	}
	cat(	textout,
		file=paste("experiments/",out_dir,"/",experiment_name,"_",i,".R",sep=""),
		append=FALSE
		)	
}

#write run_scripts

#split data frame into sets of chunk size runs
n<-seq(from=1,to=nrow(run_matrix)/chunk_size)

print(n)
filecounter=1
for(chunk in n){
	textout<-c(	"#!/bin/bash\n",
			"#BSUB -q short\n",
			paste("#BSUB -J",filecounter,"\n"),
			paste("#BSUB -oo",experiment_name,"-%J.out\n", sep=""),
			paste("#BSUB -eo",experiment_name,"-%J.err\n", sep=""),
			". /etc/profile\n",
			"module add R\n")
	chunklist<-seq(from=1,to=chunk_size)+((filecounter-1)*chunk_size)	
	for(run in chunklist){
		textout<-c(textout,paste("Rscript experiment_run.R experiments/",out_dir,"/",experiment_name,"_",run,".R ",experiment_name,"\n",sep="") )
		counter=counter+1
	}
	cat(textout,file=paste(experiment_name,"_set_",filecounter,".sh",sep=""),append=FALSE)

	filecounter=filecounter+1
}





