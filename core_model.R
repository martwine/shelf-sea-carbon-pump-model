#initalise model

source("eval_CO2_sys.R")
library("seacarb")
source("K_calcs_Johnson_OS.R")
source("default_config.R")



#########################################################################################
############## process and state variable evaluation functions ##########################
#########################################################################################

calc_as_flux<-function(pCO2,pCO2_atmos,temp,wind){
	deltapCO2<-KH_Molar_per_atmosphere("CO2",temp,S=35)*(pCO2_atmos-pCO2)*1e-6 #mol/l
	k_w<-Nightingkw("CO2",temp,wind,35)*3600*24*1.3 #1.3 scaling factor between short and long-term winds
	flux<-1e6*deltapCO2*1000*k_w #umol/m2/day
	flux
}

calc_DIC_uptake_from_NO3<-function(dNO3){
	-dNO3*redfield
}

calc_slDON_deg<-function(slDON){
	slDON*slDON_deg
}


calc_slDOC_deg<-function(slDOC){
	slDOC*slDOC_deg
}

calc_prod_slDON<-function(dNO3){
	-dNO3*nitrate_to_slDON_conv
}

calc_prod_slDOC<-function(dNO3){
	calc_prod_slDON(dNO3)*slDOM_C_TO_N
}

calc_PON_flux<-function(dNO3){
	-dNO3*(1-nitrate_to_slDON_conv)*SMLD/BMLD
}

calc_POC_flux<-function(dNO3){
	((-dNO3*redfield)-calc_prod_slDOC(dNO3))*SMLD/BMLD
}






#calculate  production during summer period
calc_TEPC_prod<-function(timestep){
	if(timestep > (BLOOM_START_DAY+BLOOM_DURATION) && timestep < mix_day){
		NH4_turnover*redfield*TEP_fraction	
	} else {
		0
	}	
}


calc_TEPC_deg<-function(TEPC){
	TEPC*TEPdeg	
}

eval_TEPC<-function(TEPC,timestep){
	TEPC+(calc_TEPC_prod(timestep)*SMLD/BMLD)-calc_TEPC_deg(TEPC)	
}

calc_PON_deg<-function(PON){
	PON*PONdeg
}

calc_POC_deg<-function(POC){
	POC*POCdeg
}

eval_slDON<-function(dNO3, slDON){
	slDON+calc_prod_slDON(dNO3)-calc_slDON_deg(slDON)
}

eval_slDOC<-function(dNO3, slDOC){
	slDOC+calc_prod_slDOC(dNO3)-calc_slDOC_deg(slDOC)
}

eval_DIC<-function(DIC,dNO3,pCO2,pCO2_atmos,temp,slDOC,POC,TEPC,depth,wind,timestep){
	remin_stuff<-ifelse(depth==SMLD,0,calc_POC_deg(POC)+calc_TEPC_deg(TEPC))
	DIC-calc_DIC_uptake_from_NO3(dNO3)+((calc_as_flux(pCO2,pCO2_atmos,temp,wind)/depth)/1000)+calc_slDOC_deg(slDOC)-calc_TEPC_prod(timestep)+remin_stuff
}

eval_PON<-function(PON,dNO3){
	PON+calc_PON_flux(dNO3)-calc_PON_deg(PON)
}

eval_POC<-function(POC,dNO3){
	POC+calc_POC_flux(dNO3)-calc_POC_deg(POC)
}

eval_BML_DIC<-function(BML_DIC, POC, TEPC){
	BML_DIC+calc_POC_deg(POC)+calc_TEPC_deg(TEPC)
}

eval_BML_NO3<-function(BML_NO3,PON){
	BML_NO3+calc_PON_deg(PON)
}

calc_mix<-function(sml_conc,bml_conc){
	(sml_conc*SMLD + bml_conc*BMLD)/COLUMN_DEPTH
}

eval_C_inventory<-function(depth, DIC, BML_DIC, slDOC,TEPC,POC){
	if(depth==SMLD && MODE==1){
		DIC*1000*SMLD
	} else if (depth==SMLD) {
		(DIC*1000*SMLD) + (BML_DIC*1000*BMLD) + (slDOC*1000*SMLD) + (TEPC*1000*BMLD) + (POC*1000*BMLD)
	} else {
		(DIC*1000 + slDOC*1000 + TEPC*1000 + POC*1000) *depth
	}
}

###################################################################################
##       Evaluate a timestep                                                      #
###################################################################################

eval_timestep<-function(timestep,current_state){
	
	# get the boundary conditions for the current timestep
	timestep_row<-box[timestep,]	
	pCO2_atmos<-timestep_row$pCO2_atmos
	dNO3<-timestep_row$dNO3
	temp<-timestep_row$temp
	wind<-timestep_row$wind

	# get the current state of the model at the end of the previous timestep
	DIC<-current_state$DIC
	slDOC<-current_state$slDOC
	slDON<-current_state$slDON
	pCO2<-current_state$pCO2
	
	# just get BML ones for the 2-box model
	if (MODE==2){
		TEPC<-current_state$TEPC
		POC<-current_state$POC
		PON<-current_state$PON
		BML_DIC<-current_state$BML_DIC
		BML_NO3<-current_state$BML_NO3
	}
	stepdata<-new.env()

	#define depths for 2-box model
	if(MODE==2){
		depth<-ifelse(timestep >= BLOOM_START_DAY && timestep < mix_day,SMLD,COLUMN_DEPTH)
	} else {depth=SMLD}

	if(MODE==2 && timestep==mix_day){
		#do the mixing
		DIC<-calc_mix(DIC,BML_DIC)
		NO3<-calc_mix(0,BML_NO3)
		slDON<-calc_mix(slDON,0)
		slDOC<-calc_mix(slDOC,0)
		pCO2<-carb(flag=15,init_TA*1e-6,DIC*1e-6)$pCO2[1]
		PON<-calc_mix(0,PON)
		POC<-calc_mix(0,POC)
		TEPC<-calc_mix(0,TEPC)
	}	
	
	stepdata$slDOC<-eval_slDOC(dNO3, slDOC)
	stepdata$slDON<-eval_slDON(dNO3, slDON)
	stepdata$airseaFlux<-calc_as_flux(pCO2,pCO2_atmos,temp,wind)
	stepdata$DIC<-eval_DIC(DIC,dNO3,pCO2,pCO2_atmos,temp,slDOC,POC,TEPC,depth,wind,timestep)
	stepdata$pCO2<-carb(flag=15,init_TA*1e-6,stepdata$DIC*1e-6)$pCO2[1]
	stepdata$deltapCO2<-pCO2_atmos-stepdata$pCO2
	if(MODE==2){
		stepdata$TEPC<-eval_TEPC(TEPC,timestep)
		stepdata$PON<-eval_PON(PON,dNO3)
		stepdata$POC<-eval_POC(POC,dNO3)
		stepdata$BML_DIC<-ifelse(depth==SMLD,eval_BML_DIC(BML_DIC,POC,TEPC),stepdata$DIC)
		stepdata$BML_NO3<-ifelse(depth==SMLD,eval_BML_NO3(BML_NO3,PON),ifelse(timestep==mix_day,NO3,BML_NO3))	
		stepdata$total_C<-eval_C_inventory(depth,stepdata$DIC,stepdata$BML_DIC,stepdata$slDOC,stepdata$TEPC,stepdata$POC)
	} else {
		stepdata$total_C<-eval_C_inventory(depth,stepdata$DIC,BML_DIC=0,stepdata$slDOC,stepdata$TEPC,stepdata$POC)
	}

	print(depth)
	print(paste("total_C_change",stepdata$total_C-current_state$total_C))	
	print(paste("air-sea flux",stepdata$airseaFlux))	
	print(stepdata$POC)
	as.data.frame(as.list(stepdata))
}


################################################################################
#     model control                                                            #
################################################################################

model_run<-function(){
	source("initialise_model.R")
	timestep = 1
	model_output<-data.frame(pCO2=init_pCO2,DIC=init_DIC,slDON=0,slDOC=0,airseaFlux=0,deltapCO2=delta_pCO2,TEPC=0,PON=0,POC=0)
	
	if(MODE==2){
		model_output$BML_DIC=init_DIC
		model_output$BML_NO3=WINTER_NITRATE
		
	}
	model_output$total_C=eval_C_inventory(depth=COLUMN_DEPTH,init_DIC,init_DIC,0,0,0)
	print(model_output)	
	while(timestep<run_length){
		current_state<-model_output[timestep,]
		model_output<-rbind(model_output,eval_timestep(timestep,current_state))
		timestep=timestep+1
	}
	model<-cbind(box,model_output)
	

	model	
}




