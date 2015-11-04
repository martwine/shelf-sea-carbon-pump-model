

source("eval_CO2_sys.R")
library("seacarb")
source("K_calcs_Johnson_OS.R")
source("default_config.R")



#########################################################################################
############## process and state variable evaluation functions ##########################
#########################################################################################

calc_as_flux<-function(pCO2,pCO2_atmos,temp,wind){
	deltapCO2<-KH_Molar_per_atmosphere("CO2",temp,S=35)*(pCO2_atmos-pCO2)*1e-6 #mol/l
	k_w<-Nightingkw("CO2",temp,wind,35)*3600*24*KW_SCALING_FACTOR # default 1.3 scaling factor between short and long-term winds
	flux<-1e6*deltapCO2*1000*k_w #umol/m2/day
	flux
}

calc_DIC_uptake_from_NO3<-function(dNO3){
	-dNO3*redfield
}

calc_slDON_deg<-function(slDON,temp){
	slDON*slDON_deg*Q10_rate_scale(temp)
}


calc_slDOC_deg<-function(slDOC,temp){
	slDOC*slDOC_deg*Q10_rate_scale(temp)
}

calc_prod_slDON<-function(dNO3){
	-dNO3*nitrate_to_slDON_conv
}

calc_prod_slDOC<-function(dNO3){
	calc_prod_slDON(dNO3)*slDOM_C_TO_N
}

calc_PON_flux<-function(dNO3,overconsumption){
	conc_change=(-dNO3*(1-nitrate_to_slDON_conv)+(overconsumption/redfield))
	amount=conc_change*1000*SMLD
	amount
}

calc_POC_flux<-function(dNO3,overconsumption){
	conc_change<-(((-dNO3*redfield)-calc_prod_slDOC(dNO3))+overconsumption)
	amount<-conc_change*1000*SMLD
	amount
}


calc_remin_overconsumption<-function(PON,slDON,POC,slDOC,temp,bottomtemp,timestep){
	#organic N that gets remineralised above redfield makes new POC at redfield. Called only in summer
	#if remineralisation is C rich, return zero (a simplification)
	if(timestep > (SPRING_START_DAY+SPRING_DURATION) && timestep < mix_day && timestep > 1e6){
		x<-(calc_slDON_deg(slDON,temp)*redfield)-calc_slDOC_deg(slDOC,temp)
		slDON_remin_C_fixation<-ifelse(x>0,x,0)
		
		if (MODE==2){
			conc_change<-slDON_remin_C_fixation

			
		} else {
			y<-(calc_PON_deg(PON,bottomtemp)*redfield)-calc_POC_deg(POC,bottomtemp)
			PON_remin_C_fixation<-ifelse(y>0,y,0)
			conc_change<-slDON_remin_C_fixation + PON_remin_C_fixation

		}
		amount=conc_change*1000*SMLD
		amount
	} else {0}
}

Q10_rate_scale<-function(temp){
	Q10^((temp-Q10_REF_TEMP)/10)
}


#calculate TEP production during summer period
calc_TEPC_prod<-function(timestep){
	if(timestep > (SPRING_START_DAY+SPRING_DURATION) && timestep < mix_day){
		NH4_turnover*redfield*TEP_fraction	
	} else {
		0
	}	
}


calc_TEPC_deg<-function(TEPC,temp){
	TEPC*TEPdeg*Q10_rate_scale(temp)	
}

eval_TEPC<-function(TEPC,temp,timestep,..){
	TEPC+(calc_TEPC_prod(timestep)*ifelse(MODE==2,SMLD/BMLD,1))-calc_TEPC_deg(TEPC,temp)	
}

calc_PON_deg<-function(PON, temp){
	deg<-PONdeg*Q10_rate_scale(temp)
	ifelse(PON-deg>0,deg,PON)
}

calc_POC_deg<-function(POC, temp){
	deg<-POCdeg*Q10_rate_scale(temp)
    ifelse(POC-deg>0,deg,POC)
}

eval_slDON<-function(dNO3, slDON,temp){
	slDON+calc_prod_slDON(dNO3)-calc_slDON_deg(slDON,temp)
}

eval_slDOC<-function(dNO3, slDOC, temp){
	slDOC+calc_prod_slDOC(dNO3)-calc_slDOC_deg(slDOC,temp)
}

eval_DIC<-function(DIC,dNO3,pCO2,pCO2_atmos,temp,bottomtemp,slDOC,POC,TEPC,depth,wind,timestep,overconsumption,resusp_DIC,...){
	#remin to DIC in SMLD in 1-box mode, in 2 box mode this happens at depth
	remin_stuff<-ifelse(depth==SMLD&&MODE==2,0,((calc_POC_deg(POC,bottomtemp))/(1000*depth))+calc_TEPC_deg(TEPC,bottomtemp)+(resusp_DIC/(1000*depth)))
   print(calc_POC_deg(POC,bottomtemp))
    print(depth)
    print(remin_stuff   )
	DIC-calc_DIC_uptake_from_NO3(dNO3)+((calc_as_flux(pCO2,pCO2_atmos,temp,wind)/depth)/1000)+calc_slDOC_deg(slDOC,temp)-calc_TEPC_prod(timestep)+remin_stuff-overconsumption/(1000*SMLD)
}


eval_PON<-function(PON,dNO3,overconsumption,bottomtemp){
	PON+calc_PON_flux(dNO3,overconsumption)-calc_PON_deg(PON,bottomtemp)
}

eval_POC<-function(POC,dNO3,overconsumption,bottomtemp){
	POC+calc_POC_flux(dNO3,overconsumption)-calc_POC_deg(POC,bottomtemp)
}

eval_BML_DIC<-function(BML_DIC, POC, TEPC, bottomtemp, resusp_DIC){
  BML_DIC+(calc_POC_deg(POC,bottomtemp)/(1000*BMLD))+calc_TEPC_deg(TEPC,bottomtemp)+resusp_DIC/(1000*BMLD)
}

eval_BML_NO3<-function(BML_NO3,PON,bottomtemp){
	BML_NO3+calc_PON_deg(PON,bottomtemp)/(1000*BMLD)
}

calc_mix<-function(sml_conc,bml_conc){
	(sml_conc*SMLD + bml_conc*BMLD)/COLUMN_DEPTH
}

eval_C_inventory<-function(depth, DIC, BML_DIC, slDOC,TEPC,POC){
	if(depth==SMLD&&MODE==2) {
		(DIC*1000*SMLD) + (BML_DIC*1000*BMLD) + (slDOC*1000*SMLD) + (TEPC*1000*BMLD) + (POC)
	} else {
		(DIC*1000 + slDOC*1000 + TEPC*1000)*depth + POC
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
	bottomtemp<-timestep_row$bottomtemp
	wind<-timestep_row$wind
	jday<-timestep_row$jday
	#include constant DIC flux from sediment (default value 0)
	if(jday> (SPRING_START_DAY+SPRING_DURATION) && jday < mix_day){
	  resusp_DIC=resusp_DIC_SUMMER
	}else{
	  resusp_DIC=resusp_DIC_WINTER
	}
	
	#add in DIC release from trawl event on TRAWL_DAY
	if(jday==TRAWL_DAY){resusp_DIC=resusp_DIC+TRAWL_DIC_RELEASE}

	# get the current state of the model at the end of the previous timestep
	DIC<-current_state$DIC
	slDOC<-current_state$slDOC
	slDON<-current_state$slDON
	pCO2<-current_state$pCO2
	TEPC<-current_state$TEPC
	POC<-current_state$POC
	PON<-current_state$PON	
	Benthic_POC<-current_state$Benthic_POC
	Benthic_PON<-current_state$Benthic_PON

	# just get BML ones for the 2-box model
	if (MODE==2){
		BML_DIC<-current_state$BML_DIC
		BML_NO3<-current_state$BML_NO3
	}
	stepdata<-new.env()

	#define depths for 2-box model
	if(MODE==2){
		depth<-ifelse(jday >= SPRING_START_DAY && jday < mix_day,SMLD,COLUMN_DEPTH)
	} else {depth=SMLD}

	if(MODE==2 && jday==mix_day){
		#do the mixing
		DIC<-calc_mix(DIC,BML_DIC)
		NO3<-calc_mix(0,BML_NO3)
		slDON<-calc_mix(slDON,0)
		slDOC<-calc_mix(slDOC,0)
		pCO2<-carb(flag=15,init_TA*1e-6,DIC*1e-6,T=temp)$pCO2[1]
		TEPC<-calc_mix(0,TEPC)
	}
	
	if(jday==SPRING_START_DAY){
		Benthic_POC<-Benthic_POC + POC
		POC<-0
		Benthic_PON<-Benthic_PON + PON
		PON<-0
	} 
		
	stepdata$remin_overconsumption<-calc_remin_overconsumption(PON,slDON,POC,slDOC,temp,bottomtemp,timestep=jday)
	stepdata$slDOC<-eval_slDOC(dNO3, slDOC, temp)
	stepdata$slDON<-eval_slDON(dNO3, slDON, temp)
	stepdata$airseaFlux<-calc_as_flux(pCO2,pCO2_atmos,temp,wind)
	stepdata$DIC<-eval_DIC(DIC,dNO3,pCO2,pCO2_atmos,temp,bottomtemp,slDOC,POC,TEPC,depth,wind,timestep=jday,stepdata$remin_overconsumption,resusp_DIC)
	stepdata$pCO2<-carb(flag=15,init_TA*1e-6,stepdata$DIC*1e-6,T=temp)$pCO2[1]
	stepdata$deltapCO2<-pCO2_atmos-stepdata$pCO2
	stepdata$TEPC<-eval_TEPC(TEPC,bottomtemp,timestep=jday)
	stepdata$PON<-eval_PON(PON,dNO3,stepdata$remin_overconsumption,bottomtemp)
	stepdata$POC<-eval_POC(POC,dNO3,stepdata$remin_overconsumption,bottomtemp)
	stepdata$Benthic_POC<-Benthic_POC
	stepdata$Benthic_PON<-Benthic_PON
	
	if(MODE==2){

		stepdata$BML_DIC<-ifelse(depth==SMLD,eval_BML_DIC(BML_DIC,POC,TEPC,bottomtemp,resusp_DIC),stepdata$DIC)
		stepdata$BML_NO3<-ifelse(depth==SMLD,eval_BML_NO3(BML_NO3,PON,bottomtemp),ifelse(timestep==mix_day,NO3,BML_NO3))	
		stepdata$total_C<-eval_C_inventory(depth,stepdata$DIC,stepdata$BML_DIC,stepdata$slDOC,stepdata$TEPC,stepdata$POC)
	} else {
		stepdata$total_C<-eval_C_inventory(depth,stepdata$DIC,BML_DIC=0,stepdata$slDOC,stepdata$TEPC,stepdata$POC)
	}
	
	#print(paste("total_C_change",stepdata$total_C-current_state$total_C))
	#print(paste("air-sea flux",stepdata$airseaFlux))
	#print(paste("TIMESTEP",timestep))
	#print(paste("POC",stepdata$POC))
	#print(paste("PON",stepdata$PON))
    #print(depth)
	as.data.frame(as.list(stepdata))
}


################################################################################
#     model control                                                            #
################################################################################

model_run<-function(){
	source("initialise_model.R")
	timestep = 1
	model_output<-data.frame(pCO2=init_pCO2,DIC=init_DIC,slDON=0,slDOC=0,airseaFlux=0,deltapCO2=delta_pCO2,TEPC=0,PON=0,POC=0,remin_overconsumption=0,Benthic_POC=0,Benthic_PON=0)
	
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




