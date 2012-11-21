#initalise model

source("eval_CO2_sys.R")
library("seacarb")
source("K_calcs_Johnson_OS.R")


#calc bottom mixed layer depth
BMLD<-COLUMN_DEPTH-SMLD

################################################################################################
###############    Boundary conditions  ########################################################


#Temperature cycle calculated from max and min values (1st March min, 1st September max)

surface_temp_cycle<-MIN_SURFACE_TEMP+(0.5*(MAX_SURFACE_TEMP-MIN_SURFACE_TEMP))*(1+sin((seq(from=-90,to=270,length.out=365)*pi)/180))

surface_temp<-c(surface_temp_cycle[(365-OFFSET):365],surface_temp_cycle[1:(365-OFFSET-1)])



#define run_length as from 1st Jan 1 year, to immediately before spring bloom the next
run_length=365+BLOOM_START_DAY


#calc jday of mixing event
mix_day<-BLOOM_START_DAY+BLOOM_DURATION+SUMMER_LENGTH

#define surface nitrate 'envelope'
nday<-c(1,BLOOM_START_DAY,BLOOM_START_DAY+BLOOM_DURATION,run_length)
nitrate<-c(WINTER_NITRATE,WINTER_NITRATE,0,0)
nitrate<-approx(nday,nitrate,n=run_length)


#get a simple seasonal cycle of atmospheric CO2
min_pCO2<-AVERAGE_pCO2-(AMPLITUDE/2)
pCO2_atmos<-(min_pCO2+(AMPLITUDE/2)*(1+cos((seq(from=0,to=360,length.out=365)*pi)/180)))

init_pCO2<-min_pCO2+AMPLITUDE


#initalise DIC from pCO2 and estimated TA 
init_DIC<-1e6*carb(flag=24,(init_pCO2-delta_pCO2),init_TA*1e-6)$DIC[1]
print(paste("init DIC:", init_DIC))

#wind speed cycle
wind_speed<-MIN_WIND+(0.5*(MAX_WIND-MIN_WIND))*(1+sin((seq(from=90,to=450,length.out=365)*pi)/180))



#initalise_data
box<-as.data.frame(nitrate)
colnames(box)<-c("day","nitrate")
box$pCO2_atmos<-c(pCO2_atmos,pCO2_atmos)[1:run_length]
box$temp<-c(surface_temp,surface_temp)[1:run_length]
box$wind<-c(wind_speed,wind_speed)[1:run_length]

#delta nitrate column - gives change between each timestep and previous
box$dNO3<-c(0,diff(box$nitrate))



#########################################################################################
############## process and state variable evaluation functions ##########################
#########################################################################################

calc_as_flux<-function(pCO2,pCO2_atmos,temp,wind){
	deltapCO2<-KH_Molar_per_atmosphere("CO2",temp,S=35)*(pCO2_atmos-pCO2)*1e-6 #mol/l
	k_w<-Nightingkw("CO2",temp,wind,35)*3600*24
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


#instantly remineralise POM in the bottom layer
calc_NO3_BML<-function(dNO3){
	-dNO3*(1-nitrate_to_slDON_conv)*SMLD/BMLD
}

calc_DIC_BML<-function(dNO3){
	-dNO3*(1-nitrate_to_slDON_conv)*redfield*SMLD/BMLD
}


eval_slDON<-function(dNO3, slDON){
	slDON+calc_prod_slDON(dNO3)-calc_slDON_deg(slDON)
}

eval_slDOC<-function(dNO3, slDOC){
	slDOC+calc_prod_slDOC(dNO3)-calc_slDOC_deg(slDOC)
}

eval_DIC<-function(DIC,dNO3,pCO2,pCO2_atmos,temp,slDOC,depth,wind){
	DIC-calc_DIC_uptake_from_NO3(dNO3)+((calc_as_flux(pCO2,pCO2_atmos,temp,wind)/depth)/1000)+calc_slDOC_deg(slDOC)
}

eval_BML_DIC<-function(dNO3, BML_DIC){
	BML_DIC+calc_DIC_BML(dNO3)
}

eval_BML_NO3<-function(dNO3, BML_NO3){
	BML_NO3+calc_NO3_BML(dNO3)
}

calc_mix<-function(sml_conc,bml_conc){
	(sml_conc*SMLD + bml_conc*BMLD)/COLUMN_DEPTH
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
	}	
	
	stepdata$slDOC<-eval_slDOC(dNO3, slDOC)
	stepdata$slDON<-eval_slDON(dNO3, slDON)
	stepdata$airseaFlux<-calc_as_flux(pCO2,pCO2_atmos,temp,wind)
	stepdata$DIC<-eval_DIC(DIC,dNO3,pCO2,pCO2_atmos,temp,slDOC,depth,wind)
	stepdata$pCO2<-carb(flag=15,init_TA*1e-6,stepdata$DIC*1e-6)$pCO2[1]
	stepdata$deltapCO2<-pCO2_atmos-pCO2
	if(MODE==2){
		stepdata$BML_DIC<-ifelse(depth==SMLD,eval_BML_DIC(dNO3, BML_DIC),stepdata$DIC)
		stepdata$BML_NO3<-ifelse(depth==SMLD,eval_BML_NO3(dNO3, BML_NO3),ifelse(timestep==mix_day,NO3,BML_NO3))	
	}


	print(stepdata$DIC)
	as.data.frame(as.list(stepdata))
}


################################################################################
#     model control                                              #
################################################################################

model_run<-function(){

	timestep = 1
	model_output<-data.frame(pCO2=init_pCO2,DIC=init_DIC,slDON=0,slDOC=0,airseaFlux=0,deltapCO2=delta_pCO2)
	if(MODE==2){
		model_output$BML_DIC=init_DIC
		model_output$BML_NO3=WINTER_NITRATE
	}
	
	while(timestep<run_length){
		current_state<-model_output[timestep,]
		model_output<-rbind(model_output,eval_timestep(timestep,current_state))
		timestep=timestep+1
	}
	model<-cbind(box,model_output)

	model	
}




