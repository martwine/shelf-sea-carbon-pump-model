#initalise model

source("eval_CO2_sys.R")
library("seacarb")
source("K_calcs_Johnson_OS.R")




###############################################################################
###### physical set-up ########################################################

#surface mixed layer depth
sMLD<-20 #m

#water column depth
column_depth<-100 #m

#calc bottom mixed layer depth
bMLD<-column_depth-sMLD

################################################################################################
###############    Boundary conditions  ########################################################


#temperature is important for respiration rates and CO2 system parameters, also air-sea flux

#min_surface_temp<-4
#max_surface_temp<-16

min_surface_temp<-10
max_surface_temp<-10



#initially let's just have a sine wave temp

surface_temp<-min_surface_temp+(0.5*(max_surface_temp-min_surface_temp))*(1+sin((seq(from=-90,to=270,length.out=365)*pi)/180))



#nitrate concentrations drive the model initally. Spring bloom is a simple linear decrease from winter nitrate to zero, where start date and duration are defined

winter_nitrate<-12 #uM
bloom_start_day<-88 #jday
bloom_duration<-15 #days

#summer length defines when the BML and SML are mixed into a single box again
summer_length<-150 #days

#define run_length as from 1st Jan 1 year, to immediately before spring bloom the next
run_length=365+bloom_start_day


#calc jday of mixing event
mix_day<-bloom_start_day+bloom_duration+summer_length

#define surface nitrate 'envelope'
nday<-c(1,bloom_start_day,bloom_start_day+bloom_duration,run_length)
nitrate<-c(winter_nitrate,winter_nitrate,0,0)
nitrate<-approx(nday,nitrate,n=run_length)


#winter atmospheric pCO2 North Sea (est from Weybourne (Phil Wilson))

winter_pCO2_atmos<-397

#get a simple seasonal cycle of atmospheric CO2
average_pCO2_atmos<-390
min_pCO2<-average_pCO2_atmos-(winter_pCO2_atmos-average_pCO2_atmos)
pCO2_atmos<-(min_pCO2+7*(1+cos((seq(from=0,to=360,length.out=365)*pi)/180)))
#pCO2_atmos<-seq(from=average_pCO2_atmos, to=average_pCO2_atmos,length.out=365) #constant pCO2 for reference run


#estimate of TA from artioli et al 2012
init_TA<-2300 #uM

#initalise DIC from pCO2 and estimated TA assuming equilibration
init_DIC<-1e6*carb(flag=24,winter_pCO2_atmos,init_TA*1e-6)$DIC[1]


#initalise_data
box<-as.data.frame(nitrate)
colnames(box)<-c("day","nitrate")
box$pCO2_atmos<-c(pCO2_atmos,pCO2_atmos)[1:run_length]
box$temp<-c(surface_temp,surface_temp)[1:run_length]
#delta nitrate column
box$dNO3<-c(0,diff(box$nitrate))

#make a second box (BML) which simply respires all the POM at redfield and mixes back in 

##########################################################################################
################  Parameter values  ######################################################
##########################################################################################

redfield<-6.6 #redfieldian C:N

# ammonium turnover rate
#rho_NH4<-0.25 #uM day-1

# C:N ratio of semi-labile DOM produced during bloom
sl_DOM_CtoN<-redfield

# Proportion of spring bloom N turned into semi-labile DON
nitrate_slDOM_conv<-0

# degradation rate of sl_DON
sl_DON_deg<-0.02 #day-1

#degradation rate of sl_DOC
sl_DOC_deg<-0.01 #day-1

#air-sea CO2 flux params
k_w<-2.9 #m day-1


#########################################################################################
############## process and state variable evaluation functions ##########################
#########################################################################################

calc_as_flux<-function(pCO2,pCO2_atmos,temp){
	delta_pCO2<-KH_Molar_per_atmosphere("CO2",temp,S=35)*(pCO2_atmos-pCO2)*1e-6 #mol/l
	flux<-1e6*delta_pCO2*1000*k_w #umol/m2/day

	
}

calc_DIC_uptake_from_NO3<-function(dNO3){
	-dNO3*redfield
}

calc_slDON_deg<-function(slDON){
	slDON*sl_DON_deg
}


calc_slDOC_deg<-function(slDOC){
	slDOC*sl_DOC_deg
}

calc_prod_slDON<-function(dNO3){
	-dNO3*nitrate_slDOM_conv
}

calc_prod_slDOC<-function(dNO3){
	calc_prod_slDON(dNO3)*sl_DOM_CtoN
}


#instantly remineralise POM in the bottom layer
calc_NO3_BML<-function(dNO3){
	-dNO3*(1-nitrate_slDOM_conv)*sMLD/bMLD
}

calc_DIC_BML<-function(dNO3){
	-dNO3*(1-nitrate_slDOM_conv)*redfield*sMLD/bMLD
}


eval_slDON<-function(dNO3, slDON){
	slDON+calc_prod_slDON(dNO3)-calc_slDON_deg(slDON)
}

eval_slDOC<-function(dNO3, slDOC){
	slDOC+calc_prod_slDOC(dNO3)-calc_slDOC_deg(slDOC)
}

eval_DIC<-function(DIC,dNO3,pCO2,pCO2_atmos,temp,slDOC,depth){
	DIC-calc_DIC_uptake_from_NO3(dNO3)+((calc_as_flux(pCO2,pCO2_atmos,temp)/depth)/1000)+calc_slDOC_deg(slDOC)
}

eval_BML_DIC<-function(dNO3, BML_DIC){
	BML_DIC+calc_DIC_BML(dNO3)
}

eval_BML_NO3<-function(dNO3, BML_NO3){
	BML_NO3+calc_NO3_BML(dNO3)
}

calc_mix<-function(sml_conc,bml_conc){
	(sml_conc*sMLD + bml_conc*bMLD)/column_depth
}


###################################################################################
##       Evaluate a timestep                                                      #
###################################################################################

eval_timestep<-function(timestep,current_state){
	timestep_row<-box[timestep,]	
	pCO2_atmos<-timestep_row$pCO2_atmos
	dNO3<-timestep_row$dNO3
	temp<-timestep_row$temp

	DIC<-current_state$DIC
	slDOC<-current_state$slDOC
	slDON<-current_state$slDON
	pCO2<-current_state$pCO2
	BML_DIC<-current_state$BML_DIC
	BML_NO3<-current_state$BML_NO3
	stepdata<-new.env()

	if(timestep<mix_day){
		stepdata$slDOC<-eval_slDOC(dNO3, slDOC)
		stepdata$slDON<-eval_slDON(dNO3, slDON)
		stepdata$airseaFlux<-calc_as_flux(pCO2,pCO2_atmos,temp)
		if(timestep<bloom_start_day){
			stepdata$DIC<-eval_DIC(DIC,dNO3,pCO2,pCO2_atmos,temp,slDOC,depth=column_depth)
		} else {
			stepdata$DIC<-eval_DIC(DIC,dNO3,pCO2,pCO2_atmos,temp,slDOC,depth=sMLD)
		}		
		stepdata$pCO2<-carb(flag=15,init_TA*1e-6,stepdata$DIC*1e-6)$pCO2[1]
		stepdata$BML_DIC<-eval_BML_DIC(dNO3, BML_DIC)
		stepdata$BML_NO3<-eval_BML_NO3(dNO3, BML_NO3)
	} else if (timestep==mix_day) {
		#do the mixing
		mixed_DIC<-calc_mix(DIC,BML_DIC)
		mixed_NO3<-calc_mix(0,BML_NO3)
		mixed_slDON<-calc_mix(slDON,0)
		mixed_slDOC<-calc_mix(slDOC,0)
		mixed_pCO2<-carb(flag=15,init_TA*1e-6,mixed_DIC*1e-6)$pCO2[1]
		
		#then run the timestep
		stepdata$slDOC<-eval_slDOC(dNO3, mixed_slDOC)
		stepdata$slDON<-eval_slDON(dNO3, mixed_slDON)
		stepdata$airseaFlux<-calc_as_flux(mixed_pCO2,pCO2_atmos,temp)
		stepdata$DIC<-eval_DIC(mixed_DIC,dNO3,mixed_pCO2,pCO2_atmos,temp,slDOC,depth=column_depth)
		stepdata$pCO2<-carb(flag=15,init_TA*1e-6,stepdata$DIC*1e-6)$pCO2[1]
		stepdata$BML_DIC<-stepdata$DIC
		stepdata$BML_NO3<-mixed_NO3			
			
			
			
	} else {
		stepdata$slDOC<-eval_slDOC(dNO3, slDOC)
		stepdata$slDON<-eval_slDON(dNO3, slDON)
		stepdata$airseaFlux<-calc_as_flux(pCO2,pCO2_atmos,temp)
		stepdata$DIC<-eval_DIC(DIC,dNO3,pCO2,pCO2_atmos,temp,slDOC,depth=column_depth)
		stepdata$pCO2<-carb(flag=15,init_TA*1e-6,stepdata$DIC*1e-6)$pCO2[1]
		stepdata$BML_DIC<-stepdata$DIC
		stepdata$BML_NO3<-BML_NO3
	}

	print(stepdata$DIC)
	as.data.frame(as.list(stepdata))
}


################################################################################
#     model control                                                            #
################################################################################

model_run<-function(run_name){

	timestep = 1
	model_output<-data.frame(pCO2=winter_pCO2_atmos,DIC=init_DIC,slDON=0,slDOC=0,airseaFlux=0,BML_DIC=init_DIC, BML_NO3=winter_nitrate)

	
	while(timestep<run_length){
		current_state<-model_output[timestep,]
		model_output<-rbind(model_output,eval_timestep(timestep,current_state))
		timestep=timestep+1
	}
	model<-cbind(box,model_output)
	write.csv(model,file=paste(run_name,".csv"))
	pdf(paste(run_name,".pdf"))
		
	
		par(mfrow=c(6,1),mar=c(3,4,0,0))
		plot(model$day,model$temp,type="l")
		plot(model$day,model$pCO2_atmos,type="l")
		plot(model$day,model$nitrate,type="l")
		plot(model$day,model$DIC,type="l")
		plot(model$day,model$slDOC,type="l")
		plot(model$day,model$slDON,type="l")	
	dev.off()
	print(cumsum(model$airseaFlux))
	model	
}

model_run("two-box ref run 1")


