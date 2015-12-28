
source("eval_CO2_sys.R")
library("seacarb")

################################################################################################
###############   Initialisation ########################################################


	#if mode is 1, set SMLD to COLUMN_DEPTH and bottom temps to top temps
	if(MODE==1)
		{
			SMLD<-COLUMN_DEPTH
			MIN_BOTTOM_TEMP<-MIN_SURFACE_TEMP
			MAX_BOTTOM_TEMP<-MAX_SURFACE_TEMP
		}
	#calc bottom mixed layer depth
	BMLD<-COLUMN_DEPTH-SMLD

    #round to integers to prevent breakage
    SUMMER_LENGTH<-round(SUMMER_LENGTH)
    SPRING_START_DAY<-round(SPRING_START_DAY)
    TRAWL_DAY<-round(TRAWL_DAY)

	#run_length=365+SPRING_START_DAY

	#multiyear run
	run_length=(365*10)+SPRING_START_DAY

	#make jday term for mixing events etc
	jday<-rep(seq(365),ceiling(run_length/365))[1:run_length]


	#get number of repeats of each of the annual forcings
	n_repeats=ceiling(run_length/365)

	#Temperature cycle calculated from max and min values (1st March min, 1st September max)

	#APPLY warming for climate scenarios
		MAX_SURFACE_TEMP=	MAX_SURFACE_TEMP+CLIMATE_WARMING

	MIN_SURFACE_TEMP=MIN_SURFACE_TEMP+CLIMATE_WARMING
	## keep things simple and control bottom temps from surface temps (see init file)
	
	MIN_BOTTOM_TEMP<-MIN_SURFACE_TEMP
	MAX_BOTTOM_TEMP<-MAX_SURFACE_TEMP-SUMMER_T_DIFF
	
	
	
	surface_temp_cycle<-MIN_SURFACE_TEMP+(0.5*(MAX_SURFACE_TEMP-MIN_SURFACE_TEMP))*(1+sin((seq(from=-90,to=270,length.out=365)*pi)/180))

	surface_temp<-c(surface_temp_cycle[(365-OFFSET):365],surface_temp_cycle[1:(365-OFFSET-1)])

	bottom_temp_cycle<-MIN_BOTTOM_TEMP+(0.5*(MAX_BOTTOM_TEMP-MIN_BOTTOM_TEMP))*(1+sin((seq(from=-90,to=270,length.out=365)*pi)/180))

	bottom_temp<-c(bottom_temp_cycle[(365-OFFSET):365],bottom_temp_cycle[1:(365-OFFSET-1)])


	#calc jday of mixing event
	mix_day<-SPRING_START_DAY+SPRING_DURATION+SUMMER_LENGTH

	#define surface nitrate 'envelope'
	#nday<-c(1,SPRING_START_DAY,SPRING_START_DAY+SPRING_DURATION,run_length)
	nday<-c(1,SPRING_START_DAY,SPRING_START_DAY+SPRING_DURATION,365)
	nitrate<-c(WINTER_NITRATE,WINTER_NITRATE,0,0)
	nitrate<-approx(nday,nitrate,n=365)
	nitrate_cycle<-rep(nitrate$y,times=n_repeats)[1:run_length]

	#get a simple seasonal cycle of atmospheric CO2
	min_pCO2<-AVERAGE_pCO2-(AMPLITUDE/2)
	pCO2_atmos<-(min_pCO2+(AMPLITUDE/2)*(1+cos((seq(from=0,to=360,length.out=365)*pi)/180)))
	
	init_pCO2<-min_pCO2+AMPLITUDE-delta_pCO2


	#initalise DIC from pCO2 and estimated TA 
	init_DIC<-1e6*carb(flag=24,init_pCO2,init_TA*1e-6,T=surface_temp[1])$DIC[1]
	print(paste("init DIC:", init_DIC))
	print("surface temps")
	print(c(MAX_SURFACE_TEMP,MIN_SURFACE_TEMP))
	
	#wind speed cycle
	wind_speed<-MIN_WIND+(0.5*(MAX_WIND-MIN_WIND))*(1+sin((seq(from=90,to=450,length.out=365)*pi)/180))
	
	#calculate orgC degradation rates as function of orgN degradation rates
	slDOC_deg<-slDOC_deg_factor*slDON_deg
	POCdeg<-PONdeg*POCdeg_factor

  #calculate release per trawl of POC - assumes density of sediments is 1 therefore trawl depth of 1cm resuspends 10 kg of material
	TRAWL_POC_RELEASE = SEDIMENT_POC_CONTENT*TRAWL_DEPTH*10
	print(TRAWL_POC_RELEASE)
	#converyt to umoles
	TRAWL_POC_RELEASE = TRAWL_POC_RELEASE*1e6/12
	print(TRAWL_POC_RELEASE)
	

	#initalise_data
	box<-data.frame(seq(run_length),nitrate_cycle)
	colnames(box)<-c("day","nitrate")
	box$jday<-jday
	box$pCO2_atmos<-rep(pCO2_atmos,times=n_repeats)[1:run_length]
	box$temp<-rep(surface_temp,times=n_repeats)[1:run_length]
	box$bottomtemp<-rep(bottom_temp,times=n_repeats)[1:run_length]
	box$wind<-rep(wind_speed,times=n_repeats)[1:run_length]
	
	#delta nitrate column - gives change between each timestep and previous
	box$dNO3<-c(0,ifelse(diff(box$nitrate)<0,diff(box$nitrate),0))


