
source("eval_CO2_sys.R")
library("seacarb")

################################################################################################
###############   Initialisation ########################################################

	#calc bottom mixed layer depth
	BMLD<-COLUMN_DEPTH-SMLD

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