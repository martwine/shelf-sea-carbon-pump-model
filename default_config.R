#################################################################
### DEFAULT CONFIG FILE FOR Shelf Sea Biogeochemistry  model ####
#################################################################


RUN_NAME<-"default_config"

out_dir<-"results"

## mode can be 1-box (surface) or 2-vertical boxes (surface and 
## depth), which are mixed for the winter

#MODE=1
MODE=2


#################################
##### WATER COLUMN DEPTHS   #####
#################################

## surface mixed layer depth
SMLD<-20

## water column depth
COLUMN_DEPTH<-100

##################################
###### BOUNDARY CONDITIONS  ######
##################################

## sine wave temperature variation defined from these values

MIN_SURFACE_TEMP<-4
MAX_SURFACE_TEMP<-16

MIN_BOTTOM_TEMP<-2
MAX_BOTTOM_TEMP<-8

OFFSET<-59 # days to 1st March (temperature min date)

## sine wave pCO2 based based on Weybourne data (from Phil Wilson)
## pCO2 max in Jan, min in July
AVERAGE_pCO2<-390 #uatm
AMPLITUDE<-14 #uatm

## sine wave wind speed seems to reasonably represent North Sea winds e.g. http://www.windfinder.com/windstats/windstatistic_humber_1_northsea.htm (1 knot = 0.5144 m/s) max Jan, min July
MIN_WIND<-7.2 #m/s
MAX_WIND<-11.3 #m/s

## model is initialised with winter nitrate and the bloom period is described by start day and length to start the model off

WINTER_NITRATE<-12 #uM
BLOOM_START_DAY<-88 #jday
BLOOM_DURATION<-15 #days

## summer length defines when the BML and SML are mixed into a single box again (they separate at the onset of the spring bloom)
SUMMER_LENGTH<-150 #days

## estimate of total alkalinity from artioli et al 2012
init_TA<-2300 #uM

## delta pCO2 defines starting disequilibium between atmos and surface ocean (+ve indicates pCO2 atmos > pCO2 ocean)
delta_pCO2<-0


##################################
######## PARAMETER VALUES ########
##################################

# ammonium turnover rate
NH4_turnover<-0.25 #uM day-1

# fraction of summer recycled production which produces TEPC
TEP_fraction<-0.2

redfield<-6.6 #C:N for 'Redfieldian' processes

## C:N ratio of semi-labile DOM produced during bloom
slDOM_C_TO_N<-redfield

## Proportion of spring bloom N turned into semi-labile DON
nitrate_to_slDON_conv<-0.5

# degradation rate of sl_DON in SML
slDON_deg<-0.02 #day-1

#degradation rate of sl_DOC in SML
slDOC_deg<-0.01 #day-1

#degradation rate of PON in BML
### need to find realistic values of these

PONdeg<-0.02 #day-1
POCdeg<-0.01 #day-1

TEPdeg<-0.005 #day-1


#"Q10" values for degradations
Q10_REF_TEMP<-15 #celcius

Q10<-1



