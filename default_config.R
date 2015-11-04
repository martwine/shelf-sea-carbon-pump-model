#################################################################
### DEFAULT CONFIG FILE FOR Shelf Sea Biogeochemistry  model ####
#################################################################


RUN_NAME<-"default_config_phys"

out_dir<-"results"

## mode can be 1-box (surface) or 2-vertical boxes (surface and 
## depth), which are mixed for the winter

#MODE=1
MODE=2

######################################
#                                    #
#  ###  #  # #   #  ### #  ###  ###  #
#  #  # #  #  # #  #    # #    #     #
#  ###  ####   #    ##  # #     ##   #
#  #    #  #   #      # # #       #  #
#  #    #  #   #   ###  #  ### ###   #
#                                    #
######################################

    ## surface mixed layer depth
    SMLD<-20

    ## water column depth
    COLUMN_DEPTH<-100


    ## sine wave temperature variation defined from these values

    MIN_SURFACE_TEMP<-6
    MAX_SURFACE_TEMP<-12


    OFFSET<-59 # days to 1st March (temperature min date)

    ## sine wave pCO2 based based on Weybourne data (from Phil Wilson)
    ## pCO2 max in Jan, min in July
    AVERAGE_pCO2<-390 #uatm
    AMPLITUDE<-14 #uatm

    ## sine wave wind speed seems to reasonably represent North Sea winds e.g. http://www.windfinder.com/windstats/windstatistic_humber_1_northsea.htm (1 knot = 0.5144 m/s) max Jan, min July
    MIN_WIND<-7.2 #m/s
    MAX_WIND<-11.3 #m/s


    ## summer length defines when the BML and SML are mixed into a single box again (they separate at the onset of the spring bloom)
    SUMMER_LENGTH<-150 #days

    ## estimate of total alkalinity from artioli et al 2012
    init_TA<-2300 #uM

    ## delta pCO2 defines starting disequilibium between atmos and surface ocean (+ve indicates pCO2 atmos > pCO2 ocean)
    delta_pCO2<-0


    # scaling factor for gas exchange rate - used to scale Nightngale
    # default value 1.3 to translate quadratic for long-term averaged winds (see Wanninkhof 1992)
    KW_SCALING_FACTOR<-1.3


    SPRING_START_DAY<-88 #jday
    SPRING_DURATION<-15 #days

########################
# Benthic fluxes       #
########################

 # fractions of particulate matter permanently buried per timestep
 BURIAL_FRAC_POC = 0.005
 BURIAL_FRAC_TEPC = 0.005

 #in umol m-2 day-1
 resusp_DIC_SUMMER = 0  
 resusp_DIC_WINTER = 0  
 
######################
# Trawl events       #
######################
 
 TRAWL_DAY = 10
 TRAWL_DIC_RELEASE = 1e6 #12g C per m2



########################
# REDFIELDIAN SETTINGS #
########################
    # i.e settings needed for a purely redfieldian / classic run

    WINTER_NITRATE<-12 #uM

    redfield<-6.6 #C:N for 'Redfieldian' processes

    #degradation rate of POC on bottom
    ### between about 1 and 50 mmol per m2 per day

    POCdeg<-1000 #umol m-2 day-1



############################
# NON-REDFIELDIAN SETTINGS #
############################

    #ratio of POC degradation to PON degradation
    # redfield for redfieldian run
    PONdeg_factor<-redfield

    ## Proportion of spring bloom N turned into semi-labile DON
    # zero for redfieldian run
    nitrate_to_slDON_conv<-0


    # ammonium turnover rate
    # zero for redfieldian run
    NH4_turnover<-0 #uM day-1

    # fraction of summer recycled production which produces TEPC
    TEP_fraction<-0.1

    ## C:N ratio of semi-labile DOM produced during bloom
    slDOM_C_TO_N<-redfield


    # degradation rate of sl_DON in SML
    slDON_deg<-0.02 #day-1

    #degradation rate of sl_DOC in SML as proportion of DON degradation
    slDOC_deg_factor<-0.5


    TEPdeg<-0.005 #day-1


    #"Q10" values for degradations
    Q10_REF_TEMP<-15 #celcius

    Q10<-2


    MIN_BOTTOM_TEMP<-5
    MAX_BOTTOM_TEMP<-10
