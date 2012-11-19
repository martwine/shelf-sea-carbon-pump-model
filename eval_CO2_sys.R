# evaluate CO2 system parameters given TA and DIC, salinity and temperature

# this turned out to be easy thanks to the seacarb package! More might end up in here later
library("seacarb")

eval_CO2_sys<-function(TA,DIC,S,T){
	carb(flag=15, var1=TA, var2=DIC, S=35, T=25, P=0, Pt=0, Sit=0, k1k2="x", kf="x", ks="d", pHscale="T", b="l10")
}


