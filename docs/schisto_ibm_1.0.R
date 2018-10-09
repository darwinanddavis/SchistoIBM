# Schisto IBM R/Netlogo 1.0

#20-6-18
# init
# early model structure for DEBIBM for schisto

# ---------------------------------------------------------------------------
# ------------------- initial Mac OS and R config ---------------------------
# ---------------------------------------------------------------------------
# if already loaded, uninstall RNetlogo and rJava
p<-c("rJava", "RNetLogo")
remove.packages(p)

# if using Mac OSX El Capitan+ and not already in JQR, download and open JGR 

# for a rJava error, run the following in terminal (src: https://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite)
# sudo ln -s $(/usr/libexec/java_home)/jre/lib/server/libjvm.dylib /usr/local/lib
# then install rJava from source
install.packages("rJava", repos = "https://cran.r-project.org/", type="source")
library(rJava)

# load JGR after downloading 
Sys.setenv(NOAWT=1)
install.packages("JGR")
library(JGR)
Sys.unsetenv("NOAWT")
JGR() # open JGR

# ---------------------------------------------------------------------------
# ------ the rest of the script needs to be run in JGR if using Mac OSX -----
# ---------------------------------------------------------------------------
# install RNetlogo from source if haven't already
# install.packages("/Users/malishev/Documents/Melbourne Uni/Programs/R code/RNetlogo/RNetLogo_1.0-4.tar.gz", repos = NULL, type="source")
install.packages("RNetLogo", repos = "https://cran.r-project.org/", type="source")

# ------------------- model setup ---------------------------
# get packages
packages <- c("adehabitatHR","rgeos","sp", "maptools", "raster","rworldmap","rgdal","dplyr")
if (require(packages)) {
  install.packages(packages,dependencies = T)
  require(packages)
}
lapply(packages,library,character.only=T)

#source DEB and heat budget models from https://github.com/darwinanddavis/MalishevBullKearney
source('DEB.R')

# set dirs
setwd("<your working dir>") # set wd
results.path<- "<dir path to store result outputs>" # set results path

# apply interpolation functions
# velfun<- approxfun(time, micro_sun[,5], rule = 2)

# *************************************************************************************
# ************************** start energetics model setup *******************************

debparams <- "<df of DEB parameters>" # optional if not setting manually
debpars=as.data.frame(read.csv(debparams,header=FALSE))$V1 # read in DEB params

# set core parameters
z=debpars[8] # zoom factor (cm)
F_m = 13290 # max spec searching rate (l/h.cm^2)
kap_X=debpars[11] # digestion efficiency of food to reserve (-)
v=debpars[13] # energy conductance (cm/h)
kap=debpars[14] # kappa, fraction of mobilised reserve to growth/maintenance (-)
kap_R=debpars[15] # reproduction efficiency (-)
p_M=debpars[16] # specific somatic maintenance (J/cm3)
k_J=debpars[18] # maturity maint rate coefficient (1/h)
E_G=debpars[19] # specific cost for growth (J/cm3)
E_Hb=debpars[20] # maturity at birth (J)
E_Hp=debpars[21] # maturity at puberty (J)
h_a=debpars[22]*10^-1 # Weibull aging acceleration (1/h^2)
s_G=debpars[23] # Gompertz stress coefficient (-)

# set thermal respose curve paramters
T_REF = debpars[1]-273
TA = debpars[2] # Arrhenius temperature (K)
TAL = debpars[5] # low Arrhenius temperature (K)
TAH = debpars[6] # high Arrhenius temperature (K)
TL = debpars[3] # low temp boundary (K)
TH = debpars[4] # hight temp boundary (K)

# set auxiliary parameters
del_M=debpars[9] # shape coefficient (-) 
E_0=debpars[24] # energy of an egg (J)
mh = 1 # survivorship of hatchling in first year
mu_E = 585000 # molar Gibbs energy (chemical potential) of reserve (J/mol)
E_sm=186.03*6

# set initial state
E_pres_init = (debpars[16]*debpars[8]/debpars[14])/(debpars[13]) # initial reserve
E_m <- E_pres_init
E_H_init = debpars[21] + 5

#### change inital size here by multiplying by < 0.85 ####
V_pres_init = (debpars[26] ^ 3) * 0.85 
d_V<-0.3
mass <- V_pres_init + V_pres_init*E_pres_init/mu_E/d_V*23.9

# ************************** end energetics model setup *******************************
# *************************************************************************************

# *************************************************************************************
# ************************** start netlogo model setup *******************************

nl.path<- "<dir path to Netlogo program>" 
ver <-"<version number of Netlogo>" # type in version of Netlogo e.g. "6.0.1"
# if error, try adding "/app" to end of dir path for running in Windows and El Capitan for OSX
# nl.path<-"<dir path to Netlogo program>/app" 
NLStart(nl.path)
NLStart(nl.path, nl.jarname = paste0("netlogo-",ver,".jar"))
model.path<- "<dir path to Netlogo model>"
NLLoadModel(model.path)

# ************************** setup netlogo model
### refer to https://github.com/darwinanddavis/SchistoIBM/blob/master/schisto_ibm.html ### 

### environment params ### 
# resources
if(density=="high"){ # high/low food density in env
  NL_food<-100000L # Food patches
}else{ 
  NL_food<-1000L
}
# dRdt

# miracidia
# dMdt

# cercariae 
# dCdt

### individual params ###
NL_days<-100 # No. of days to simulate
VTMIN<- 26 # lower activity thermal limit
VTMAX<- 35 # upper activity thermal limit

# update initial conditions for DEB model 
Es_pres_init<-(E_sm*gutfull)*V_pres_init
acthr<-1
# Tb_init<-20
step = 1/24
debout<-DEB(step = step, z = z, del_M = del_M, F_m = F_m * 
    step, kap_X = kap_X, v = v * step, kap = kap, p_M = p_M * 
    step, E_G = E_G, kap_R = kap_R, k_J = k_J * step, E_Hb = E_Hb, 
    E_Hj = E_Hb, E_Hp = E_Hp, h_a = h_a/(step^2), s_G = s_G, 
    T_REF = T_REF, TA = TA, TAL = TAL, TAH = TAH, TL = TL, 
    TH = TH, E_0 = E_0, E_pres=E_pres_init, V_pres=V_pres_init, E_H_pres=E_H_init, acthr = acthr, breeding = 1, Es_pres = Es_pres_init, E_sm = E_sm)

Es_pres_init<-(E_sm*gutfull)*V_pres_init
X_food<-3000
V_pres<-debout[2]
wetgonad<-debout[19]
wetstorage<-debout[20]
wetfood<-debout[21]
ctminthresh<-120000
# Tairfun<-Tairfun_shd
# Tc_init<-Tairfun(1)+0.1 # Initial core temperature

NL_T_b_min <- VTMIN         # Min foraging T_b
NL_T_b_max <- VTMAX        # Max foraging T_b
NL_reserve_max <- E_m    # Maximum reserve level
NL_reserve <- E_m        # Initial reserve density
NL_maint <- round(p_M, 3)       # Maintenance cost

# ************************** end netlogo model setup **********************************
# *************************************************************************************

# *************************************************************************************
# **************************** start netlogo sim  *************************************

sc<-1 # set no. of desired simualations---for automating writing of each sim results to file. N = N runs
for (i in 1:sc){ # ********************** start netlogo sim loop

NLCommand("set Shade-patches",NL_shade,"set Food-patches",NL_food,"set No.-of-days",NL_days,"set T_opt_lower precision", NL_T_b_min, "2","set T_opt_upper precision", NL_T_b_max, "2",
"set reserve-level", NL_reserve, "set Maximum-reserve", NL_max_reserve, "set Maintenance-cost", NL_maint,
"set Movement-cost precision", NL_move, "3", "set gutthresh", NL_gutthresh, 'set gutfull', gutfull, 'set V_pres precision', V_pres, "5", 'set wetstorage precision', wetstorage, "5", 
'set wetfood precision', wetfood, "5", 'set wetgonad precision', wetgonad, "5","setup")

NL_ticks<-NL_days / (2 / 60 / 24) # No. of NL ticks (measurement of days)
NL_T_opt_l<-NLReport("[T_opt_lower] of turtle 0")
NL_T_opt_u<-NLReport("[T_opt_upper] of turtle 0")

# data frame setup for homerange polygon
turtles<-data.frame() # make an empty data frame
NLReport("[X] of turtle 0"); NLReport("[Y] of turtle 0")
who<-NLReport("[who] of turtle 0")

debcall<-0 # check for first call to DEB
stepcount<-0 # DEB model step count

for (i in 1:NL_ticks){ # ********************** start netlogo sim loop
stepcount<-stepcount+1
NLDoCommand(1, "go")

######### Example of reporting NL variables
shade <- NLGetAgentSet("in-shade?","turtles", as.data.frame=T); shade<-as.numeric(shade) # returns an agentset of whether turtle is currently on shade patch

# **************************** start deb sim  *************************************

if(stepcount==1) { # run DEB loop every time step (2 mins)
stepcount<-0

# report activity state
actstate<-NLReport("[activity-state] of turtle 0")
 # Reports true if turtle is in food 
actfeed<-NLGetAgentSet("in-food?","turtles", as.data.frame=T); actfeed<-as.numeric(actfeed)
 
n<-1 # time steps
step<-2/1440 # step size (2 mins). For hourly: 1/24
# update direct movement cost
if(actstate == "S"){
	NLCommand("set Movement-cost", NL_move)
	}else{
		NLCommand("set Movement-cost", 1e-09)
		} 
# if within activity range, it's daytime, and gut below threshold 
if(Tbs$Tc>=VTMIN & Tbs$Tc<=VTMAX & Zen!=90 & gutfull<=NL_gutthresh){ 
  acthr=1 # activity state = 1 
if(actfeed==1){ # if in food patch
	X_food<-NLReport("[energy-gain] of turtle 0") # report joules intake
	}
	}else{
		X_food = 0 
		acthr=0
		}

if(debcall==0){ # calculate deb output 
	# initialise DEB
	debout<-matrix(data = 0, nrow = n, ncol = 26)
	deb.names<-c("E_pres","V_pres","E_H_pres","q_pres","hs_pres","surviv_pres","Es_pres","cumrepro","cumbatch","p_B_past","O2FLUX","CO2FLUX","MLO2","GH2OMET","DEBQMET","DRYFOOD","FAECES","NWASTE","wetgonad","wetstorage","wetfood","wetmass","gutfreemass","gutfull","fecundity","clutches")
	colnames(debout)<-deb.names
	# initial conditions
	debout<-DEB(E_pres=E_pres_init, V_pres=V_pres_init, E_H_pres=E_H_init, acthr = acthr, Tb = Tb_init, breeding = 1, Es_pres = Es_pres_init, E_sm = E_sm, step = step, z, del_M = del_M, F_m = F_m * 
    step, kap_X = kap_X, v = v * step, kap = kap, p_M = p_M * 
    step, E_G = E_G, kap_R = kap_R, k_J = k_J * step, E_Hb = E_Hb, 
    E_Hj = E_Hb, E_Hp = E_Hp, h_a = h_a/(step^2), s_G = s_G, 
    T_REF = T_REF, TA = TA, TAL = TAL, TAH = TAH, TL = TL, 
    TH = TH, E_0 = E_0)
	debcall<-1
	}else{
		debout<-DEB(step = step, z = z, del_M = del_M, F_m = F_m * 
    step, kap_X = kap_X, v = v * step, kap = kap, p_M = p_M * 
    step, E_G = E_G, kap_R = kap_R, k_J = k_J * step, E_Hb = E_Hb, 
    E_Hj = E_Hb, E_Hp = E_Hp, h_a = h_a/(step^2), s_G = s_G, 
    T_REF = T_REF, TA = TA, TAL = TAL, TAH = TAH, TL = TL, 
    TH = TH, E_0 = E_0, 
		  X=X_food,acthr = acthr, Tb = Tbs$Tc, breeding = 1, E_sm = E_sm, E_pres=debout[1],V_pres=debout[2],E_H_pres=debout[3],q_pres=debout[4],hs_pres=debout[5],surviv_pres=debout[6],Es_pres=debout[7],cumrepro=debout[8],cumbatch=debout[9],p_B_past=debout[10])
		}
mass<-debout[22]
gutfull<-debout[24]
NL_reserve<-debout[1]
V_pres<-debout[2]
wetgonad<-debout[19]
wetstorage<-debout[20]
wetfood<-debout[21] 

#update NL wetmass properties 
NLCommand("set V_pres precision", V_pres, "5")
NLDoCommand("plot xcor ycor")
NLCommand("set wetgonad precision", wetgonad, "5")
NLDoCommand("plot xcor ycor")
NLCommand("set wetstorage precision", wetstorage, "5")
NLDoCommand("plot xcor ycor")
NLCommand("set wetfood precision", wetfood, "5")
NLDoCommand("plot xcor ycor") 
} # **************************** end deb sim loop

NLCommand("set reserve-level", NL_reserve) # update reserve
NLCommand("set gutfull", debout[24])# update gut level

# **************************** end deb sim  *************************************

# generate results, with V_pres, wetgonad, wetstorage, and wetfood from debout
if(i==1){
	results<-cbind(tick,Tb,rate,shade,V_pres,wetgonad,wetstorage,wetfood,NL_reserve) 
	}else{
		results<-rbind(results,c(tick,Tb,rate,shade,V_pres,wetgonad,wetstorage,wetfood,NL_reserve))
		}
results<-as.data.frame(results)

# generate data frames for homerange polygon
if (tick == NL_ticks - 1){
	X<-NLReport("[X] of turtle 0"); head(X)
	Y<-NLReport("[Y] of turtle 0"); head(Y)
	turtles<-data.frame(X,Y)
	who1<-rep(who,NL_ticks); who # who1<-rep(who,NL_ticks - 1); who 
	turtledays<-rep(1:NL_days,length.out=NL_ticks,each=720) 
	turtle<-data.frame(ID = who1,days=turtledays)
	turtles<-cbind(turtles,turtle)
	}

} # *************** end netlogo sim loop

# get hr data
spdf<-SpatialPointsDataFrame(turtles[1:2], turtles[3]) # creates a spatial points data frame (adehabitatHR package)
homerange<-mcp(spdf,percent=95)

# writing new results
if (exists("results")){  #if results exist
	sc<-sc-1 
	nam <- paste("results", sc, sep = "") # generate new name with added sc count
	rass<-assign(nam,results) #assign new name to results. call 'results1, results2 ... resultsN'
	namh <- paste("turtles", sc, sep = "")  #generate new name with added sc count
	rassh<-assign(namh,turtles) #assign new name to results. call 'results1, results2 ... resultsN'
	nams <- paste("spdf", sc, sep = "") 
	rasss<-assign(nams,spdf) 
	namhr <- paste("homerange", sc, sep = "")  
	rasshr<-assign(namhr,homerange) 

	fh<-results.path; fh
	for (i in rass){
		# export all results
		write.table(results,file=paste(fh,nam,".R",sep=""))
		}
	for (i in rassh){
		# export turtle location data
		write.table(turtles,file=paste(fh,namh,".R",sep=""))
		}
		#export NL plots
		month<-"sep"
		#spatial plot
		sfh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_move","",sep="");sfh
		NLCommand(paste("export-plot \"Spatial coordinates of transition between activity states\" \"",results.path,sfh,".csv\"",sep=""))
		#temp plot 
		tfh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_temp",sep="")
		NLCommand(paste("export-plot \"Body temperature (T_b)\" \"",results.path,tfh,".csv\"",sep=""))
		#activity budget
		afh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_act","",sep="");afh
		NLCommand(paste("export-plot \"Global time budget\" \"",results.path,afh,".csv\"",sep=""))
		#text output
		xfh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_txt",sep="");xfh
		NLCommand(paste("export-output \"",results.path,xfh,".csv\"",sep=""))
		#gut level
		gfh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_gut","",sep="");gfh
		NLCommand(paste("export-plot \"Gutfull\" \"",results.path,gfh,".csv\"",sep=""))
		#wet mass 
		mfh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_wetmass","",sep="");mfh
		NLCommand(paste("export-plot \"Total wetmass plot\" \"",results.path,mfh,".csv\"",sep=""))
		#movement cost (loco) 
		lfh<-paste(month,NL_days,round(mass,0),NL_shade,as.integer(NL_food*10),"_",sc,"_loco","",sep="");lfh
		NLCommand(paste("export-plot \"Movement costs\" \"",results.path,lfh,".csv\"",sep=""))
	}
} # ********************** end netlogo sim loop

# example files of results output in results.path
list.files(results.path)

#*************************** end netlogo sim *******************************
#***************************************************************************


