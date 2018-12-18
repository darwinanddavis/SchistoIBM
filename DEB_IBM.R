## DEB IBM 
# check bottom of page for diagnostics for running netlogo in rstudio

# version 

#16-11-18
# changed resource wave eq to r = r + alpha * pars["r"] * sin(2 * pi * t/p)  
# new damage density deb params (IndividualModel_IBM2.c, ILL_shrink_damageA5.Rda)
# added alpha and periodicity (p) param space to resource dynamics  

# 11-12-18
# changed overdamped and periodic food dynamics to cyclical resource dynamics
# fixed cyclical resource dynamics

# 28-11-18
# added overdamped and periodic food dynamics to nl loop 

# 27-11-18
# fixed NAs in 'create snails' command (Env_G param)  
# list to check if 'create snails' command isn't producing NAs from Env_G
# set pop density outputs in NL loop to integer to pass into Env_G and rbinom func

# 23-11-18
# all user inputs at beginning of doc  

# 22-11-19 
# added debfunction.txt and pars.txt for defining params  

# 19-11-18  
# added  "DEB_INF_GUTS_IBM_1.1.nlogo" as test model    

# TO DO

# set netlogo in windows to /app folder and /java for mac
# remove .o .so and .dll from github 

# change r and food 
  # - plot food against r as heatmap
# define 
	# - starvation-related hazard rate
	# - pars and debfunction
# find papers on periodicity in resource loads in pops (Nisbet, Gurney, daphnia, etc)
# E.g. Nisbet Gurney Population dynamics in a periodically varying environment
# verify what volume of food density is reasonable  (F)
# - sin wave of resource change 
# - step function of resources (on/off season) (most unrealistic)   
# - non-regenerative detritus (event-based)   
# R = algae supply rate, don't vary K

# heat map of where peaks or resources and peaks of cercariae occur    

# fix days parameter in NL   
# output NL plots to R 

# OUTPUTS
# survival and shell length results from DEBstep
# plot of rP vs. P/V (parasite biomass outcome)  

# --------------- --------------- --------------- ---------------
# --------------- --------------- --------------- ---------------

### Files required  
# "DEB_INF_GUTS_IBM.nlogo"
# "ILL_shrink_damageA5.Rda"
# "IndividualModel_IBM2.c"
# "IndividualModel_IBM.so" # Mac OSX. generated from C
# "IndividualModel_IBM.o" # Mac OSX. generated from C  
# "IndividualModel_IBM.dll" # Windows. generated from C  

test.java <- 1 # 1 = run java diagnostics  
mac <- 1 # mac or windows system? 1 = mac, 0 = windows 
gui <- 1 # display the gui? 1 = yes, 0 = no
pck <- 0 # if not already, install rnetlogo and rjava from source? 1 = yes, 0 = already installed 
save_to_file <- 0 # 1 = save plots to local dir, 0 = plot in current R session
 
# run java test  
if(test.java==1){
  require(RCurl)
  script <- getURL("https://raw.githubusercontent.com/darwinanddavis/SchistoIBM/master/mac/java_test.R", ssl.verifypeer = FALSE)
  eval(parse(text = script))
}

#### :three: [GCC compiler in R (unconfirmed)](https://stackoverflow.com/questions/1616983/building-r-packages-using-alternate-gcc)
#### [Running Netlogo 6.0.+](https://github.com/NetLogo/NetLogo/issues/1282)

# Search "@netlogo" for netlogo code in file  

# ------------------------------------------------------------------------------------------------------------- 

# ------------------------------------------------------- 
# ------------ Running NetLogo in Mac -------------------

# if using Mac OSX El Capitan+ and not already in JQR, download and open JGR 
if(mac==1){
  install.packages('JGR',,'http://www.rforge.net/')
  library(JGR)
  JGR::JGR()
}

# ------------------------------------------------------------------
# ------------------------- JGR onwards ----------------------------
# ------------------------------------------------------------------

####################################end  set user inputs ####################################end 
# isolate sensitive data:
# "ILL_shrink_damageA5.Rda"
seninf <- 1 # 1 = keep files local; 0 = make files public  
if(seninf==1){pp <- "/Users/malishev/Documents/Emory/research/schisto_ibm/SchistoIBM_/"}else{pp <-""}

# set paths
mac <- 1 # mac or windows system? 1 = mac, 0 = windows 
gui <- 1 # display the gui? 1 = yes, 0 = no
pck <- 0 # if not already, install rnetlogo and rjava from source? 1 = yes, 0 = already installed 
save_to_file <- 0 # 1 = save plots to local dir, 0 = plot in current R session

# set dir paths  
wd <- "/Users/malishev/Documents/Emory/research/schisto_ibm/SchistoIBM" # set working directory  
ver_nl <-"6.0.4" # type in Netlogo version. found in local dir. 
ver_gcc <-"4.6.3" # NULL # type in gcc version (if known). leave as "NULL" if unknown   
nl.path <- "/Users/malishev/Documents/Melbourne Uni/Programs/" # set path to Netlogo program location

# -------------------------------------------- -------------------------------------------- -------------------------------------------- 
# ------------------------------ define starting conditions for simulation model @netlogo ----------------------------------------------
# -------------------------------------------- -------------------------------------------- --------------------------------------------  
resources <- "cyclical" # set resources: "cyclical" or "event"
n.ticks <- 150 # set number of simulation ticks
day <- 1 # number of days to run simulation  
cs <- list() # diagnostics list for checking NAs in 'create snails' command  

####################################end set user inputs ####################################end 

# set model paths
setwd(wd)
nl.model <- list.files(pattern="*.nlogo")[1] # Netlogo model
nl.path <- paste0(nl.path,"NetLogo ",ver_nl,"/Java/"); nl.path
model.path <- paste0(wd,"/"); model.path # set path to Netlogo model  

# if already loaded, uninstall RNetlogo and rJava
if(pck==1){
  p<-c("rJava", "RNetLogo"); remove.packages(p)
  # then install rJava and RNetLogo from source
  install.packages("rJava", repos = "https://cran.r-project.org/", type="source"); library(rJava)
  install.packages("RNetLogo", repos = "https://cran.r-project.org/", type="source"); library(RNetLogo)
}
library(rJava); library(RNetLogo) 

# check pck versions
installed.packages()["RNetLogo","Version"] 
installed.packages()["rJava","Version"]

# check rJava version  
.jinit()
.jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
# get latest Java/Oracle version: https://www.oracle.com/technetwork/java/javase/downloads/index-jsp-138363.html  

### install relevant packages   
packages <- c("Matrix","deSolve","mvtnorm","LaplacesDemon","coda","adaptMCMC","ggplot2", "RCurl","RColorBrewer")

if(require(packages)){
  install.packages(packages,dependencies = T)
  require(packages)
}
lapply(packages,library,character.only=T)

#####################################################
# load plot function 
require(RCurl)
script <- getURL("https://raw.githubusercontent.com/darwinanddavis/plot_it/master/plot_it.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

require(RColorBrewer)
display.brewer.all()
# Set global plotting parameters
cat("plot_it( \n0 for presentation, 1 for manuscript, \nset colour for background, \nset colour palette 1. use 'display.brewer.all()', \nset colour palette 2. use 'display.brewer.all()', \nset alpha for colour transperancy, \nset font style \n)")

plot_it(0,"blue","Spectral","Greens",1,"mono") # set plot function params       
plot_it_gg("white") # same as above for ggplot   
#####################################################

### Install rtools and gcc for using C code and coda package 
#### https://cran.r-project.org/bin/macosx/tools/

# make below paste function from gcc
if(mac==1){ #### Mac OSX
  rtools <- "/usr/local/clang6/bin"
  gcc <- paste0("usr/local/clang6/gcc-",ver_gcc,"/bin")
  # Mac OSX
  }else{ #### Windows 
  rtools <- "C:\\Rtools\\bin"
  gcc <- paste0("C:\\Rtools\\gcc-",ver_gcc,"\\bin")
}

#### point to path on comp to access rtools and gcc for C compiler  
path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
new_path <- c(rtools, gcc, path)
new_path <- new_path[!duplicated(tolower(new_path))]
Sys.setenv(PATH = paste(new_path, collapse = ";"))

if(mac==1){
    # dyn.unload("IndividualModel_IBM2.so") # unLoad .so (Mac OSX
    system("R CMD SHLIB IndividualModel_IBM2.c") # generates .o and .so files 
    dyn.load("IndividualModel_IBM2.so") # Load .so (Mac OSX)
}else{
  # compile my model from C definition
  dyn.unload("IndividualModel_IBM2.dll") # unload dll (Windows only)
  system("R CMD SHLIB IndividualModel_IBM2.c")
  dyn.load("IndividualModel_IBM2.dll") # Load dll (Windows only)
}

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# point where windows/mac merge  

# load DEB starvation model parameters and create mcmc (and convert mcmc chain to coda format)
samps = readRDS(paste0(pp,"ILL_shrink_damageA5.Rda"))
# lpost = samps$log.p #
# samps = convert.to.coda(samps)
# samps = cbind(samps, lpost)
# samps <- as.mcmc(samps[, c("iM", "k", "M", "EM", "Fh", "muD", "DR", "fe", "yRP",
#                            "ph", "yPE", "iPM", "eh", "mP", "alpha", "yEF", "LM",
#                            "kd", "z", "kk", "hb", "theta", "mR", "yVE", "yEF2",
#                            "sd.LI1", "sd.LU1", "sd.EI1", "sd.EU1", "sd.W1", 
#                            "sd.LI2", "sd.LU2", "sd.EI2", "sd.EU2", "sd.W2",
#                            "lpost")])
samps <- as.mcmc(samps[, c("iM", "k", "M", "EM", "Fh", "muD", "DR", "fe", "yRP",
                           "ph", "yPE", "iPM", "eh", "mP", "alpha", "yEF", "LM",
                           "kd", "z", "kk", "hb", "theta", "mR", "yVE", "yEF2",
                           "sd.L", "sd.E", "sd.W", "lpost")])

### summarise and plot estimated params
svar <- "M" # select variable 
sampsvar <- samps[,svar] # pull from mcmc
summary(sampsvar) # get mean, sd, se, and quantiles for each input variable  

den <- density(sampsvar) # get AUC
densplot(sampsvar, show.obs = F,type="n") # density estimate of each variable
polygon(den, col=adjustcolor(colv,alpha=0.5),border=colv) # fill AUC 

plot(sampsvar,trace=T,density=T,col=colv) # traceplot (below) and density plot (above)
# intensive  
traceplot(sampsvar,smooth=T,type="l",lwd=0.3,xlim=c(0,length(sampsvar)),col=colv[2],xlab=paste0("Iterations"),ylab=paste0("Sampled values"),main=paste0("Sampled values over iterations for ",svar)) # iterations vs sampled valued per variable 

### save plots to PDF
traceplot <- 0 # include traceplot? intensive!!!

if(save_to_file==1){
	par(mfrow=c(1,1))
	plotlist <- list()
	pdf("mcmc_vars.pdf",onefile = T,paper="a4")
	for(i in colnames(samps)){
		par(bty="n", las = 1)
		if(traceplot==1){
		traceplot(sampsvar,smooth=T,type="l",xlim=c(0,length(sampsvar)),col=colv[2],xlab=paste0("Iterations"),ylab=paste0("Sampled values"),main=paste0("Sampled values over iterations for ",svar)) # iterations vs sampled valued per variable
		}
	svar <- i # select variable 
	sampsvar <- samps[,svar] # pull from mcmc
	den <- density(sampsvar) # get AUC
	densplot(sampsvar, show.obs = F,type="n",main=paste0("Density estimate of ",i)) # density estimate of each variable
	polygon(den, col=adjustcolor(colv,alpha=0.5),border=colv) # fill AUC 
}  
dev.off()
}

# get the best fit DEB parameters to match the data (using mcmc)
read.csv("pars.txt",header=T,sep="/",fill=T,flush=T,strip.white=T,row.names=NULL)
pars = as.vector(data.frame(samps)[max(which(data.frame(samps)$lpost >= max(data.frame(samps)$lpost) -0.001)),])
pars["Fh"] = 0.25
pars["ENV"] = 500 # Units: L
pars["r"] = 1   # Units: day-1
pars["step"] = 1  # Units: day
pars["epsilon"] = 20 # Units: L host-1, day-1 (Rounded estimate from Civitello and Rohr)
pars["sigma"] = 0.5 
pars["m_M"] = 1   # Units: day-1
pars["m_Z"] = 1   # Units: day-1
pars["M_in"] = 10
pars["K"] = 1
######

### -------------------- solve deb state for each time step --------------------
# display list of param definitions
read.csv("debfunction.txt",header=T,sep="/",fill=T,flush=T,strip.white=T,row.names=NULL)
DEB = function(step, Food, L, e, D, RH, P, RP, DAM, HAZ, iM, k, M, EM, 
               Fh, muD, DR, yRP, ph, yPE, iPM, eh, mP, alpha, yEF,
               LM, kd, z, kk, hb, theta, mR, yVE, ENV, Lp){
  # starting conditions 
      initials = c(Food=Food, L=L, e=e, D=D, RH=RH, P=P, RP=RP, DAM=DAM, HAZ=HAZ)
      # deb parameters
      parameters = c(iM, k, M, EM, Fh, muD, DR, yRP, ph, yPE, iPM,
                     eh, mP, alpha, yEF, LM, kd, z, kk, hb, theta, mR, yVE, ENV, Lp)
      # estimate starting deb conditions using fitted params by solving ode's
      ## return survival and host shell length  
      DEBstep <- lsoda(initials, c(0,step), func = "derivs", dllname = "IndividualModel_IBM2", 
                             initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=500000,
                             as.numeric(parameters),  rtol=1e-6, atol=1e-6, hmax=1)
      DEBstep[2, 2:12] # 12 = survival
}

# deb output for each timestep 
result = DEB(step=1, Food=5, L=10, e=0.9, D=as.numeric(pars["DR"]), RH=0, P=0, RP=0, DAM=0, HAZ=0, iM=pars["iM"], k=pars["k"], M=pars["M"], EM=pars["EM"], 
    Fh=pars["Fh"], muD=pars["muD"], DR=pars["DR"], yRP=pars["yRP"], ph=pars["ph"], yPE=pars["yPE"], iPM=pars["iPM"], eh=pars["eh"],
    mP=pars["mP"], alpha=pars["alpha"], yEF=pars["yEF"], LM=pars["LM"], kd=pars["kd"], z=pars["z"], kk=pars["kk"], hb=pars["hb"],
    theta=pars["theta"], mR=pars["mR"], yVE=pars["yVE"], ENV=pars["ENV"], Lp=10)

### Exposure submodel
# pass the deb state vars into infection model 
Infection = function(snail.stats, miracidia, parameters){
  # Parameters
  epsilon = as.numeric(parameters["epsilon"])
  sigma = as.numeric(parameters["sigma"])
  ENV = as.numeric(parameters["ENV"])
  m_M = as.numeric(parameters["m_M"])
  step = as.numeric(parameters["step"])
  
  # Later calculations depend on exposure probabilities
  exp.rates = epsilon/ENV*(snail.stats[,"L"]>0) # This is just to get uniform exposure rates
  sum.exp.rates = sum(exp.rates)
  
  # Probabilities for fate of miracidia
  ## Still in water
  P.left.in.water = exp(-(m_M+sum(exp.rates))*step)
  ## Infect a snail
  P.infects.this.snail = (1 - P.left.in.water)*(sigma*exp.rates/(m_M+sum.exp.rates)) 
  ## Die in water or fail to infect
  P.dead = (1 - P.left.in.water)*(m_M/(m_M+sum.exp.rates)) + sum((1 - P.left.in.water)*((1-sigma)*exp.rates/(m_M+sum.exp.rates)))
  
  prob.vector = c(P.infects.this.snail, P.left.in.water, P.dead)
  
  # Multinomial outcome from number of miracidia in env based on their survival probability
  rmultinom(n=1, size=miracidia, prob=prob.vector)
  #sum(P.left.in.water, P.invades.this.snail, P.dead)
}

### update all the snails @netlogo
update.snails = function(who, new.L, new.e, new.D, new.RH, new.P, new.RP, new.DAM, new.HAZ, new.LG){
  paste("ask snail", who, 
        "[set L", new.L,
        "set ee", new.e,
        "set D", new.D,
        "set RH", new.RH,
        "set P", new.P,        
        "set RPP", new.RP,       
        "set DAM", new.DAM,
        "set HAZ", new.HAZ,
		"set LG", new.LG, # new max length
        "]")
}

#Example update
#paste(mapply(update.snails, who=snail.stats[,"who"], new.L=L, new.e=e, new.D=D, new.RH=RH, new.P=P, new.RP=RP, new.DAM=DAM, new.HAZ=HAZ), collapse=" ")

geterrmessage() # check if there were any error messages  


################################################################################	
######################### load Netlogo ###################################### @netlogo
################################################################################
### Start netlogo and the model

# working NLStart in RStudio. works with gui=F (2018/09/24)
if(gui==0){
  	NLStart(nl.path,gui=F,nl.jarname = paste0("netlogo-",ver_nl,".jar")) # open netlogo without a gui  
	}else{
	NLStart(nl.path,nl.jarname = paste0("netlogo-",ver_nl,".jar")) # open netlogo
}
NLLoadModel(paste0(model.path,nl.model),nl.obj=NULL) # load model  
# if java.lang error persists, try copying all .jar files from the 'Java' folder where Netlogo is installed into the main Netlogo folder   	

# set type of resource dynamics @netlogo
set_resources<-function(resources){ # set resource input in env  
  if (resources == "cyclical"){
    NLCommand("set resources \"cyclical\" ") 
  }else{
    NLCommand("set resources \"event\" ") 
  }
}
set_resources(resources) # set resources: "cyclical" or "event"  @netlogo

################################################################################	
######################### start Netlogo sim ####################################
################################################################################
testrun <-0 # do a quick testrun to see plots
 ifelse(testrun==1,n.ticks<-5,n.ticks<-100)

#pdf("/Users/malishev/Documents/Emory/research/schisto_ibm/plots/resouce_dynamics_test.pdf",onefile=T,paper="a4")

#cyclical resource param space
alpha_pars <- c(0,0.25,0.5,1) # alphas
p_pars <- c(10,20,50,100) # ps
alpha_list <- list()
p_list <- list()
env_z_list <- list() # cerc list  
env_f_list <- list() # food list  
master <- list() # master list for cerc density (Env_Z) 
food_master <- list() # master list for food dynamics (Env_F) 

#define plot window
graphics.off() 
plot.matrix <- matrix(c(length(alpha_pars),length(p_pars)))
par(mfrow=plot.matrix)

for(alpha in alpha_pars){ # loop through alphas 
	for(p in p_pars){ # loop through periodicity (p)

# ---------------------------- start sim model ----------------------------
NLCommand("setup")
Env_G = integer() # make sure Env_G is an integer  
for(t in 1:n.ticks){ # @netlogo
  snail.stats = NLGetAgentSet(c("who", "L", "ee", "D", "RH", "P", "RPP", "DAM", "HAZ", "LG"), "snails")
  N.snails = length(snail.stats[,"L"])
  environment = as.numeric(NLGetAgentSet(c("F", "M", "Z", "G"), "patches"))
  
  # Infect snails
  Infection.step = as.vector(Infection(snail.stats, environment[2], pars)) # Who gets infected
  snail.stats[which(Infection.step[1:N.snails] > 0),"P"] = snail.stats[which(Infection.step[1:N.snails] > 0),"P"] + 2.85e-5 # add biomass of one miracidia
  
  # Update DEBS, HAZ=0 so survival probs are calculated for the current day
  snail.update = t(mapply(DEB, L=snail.stats[,2], e=snail.stats[,3], D=snail.stats[,4], RH=snail.stats[,5],
                          P=snail.stats[,6], RP=snail.stats[,7], DAM=snail.stats[,8], Lp=snail.stats[,10], 
                          MoreArgs = list(step=1, HAZ=0, Food=environment[1], iM=pars["iM"], k=pars["k"], M=pars["M"], EM=pars["EM"], Fh=pars["Fh"], muD=pars["muD"],
                                          DR=pars["DR"], yRP=pars["yRP"], ph=pars["ph"], yPE=pars["yPE"], iPM=pars["iPM"], eh=pars["eh"],
                                          mP=pars["mP"], alpha=pars["alpha"], yEF=pars["yEF"], LM=pars["LM"], kd=pars["kd"], z=pars["z"], 
                                          kk=pars["kk"], hb=pars["hb"], theta=pars["theta"], mR=pars["mR"], yVE=pars["yVE"], ENV=pars["ENV"])))
  
  L = snail.update[,"L"] # host structural length
  e = snail.update[,"e"] # host scaled reserve density    
  D = snail.update[,"D"] # host development 
  RH = snail.update[,"RH"] # host energy to reproduction buffer  
  DAM = snail.update[,"DAM"] # host damage from starvation  
  HAZ = snail.update[,"HAZ"] # host hazard rate from starvation   
  LG = snail.update[,"LG"] # host shell length  
  P = snail.update[,"P"] # parasite mass (sum within host)
  RP = snail.update[,"RP"] # parasite reproductive buffer  
  ingestion = sum(environment[1] - snail.update[,"Food"]) # food intake by host from environment 
  
  Eggs = floor(RH/0.015)  # Figure out how many (whole) eggs are released  
  RH = RH %% 0.015        # Remove released cercariae from the buffer
  Cercs = floor(RP/4e-5)  # Figure out how many (whole) cercs are released
  RP = RP %% 4e-5         # Remove released cercariae from buffer
 # set pop density outputs to integer to pass into Env_G and rbinom func
  Eggs = as.integer(Eggs); Cercs = as.integer(Cercs)
  
  # Update environment 
  Env_M = as.numeric(Infection.step[N.snails + 1] + pars["M_in"])
  Env_Z = as.numeric(environment[3]*exp(-pars["m_Z"]*pars["step"]) + sum(Cercs)/pars["ENV"])
  Env_G = as.integer(Env_G)
  Env_G[day] = max(0, sum(Eggs),na.rm=T)  
  Env_G[is.na(Env_G)] <- 0 # turn NAs to 0 to feed into rbinom function  
  # define food dynamics @netlogo
  if(resources == "cyclical"){
	alpha <- alpha # amplitude of resources
	p <- p  # periodicity (time range of resource cycles)  
	r_t <- pars["r"] + alpha * pars["r"] * sin(2 * pi * t/p) 
	
	#Env_F = Env_F + (alpha * sin(2 * pi * t/p)) # core resource dynamics eq
	#Env_F = max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-pars["r"]*pars["step"])) - ingestion)) # Analytical soln to logistic - ingestion (alphas [1,100])
	#Env_F = max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-pars["r"]*pars["step"])) + (alpha * sin(2 * pi * t/p)) - ingestion)) # Analytical soln to logistic - ingestion with external resource addition (alphas [1,100])
	#Env_F = max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-r_t*pars["step"])) - ingestion)) # Analytical soln to logistic - ingestion with resource growth wave (alphas [1,100])
	#Env_F = max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-pars["r"]*pars["step"])) + (alpha * pars["r"] * sin(2 * pi * t/p)) - ingestion)) # Analytical soln to logistic - ingestion with equilibrium resource growth wave (alphas [0,1])
	Env_F = max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-r_t*pars["step"])) - ingestion)) # Analytical soln to logistic - ingestion with equilibrium resource growth wave (alphas [0,1])
	} # end food dynamics  

  # Command back to NL @netlogo
  NLCommand("ask patch 0 0 [set F", Env_F, "set M", Env_M, "set Z", Env_Z, "set G", Env_G[day], "]")
  snail.commands = paste(mapply(update.snails, who=snail.stats[,"who"], new.L=L, new.e=e, new.D=D, new.RH=RH, new.P=P, new.RP=RP, new.DAM=DAM, new.HAZ=HAZ, new.LG=LG), collapse=" ")
  NLCommand(snail.commands) 
  if(day > 10){
	 NLCommand("create-snails ", rbinom(n=1, size=Env_G[day - 10], prob=0.5), "[set L 0.75 set ee 0.9 set D 0 set RH 0 set P 0 set RPP 0 set DAM 0 set HAZ 0 set LG 0.75]")
	} # end create snails
  NLCommand("go")
  #cs[t] <- rbinom(n=1, size=Env_G[day - 10], prob=0.5) # list to check 'create snails' output doesn't produce NAs
  day = day + 1 
 if(testrun==1){
	  env_z_list[t] <- Env_F + p # use to test plot outputs quickly (plots food + p value to show amplitude)  
	}else{
		env_z_list[t] <- Env_Z # get cercariae density 
		env_f_list[t] <- Env_F # get food growth  
		} # end testrun
} # end p_pars list
  env_z_list<-as.numeric(env_z_list)  
  env_f_list<-as.numeric(env_f_list)
  master[[length(master)+1]] <- env_z_list # store cerc in master list
  food_master[[length(food_master)+1]] <- env_f_list # store food in master list 
  # plot outputs 
  plot(env_z_list,type="l",las=1,bty="n",ylim=c(0,do.call(max,master)),col=round(do.call(max,master)),
	main=paste0("amplitude = ",alpha, "; periodicity = ", p),ylab="Cercariae density",xlab="Days") 
  text(which(env_z_list==max(env_z_list)),max(env_z_list),paste0("a= ",alpha," \n p= ",p),#col=max(env_z_list),
  )
#dev.off()
	} # ------------------ end p_pars 
} # --------------------------------- end alphas



ee<- env_f_list
plot(ee,type="l",las=1,bty="n")
abline(h=c(1),lty=2)

for(m in food_master[1:18]){
	plot(m,type="l",las=1,bty="n",ylim=c(0,do.call(max,master)),col=round(max(m)),
	#main = bquote("amplitude = " ~ .(alpha_pars[i])))
	main=paste0("FOOD amplitude = ",alpha, "; periodicity = ", p),xlab="Time",axis=T)		
	text(which(m==max(m)),max(m),paste0("a= ",alpha," \np= ",p)#col=round(max(m)),
	)
	#points(f,type="l",las=1,bty="n",ylim=c(0,do.call(max,food_master)),col=round(max(f)),add=T)
	}# master  



# plot master 
# define plot window
graphics.off() 
plot.matrix <- matrix(c(length(alpha_pars),length(p_pars)))
par(mfrow=plot.matrix)


# i <- 1 #cycle through alpha_pars values. uncheck bquote in below plot code 
for(m in food_master){
	plot(m,type="l",las=1,bty="n",ylim=c(0,do.call(max,master)),col=round(max(m)),
	#main = bquote("amplitude = " ~ .(alpha_pars[i])))
	main=paste0("FOOD amplitude = ",alpha, "; periodicity = ", p),xlab="Time",axis=T)		
	text(which(m==max(m)),max(m),paste0("a= ",alpha," \np= ",p)#col=round(max(m)),
	)
	#points(f,type="l",las=1,bty="n",ylim=c(0,do.call(max,food_master)),col=round(max(f)),add=T)
	}# master  

# plot master with ggplot 
y_m <- melt(master);y_m
require(reshape2);require(ggplot2); require(ggthemes)
ggplot() +
  geom_line(data = y_m, aes(x = rep.int(1:n.ticks,max(L1)) , y = value, group = L1,
	colour=factor(L1)),
	linetype=y_m$L1) +
  theme_tufte()
+ # geom_text(x=,y=,label = max(value),check_overlap = TUE)

### plot param space  
packages <- c("Interpol.T","lubridate","ggExtra","tidyr","ggthemes","ggplot2")
install.packages(packages,dependencies = T)

# get tbl_df tibble of (ORIGINAL)
#   stationid   day  hour month  year  temp
 #1 T0001         1     0 Jan    2004  -1.7
# 2 T0001         1     1 Jan    2004  -1.8
# 3 T0001         1     2 Jan    2004  -1.8 

#   Station  alpha  p     r   env_list  density
#   <fct>     <int> <int> <ord> <dbl> <dbl>
 #1 T0001      0     10  0.1   cerc  -1.7
# 2 T0001      0.25  20  0.25  food  -1.8
# 3 T0001      0.5   50  1     infec -1.8

master <- tbl_df(master)

p <-ggplot(master,aes(day,hour,fill=temp))+
  geom_tile(color= "white",size=0.1) +
  scale_fill_viridis(name="Hrly Temps C",option ="magma")
p <-p + facet_grid(year~month)
p <-p + scale_y_continuous(trans = "reverse", breaks = unique(df$hour))
p <-p + scale_x_continuous(breaks =c(1,10,20,31)) 
p <- p + theme_tufte() + 
  theme(legend.title=element_text(size=8)) +
  labs(title= paste("Hourly Temps - Station",statno), x="Day", y="Hour Commencing") +
  theme(legend.position = "bottom")
p




# NLQuit()



# ----------------------------------- END SIMULATION ------------------------------------------------- 
# --------------------------------------------------------------------------------------------------- 
# ---------------------------------------------------------------------------------------------------  

### ---------------------------- Netlogo diagnostics attempted ---------------------------------     
# 27-11-18

# TS step for rbinom producing NAs (https://stackoverflow.com/questions/13264001/r-rbinom-na-and-matrices-how-to-ignore-na-yet-retain-them)
NLCommand("create-snails ", rbinom(n=1, size=as.vector(which(!is.na(Env_G)))[day - 10], prob=0.5), "[set L 0.75 set ee 0.9 set D 0 set RH 0 set P 0 set RPP 0 set DAM 0 set HAZ 0 set LG 0.75]")

#25-9-18
# 1. have tried this
Sys.setenv(NOAWT=1) 
Sys.unsetenv("NOAWT") 

# 2. have attempted this with NL v. 5.3.0 and 6.0.4
ver_nl <-"6.0.4" # type in Netlogo version
nl.path <- "/Users/malishev/Documents/Melbourne Uni/Programs/" ; nl.path # set path to Netlogo program location

nl.path <- paste0(nl.path,"NetLogo ",ver_nl,"/"); nl.path # opt 1
nl.path <- paste0(nl.path,"NetLogo ",ver_nl,"/Java/"); nl.path # opt 2 
nl.path <- paste0(nl.path,"NetLogo ",ver_nl,"/Java"); nl.path # opt 3 

model.path <- paste0(wd,"/"); model.path # set path to Netlogo model  
nl.model <- "DEB_INF_IBM_almost_working2.nlogo" # name of Netlogo model
NLStart(nl.path,nl.jarname = paste0("netlogo-",ver_nl,".jar")) # open netlogo

### ------------------------------- Alternative food dynamics ------------------
	### periodic food dynamics 
	# eq 1a in Abrams2004)
	# Env_F = I in Abrams2004 
	# params to change?
	# r = Q (variation in resource growth)
	# gamma (sin * gamma) = resource amplitude 	
	### Env_F <- Env_F * (1 + sin(2 * pi * pars["step"] / pars["r"])) - ingestion

	### overdamped food cycle (from Nisbet et al. 1976)
	#pars["r"] * pars["step"] > pi/2 # oscillatory food cycle 
	# undefined params
	## a 
	## f
	## phi
	## b
	### pars["r"] <- pars["r"]*(1 + a * cos(2 * pi * f * pars["step"] + phi)) # eq 2 from Nisbet et al. 1976 	
	### pars["K"] <- pars["K"]*(1 + b * cos * 2 * pi * f * pars["step"])

# plotting notes

# insert posthoc titles into plot 
# https://stackoverflow.com/questions/21333083/how-to-use-an-atomic-vector-as-a-string-for-a-graph-title-in-r


# ----- resource dynamics test
#graphics.off()
par(mfrow=c(1,3));lplot <- function(...) plot(..., type="l",ylim=c(-1,1))

RR1<-list();RR2<-list();RR3<-list()
alpha<-0.8
p<-10
pars["r"] <- 1
ingestion <-0

for(t in 1:100){
	  r_t <- pars["r"] + alpha * pars["r"] * sin(2 * pi * t/p) 
	# non-cyclical
	RR1[t]<-max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-pars["r"]*pars["step"])) - ingestion)) # Analytical soln to logistic - ingestion (alphas [1,100])
	# cyclical
	RR2[t]<-max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-r_t*pars["step"])) - ingestion))
	# external addition of cyclical 
	RR3[t]<-max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-pars["r"]*pars["step"])) + r_t - ingestion)) # Analytical soln to logistic - ingestion with equilibrium resource growth wave (alphas [0,1])
	lplot(as.numeric(RR1)); lplot(as.numeric(RR2)); lplot(as.numeric(RR3))
	}
# --------- end test
