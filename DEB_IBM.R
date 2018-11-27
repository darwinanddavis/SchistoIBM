## DEB IBM 
# check bottom of page for diagnostics for running netlogo in rstudio

# version  

# 27-11-18
# list to check if 'create snails' command isn't producing NAs from Env_G
# set pop density outputs in NL loop to integer to pass into Env_G and rbinom func

# 23-11-18
# all user inputs at beginning of doc  

# 22-11-19 
# added debfunction.txt and pars.txt for defining params  

# 19-11-18  
# added  "DEB_INF_GUTS_IBM_1.1.nlogo" as test model    

# TO DO
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
# "FullStarve_shrink_production2.Rda"
# "IndividualModel_IBM.c"
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
# "FullStarve_Shrink_adaptMCMC_original_DAM_run2.Rda"
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

# define starting conditions for simulation model @netlogo  
resources <- "cyclical" # set resources: "cyclical" or "event"
n.ticks <- 50 # set number of simulation ticks
day <- 1 # number of days to run simulation    

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
    # dyn.unload("IndividualModel_IBM.so") # unLoad .so (Mac OSX
    system("R CMD SHLIB IndividualModel_IBM.c") # generates .o and .so files 
    dyn.load("IndividualModel_IBM.so") # Load .so (Mac OSX)
}else{
  # compile my model from C definition
  dyn.unload("IndividualModel_IBM.dll") # unload dll (Windows only)
  system("R CMD SHLIB IndividualModel_IBM.c")
  dyn.load("IndividualModel_IBM.dll") # Load dll (Windows only)
}

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# point where windows/mac merge  

# load DEB starvation model parameters and create mcmc
samps = readRDS(paste0(pp,"FullStarve_Shrink_adaptMCMC_original_DAM_run2.Rda"))
lpost = samps$log.p #
samps = convert.to.coda(samps) # convert mcmc chain to coda format
samps = cbind(samps, lpost)
samps <- as.mcmc(samps[, c("iM", "k", "M", "EM", "Fh", "muD", "DR", "fe", "yRP",
                           "ph", "yPE", "iPM", "eh", "mP", "alpha", "yEF", "LM",
                           "kd", "z", "kk", "hb", "theta", "mR", "yVE", "yEF2",
                           "sd.LI1", "sd.LU1", "sd.EI1", "sd.EU1", "sd.W1", 
                           "sd.LI2", "sd.LU2", "sd.EI2", "sd.EU2", "sd.W2",
                           "lpost")])

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

## solve deb state for each time step 
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
      DEBstep <- lsoda(initials, c(0,step), func = "derivs", dllname = "IndividualModel_IBM", 
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

# ---------------------------- start sim model ----------------------------
NLCommand("setup")
Env_G = integer() # make sure Env_G is an integer  
cs <- list()
for(t in 1:n.ticks){ # @netlogo
  snail.stats = NLGetAgentSet(c("who", "L", "ee", "D", "RH", "P", "RPP", "DAM", "HAZ","LG"), "snails")
  N.snails = length(snail.stats[,"L"])
  environment = as.numeric(NLGetAgentSet(c("F", "M", "Z", "G"), "patches"))
  
  # Infect snails
  Infection.step = as.vector(Infection(snail.stats, environment[2], pars)) # Who gets infected
  snail.stats[which(Infection.step[1:N.snails] > 0),"P"] = snail.stats[which(Infection.step[1:N.snails] > 0),"P"] + 2.85e-5
  
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
  Eggs = as.integer(Eggs); Cercs = as.integer(Cercs); Env_G = as.integer(Env_G)
  
  # Update environment
  Env_F = max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-pars["r"]*pars["step"])) - ingestion)) # Analytical soln to logistic - ingestion
  Env_M = as.numeric(Infection.step[N.snails + 1] + pars["M_in"])
  Env_Z = as.numeric(environment[3]*exp(-pars["m_Z"]*pars["step"]) + sum(Cercs)/pars["ENV"])
  Env_G[day] = max(0, sum(Eggs))  
  
  # define food dynamics @netlogo
  NLCommand(" ask patches [set F F ]") 
  
  # Command back to NL @netlogo
  NLCommand("ask patch 0 0 [set F", Env_F, "set M", Env_M, "set Z", Env_Z, "set G", Env_G[day], "]")
  snail.commands = paste(mapply(update.snails, who=snail.stats[,"who"], new.L=L, new.e=e, new.D=D, new.RH=RH, new.P=P, new.RP=RP, new.DAM=DAM, new.HAZ=HAZ, new.LG=LG), collapse=" ")
  NLCommand(snail.commands) 
  if(day > 10){
	 NLCommand("create-snails ", rbinom(n=1, size=Env_G[day - 10], prob=0.5), "[set L 0.75 set ee 0.9 set D 0 set RH 0 set P 0 set RPP 0 set DAM 0 set HAZ 0 set LG 0.75]")
}
  NLCommand("go")
  cs[t] <- rbinom(n=1, size=Env_G[day - 10], prob=0.5) # list to check 'create snails' output doesn't produce NAs
  day = day + 1
}

NLQuit()

# TS step for rbinom producing NAs 
NLCommand("create-snails ", rbinom(n=1, size=as.vector(which(!is.na(Env_G)))[day - 10], prob=0.5), "[set L 0.75 set ee 0.9 set D 0 set RH 0 set P 0 set RPP 0 set DAM 0 set HAZ 0 set LG 0.75]")

### ---------------------------- Netlogo diagnostics attempted ---------------------------------     
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

