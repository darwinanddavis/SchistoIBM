## DEB IBM 
# check bottom of page for diagnostics for running netlogo in rstudio

# version 1.3 (post sicb ver with molluscicide)
# DEB_INF_GUTS_IBM_1.1.nlogo
# IndividualModel_IBM3.c (generates .o and .so files)  
# FullStarve_shrink_dilute_damage3.Rda

# version 1.2 (25-1-19)
# DEB_INF_GUTS_IBM_1.1.nlogo
# IndividualModel_IBM3.c (generates .o and .so files)  
# FullStarve_shrink_dilute_damage3.Rda

# version 1.1
# DEB_INF_GUTS_IBM_1.1.nlogo
# IndividualModel_IBM2.c
# ILL_shrink_damageA5.Rda

# 6-6-19
# added rep loop in model input for me_days and me_im 

# 31-5-19
# updated model input for detr_impact and me_days 
# changed detritus impact to 50%
# final nih sims 

# 30-5-19
# detr doesn't rebound after detr_impact

# 29-5-19
# removed Env_G[day] <- max(0, sum(Eggs),na.rm=T) after killing off eggs in env from me event
# added me_days for choosing independent me_days (multiple me events per sim)
# added detr_impact and detr_impact_days 

# 25-5-19 (v. 1.4)
# molluscicide event and impact
# set hb to hb_pars[1] (0.001) on non snail control days 
# added loop to automate save of sim outputs/plots 

# 24-5-19
# added infected host length, summed host eggs, and mean host eggs

# 14-5-19
# mollusciciding kills eggs for all previous 10 days based on me_im [prob = exp(-me_im)]
# fixed empty agent set error with Env_G[day] <- 0 for no hosts (error: Error in NLCommand("ask patch 0 0 [set F", Env_F, "set M", Env_M, "set Z", )
# turned NULLs into NAs to get numeric values for hl_list and pmass_list (14-5-19. error: cannot coerce double)

# 13-5-19
# initial me and me_im sims 

# 2-5-19
# added rep number for sicb fig 3
# made sim input section user-defined 

# 1-5-19
# fixed create-snails error for me days 

# 22-4-19 (v. 1.3)
# NL agentset error (v. 1.3)
## if snails > 0
##   also kill eggs with me event
## if no snails
##   no miracidia loss due to infection
##   No cerceriae shedding
##   no snail eggs
##   food grows without feeding
# added molluscicide impact (me_im and me_im_pars)
# added legend to plots (sim_type). need to fix host biomass output

# 18-4-19 (v. 1.2)
# fixed snail control 

# 15-4-19
# added chi param and host biomass, mean host length, and summed parasite mass to master output  

# 9-4-19
# added molluscicide events
# changed yRP to 0.824 to match 2 hour hangover period  

# 2-4-19 
# changed hb pars to match sicb paper
# added file name (global_output_fh) to before sim

# 25-3-19
# added 0.04 to hb param range 

# 5-3-19
# updated hb rates for sicb paper 

# 11-2-19
# added harwell script

# 4-2-19
# added multi-panel plots for sim outputs 

# 25-1-19 (v.1.2)
# MCMC DEB params (full starve model Rda)
# new detritus input
# new food growth in C file
# new mcmc deb estimates

# 23-1-19
# fixed individual snail deb update when using detritus: Food=environment[1]*(snail.update[,2]^2)/sum(snail.update[,2]^2), 
# added detritus supply as food source  
# added detritus param: pars["d_Z"] = 0.2 # detritus mg L-1 day-1
# removed resource cycle toggle from within sim loop to define param section (now alpha and rho of 0 = no cycle)

# 22-1-19
# added host length and parasite biomass to sim outputs 

# 15-1-19
# set LHS parameter space 

# 31-12-18
# added molluscide events (me_pars) for 95% host adult and eggs in env mortality rate (per day)

# 21-12-18
# added adult, juv, and infected host pop that sheds to sim results output
# fixed install packages section for windows 

# 20-12-18
# added host length and parasite biomass to model outputs 
# removed .so .o and .dll files from github and added to .gitignore and sensitive files dir 

#18-12-18
# changed p to rho
# added rgrowth (pars rg) (rg)
# added master lists for total and infected hosts 
# changed resource wave eq to r = r + alpha * r * sin(2 * pi * t/rho)  

#16-12-18
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

###### TO DO ######

# for new snail deb feeding 
## in new C file
### redefine dFdf (dfdt = -iM * f * sum(L^2) + rF(1-K/F) + Det)  


# LHC sampling of r, alpha, and rho for different host sizes (for cyclical food)  
# get ratio of infected hosts and shedding infected hosts     

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

#############################################################################################
#################################### Mac OSX test ###########################################
#############################################################################################

# Search "@netlogo" for netlogo code in file  

### Files required  
# "DEB_IBM.R"
# "DEB_INF_GUTS_IBM.nlogo"
# "FullStarve_shrink_dilute_damage3.Rda"
# "IndividualModel_IBM3.c"
# "IndividualModel_IBM3.so" # Mac OSX. generated from C
# "IndividualModel_IBM3.o" # Mac OSX. generated from C  
# "IndividualModel_IBM.dll" # Windows. generated from C  

test.java <- 0 # 1 = run java diagnostics  

# run java test  
install.packages("RCurl")
if(test.java==1){
  require(RCurl)
  script <- getURL("https://raw.githubusercontent.com/darwinanddavis/SchistoIBM/master/mac/java_test.R", ssl.verifypeer = FALSE)
  eval(parse(text = script))
  # check rJava version  
  .jinit()
  .jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
  # get latest Java/Oracle version: https://www.oracle.com/technetwork/java/javase/downloads/index-jsp-138363.html  
  
}
# :three: [GCC compiler in R (unconfirmed)](https://stackoverflow.com/questions/1616983/building-r-packages-using-alternate-gcc)
# [Running Netlogo 6.0.+](https://github.com/NetLogo/NetLogo/issues/1282)

#################################  Running NetLogo in Mac ##################################
# if using Mac OSX El Capitan+ and not already in JQR, download and open JGR 
mac <- 0

if(mac==1){
  install.packages('JGR',,'http://www.rforge.net/')
  library(JGR)
  JGR::JGR()
}

#############################################################################################
############################# Windows or JGR onwards ########################################
#############################################################################################

####################################  set user inputs ####################################### 
# isolate sensitive data:
# "FullStarve_shrink_dilute_damage3.Rda"

# set user outputs
snab <- 1 # 1 = use remote access (snab comp), 0 = run model on your comp 
mac <- 0 # mac or windows system? 1 = mac, 0 = windows 
gui <- 0 # display the gui? 1 = yes, 0 = no
pck <- 0 # if not already, install rnetlogo and rjava from source? 1 = yes, 0 = already installed 
save_to_file <- 0 # 1 = save simulation outputs to local dir, 0 = plot in current R session
mcmcplot <- 0 # 1 = save mcmc plots to dir
traceplot <- 0 # 1 = include traceplot in mcmc plots? intensive!!!

# set dir paths (for "/" for both Windows and Mac)
if(snab==1){
  # set dir paths (for "/" for both Windows and Mac)
  wd <- "R:/CivitelloLab/matt/schisto_ibm" # set working directory  
  ver_nl <-"6.0.4"# type in Netlogo version. found in local dir. 
  ver_gcc <-"4.6.3" # NULL # type in gcc version (if known). leave as "NULL" if unknown   
  nl.path <- "C:/Program Files" # set path to Netlogo program location
}else{
  wd <- "/Users/malishev/Documents/Emory/research/schisto_ibm/SchistoIBM" # set working directory  
  ver_nl <-"6.0.4" # type in Netlogo version. found in local dir. 
  ver_gcc <-"4.6.3" # NULL # type in gcc version (if known). leave as "NULL" if unknown   
  nl.path <- "/Users/malishev/Documents/Melbourne Uni/Programs" # set path to Netlogo program location
}

# define starting conditions for simulation model @netlogo
n.ticks <- 150 # set number of days to simulate
day <- 1 # number of days to run simulation  
resources <- "event" # set resources: "cyclical" or "event"
resource_type <- "detritus" # set resource type as "algae" or "detritus"

####################################  set model paths #######################################
# load files
deb_samps <- "FullStarve_shrink_dilute_damage3.Rda"
deb_compile <- "IndividualModel_IBM3"

setwd(wd)
nl.model <- list.files(pattern="*.nlogo") ;nl.model # Netlogo model
if(mac==1){
  nl.path <- paste0(nl.path,"/NetLogo ",ver_nl,"/Java/"); cat("Mac path:",nl.path)
}else{
  nl.path <- paste0(nl.path,"/NetLogo ",ver_nl,"/app/"); cat("Windows path:",nl.path)
}
model.path <- paste0(wd,"/"); model.path # set path to Netlogo model   

####################################  load packages #######################################
# if already loaded, uninstall RNetlogo and rJava
if(pck==1){
  p<-c("rJava", "RNetLogo"); remove.packages(p)
  # then install rJava and RNetLogo from source
  if(mac==1){
    install.packages("rJava", repos = "https://cran.r-project.org/", type="source"); library(rJava)
    install.packages("RNetLogo", repos = "https://cran.r-project.org/", type="source"); library(RNetLogo)
  }
}
require(RNetLogo)
require(rJava)

# check pck versions
installed.packages()["RNetLogo","Version"] 
installed.packages()["rJava","Version"]

# install relevant packages   
packages <- c("Matrix","deSolve","mvtnorm","LaplacesDemon","coda","adaptMCMC","sp","RNetLogo","ggplot2","RCurl","RColorBrewer","Interpol.T","lubridate","ggExtra","tidyr","ggthemes","reshape2","pse","sensitivity","beepr")  
if(require(packages)){
  install.packages(packages,dependencies = T)
}
# load annoying packages manually because they're stubborn 
if(mac==0){
  install.packages("RNetLogo")
  install.packages("RCurl")
  install.packages("Interpol.T")
  install.packages("lubridate")
  install.packages("tidyr")
  install.packages("ggthemes")
  install.packages("ggExtra")
  install.packages("beepr")
}
ppp <- lapply(packages,require,character.only=T)
if(any(ppp==F)){cbind(packages,ppp);cat("\n\n\n ---> Check packages are loaded properly <--- \n\n\n")}

cs <- list() # diagnostics list for checking NAs in create snails command  

# load plot function 
script <- getURL("https://raw.githubusercontent.com/darwinanddavis/plot_it/master/plot_it.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
display.brewer.all()
# Set global plotting parameters
cat("plot_it( \n0 for presentation, 1 for manuscript, \nset colour for background, \nset colour palette 1. use 'display.brewer.all()', \nset colour palette 2. use 'display.brewer.all()', \nset alpha for colour transperancy, \nset font style \n)")
plot_it(0,"blue","YlOrRd","Greens",1,"mono") # set plot function params       
plot_it_gg("white") # same as above for ggplot     

# load harwell script
script <- getURL("https://raw.githubusercontent.com/darwinanddavis/harwell/master/harwell.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
################################  compile packages and load files ###################################

### Install rtools and gcc for using C code and coda package 
#### https://cran.r-project.org/bin/macosx/tools/

# define paths for gcc compiler 
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
  # dyn.unload("IndividualModel_IBM3.so") # unLoad .so (Mac OSX
  system(paste0("R CMD SHLIB ",deb_compile,".c")) # generates .o and .so files 
  dyn.load(paste0(deb_compile,".so")) # Load .so (Mac OSX)
}else{
  # compile model from C definition
  # dyn.unload(paste0(deb_compile,".dll")) # unload dll (Windows only)
  system(paste0("R CMD SHLIB ",deb_compile,".c"))
  dyn.load(paste0(deb_compile,".dll"))# Load dll (Windows only)
}

####################################  load deb params #######################################

# load DEB starvation model parameters and create mcmc (and convert mcmc chain to coda format)
samps = readRDS(deb_samps)
samps <- as.mcmc(samps[, c("iM", "k", "M", "EM", "Fh", "muD", "DR", "fe", "yRP",
                           "ph", "yPE", "iPM", "eh", "mP", "alpha", "yEF", "LM",
                           "kd", "z", "kk", "hb", "theta", "mR", "yVE", "yEF2",
                           "sd.LI1", "sd.LU1", "sd.EI1", "sd.EU1", "sd.W1",  "sd.LI2",
                           "sd.LU2", "sd.EI2", "sd.EU2", "sd.W2", "gammaH", "gammaP", "lpost")])


# ---------- summarise and plot estimated params
svar <- "M" # select variable 
sampsvar <- samps[,svar] # pull from mcmc
summary(sampsvar) # get mean, sd, se, and quantiles for each input variable  

den <- density(sampsvar) # get AUC
densplot(sampsvar, show.obs = F,type="n") # density estimate of each variable
polygon(den, col=adjustcolor(colv,alpha=0.5),border=colv) # fill AUC 

plot(sampsvar,trace=T,density=T,col=colv) # traceplot (below) and density plot (above)
# intensive  
traceplot(sampsvar,smooth=T,type="l",lwd=0.3,xlim=c(0,length(sampsvar)),col=colv[2],xlab=paste0("Iterations"),ylab=paste0("Sampled values"),main=paste0("Sampled values over iterations for ",svar)) # iterations vs sampled valued per variable 

if(mcmcplot==1){
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
} # end mcmcplot
# ----------

# get the best fit DEB parameters to match the data (using mcmc)
read.csv("pars.txt",header=T,sep="/",fill=T,flush=T,strip.white=T,row.names=NULL)
pars = as.vector(data.frame(samps)[max(which(data.frame(samps)$lpost >= max(data.frame(samps)$lpost) -0.001)),])
pars["Fh"] = 2 # f_scaled (for v.1.1)
pars["ENV"] = 500 # Units: L
pars["r"] = 1   # Units: day-1
pars["step"] = 1  # Units: day
pars["epsilon"] = 20 # Units: L host-1, day-1 (Rounded estimate from Civitello and Rohr)
pars["sigma"] = 0.5 
pars["m_M"] = 1   # Units: day-1
pars["m_Z"] = 1   # Units: day-1
pars["M_in"] = 10
pars["K"] = 5
pars["Det"] = 0.1 # Units mg C/L-1 d-1 (detritus) 
pars["yRP"] = 0.824 # new  0.0471 * (1/0.4) * 7 [expected shedding output per week * (snails shed 40% of their total cercariae during 9-11 AM.) * 7 days] 

####################################  solve deb eqs #######################################
# display list of param definitions
read.csv("debfunction.txt",header=T,sep="/",fill=T,flush=T,strip.white=T,row.names=NULL)
DEB = function(step, Food, L, e, D, RH, P, RP, DAM, HAZ, iM, k, M, EM, 
               Fh, muD, DR, yRP, ph, yPE, iPM, eh, mP, alpha, yEF,
               LM, kd, z, kk, hb, theta, mR, yVE, ENV, Lp, SAtotal, r, K, Det){
  # starting conditions 
  initials = c(Food=Food, L=L, e=e, D=D, RH=RH, P=P, RP=RP, DAM=DAM, HAZ=HAZ)
  # deb parameters
  parameters = c(iM, k, M, EM, Fh, muD, DR, yRP, ph, yPE, iPM,
                 eh, mP, alpha, yEF, LM, kd, z, kk, hb, theta, mR, yVE, ENV, Lp, SAtotal, r, K, Det)
  # estimate starting deb conditions using fitted params by solving ode's
  ## return survival and host shell length  
  DEBstep <- lsoda(initials, c(0,step), func = "derivs", dllname = deb_compile, 
                   initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=500000,
                   as.numeric(parameters),  rtol=1e-6, atol=1e-6, hmax=1)
  DEBstep[2, 2:12] # 12 = survival
} # end deb model

### deb output for each timestep 
result = DEB(step=1, Food=5, L=10, e=0.9, D=as.numeric(pars["DR"]), RH=0, P=0, RP=0, DAM=0, HAZ=0, iM=pars["iM"], k=pars["k"], M=pars["M"], EM=pars["EM"], 
             Fh=pars["Fh"], muD=pars["muD"], DR=pars["DR"], yRP=pars["yRP"], ph=pars["ph"], yPE=pars["yPE"], iPM=pars["iPM"], eh=pars["eh"],
             mP=pars["mP"], alpha=pars["alpha"], yEF=pars["yEF"], LM=pars["LM"], kd=pars["kd"], z=pars["z"], kk=pars["kk"], hb=pars["hb"],
             theta=pars["theta"], mR=pars["mR"], yVE=pars["yVE"], ENV=pars["ENV"], Lp=10,SAtotal=7007.822, r=pars["r"], K=pars["K"], Det=pars["Det"]) 

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
} # end infection model 

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
} # end host update

#Example update
#paste(mapply(update.snails, who=snail.stats[,"who"], new.L=L, new.e=e, new.D=D, new.RH=RH, new.P=P, new.RP=RP, new.DAM=DAM, new.HAZ=HAZ), collapse=" ")

geterrmessage() # check if there were any error messages  

###########################################################################################
####################################  load netlogo ######################################## 
###########################################################################################
# @netlogo

# working NLStart in RStudio. works with gui=F (2018/09/24)
if(gui==0){
  NLStart(nl.path,gui=F,nl.jarname = paste0("netlogo-",ver_nl,".jar")) # open netlogo without a gui  
}else{
  NLStart(nl.path,nl.jarname = paste0("netlogo-",ver_nl,".jar")) # open netlogo
}

NLLoadModel(paste0(model.path,nl.model),nl.obj=NULL) # load model  
# if java.lang error persists on Mac, try copying all .jar files from the 'Java' folder where Netlogo is installed into the main Netlogo folder   	

################################################################################################
####################################  start netlogo sim ######################################## 
################################################################################################

# WHO sims
# every 30 days 
# 0.99, 0.95, 0.90, 0.75, 0.5
# 3 reps

#grant proposal sims
# rg = 0.25, det = 0.1 
# infected length, egg density, mean eggs 


# me sims (14-5-19)
# season = 150 days 
# every 2 weeks at different me_im (4 sims)  # me_event 1
# every month at different me_im (4 sims)   # me_event 2
# every 2 months at different me_im (4 sims)   # me_event 3
# 0.25 for algae and detritus 
# 4 reps each 

################################# 14-15-19 sims

# need to run me_event 1,2,3 and me_im_event 1,2,3,4

me_event <- 1 # 1 = on day 30, every month, skip a month 
me_im_event <- 1 # 1 = 0.69, 2 = 1.39, 3 = 2.30, 4 = 3.0, 5 = 4.60 for me_im
resource_type="algae" # detritus # set resource type
detr_impact = 0 # run detritus impact?
snail_control_days = 1 # run multiple molluscicide days per sim
rep_num <- 5 # number of replications



# start save multiple sims to dir ---------------------------------------
for(rn in 1:rep_num){  
  # # for algae and det 
  # for(rn in 1:3){ # reps
  for(me_im_event in 1:5){
    for(me_event in 1:3){ # run all me event scenarios and save to file.
      
      resources="event" # set resource cycles
      snail_control = 1 # run molluscicide sims?
      replication = 0
      sicb_sims = 0 # run sicb sims? # 2-5-19
      fig2 <- 0 # run fig 2 in sicb? 0 = fig 3
      
      # set type of resource input @netlogo
      set_resource_type<-function(resource_type){ # set resource input in env  
        if(resource_type == "detritus"){NLCommand("set resource_type \"detritus\" ")}else{NLCommand("set resource_type \"algae\" ")}}
      set_resource_type(resource_type) # set resource type: "detritus" or "algae"  @netlogo
      
      # set type of resource dynamics @netlogo
      set_resources<-function(resources){ # set resource input in env  
        if (resources == "cyclical"){NLCommand("set resources \"cyclical\" ")}else{NLCommand("set resources \"event\" ")}}
      set_resources(resources) # set resources: "cyclical" or "event"  @netlogo
      cat("\nResource type = ",resource_type,"\nResources = ",resources)
      
      # OG scenario
      # Fh = c(0.5,1,2,5)
      # K = c(1,2,5,10) 
      
      if(save_to_file==1){pdf(paste0(wd,"/master_sim.pdf"),onefile=T,paper="a4")}
      
      if(sicb_sims == 1){
        # sicb sims (31-3-19)
        # detr and algae values for reasonable cerc outputs
        # FIG 2
        # detr_pars <- c(0.1,0.25, 0.5) # detritus input (mg L^-1 day^-1)
        # rg_pars <- c(0.1,0.25,0.5)
        # hb_pars <- 0.001
        # FIG 3
        # detr_pars <- c(0.1,0.25, 0.5) # detritus input (mg L^-1 day^-1)
        # rg_pars <- c(0.1,0.25,0.5)
        # hb_pars <- 0.001
        # hb_pars = c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1)
        
        if(fig2 == 1){ # sicb fig 2
          
          # algae params
          rg_pars <- 1 # resource growth rates (rs)
          
          # detritus params
          detr_pars <- 0.5 # detritus input (mg L^-1 day^-1)
          
          # mortality params 
          hb_pars <- 0.001
          
        }else{ # sicb fig 3
          
          replication = 1 # run reps for fig 3
          
          # algae params
          rg_pars <- c(0.1,0.25,0.5) # resource growth rates (rs)
          
          # detritus params
          detr_pars <- c(0.1,0.25,0.5) # detritus input (mg L^-1 day^-1)
          
          # mortality params 
          hb_pars <- c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1) ; hb_pars
          
        }
        
        me_pars <- 0 # set molluscicide events to 0 
        if(resource_type=="detritus"){
          detr_pars <- detr_pars; alpha_pars <- 0; rho_pars <- 1; rg_pars <- 0;cat("detritus input = ",detr_pars,"\nrg = ",rg_pars,"\nhb = ", hb_pars)
        }else{detr_pars <- 0;cat("detritus input = ", detr_pars,"\nrg = ",rg_pars)}
        # set resource to cycle or be constant
        if(resource_type=="algae"){
          if(resources=="cyclical"){
            rg_pars <- rg_pars # resource growth rates (rs)
            alpha_pars <- c(0,0.25,0.5,0.75,1) # amplitude of resources (alphas)
            rho_pars <- c(1,seq(10,n.ticks,10)) # periodicity of resources (rhos)
            cat("alphas = ",alpha_pars,"\nrhos = ",rho_pars,"\nrgs = ",rg_pars)
          }else{alpha_pars <- 0; rho_pars <- 1; rg_pars <- rg_pars;cat("\nalphas = ",alpha_pars,"\nrhos = ",rho_pars,"\nrg = ",rg_pars,"\nhb = ", hb_pars)}
        }
        
      }else{
        
        # algae params
        rg_pars <- 0.25 # resource growth rates (rs)
        
        # detritus params
        detr_pars <- 0.25 # detritus input (mg L^-1 day^-1)
        
        # mortality params 
        hb_pars <- c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1); hb_pars
        
      } # end sicb sims
      
      # # define param sample space with LHS
      # require(sp)
      # require(pse)
      # require(lhs)
      # lhsmodel <- function(params){
      #   params <- factors_space[[2]]*factors_space[[3]]*factors_space[[4]]
      # }
      # factors <- c("alpha","rho","rg","me") # name of params
      # factors_space <- list(alpha_pars,rho_pars,rg_pars,me_pars)
      # q <- rep("qnorm",length(factors)) # apply the dist to be used
      # q.arg <- list(list(alpha_pars),list(rho_pars),list(rg_pars),list(me_pars)) # inputs for dist q
      # # list(list(mean=1.7, sd=0.3), list(mean=40, sd=1),list(min=1, max=50) )
      # N <- prod(as.numeric(summary(factors_space)[,1]))
      # lhs_model <- LHS(model=lhsmodel,factors=factors,N=N,q=q,q.arg=q.arg,nboot=100)
      # lhs_data <- get.data(lhs_model) # param space from LHS
      # lhs_results <- get.results(lhs_model)
      # get.N(lhs_model) # get the number of output points in hypercube
      # lhs_data
      
      Env_G = numeric() # create empty environment vector 
      day <- 1 # reset sim days 
      
      # individual outputs
      cerc_list <- list() # cercariae   
      food_list <- list() # food in env 
      juv_list <- list() # juvenile hosts
      adult_list <- list() # adult hosts 
      infec_list <- list() # infected hosts
      infec_shed_list <- list() # infected shedding hosts
      hl_list <- list() # host length
      pmass_list <- list() # parasite biomass 
      host_biomass_list <- list() # host biomass list 
      egg_list <- list() # egg density in env
      egg_mean_list <- list() # mean host eggs list
      infec_shed_length_list <- list() # length of infected shedding hosts
      
      # master outputs
      cerc_master <- list() # master list for cerc density (Env_Z) 
      food_master <- list() # master list for food dynamics (Env_F) 
      juv_master <- list() # master list for total host pop () 
      adult_master <- list() # master list for total host pop () 
      infec_master <- list() # master list for infected host pop () 
      infec_shed_master <- list() # master list for infected shedding host pop
      hl_master <- list() # master list for host length
      pmass_master <- list() # master list for parasite biomass 
      host_biomass_master <- list() # master list for host biomass
      egg_master <- list() # master list for egg density in env
      egg_mean_master <- list() # master list for mean host eggs 
      infec_shed_length_master <- list() # master list for length of infected shedding hosts
      
      # define plot window
      plot.matrix <- matrix(c(length(rg_pars),length(detr_pars)))
      par(mfrow=plot.matrix)
      
      if(snail_control == 1){
        
        me_im_pars <-  c(0.69, 1.39, 2.3, 3, 4.6) # 2.3 = 90% snail mortality from molluscicide event (per day)
        
        # if(me_event==1){me_pars = 30;me_day = "30"}
        # if(me_event==2){me_pars = seq(30,120,30); me_day = "monthly"}
        # if(me_event==3){me_pars = c(30,90,120); me_day = "skip"}
        
        me_pars = 1
        if(me_event==1){me_days = 30;me_fh = "30"}
        if(me_event==2){me_days = seq(30,120,30); me_fh = "monthly"}
        if(me_event==3){me_days = c(30,90,120); me_fh = "skip"}
        
        if(me_im_event==1){me_im_pars = me_im_pars[1]}
        if(me_im_event==2){me_im_pars = me_im_pars[2]}
        if(me_im_event==3){me_im_pars = me_im_pars[3]}
        if(me_im_event==4){me_im_pars = me_im_pars[4]}
        if(me_im_event==5){me_im_pars = me_im_pars[5]}
        
        # molluscicide events (me)
        # me_pars <- c(7, 14, 28, 56) # every x^2 days
        # me2
        # The dry seasons in Kenya are generally from mid-June to October, 
        # and from late-December to mid-March
        # me_pars <- seq(10,n.ticks,10) #  When transmission is seasonal, it is recommended to carry out regular applications of molluscicide (WHO 2017)
        # me3
        # me_pars <- XX # when infected snails appear 
        # me_pars = 1000000
        
        hb_pars <- 0.001
        if(resource_type=="algae"){detr_pars <- 0; algae <- rg_pars}else{detr_pars <- detr_pars; rg_pars <- 0}
        cat("\nalgae:",rg_pars,"\ndetritus:",detr_pars,"\nrho:",0,"\nalpha:",alpha_pars,"\nmortality (if not mollusciciding):",hb_pars,"\nmolluscicide days:",me_pars, "\nmolluscicide impact: ",me_im_pars)
      }else{
        me_pars <- n.ticks+1
        me_im_pars <- 0
        hb_pars <- hb_pars
        cat("Snail control will occur every ",max(me_pars)/length(me_pars)," days \n Mortality is ",hb_pars) 
      }
      
      cat("Mollusciding on day", me_pars)
      
      # rg_pars <- 0
      # detr_pars <- 0.25
      # hb_pars <- 0.001
      # me_im_pars
      # alpha_pars <- 0
      # rho_pars <- 1
      me_im_pars
      
      # file name to save results
      # global_output_fh = paste0(wd,"/nih_grant_sims/global_output_",resource_type,"_",resources,"_me",me_event,"_meim",me_im_event,"_rep",rn,".R"); global_output_fh
      
      if(resource_type=="algae"){detr =0; rg = 0.25;detr_impact=0;detr_impact_days=n.ticks+1}else{detr=0.25;rg=0}
      if(detr_impact==1){
        detr_impact_days = 30 # set detritus impact days 
      }else{
        detr_impact_days = n.ticks+1 # no detritus impact 
      }
      
      # if(snail_control_days==1){
      #   me_pars = n.ticks+1 # for using me_pars loop
      #   me_im = me_im_pars
      #   me_days = c(30,60,90,120)
      #   fhhh = "30120"
      # }else{
      #   hb_pars <- 0.001
      #   me_pars = me_pars # for using me_pars loop
      #   me_im = me_im_pars
      #   me_days =  n.ticks+1
      #   fhhh = "" 
      # }
      
      
      fhh = paste0(resource_type,"_meday_",me_fh,"_meim",me_im_event,"_rep",rn);fhh
      # fhh = paste0("global_output_",resource_type,"_det30_new50");fhh
      
      if(resource_type=="algae"){detr_pars <- 0; algae <- rg_pars}else{detr_pars <- detr_pars; rg_pars <- 0}
      if(length(me_days)== 1 | 2 | 3){me_days = rep(me_days,4)} # set vector same length as in sim model
      if(length(detr_impact_days)==1){detr_impact_days = rep(detr_impact_days,2)} # set vector same length as in sim model
      cat("\nalgae:",rg_pars,"\ndetritus:",detr_pars,"\nrho:",0,"\nalpha:",alpha_pars,"\nmortality (if not mollusciciding):",hb_pars,"\nmolluscicide days:",me_pars, "\nmolluscicide impact: ",me_im_pars)
      cat("\nM days = ",me_day); cat("\nM impact = ",me_im); cat("\nDet impact days = ",detr_impact_days) 
      global_output_fh = paste0(wd,"/me_detr_sims/",fhh,".R"); global_output_fh
      
      
      ####################################  start netlogo sim ######################################## 
      for(hb in hb_pars){
        for(detr in detr_pars){ # loop through detritus inputs
          for(alpha in alpha_pars){ # loop through alphas (amplitude in food cycle)
            for(rho in rho_pars){ # loop through rhos (periodicity of food cycle)
              for(rg in rg_pars){ # loop through rgs (food growth rates)
                for(me in me_pars){ # loop through mes (molluscicide events)
                  for(me_im in me_im_pars){ # loop through mes (molluscicide events)
                    
                    NLCommand("setup")
                    day <- 1 # reset days 
                    Env_G[is.na(Env_G)] <- 0 # turn NAs to 0 to feed into rbinom function 
                    for(t in 1:n.ticks){ # start nl sim  @netlogo
                      # define food dynamics for cyclical algal (logistic food growth equation) or detritus food sources
                      alpha <- alpha # amplitude of resources
                      rho <- rho  # periodicity (time range of resource cycles)  
                      rg <- rg # resource growth rate 
                      rg_t <- rg + alpha * rg * sin(2 * pi * t/rho) # equilibrium cyclical resource dynamics (19-12-18)
                      pars["r"] <- rg_t # set resource growth rate 
                      pars["Det"] <- detr # Units mg C/L-1 d-1 (detritus)
                      
                      # @detr_impact  
                      if(detr_impact==1){ # detritus impact on select days
                        if(day == detr_impact_days[1] | day == detr_impact_days[2]){ 
                          detr_95 <- detr * 0.5 # reduce detr supply rate to 25%
                          pars["Det"] <- detr_95 # set det to detr95
                          detr <- pars["Det"] # replace exisiting det for forthcoming days 
                        }
                      }
                      # else{ # no detritus impact
                      #     pars["Det"] <- detr # Units mg C/L-1 d-1 (detritus)
                      #     }
                      
                      # snail control (me events)
                      if(snail_control==1){
                        # if(day==me){ # for me_pars loop (separate sims per me_event)
                        # if(day %% me_day1==0){ # @me_day (v. 1.4) (multiple me_events per sim)
                        if(day == me_days[1] | day == me_days[2] | day == me_days[3] | day == me_days[4]){ # @me_day (v. 1.4) (select multiple me_events per sim)
                          hb <- me_im
                          cat("\n---------- DEAD SNAILS \ndays = ", day,"\nnumber of snails = ", length(snail.stats$L))
                        }else{hb <- hb_pars[1]} # hb on non-snail control days 
                      }else{hb <- hb_pars[1]} # hb when not using snail control 
                      
                      # set environment variables 
                      Env_G[is.na(Env_G)] <- 0 # turn NAs to 0 to feed into rbinom function 
                      environment = as.numeric(NLGetAgentSet(c("F", "M", "Z", "G"), "patches")) # calc food, free miracidia, cercariae released, and eggs, per patch
                      
                      # ----------------------------------- if there are hosts 
                      if(NLReport("count snails") > 0){
                        # set host variables
                        snail.stats = NLGetAgentSet(c("who", "L", "ee", "D", "RH", "P", "RPP", "DAM", "HAZ", "LG"), "snails")
                        N.snails = length(snail.stats[,"L"])
                        
                        # Infect snails
                        Infection.step = as.vector(Infection(snail.stats, environment[2], pars)) # Who gets infected
                        snail.stats[which(Infection.step[1:N.snails] > 0),"P"] = snail.stats[which(Infection.step[1:N.snails] > 0),"P"] + 2.85e-5 # add biomass of one miracidia
                        
                        # Update DEBS, HAZ=0 so survival probs are calculated for the current day
                        snail.update = t(mapply(DEB, L=snail.stats[,2], e=snail.stats[,3], D=snail.stats[,4], RH=snail.stats[,5],
                                                P=snail.stats[,6], RP=snail.stats[,7], DAM=snail.stats[,8], Lp=snail.stats[,10],# Food=environment[1]*(snail.stats[,2]^2)/sum(snail.stats[,2]^2), # update food availability per snail 
                                                MoreArgs = list(step=1, HAZ=0, Food=environment[1],# constant food available (23-1-19)
                                                                iM=pars["iM"], k=pars["k"], M=pars["M"], EM=pars["EM"], Fh=pars["Fh"], 
                                                                muD=pars["muD"],
                                                                DR=pars["DR"], yRP=pars["yRP"], ph=pars["ph"], yPE=pars["yPE"], iPM=pars["iPM"], eh=pars["eh"],
                                                                mP=pars["mP"], alpha=pars["alpha"], yEF=pars["yEF"], LM=pars["LM"], kd=pars["kd"], z=pars["z"], 
                                                                kk=pars["kk"], 
                                                                hb=hb,
                                                                theta=pars["theta"], mR=pars["mR"], yVE=pars["yVE"], SAtotal= sum(snail.stats[,2]^2), 
                                                                ENV=pars["ENV"], r=pars["r"], K=pars["K"], 
                                                                Det=pars["Det"]))) # detritus (Det) defined in C file
                        
                        L = snail.update[,"L"] # host structural length
                        e = snail.update[,"e"] # host scaled reserve density    
                        D = snail.update[,"D"] # host development 
                        RH = snail.update[,"RH"] # host energy to reproduction buffer  
                        DAM = snail.update[,"DAM"] # host damage from starvation  
                        HAZ = snail.update[,"HAZ"] # host hazard rate from starvation   
                        LG = snail.update[,"LG"] # host shell length  
                        P = snail.update[,"P"] # parasite mass (sum within host)
                        RP = snail.update[,"RP"] # parasite reproductive buffer  
                        
                        # ingestion = environment[1] - sum(snail.update[,"Food"]) # food intake by host from environment (for v.1.1)
                        chi <- pars["M"]/(1 + pars["EM"]) # length to volume conversion factor for getting biomass
                        host_biomass <- sum(chi * L^3) # get total host biomass
                        
                        Eggs = floor(RH/0.015)  # Figure out how many (whole) eggs are released  
                        # if(day==me){Eggs <- Eggs[1:round(0.1*length(Eggs))]} # kill off 90% of snail eggs in water with molluscicide event  
                        RH = RH %% 0.015        # Remove released cercariae from the buffer
                        Cercs = floor(RP/4e-5)  # Figure out how many (whole) cercs are released
                        RP = RP %% 4e-5         # Remove released cercariae from buffer
                        Eggs = as.integer(Eggs); Cercs = as.integer(Cercs) 
                        eggs_mean <- mean(Eggs[Eggs > 0]) # get mean eggs in env
                        eggs_mean <- as.numeric(eggs_mean)
                        
                        # Update environment 
                        Env_M = as.numeric(Infection.step[N.snails + 1] + pars["M_in"]) # total miracidia density 
                        Env_Z = as.numeric(environment[3]*exp(-pars["m_Z"]*pars["step"]) + sum(Cercs)/pars["ENV"]) # total cerc density
                        Env_G = as.integer(Env_G) # set pop density outputs to integer to pass into Env_G and rbinom func
                        
                        # @snail_control
                        # kill off 90% of snail eggs in water with molluscicide event 
                        # ifelse(day==me,Env_G[day] <- max(0, 0.1*sum(Eggs),na.rm=T),Env_G[day] <- max(0, sum(Eggs),na.rm=T)) # for me_pars loop (separate sims per me_event)
                        # ifelse(day %% me_day1==0,Env_G[day] <- max(0, 0.1*sum(Eggs),na.rm=T),Env_G[day] <- max(0, sum(Eggs),na.rm=T)) # @me_day (v. 1.4) (multiple me_events per sim)
                        ifelse(day == me_days[1] | day == me_days[2] | day == me_days[3] | day == me_days[4],Env_G[day] <- max(0, 0.1*sum(Eggs),na.rm=T),Env_G[day] <- max(0, sum(Eggs),na.rm=T)) # @me_day (v. 1.4)
                        
                        # Env_G[day] <- max(0, sum(Eggs),na.rm=T)
                        egg <- sum(Eggs) # get summed host eggs 
                        
                        Env_G[is.na(Env_G)] <- 0 # turn NAs to 0 to feed into rbinom function  
                        # Env_F = max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-pars["r"]*pars["step"])) - ingestion)) # Analytical soln to logistic - ingestion (alphas [1,100]) (original r growth equation)
                        # Env_F = max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-rg_t*pars["step"])) - ingestion)) # Analytical soln to logistic - ingestion with equilibrium resource growth wave (rg_t) (alphas [0,1]) (for v.1.1)     
                        # F = F * exp(- r + alpha * r * sin(2 * pi * t/rho) * s) * (1 - F/K) - f(i_{M} * sum(L^2) # v. 1.2 algae and detritus with cyclical algal growth
                        # r_t <- pars["r"] + alpha * pars["r"] * sin(2 * pi * t/rho) # equilibrium resource dynamics (static)
                        Env_F = max(0.001, snail.update[1,"Food"]) # algal or detritus food sources (for v.1.2)
                        
                        # @detr_impact  
                        if(detr_impact==1){ # detritus impact on select days 
                          if(day == detr_impact_days[1] | day == detr_impact_days[2]){
                            Env_F_95 = max(0.001, snail.update[1,"Food"] * 0.5) # reduce detritusto 25% (v.1.4)
                            Env_F = Env_F_95
                          }
                        }
                        # else{ # no detritus impact (normal conditions)
                        #   Env_F = max(0.001, snail.update[1,"Food"]) # algal or detritus food sources (for v.1.2)
                        #   }
                        
                        # Command back to NL @netlogo
                        snail.commands = paste(mapply(update.snails, who=snail.stats[,"who"], new.L=L, new.e=e, new.D=D, new.RH=RH, new.P=P, new.RP=RP, new.DAM=DAM, new.HAZ=HAZ, new.LG=LG), collapse=" ")
                        NLCommand(snail.commands) 
                        
                        hl_list[t] <- mean(L) # get average host lengths per model step 
                        pmass_list[t] <- sum(P) # get total parasite mass per model step
                        host_biomass_list <- host_biomass # get host mass 
                        egg_mean_list[t] <- eggs_mean # get mean eggs in env
                        
                      }else{ # ----------------------------------- if there are no hosts 
                        
                        Env_M = as.numeric(environment[2]*exp(-pars["m_M"]*pars["step"]) + pars["M_in"]) # total miracidia density 
                        Env_Z = as.numeric(environment[3]*exp(-pars["m_Z"]*pars["step"])) # total cerc density
                        Env_G = as.integer(Env_G) # set pop density outputs to integer to pass into Env_G and rbinom func
                        Env_G[day] <- 0
                        Env_F = ifelse(pars["Det"] == 0, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-pars["r"]*pars["step"]))), as.numeric(environment[1] + pars["Det"]))
                        
                        # @detr_impact  
                        if(detr_impact==1){ # detritus impact on select days 
                          if(day == detr_impact_days[1] | day == detr_impact_days[2]){
                            Env_F = as.numeric(environment[1] + pars["Det"]) * 0.5
                          }
                        }
                        # else{ # no detritus impact (normal conditions)
                        #     Env_F = ifelse(pars["Det"] == 0, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-pars["r"]*pars["step"]))), as.numeric(environment[1] + pars["Det"]))
                        #   }
                        
                        egg <- 0
                      } # end no hosts
                      
                      # update patch variables 
                      NLCommand("ask patch 0 0 [set F", Env_F, "set M", Env_M, "set Z", Env_Z, "set G", Env_G[day], "]")
                      
                      # kill snail eggs with molluscicide event  
                      if(day > 10){
                        if(snail_control==1){ 
                          # if(day==me){ # for me_pars loop (separate sims per me_event)
                          # if(day %% me_day1==0){ # @me_day (v. 1.4) (multiple me_events per sim)
                          if(day == me_days[1] | day == me_days[2] | day == me_days[3] | day == me_days[4]){ # @me_day (v. 1.4) (multiple me_events per sim)
                            Env_G[(day - 10):(day - 1)] <- rbinom(n=10,size=Env_G[(day - 10):(day - 1)],prob=exp(-me_im)) # mollusciciding kills eggs for all previous 10 days 
                            create_snails <- rbinom(n=1, size=Env_G[day - 10], prob=0.5) # keep normal prob for me day  
                          }else{create_snails <- rbinom(n=1, size=Env_G[day - 10], prob=0.5)}
                        }else{create_snails <- rbinom(n=1, size=Env_G[day - 10], prob=0.5)}
                        NLCommand("create-snails ", create_snails, "[set L 0.75 set ee 0.9 set D 0 set RH 0 set P 0 set RPP 0 set DAM 0 set HAZ 0 set LG 0.75]")
                      } # end create snails
                      
                      NLCommand("go") # run @netlogo sim steps
                      #cs[t] <- rbinom(n=1, size=Env_G[day - 10], prob=0.5) # list to check 'create snails' output doesn't produce NAs
                      day = day + 1 
                      Env_G[is.na(Env_G)] <- 0 # turn NAs to 0 to feed into rbinom function 
                      # results outputs
                      cerc_list[t] <- Env_Z # get cercariae density 
                      food_list[t] <- Env_F # get food growth
                      juv_list[t] <- length(which(snail.stats$RH==0)) # get juvenile hosts
                      adult_list[t] <- length(which(snail.stats$RH>0)) # get adult hosts
                      infec_list[t] <- length(which(snail.stats$P>0)) # get just infected hosts
                      infec_shed_list[t] <- length(which(snail.stats$RP>0)) # get infected hosts that are shedding
                      # get length of infected hosts
                      ish <- subset(snail.stats,RP>0) # get infected hosts
                      infec_shed_length_list[t] <- mean(ish$L) # get their mean length
                      
                      egg_list[t] <- egg # save to host egg list 
                      egg_mean_list[t] <- eggs_mean # get mean eggs in env 
                      
                      
                    } # --------------------------------------- end nl sim
                    
                    # turn NULLs into NAs to get numeric values below (14-5-19. error: cannot coerce double)
                    hl_list <- lapply(hl_list, function(x) ifelse(x == "NULL", NA, x))
                    pmass_list <- lapply(pmass_list, function(x) ifelse(x == "NULL", NA, x))
                    egg_mean_list <- lapply(egg_mean_list, function(x) ifelse(x == "NaN", NA, x))
                    infec_shed_length_list <- lapply(infec_shed_length_list, function(x) ifelse(x == "NULL", NA, x))
                    
                    # save individual outputs 
                    cerc_list <- as.numeric(cerc_list) 
                    food_list <- as.numeric(food_list)
                    juv_list <- as.numeric(juv_list)
                    adult_list <- as.numeric(adult_list)
                    infec_list <- as.numeric(infec_list)
                    infec_shed_list <- as.numeric(infec_shed_list)
                    hl_list <- as.numeric(hl_list)
                    pmass_list <- as.numeric(pmass_list)
                    host_biomass_list <- as.numeric(host_biomass_list)
                    egg_list <- as.numeric(egg_list)
                    egg_mean_list <- as.numeric(egg_mean_list)
                    infec_shed_length_list <- as.numeric(infec_shed_length_list)
                    
                    # save master outputs 
                    cerc_master[[length(cerc_master)+1]] <- cerc_list # cerc master list
                    food_master[[length(food_master)+1]] <- food_list # food master list
                    juv_master[[length(juv_master)+1]] <- juv_list # juv pop master list
                    adult_master[[length(adult_master)+1]] <- adult_list # adult pop master list
                    infec_master[[length(infec_master)+1]] <- infec_list # infected host pop master list
                    infec_shed_master[[length(infec_shed_master)+1]] <- infec_shed_list # infected shedding host pop master list
                    hl_master[[length(hl_master)+1]] <- hl_list # host length master
                    pmass_master[[length(pmass_master)+1]] <- pmass_list # parasite mass master
                    host_biomass_master[[length(host_biomass_master)+1]] <- host_biomass_list # host biomass master
                    egg_master[[length(egg_master)+1]] <- egg_list # summed egg master
                    egg_mean_master[[length(egg_mean_master)+1]] <- egg_mean_list # mean egg master
                    infec_shed_length_master[[length(infec_shed_length_master)+1]] <- infec_shed_length_list # mean infected shedding host length master
                    
                    # day_master[[length(day_master)+1]] <- day_list
                    ### plot outputs 
                    #   plot(cerc_list,type="l",las=1,bty="n",ylim=c(0,do.call(max,cerc_master)),col=round(do.call(max,cerc_master)),
                    # 	main=paste0("alpha = ",alpha, "; rho = ", rho, "; r = ", rg),ylab="Cercariae density",xlab="Days") 
                    #   paste0(expression("alpha = ",alpha, "; rho = ", rho, "; r = ", rg)) 
                    #   text(which(cerc_list==max(cerc_list)),max(cerc_list),paste0("a= ",alpha," \n p= ",rho)#,col=max(cerc_list),
                    #        )
                    #abline(h=which(cerc_list==max(cerc_list)),type=3,col=round(do.call(max,cerc_master))) # draw line at max peak
                    if(save_to_file==1){dev.off()}
                  } # ----end me_im
                } # --------------- end mes
              } # ------------------------------ end rgs
            } # --------------------------------------------- end rhos
          } # ----------------------------------------------------------- end alphas
        } # ------------------------------------------------------------------------- end detritus
      } # end hb pars
      # ####################################  end netlogo sim ######################################## 
      
      # results output 
      # save sim results to dir 
      global_output <- list(cerc_master,food_master,juv_master, adult_master,infec_master,infec_shed_master,hl_master,pmass_master,host_biomass_master, egg_master, egg_mean_master,infec_shed_length_master) 
      str(global_output)
      
      # turn NaN and NA into 0s
      global_output <- rapply(global_output, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
      global_output <- rapply(global_output, f=function(x) ifelse(is.na(x),0,x), how="replace" )
      
      # global_output <- cerc_master # save just cercs
      saveRDS(global_output,global_output_fh) # save to dir 
      cat("Output saved in ", global_output_fh)
      # read in saved sim results
      cat("order = cerc, food, juv, adult, infected, infected shedding, mean host length, mean parasite mass, summed host biomass", "mean host eggs", "infected host length")
      
    } # end me_im_event sim runs
  } # end me_im_event sims
  
} # end rep_num loop 

# ------------------------- plot individual outputs -------------------------

mm_ = readRDS(global_output_fh)
layout(matrix(c(1:16),4,4,byrow=T))
# plot master
mm <- mm_[[2]] # 1=cercs, 2=food
y_m <- melt(mm);y_m
ggplot() +
  geom_point(data = y_m, aes(x = rep.int(1:n.ticks,max(L1)) , y = value, group = L1, colour=factor(L1)), ) +
  geom_line(data = y_m, aes(x = rep.int(1:n.ticks,max(L1)) , y = value, group = L1, colour=factor(L1)), ) +
  #linetype=y_m$L1) +
  theme_tufte() 
# +  geom_text(x=,y=,label = max(value),check_overlap = TUE)

# ### for just cerc results ("global_output_algae_event_hb_sicb_cercs_rep1.R")
# for(pp in seq(10,110,8)){
#   plot(mm_[[pp]],type="l")
# }


#------------------------- plot all sim results in one window -------------------------
sim_type <- "detritus" # "detritus" # "algae" "me"

require(gridExtra)
gspl <- list()
K_pars = c(0.5, 1, 2, 3) 
Fh_pars = c(0.5, 1, 1.5, 2)
ttl_list <- c("cerc","food", "juv", "adult", "infec", "infec (shed)", "host L", "parasite mass","host biomass","sum host eggs","mean host eggs", "infected host length")
if(sim_type=="me"){legend_pars <- me_pars}; if(sim_type=="algae"){legend_pars <- rg_pars}; if(sim_type=="detritus"){legend_pars <- detr_pars};  
if(snail_control==0){me_pars <- n.ticks+1}

# choose sim to plot
global_sim_plot <- mm_
graphics.off()

for(g in 1:length(global_sim_plot)){
  par(bty="n", las = 1)
  mm <- global_sim_plot[[g]]
  y_m <- melt(mm);y_m
  gspl[[g]] <- ggplot() +
    # geom_point(data = y_m, aes(x = rep.int(1:n.ticks,max(L1)) , y = value, group = L1, colour=factor(L1)), ) +
    geom_line(data = y_m, aes(x = rep.int(1:length(mm_[[1]][[1]]),max(L1)) , y = value, group = L1, colour=factor(L1)) ) +
    # scale_color_manual(values = viridis(length(mm))) +
    #linetype=y_m$L1) +
    theme_tufte() +
    theme(legend.position = "right") +
    scale_color_manual(name=paste0(sim_type," days"), # stupid ggplot legend @ggplotlegend  # http://r-statistics.co/Complete-Ggplot2-Tutorial-Part2-Customizing-Theme-With-R-Code.html
                       labels = legend_pars,
                       values = unique(y_m$L1)) +
    geom_vline(xintercept = me_pars, linetype="dotted", # add me event lines  
               color = 1:5, size=0.5) +
    labs(title=ttl_list[g],x="",y="") +
    #labs(title=paste0("Fh",ttl_list1[g],"_K",ttl_list2[g]),x="",y="") +
    if(g==length(global_sim_plot)-1){
      # theme(legend.title=element_text(size=0.2), 
      #       legend.text=element_text(size=0.2)) +
      #   theme(legend.position = "right")
      labs(x="Time") 
    }else{
      theme(legend.position="none")
    } 
}
# +  geom_text(x=,y=,label = max(value),check_overlap = TUE)

ttl_list
do.call(grid.arrange,gspl[c(1:4,6,7,10,11,12)]) # plot in one window 

# setwd("R:/CivitelloLab/matt/schisto_ibm/nih_grant_sims")
# pdf(paste0(getwd(),"/",resource_type,"_me",me_event,"_meim",me_im_event,"_rep",rn,".pdf"),onefile=T,width=15,height=10,paper="a4r")
# do.call(grid.arrange,gspl[c(1:4,6,10,11,12)]) # plot in one window
# # gspl[7] # plot individual outputs, 1 = cercs
# dev.off()

fhh

# setwd("R:/CivitelloLab/matt/schisto_ibm/nih_grant_sims")
# pdf(paste0(getwd(),"/final/",fhh,".pdf"),onefile=T,width=15,height=10,paper="a4r")
mains <- c("cercs","food","juv","adult","","infec (shed)",rep("",3),"sum host eggs","mean host eggs","infec host length")
par(mfrow=c(3,3),las=1,bty="n")
for(i in c(1:4,6,10,11,12)){
  plot(mm_[[i]][[1]],type="l",main=mains[i]) # plot cercs and infected host length
  # abline(v=seq(me_day,150,me_day),lty=3)
  # abline(v=c(me_days),lty=3)
}
# dev.off()


# NLQuit()


#################################################################################################
##########################################  end body ############################################ 
#################################################################################################
