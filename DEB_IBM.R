## DEB IBM 
# check bottom of page for diagnostics for running netlogo in rstudio

# --------------- --------------- --------------- ---------------
# --------------- --------------- --------------- ---------------

### Files required  
# "DEB_INF_GUTS_IBM.nlogo"
# "FullStarve_shrink_production2.Rda"
# "IndividualModel_IBM.c"
# "IndividualModel_IBM.so" # Mac OSX. generated from C
# "IndividualModel_IBM.o" # Mac OSX. generated from C  
# "IndividualModel_IBM.dll" # Windows. generated from C  

### Troubleshooting  
#### :one: [Installing compiler toolchain for Mac OSX](https://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-macos/)  
#### :two: if rJava error, run the following in terminal (src: https://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite)
  ##### and #http://paulklemm.com/blog/2015-02-20-run-rjava-with-rstudio-under-osx-10-dot-10/
  ##### sudo ln -s $(/usr/libexec/java_home)/jre/lib/server/libjvm.dylib /usr/local/lib 

# java troubleshooting. tests are sequential. 
# https://stackoverflow.com/questions/14915898/rnetlogo-function-nlstart-fails-to-launch-gui

test.java <- 0 # 1 = run java diagnostics  
mac <- 1 # 1 = Mac, 0 = Windows
gui <- 1 # 1 = run Netlogo with a GUI  

if(test.java==1){
  # test 1  
  # test java is working
  .jinit() 
  .jnew( "java/awt/Point", 10L, 10L )
  f <- .jnew("java/awt/Frame","Hello")
  .jcall(f,"setVisible",TRUE)
  
  # test 2
  component <- .jnull()
  component <- .jcast(component, new.class = "java/awt/Component")
  message <- .jnew("java/lang/String","This is a JOptionPane test from rJava.")
  message <- .jcast(message, new.class = "java/lang/Object")
  title <- .jnew("java/lang/String","Test")
  type <- .jnew("java/lang/Integer", as.integer(2))
  f <- .jnew("javax/swing/JOptionPane")
  .jcall(f,"showMessageDialog", component, message, title, .jsimplify(type))
  
  # test 3
  .jcall("java/lang/System", "S", "getProperty", "java.vm.version")
  .jcall("java/lang/System", "S", "getProperty", "java.vm.name")
  .jcall("java/lang/System", "S", "getProperty", "java.vm.info")
  .jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
  .jcall("java/lang/System", "S", "getProperty", "sun.arch.data.model")
  
  # test 4
  .jcall("java/lang/System", "S", "getProperty", "java.awt.headless")
  Sys.getenv("NOAWT")
}

#### :three: [GCC compiler in R (unconfirmed)](https://stackoverflow.com/questions/1616983/building-r-packages-using-alternate-gcc)
#### [Running Netlogo 6.0.+](https://github.com/NetLogo/NetLogo/issues/1282)

# Search "@netlogo" for netlogo code in file  

# ------------------------------------------------------------------------------------------------------------- 

# ------------------------------------------------------- 
# ------------ Running NetLogo in Mac -------------------

# if using Mac OSX El Capitan+ and not already in JQR, download and open JGR 
if(mac==1){
  # load JGR after downloading 
  Sys.setenv(NOAWT=1)
  install.packages("JGR")
  library(JGR)
  Sys.unsetenv("NOAWT")
  JGR() # open JGR
}

# ------------------------------------------------------------------
# ------------------------- JGR onwards ----------------------------
# ------------------------------------------------------------------

#### set wd
mac <- 1
gui <- 1
wd <- "/Users/malishev/Documents/Emory/research/schisto_ibm/DEB_IBM"
setwd(wd)

### load Netlogo model  
installed.packages()["RNetLogo","Version"] # check rnetlogo version  
installed.packages()["rJava","Version"] # check rnetlogo version  

ver <-"6.0.4" # type in Netlogo version
nl.path <- "/Users/malishev/Documents/Melbourne Uni/Programs/" ; nl.path # set path to Netlogo program location
nl.path <- paste0(nl.path,"NetLogo ",ver,"/Java/"); nl.path

model.path <- paste0(wd,"/"); model.path # set path to Netlogo model  
#nl.model <- "DEB_INF_IBM_almost_working2.nlogo" # name of Netlogo model
nl.model <- "DEB_INF_GUTS_IBM.nlogo"


# if already loaded, uninstall RNetlogo and rJava
p<-c("rJava", "RNetLogo"); remove.packages(p)

# then install rJava and RNetLogo from source
install.packages("rJava", repos = "https://cran.r-project.org/", type="source"); library(rJava)
install.packages("RNetLogo", repos = "https://cran.r-project.org/", type="source"); library(RNetLogo)

# check rJava version  
.jinit()
.jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
# get latest Java/Oracle version: https://www.oracle.com/technetwork/java/javase/downloads/index-jsp-138363.html  

### install relevant packages   
packages <- c("Matrix","deSolve","mvtnorm","LaplacesDemon","coda","adaptMCMC")

if(require(packages)){
  install.packages(packages,dependencies = T)
  require(packages)
}
lapply(packages,library,character.only=T)

### Install rtools and gcc for using C code and coda package 
#### https://cran.r-project.org/bin/macosx/tools/

if(mac==1){ #### Mac OSX
  rtools <- "/usr/local/clang6/bin"
  gcc <- "usr/local/clang6/gcc-4.6.3/bin"
  # Mac OSX
  }else{ #### Windows 
  rtools <- "C:\\Rtools\\bin"
  gcc <- "C:\\Rtools\\gcc-4.6.3\\bin"
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

# load DEB starvation model parameters and create mcmc
#samps = as.mcmc(readRDS("FullStarve_shrink_production2.Rda"))
#samps = readRDS("FullStarve_shrink_adaptMCMC_original_DAM_run2.Rda")

samps = readRDS("FullStarve_Shrink_adaptMCMC_original_DAM_run2.Rda")
lpost = samps$log.p
samps = convert.to.coda(samps)
samps = cbind(samps, lpost)
samps <- as.mcmc(samps[, c("iM", "k", "M", "EM", "Fh", "muD", "DR", "fe", "yRP",
                           "ph", "yPE", "iPM", "eh", "mP", "alpha", "yEF", "LM",
                           "kd", "z", "kk", "hb", "theta", "mR", "yVE", "yEF2",
                           "sd.LI1", "sd.LU1", "sd.EI1", "sd.EU1", "sd.W1", 
                           "sd.LI2", "sd.LU2", "sd.EI2", "sd.EU2", "sd.W2",
                           "lpost")])
summary(samps)

# get the best fit DEB parameters to match the data (using mcmc)
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

DEB = function(step, Food, L, e, D, RH, P, RP, DAM, HAZ, iM, k, M, EM, 
               Fh, muD, DR, yRP, ph, yPE, iPM, eh, mP, alpha, yEF,
               LM, kd, z, kk, hb, theta, mR, yVE, ENV, Lp){
  
      initials = c(Food=Food, L=L, e=e, D=D, RH=RH, P=P, RP=RP, DAM=DAM, HAZ=HAZ)
      parameters = c(iM, k, M, EM, Fh, muD, DR, yRP, ph, yPE, iPM,
                     eh, mP, alpha, yEF, LM, kd, z, kk, hb, theta, mR, yVE, ENV, Lp)
    
      DEBstep <- lsoda(initials, c(0,step), func = "derivs", dllname = "IndividualModel_IBM", 
                             initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=500000,
                             as.numeric(parameters),  rtol=1e-6, atol=1e-6, hmax=1)
      DEBstep[2, 2:12] # 12 = survival
  
}

result = DEB(step=1, Food=5, L=10, e=0.9, D=as.numeric(pars["DR"]), RH=0, P=0, RP=0, DAM=0, HAZ=0, iM=pars["iM"], k=pars["k"], M=pars["M"], EM=pars["EM"], 
    Fh=pars["Fh"], muD=pars["muD"], DR=pars["DR"], yRP=pars["yRP"], ph=pars["ph"], yPE=pars["yPE"], iPM=pars["iPM"], eh=pars["eh"],
    mP=pars["mP"], alpha=pars["alpha"], yEF=pars["yEF"], LM=pars["LM"], kd=pars["kd"], z=pars["z"], kk=pars["kk"], hb=pars["hb"],
    theta=pars["theta"], mR=pars["mR"], yVE=pars["yVE"], ENV=pars["ENV"], Lp=10)

### Exposure submodel
Infection = function(snail.stats, miracidia, parameters){
  # Parameters
  epsilon = as.numeric(parameters["epsilon"])
  sigma = as.numeric(parameters["sigma"])
  ENV = as.numeric(parameters["ENV"])
  m_M = as.numeric(parameters["m_M"])
  step = as.numeric(parameters["step"])
  
  # Later calculations depend on exposure probabilities
  exp.rates =epsilon/ENV*(snail.stats[,"L"]>0) # This is just to get uniform exposure rates
  sum.exp.rates = sum(exp.rates)
  
  # Probabilities for fate of miracidia
  P.left.in.water = exp(-(m_M+sum(exp.rates))*step)                             # Still in water
  P.infects.this.snail = (1 - P.left.in.water)*(sigma*exp.rates/(m_M+sum.exp.rates))  # Infect a snail
  # Die in water or fail to infect
  P.dead = (1 - P.left.in.water)*(m_M/(m_M+sum.exp.rates)) + sum((1 - P.left.in.water)*((1-sigma)*exp.rates/(m_M+sum.exp.rates)))                      # die
  
  prob.vector = c(P.infects.this.snail, P.left.in.water, P.dead)
  
  # Multinomial outcome
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
  	NLStart(nl.path,gui=F,nl.jarname = paste0("netlogo-",ver,".jar")) # open netlogo without a gui  
	}else{
	NLStart(nl.path,nl.jarname = paste0("netlogo-",ver,".jar")) # open netlogo
}
NLLoadModel(paste0(model.path,nl.model),nl.obj=NULL) # load model  
# if java.lang error persists, try copying all .jar files from the 'Java' folder where Netlogo is installed into the main Netlogo folder   	
	
################################################################################	
######################### start Netlogo sim ####################################
################################################################################
NLCommand("setup")
day = 1
n.ticks = 100
Env_G = numeric()
for(t in 1:n.ticks){
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
  
  L = snail.update[,"L"]
  e = snail.update[,"e"]
  D = snail.update[,"D"]
  RH = snail.update[,"RH"]
  P = snail.update[,"P"]
  RP = snail.update[,"RP"]
  DAM = snail.update[,"DAM"]
  HAZ = snail.update[,"HAZ"]
  LG = snail.update[,"LG"]
  ingestion = sum(environment[1] - snail.update[,"Food"])
  
  Eggs = floor(RH/0.015)  # Figure out how many (whole) eggs are released
  RH = RH %% 0.015        # Remove released cercariae from the buffer
  Cercs = floor(RP/4e-5)  # Figure out how many (whole) cercs are released
  RP = RP %% 4e-5         # Remove released cercariae from buffer
  
  # Update environment
  Env_F = max(0.001, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-pars["r"]*pars["step"])) - ingestion)) # Analytical soln to logistic - ingestion
  Env_M = as.numeric(Infection.step[N.snails + 1] + pars["M_in"])
  Env_Z = as.numeric(environment[3]*exp(-pars["m_Z"]*pars["step"]) + sum(Cercs)/pars["ENV"])
  Env_G[day] = max(0, sum(Eggs))
  
  # Command back to NL
  NLCommand("ask patch 0 0 [set F", Env_F, "set M", Env_M, "set Z", Env_Z, "set G", Env_G[day], "]")
  snail.commands = paste(mapply(update.snails, who=snail.stats[,"who"], new.L=L, new.e=e, new.D=D, new.RH=RH, new.P=P, new.RP=RP, new.DAM=DAM, new.HAZ=HAZ, new.LG=LG), collapse=" ")
  NLCommand(snail.commands)
  if(day > 10){
    NLCommand("create-snails ", rbinom(n=1, size=Env_G[day - 10], prob=0.5), "[set L 0.75 set ee 0.9 set D 0 set RH 0 set P 0 set RPP 0 set DAM 0 set HAZ 0 set LG 0.75]")}
  NLCommand("go")
  
  day = day + 1
}


NLQuit()


### ---------------------------- Netlogo diagnostics attempted ---------------------------------     
#25-9-18
# 1. have tried this
Sys.setenv(NOAWT=1) 
Sys.unsetenv("NOAWT") 

# 2. have attempted this with NL v. 5.3.0 and 6.0.4
ver <-"6.0.4" # type in Netlogo version
nl.path <- "/Users/malishev/Documents/Melbourne Uni/Programs/" ; nl.path # set path to Netlogo program location

nl.path <- paste0(nl.path,"NetLogo ",ver,"/"); nl.path # opt 1
nl.path <- paste0(nl.path,"NetLogo ",ver,"/Java/"); nl.path # opt 2 
nl.path <- paste0(nl.path,"NetLogo ",ver,"/Java"); nl.path # opt 3 

model.path <- paste0(wd,"/"); model.path # set path to Netlogo model  
nl.model <- "DEB_INF_IBM_almost_working2.nlogo" # name of Netlogo model
NLStart(nl.path,nl.jarname = paste0("netlogo-",ver,".jar")) # open netlogo

