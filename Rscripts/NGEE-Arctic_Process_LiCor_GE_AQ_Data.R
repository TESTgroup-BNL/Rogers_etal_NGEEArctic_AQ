####################################################################################################
#          
#   Farquhar Model Fitting Code: A-Q curve analysis
#   Author: Shawn P. Serbin
#
#
#   Estimate aQY & Amax (light and CO2 sat) using Farquhar equation optimization
#
#       --- Uses DEoptim global parameter optimization
#       --- Outputs Rd, aQY, and Amax (when a full curve is availible) estimated parameters 
#               in a single output .csv file with assocaited information
#       --- Output A-Q fit diagnostic figures and DEoptim parameter trace plots
#
#
#   Requirements:
#   DEoptim package
#   install.packages('DEoptim')
#
#
#  	--- Last updated: 11.14.2018 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

list.of.packages <- c("httr","RCurl","readxl","tools","DEoptim","dplyr")  # packages needed for script
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

# define function to grab files from GitHub
#devtools::source_gist("gist.github.com/christophergandrud/4466237")
source_GitHubData <-function(url, sep = ",", header = TRUE) {
  require(httr)
  request <- GET(url)
  stop_for_status(request)
  handle <- textConnection(content(request, as = 'text'))
  on.exit(close(handle))
  read.table(handle, sep = sep, header = header, quote = "")
}

# load libraries
library(readxl)
library(tools)
library(DEoptim)
library(dplyr)
#--------------------------------------------------------------------------------------------------#


#---------------- *User defined settings.* --------------------------------------------------------#
### Use GitHub data source?
use_GitHub <- TRUE

### Location of R scripts.  Needed for Farquhar model optimization. Contains functions.
if (use_GitHub) {
  r.functions <- RCurl::getURL("https://raw.githubusercontent.com/TESTgroup-BNL/Rogers_etal_NGEEArctic_AQ/master/Rscripts/photo.processing.functions.R", 
                ssl.verifypeer = FALSE)
  eval(parse(text = r.functions))
} else {
  r.functions <- file.path('')  # Path to photo.processing.functions.R
}

### Define directory containing the input LiCor 6400 data file to process 
in.dir <- file.path('Rogers_etal_NGEEArctic_AQ/input_data/')  # example, or see below when using GitHub as the source
dataset <- 'NGEE-Arctic_Light_Response_Barrow_2016_Raw_Data.csv'  # using a specific data file on the local machine

### Define input data file name.
if (use_GitHub) {
  githubURL <- "https://raw.githubusercontent.com/TESTgroup-BNL/Rogers_etal_NGEEArctic_LightResponse/master/input_data/NGEE-Arctic_Light_Response_Barrow_2016_Raw_Data.csv"
  dataset <- source_GitHubData(githubURL)
  ge.data <- dataset
} else {
  dataset <- dataset  # using a specific data file on the local machine
  
  ### Define input file to be processed
  extension <- file_ext(dataset)
  if ("csv" %in% extension) {
    print("*.csv")
    ge.data <- read.table(file.path(in.dir,dataset), header=T,sep=",",quote = "")
  } else if ("xls" %in% extension | "xlsx" %in% extension) {
    print("*.xl file")
    ge.data <- read_excel(path = file.path(in.dir,dataset), sheet = 1)
    ge.data <- data.frame(ge.data)
  }
}

# clean up
if (!is.null(ge.data$Rep)){
  ge.data$Rep[ge.data$Rep==-9999] <- 1
}
ge.data[ge.data==-9999] <- NA
#ge.data[ge.data==""] <- NA
if (is.null(ge.data$deltaT) && !is.null(ge.data$Tl_minus_Ta)) {
  ge.data$deltaT <- ge.data$Tl_minus_Ta
} else if (is.null(ge.data$deltaT) && is.null(ge.data$Tl_minus_Ta)) {
  # Calculate deltaT
  ge.data$deltaT <- ge.data$Tleaf - ge.data$Tair
}
if (is.null(ge.data$Ci_Ca) && !is.null(ge.data$Ci_over_Ca)) {
  ge.data$Ci_Ca <- ge.data$Ci_over_Ca
}
summary(ge.data)    ## Summary of dataset


#--------------------------------------------------------------------------------------------------#
### Define main output directory 
out.dir <- file.path(path.expand('~/scratch/Rogers_etal_NGEEArctic_LightResponse/Rogers_etal_NGEEArctic_LightResponseCurve_Fits/'))
unlink(out.dir,recursive=T) # delete old output if rerunning
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
setwd(file.path(out.dir)) # set working directory
getwd()  # check wd
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# *********************************** QA/QC Options ***********************************
###
aQY.lower.cutoff <- 15      # lower light level to include in AQ fit
aQY.cutoff <- 80            # Upper AQ fit cutoff umols PAR  Options: E.g. 150, 100, 75
Amax.cutoff <- 450          # umols PAR
high_light.cutoff <- 2200   # umols PAR
aQY.min.n <- 3              # minimum number of obs for calculating aparent QY.  Best is 3.
Amax.min.n <- 2             # minimum number of obs for calculating Amax.  Best is 3.

### Sample QC checks
Cond.cutoff <- 0.05       ## Throw out observations with Cond < cutoff. E.g. 0.08
Ci.cutoff <- c(0,2000)    ## Throw out observations with Ci out of bounds
Tleaf.cutoff <- 0.9999    ## How much Tleaf variation to accept in curve data. E.g. 1 
# would allow variation of 1 degree around the mean Tleaf
# *********************************** QA/QC Options ***********************************


# ********************************** DEoptim Options **********************************
### DEoptim options. Controls the parameter search space for each DEoptim iteration.
lower.bound <- c(1,0.01)          ## Lower bound of Amax, & theta2
#lower.bound <- c(1,0.75)          ## Lower bound of Amax, & theta2
upper.bound <- c(100,0.9999)      ## Upper bound of Amax, & theta2
max.iters <- 1000                 ## Max iterations
RMSE.min <- 0.05                  ## Min value of RMSE to be reached (VTR) during 
NP <- 100                         ## Number of population members. For many problems it is best to set 
                                  ## NP to be at least 10 times the length of the parameter vector.
DWF <- 0.8                        ## Differential weighting factor from interval [0,2]. Default to 0.8.
CR <- 0.6                         ## Crossover probability from interval [0,1]. Default to 0.5.
# DEoptim optimization. DEoptim will stop when either 
# 1) RMSE=RMSE.min or 2) max iters are reached.
# Should be sufficiently small, e.g. 0.5. Set to 
# a low number to ensure convergance (e.g. 0.0001)
# or if unsure. This option can speed up convergence
#see ?DEoptim
# ********************************** DEoptim Options **********************************


# ************************************* References ************************************
# Thornley JHM. 1976. Mathematical models in plant physiology. London: Academic Press.

# Bernacchi, C.J., Pimentel, C., & Long, S.P. (2003). In vivo temperature response 
# functions of parameters required to model RuBP-limited photosynthesis. 
# Plant Cell and Environment, 26, 1419-1430

# Posada, J.M., Lechowicz, M.J., & Kitajima, K. (2009). Optimal photosynthetic use of 
# light by tropical tree crowns achieved by adjustment of individual leaf angles and 
# nitrogen content. Annals of Botany, 103, 795-805

# Posada, J.M., Sievanen, R., Messier, C., Perttunen, J., Nikinmaa, E., & Lechowicz, M.J. 
# (2012). Contributions of leaf photosynthetic capacity, leaf angle and self-shading to 
# the maximization of net photosynthesis in Acer saccharum: a modelling assessment. 
# Annals of Botany, 110, 731-741
# ************************************* References ************************************
#--------------------------------------------------------------------------------------------------#


#==================================================================================================#
# End of user defined options. Start of program.
#==================================================================================================#


###################################### START SCRIPT ################################################
print('**************************************************')
print('**************** STARTING SCRIPT *****************')
print('**************************************************')
print(' ')
####################################################################################################


#--------------------------------------------------------------------------------------------------#
# Load utils and Farquhar functions for data QA/QC and optimization
if (use_GitHub) {
  print("*** Loading helper functions from GitHub source ***")
} else {
  source(file.path(r.functions,'photo.processing.functions.R'))
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Setup output dirs
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
aq.fig.dir <- file.path(out.dir,"AQ_Diagnostics")
if (! file.exists(aq.fig.dir)) dir.create(aq.fig.dir)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Apply data QA/QC functions
ge.data.qc <- data.qc(data=ge.data,out.dir=out.dir,Cond.cutoff=Cond.cutoff,Ci.cutoff=Ci.cutoff,
                      Tleaf.cutoff=Tleaf.cutoff)
ge.data.qc[is.na(ge.data.qc)] <- NA
out.qc.ge.data <- data.frame(ge.data.qc$Sample.Info,ge.data.qc$GE.data)
write.csv(out.qc.ge.data,file=paste(out.dir,"/","QC_GE_Data.csv",sep=""),row.names=FALSE)
rm(out.qc.ge.data,ge.data)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Get sample info and summary stats
samples <- unique(ge.data.qc$Sample.Info)

### Get data names
data.names <- names(ge.data.qc$GE.data)
remove.nms <- match(c("PRESS","PARI","PHOTO"),toupper(data.names))
data.names <- data.names[-remove.nms]
data.names # What to summarize

index <- within(ge.data.qc$GE.data, indx <- as.numeric(interaction(ge.data.qc$Sample.Info, 
                                                                   drop=TRUE,lex.order=TRUE)))
samples <- data.frame(samples,Index=unique(index$indx))
samples <- samples[order(samples$Index),]
row.names(samples) <- seq(len=nrow(samples))
samples <- samples[,-match("Index",names(samples))]

### Obs stats
means <- index %>% group_by(indx) %>% summarise_all(funs(mean))
sdevs <- index %>% group_by(indx) %>% summarise_all(funs(sd))
names(means) <- paste0("Mean_",names(means))
means <- means[paste0("Mean_",data.names)]
sdevs <- sdevs[data.names]
CVs <- (sdevs/means)*100

# Create output dataframe of diagnostic info
out.diagnostics <- data.frame(Leaf_Area=means$Mean_Area,Mean_Tair=means$Mean_Tair,Mean_Tleaf=means$Mean_Tleaf,
                              Tleaf.CV=CVs$Tleaf,Mean_deltaT=means$Mean_deltaT,Mean_RH_R=means$Mean_RH_R,
                              Mean_RH_S=means$Mean_RH_S,Mean_VpdA=means$Mean_VpdA,Mean_VpdL=means$Mean_VpdL,
                              Mean_PARo=means$Mean_PARo,Mean_Ci_Ca=means$Mean_Ci_Ca,Ci_Ca.CV=CVs$Ci_Ca,
                              Mean_CO2R=means$Mean_CO2R,Mean_CO2S=means$Mean_CO2S,CO2S.CV=CVs$CO2S,
                              Mean_Cond=means$Mean_Cond,Cond.CV=CVs$Cond,Mean_Ci=means$Mean_Ci)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Fit Farquhar model to data
data <- data.frame(ge.data.qc$Sample.Info,ge.data.qc$GE.data)

# Setup outputs
mean_Ci <- array(data=NA,dim=dim(samples)[1])
Amb.Photo <- array(data=NA,dim=dim(samples)[1])
Rd <- array(data=NA,dim=dim(samples)[1])
aQY <- array(data=NA,dim=dim(samples)[1])
LCP <- array(data=NA,dim=dim(samples)[1])
Amax <- array(data=NA,dim=dim(samples)[1])
convex <- array(data=NA,dim=dim(samples)[1])
Anet.1500 <- array(data=NA,dim=dim(samples)[1])
Anet.1800 <- array(data=NA,dim=dim(samples)[1])
Anet.2000 <- array(data=NA,dim=dim(samples)[1])
LUE.max <- array(data=NA,dim=dim(samples)[1])
ILUE.max <- array(data=NA,dim=dim(samples)[1])
RMSE.DEoptim <- array(data=NA,dim=dim(samples)[1])
RMSE.photo <- array(data=NA,dim=dim(samples)[1])

# Main outer loop
for (i in seq_along(1:dim(samples)[1])) {
  sub.data <- merge(data,samples[i,],by=names(samples[i,]))
  sub.data = sub.data[order(sub.data$PARi),]
  print("--- Processing Sample: ")
  print(paste(as.vector(unlist(samples[i,])),collapse=" "))
  
  # Determine what level of GE processing
  chk2 <- length(which(sub.data$PARi<aQY.cutoff))
  chk3 <- length(which(sub.data$PARi>Amax.cutoff))
  
  # Find ambient photo
  keep <- which(sub.data$PARi >= 1460)
  if (length(keep)>0) {
    Amb.Photo[i] <- mean(sub.data$Photo[keep],na.rm=TRUE)
  } else {
    Amb.Photo[i] <- -9999
  }
  
  # Check if there are enough obs. Remove curves with to few observations
  if (chk2<aQY.min.n){
    mean_Ci[i] <- -9999
    Rd[i] <- -9999
    aQY[i] <- -9999
    Amax[i] <- -9999
    convex[i] <- -9999
    Anet.1500[i] <- -9999
    Anet.1800[i] <- -9999
    Anet.2000[i] <- -9999
    LCP[i] <- -9999
    LUE.max[i] <- -9999
    ILUE.max[i] <- -9999
    RMSE.DEoptim[i] <- -9999
    RMSE.photo[i] <- -9999
  } else {
    
    ## calculate mean CI within init part of the curve
    sub_data_PAR <- subset(sub.data, PARi>aQY.lower.cutoff & PARi < aQY.cutoff)
    mean_Ci[i] <- mean(sub_data_PAR$Ci)
    
    # Extract photo and PAR data
    Photo <- sub.data[,match("PHOTO",toupper(names(sub.data)))] # Find photo data
    PARi <- sub.data[,match("PARI",toupper(names(sub.data)))] # Find photo data
    
    # Find quantum efficiency first
    reg <- lm(Photo[which(PARi>aQY.lower.cutoff & PARi<=aQY.cutoff)]~PARi[which(PARi>aQY.lower.cutoff & PARi<=aQY.cutoff)],sub.data)
    aQY[i]  <- coef(reg)[[2]]
    Rd[i]  <- coef(reg)[[1]]
    if (Rd[i] <0) {
      Rd[i]  <- Rd[i] *-1
    }
    x <- seq(-200,400,0.1)
    LCP[i]  <- x[which(round((aQY[i]*x-Rd[i]),2)==0)[1]]
    
    # Calculate APAR
    IaQY <- aQY[i]*sub.data$PARi
    
    f.model <- NA # Type variable for defining plots
    
    # Full curve analysis
    if (chk2>=aQY.min.n & chk3>=Amax.min.n) {  
      f.model <- 1
      fit <- DEoptim(A.Q,lower=lower.bound, upper=upper.bound, DEoptim.control(NP=NP,
                                                               F=DWF,CR=CR, itermax=max.iters, VTR=RMSE.min, strategy=2, trace=FALSE))
      tempout <- data.frame(fit$optim)
      
      # Stats
      Amax[i] <- tempout[1,1]
      convex[i] <- tempout[2,1]
      RMSE.DEoptim[i] <- fit$optim$bestval
      mod.Photo <- (((IaQY+Amax[i])-sqrt((IaQY+Amax[i])^2-(4*convex[i]*IaQY*Amax[i])))/(2*convex[i]))-Rd[i]
      residuals <- mod.Photo-Photo
      RMSE.photo[i] <- sqrt(mean((residuals)^2))
      
      # Calculate LUE.max --> Posada et al., (2012)
      mod.PARi <- seq(0,2000,1)
      mod.IaQY <- aQY[i]*mod.PARi
      mod.Photo <- (((mod.IaQY+Amax[i])-sqrt((mod.IaQY+Amax[i])^2-(4*convex[i]*mod.IaQY*Amax[i])))/(2*convex[i]))-Rd[i]
      LUE <- mod.Photo/mod.PARi
      LUE.max[i] <- max(LUE)
      ILUE.max[i] <- mod.PARi[which.max(LUE)] # PARi of LUE.max
      
      # Calculate modeled Anet.1500
      mod.PARi <- 1500
      mod.IaQY <- aQY[i]*mod.PARi
      Anet.1500[i] <- (((mod.IaQY+Amax[i])-sqrt((mod.IaQY+Amax[i])^2-(4*convex[i]*mod.IaQY*Amax[i])))/(2*convex[i]))-Rd[i]
      rm(mod.PARi,mod.IaQY)
      
      # Calculate modeled Anet.1800
      mod.PARi <- 1800
      mod.IaQY <- aQY[i]*mod.PARi
      Anet.1800[i] <- (((mod.IaQY+Amax[i])-sqrt((mod.IaQY+Amax[i])^2-(4*convex[i]*mod.IaQY*Amax[i])))/(2*convex[i]))-Rd[i]
      rm(mod.PARi,mod.IaQY)
      
      # Calculate modeled Anet.2000
      mod.PARi <- 2000
      mod.IaQY <- aQY[i]*mod.PARi
      Anet.2000[i] <- (((mod.IaQY+Amax[i])-sqrt((mod.IaQY+Amax[i])^2-(4*convex[i]*mod.IaQY*Amax[i])))/(2*convex[i]))-Rd[i]
      rm(mod.PARi,mod.IaQY)
      
      ### Display info to console 
      print(paste("aQY: ",round(aQY[i],4)))
      print(paste("LCP: ",round(LCP[i],4)))
      print(paste("Amax: ",round(Amax[i],2)))
      print(paste("convex: ",round(convex[i],2)))
      print(paste("Anet.2000: ",round(Anet.2000[i],2)))
      print(paste("LUE.max: ",round(LUE.max[i],4)))
      print(paste("Rd: ",round(Rd[i],4)))  
      print(paste("DEoptim RMSE: ",round(RMSE.DEoptim[i],2)))
      print(paste("Photo RMSE: ",round(RMSE.photo[i],2)))
      print("-----------------------------------------------------------")
      flush.console()
      
    } else {
      f.model <- 2
      Amax[i] <- -9999
      convex[i] <- -9999
      Anet.1500[i] <- -9999
      Anet.1800[i] <- -9999
      Anet.2000[i] <- -9999
      LUE.max[i] <- -9999
      ILUE.max[i] <- -9999
      RMSE.DEoptim[i] <- -9999
      
      mod.Photo <- aQY[i]*PARi[PARi<=100] - Rd[i]
      residuals <- mod.Photo-Photo[which(PARi<=100)]
      RMSE.photo[i] <- sqrt(mean((residuals)^2))
      
    } # End if/else processing
    
    ### Generate fit diagnostics
    if (dim(sub.data)[1]>=aQY.min.n){
      sample.name <- paste(unlist(samples[i,]),collapse="_")
      sample.name <- gsub(pattern="/","-",sample.name)
      params <- data.frame(aQY=aQY[i],LCP=LCP[i],Amax=Amax[i],convex=convex[i],LUE.max[i],Rd=Rd[i],RMSE.photo[i])
      plot.ge.fit(type="A-Q",sub.data,fit,params,aq.fig.dir,sample.name,f.model)
      rm(f.model)
    }
  } # End if/else obs check
  print("** Removing temp objects ***")
  rm(chk2,chk3, reg, IaQY, f.model, fit, tempout, mod.Photo, residuals, mod.PARi, mod.IaQY,
     mod.Photo, LUE)
} # End for loop
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Generate output
output.data <- data.frame(samples,out.diagnostics,mean_Ci_init_slope=mean_Ci,Amb.Photo=Amb.Photo,
                          aQY=aQY,Amax.g=Amax,Convex=convex,
                          Anet.1500=Anet.1500,Anet.1800=Anet.1800,Anet.2000=Anet.2000,Rd=Rd,
                          LCP=LCP,LUE.max=LUE.max,ILUE.max=ILUE.max,
                          RMSE.DEoptim=RMSE.DEoptim,RMSE.photo=RMSE.photo)
if (use_GitHub) {
  output.dataset <- "Example_fitted_A-Q_response_curves.csv"
} else {
  output.dataset <- paste(strsplit(dataset,".csv"),".processed.csv",sep="")
}

write.csv(output.data,file=paste(out.dir,"/",output.dataset,sep=""),row.names=FALSE)

rm(list=ls(all=TRUE))   # clear workspace
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF