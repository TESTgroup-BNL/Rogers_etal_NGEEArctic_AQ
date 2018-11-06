####################################################################################################
#
#  	--- Last updated:  11.06.2018 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

list.of.packages <- c("httr","RCurl","readxl","tools","dplyr")  # packages needed for script
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

library(readxl)
library(tools)
library(dplyr)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Use GitHub data source?
use_GitHub <- TRUE

### Location of R scripts.
if (use_GitHub) {
  r.functions <- RCurl::getURL("https://raw.githubusercontent.com/TESTgroup-BNL/Rogers_etal_NGEEArctic_AQ/master/Rscripts/photo.processing.functions.R", 
                               ssl.verifypeer = FALSE)
  eval(parse(text = r.functions))
} else {
  r.functions <- file.path('')  # Path to photo.processing.functions.R
}

### Define directory containing the input LiCor 6400 data file to process 
in.dir <- file.path('Rogers_etal_NGEEArctic_AQ/input_data/')  # example, or see below when using GitHub as the source
dataset <- 'NGEE-Arctic_2016_AQ_Raw_Data.csv'  # using a specific data file on the local machine

### Define input data file name.
if (use_GitHub) {
  githubURL <- "https://raw.githubusercontent.com/TESTgroup-BNL/Rogers_etal_NGEEArctic_AQ/master/input_data/NGEE-Arctic_2016_AQ_Raw_Data.csv"
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

# *********************************** QA/QC Options ***********************************
### Sample QC checks
Cond.cutoff <- 0.05       ## Throw out observations with Cond < cutoff. E.g. 0.08
Ci.cutoff <- c(0,2000)    ## Throw out observations with Ci out of bounds
Tleaf.cutoff <- 0.9999    ## How much Tleaf variation to accept in curve data. E.g. 1 
# would allow variation of 1 degree around the mean Tleaf
# *********************************** QA/QC Options ***********************************

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


#--------------------------------------------------------------------------------------------------#
### Main output directory 
out.dir <- file.path(path.expand('~/scratch/Rogers_etal_NGEEArctic_AQ_Fits/Mean_Ci/2016/'))
unlink(out.dir,recursive=T) # delete old output if rerunning
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
setwd(file.path(out.dir)) # set working directory
getwd()  # check wd
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Apply data QA/QC functions
#data <- ge.data  #For debugging
ge.data.qc <- data.qc(data=ge.data,out.dir=out.dir,Cond.cutoff=Cond.cutoff,Ci.cutoff=Ci.cutoff,
                      Tleaf.cutoff=Tleaf.cutoff)
ge.data.qc[is.na(ge.data.qc)] <- NA
out.qc.ge.data <- data.frame(ge.data.qc$Sample.Info,ge.data.qc$GE.data)
write.csv(out.qc.ge.data,file=paste(out.dir,"/","QC_GE_Data.csv",sep=""),row.names=FALSE)
rm(out.qc.ge.data,ge.data)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Calculate mean Ci
# Get sample info and summary stats
samples <- unique(ge.data.qc$Sample.Info)

data.names <- names(ge.data.qc$GE.data)
remove.nms <- match(c("PRESS","PARI","PHOTO"),toupper(data.names))
data.names <- data.names[-remove.nms]
data.names # What to summarize

index <- within(ge.data.qc$GE.data, indx <- as.numeric(interaction(ge.data.qc$Sample.Info, 
                                                                   drop=TRUE,lex.order=TRUE)))

remove <- which(index$PARi > 100)
index2 <- index[-remove,]

remove <- which(samples$Sample_Barcode==1824 | samples$Sample_Barcode==1958 | samples$Sample_Barcode==1877 | samples$Sample_Barcode==2184)
samples <- samples[-remove,]

samples <- data.frame(samples,Index=unique(index2$indx))
samples <- samples[order(samples$Index),]
row.names(samples) <- seq(len=nrow(samples))
samples <- samples[,-match("Index",names(samples))]

means <- index2 %>% group_by(indx) %>% summarise_all(funs(mean))
names(means) <- paste0("Mean_",names(means))
means <- means[paste0("Mean_",data.names)]

output_data <- data.frame(samples,Mean_Tleaf=means$Mean_Tleaf,Mean_CiCa_lt_100umols=means$Mean_Ci_Ca,
                          Mean_Ci_lt_100umols=means$Mean_Ci)

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Outputs
write.csv(output_data,file=file.path(out.dir,"Mean_Ci_stats.csv"),row.names=FALSE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## EOF