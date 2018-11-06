####################################################################################################
#
#  	--- Last updated:  11.06.2018 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

library(readxl)
library(tools)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Location of R scripts.  Needed for Farquhar model optimization. Contains functions.
r.functions <- '/Volumes/TEST/Manuscripts/Rogers_Arctic_AQ/R_Scripts/Photo/'  # Path to photo.processing.functions.R
source(paste0(r.functions,'/','photo.processing.functions.R'))

# *********************************** QA/QC Options ***********************************
### Sample QC checks
Cond.cutoff <- 0.05       ## Throw out observations with Cond < cutoff. E.g. 0.08
Ci.cutoff <- c(0,2000)    ## Throw out observations with Ci out of bounds
Tleaf.cutoff <- 0.9999    ## How much Tleaf variation to accept in curve data. E.g. 1 
# would allow variation of 1 degree around the mean Tleaf
# *********************************** QA/QC Options ***********************************

### Input LI6400 dataset.  First define location of file (i.e. directory). 
in.dir <- '/Volumes/TEST/Projects/NGEE-Arctic/Data/Barrow/Gas_Exchange/AQ/field data/Final_Data_Compilation/'

#dataset <- 'NGEE-Arctic_2014_AQ_Raw_Data_AR_QA2.xlsx'
dataset <- 'NGEE-Arctic_2016_AQ_Raw_Data_AR_QA3.xlsx'

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
#out.dir <- file.path('/Volumes/TEST/Manuscripts/Rogers_Arctic_AQ/R_Ouput/Mean_Ci/2014/')
out.dir <- file.path('/Volumes/TEST/Manuscripts/Rogers_Arctic_AQ/R_Ouput/Mean_Ci/2016/')
unlink(out.dir,recursive=T) # delete old output if rerunning
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
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
#ge.data.qc$Sample.Info <- ge.data.qc$Sample.Info[,-2]

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

#which(!index$indx %in% unique(index2$indx))
remove <- which(samples$Sample_Barcode==1824 | samples$Sample_Barcode==1958 | samples$Sample_Barcode==1877 | samples$Sample_Barcode==2184)
samples <- samples[-remove,]

samples <- data.frame(samples,Index=unique(index2$indx))
samples <- samples[order(samples$Index),]
row.names(samples) <- seq(len=nrow(samples))
samples <- samples[,-match("Index",names(samples))]

means <- aggregate(.~index2$indx,data=index2,mean)
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