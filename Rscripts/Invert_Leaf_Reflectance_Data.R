####################################################################################################
#          
#   Invert leaf reflectance data to estimate transmittance and calculate absorption
#   
#   Approach: BayesianTools DEzs MCMC algorithm  (http://dream.r-forge.r-project.org/)
#
#   Author: Shawn P. Serbin
#
#
#
#
#
#    --- Last updated:  11.06.2018 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#---------------- Load required libraries ---------------------------------------------------------#
list.of.packages <- c("devtools","readxl")  # packages needed for script
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

# install required PEcAn packages (pecanproject.org)
ok <- require(PEcAn.logger) ; if (! ok) {
  devtools::install_github("PEcAnProject/pecan", subdir="base/logger", ref = "develop") # use development version of PEcAn
} else{
  print("*** Package found: PEcAn.logger ***")
}

ok <- require(PEcAn.remote) ; if (! ok) {
  devtools::install_github("PEcAnProject/pecan", subdir="base/remote", ref = "develop") # use development version of PEcAn
} else {
  print("*** Package found: PEcAn.remote ***")
}

ok <- require(PEcAn.assim.batch) ; if (! ok) {
  devtools::install_github("PEcAnProject/pecan", subdir="modules/assim.batch", ref = "develop") # use development version of PEcAn
} else {
  print("*** Package found: PEcAn.assim.batch ***")
}

ok <- require(PEcAnRTM) ; if (! ok) {
  devtools::install_github("PEcAnProject/pecan", subdir="modules/rtm", ref = "develop") # use development version of PEcAn
} else {
  print("*** Package found: PEcAnRTM ***")
}

library(PEcAnRTM)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Use GitHub data source?
use_GitHub <- TRUE
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Import leaf spectra and leaf trait / species info
if (use_GitHub) {
  spectra_file
} else {
  
}


#--------------------------------------------------------------------------------------------------#

