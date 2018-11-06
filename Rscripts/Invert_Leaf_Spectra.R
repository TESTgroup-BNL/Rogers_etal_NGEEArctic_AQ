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
#    --- Last updated:  10.10.2017 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#---------------- Load required libraries ---------------------------------------------------------#
ok <- require(devtools) ; if (! ok) {
  install.packages('devtools', repos="http://cran.rstudio.com/")
} else {
  print("*** Package found: devtools ***")
}

ok <- require(readxl) ; if (! ok) {
  install.packages('readxl', repos="http://cran.rstudio.com/")
} else {
  print("*** Package found: readxl ***")
}

ok <- require(PEcAn.logger) ; if (! ok) {
  library(devtools)
  install_github("PEcAnProject/pecan", subdir="base/logger", ref = "develop") # use development version of PEcAn
} else{
  print("*** Package found: PEcAn.logger ***")
}


ok <- require(PEcAn.assim.batch) ; if (! ok) {
  library(devtools)
  install_github("PEcAnProject/pecan", subdir="modules/assim.batch", ref = "develop") # use development version of PEcAn
} else {
  print("*** Package found: PEcAn.assim.batch ***")
}

ok <- require(PEcAnRTM) ; if (! ok) {
  library(devtools)
  #install_github("PEcAnProject/pecan", subdir="modules/rtm", ref = "develop") # use development version of PEcAn
  install_github('ashiklom/pecan/modules/rtm@rtm-no-neff')
  
} else {
  print("*** Package found: PEcAnRTM ***")
}

library(PEcAnRTM)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Import leaf spectra and leaf trait / species info
in_dir <- '/Volumes/TEST/Projects/NGEE-Arctic/Data/Barrow/Spectra/SVC_HR-1024i/2016_Data/'
spectra_file <- 'NGEE-Arctic_Barrow_2016_HR1024i_Leaf_Spectral_Reflectance.xlsx'
refl_data <- read_excel(path=file.path(in_dir,spectra_file), sheet = 2)


# need to update this part to read in spreadsheet that translates BNL barcode to species name

sample_info_dir <- '~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/Serbin_Arctic_Spectra-Trait/Data/Raw_trait_data/'
sample_info_file <- 'Barrow2016_samples_LMA.xlsx'
sample_info <- read_excel(path=file.path(sample_info_dir,sample_info_file), sheet = 1)
sample_info[sample_info==-9999] <- NA     
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Output directory
manuscript_dir <- "/Volumes/TEST/Manuscripts/Rogers_Arctic_AQ/"
out.dir <- file.path(manuscript_dir,'R_Ouput', 'Inverted_Leaf_Spectra/Range_400_700nm')
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Subset spectral data to 
Start.wave <- 400 # 480
End.wave <- 2500

abs.Start.wave <- 400  # start abs calc wavelength
abs.End.wave <- 700    # end abs calc wavelength

keep_bands <- names(refl_data)[match(paste0("Wave_",seq(Start.wave,End.wave,1)),names(refl_data))]
sub_refl_data <- data.frame(droplevels(refl_data[,keep_bands]))*0.01

if (exists("trans_data")) {
  keep_bands <- names(trans_data)[match(paste0("Wave_",seq(Start.wave,End.wave,1)),names(trans_data))]
  sub_trans_data <- data.frame(droplevels(trans_data[,keep_bands]))*0.01
}

waves <- seq(Start.wave,End.wave,1)
prospect_waves <- seq(400,2500,1)
keep_waves <- prospect_waves %in% waves
#plot(waves,final.spec[1,], type="l", col="black", lwd=3)

# grab spec info for both reflectance and transmittance inputs (if availible)

if (exists("trans_data")) {
  refl_spec_info <- refl_data[,!(names(refl_data) %in% paste0("Wave_",seq(350,2500,1))) ]
  trans_spec_info <- trans_data[,!(names(trans_data) %in% paste0("Wave_",seq(350,2500,1))) ]
} else {
  refl_spec_info <- refl_data[,!(names(refl_data) %in% paste0("Wave_",seq(350,2500,1))) ]
}


## Quick diagnostics plots
plot(seq(Start.wave,End.wave,1),sub_refl_data[1,], type= "l")

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Create spec file + info to be inverted
refl_spec_info2 <- refl_spec_info[,c(1,2,3,5,6,7)]
sample_info2 <- sample_info[,c(1,4,6,7)]  # just need to keep this
output_sample_info <- merge(refl_spec_info2, sample_info2, by="Sample_Barcode")
names_output_sample_info <- names(output_sample_info)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Run inversions - whole directory
output.LRT <- list(Spec.Info=array(NA,dim=c(dim(output_sample_info)[1],dim(output_sample_info)[2])),
                   obs.Reflectance=array(NA,dim=c(dim(sub_refl_data)[1],length(waves))),
                   mod.Reflectance=array(NA,dim=c(dim(sub_refl_data)[1],length(prospect_waves))),
                   mod.Transmittance=array(NA,dim=c(dim(sub_refl_data)[1],length(prospect_waves))))

mod.params <- array(NA,dim=c(dim(sub_refl_data)[1],15)) # P5 param output
inv.samples <- NA
# names: N.mu, N.q25, N.q975, Cab.mu, Cab.q25, Cab.q75, Car.mu, Car.q25, Car.q75,
# Cw.mu, Cw.q25, Cw.q75, Cm.mu, Cm.q25, Cm.q75, gelman.diag
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Run inversion
prospect_ver <- 5

# setup prospect error envelope list
p.refl.stats <- list(lower = array(data=NA,c(dim(sub_refl_data)[1],2101)),
                     upper = array(data=NA,c(dim(sub_refl_data)[1],2101)))

# setup inversion model
model <- function(x) prospect(x, 5)[min(which(prospect_waves %in% waves, arr.ind = TRUE)):2101,1] # subset

# setup prior
prior <- prospect_bt_prior(5)
#prior$lower <- c(1,5,1,0,0,0) # N, Cab, Car, Cw, Cm, resid
#prior$upper <- c(1,Inf,30,0,0,0) # N, Cab, Car, Cw, Cm, resid

# Plot title variable
title_var <- "BNL_Barcode"

# Set progress bar
# ?see txtProgressBar
print("Starting Inversion:")
print(paste0("Inverting: ", dim(sub_refl_data)[1]))
print(" ")
pb <- txtProgressBar(min = 0, max = dim(sub_refl_data)[1], width= 50,style = 3)
system.time(for (i in seq_along(1:dim(sub_refl_data)[1]) ) {
#system.time(for (i in seq_along(1:3) ) {
  print(" ")
  print(paste0("Inverting: ",unlist(refl_spec_info[i,title_var])))
  samples <- invert_bt(observed = t(sub_refl_data[i,]), model = model, prior = prior,
                       custom_settings = list(init = list(iterations = 2000),
                                              loop = list(iterations = 1000),
                                              other = list(max_iter = 1e6, min_samp = 5000, threshold = 1.1)))
  samples_burned <- PEcAn.assim.batch::autoburnin(BayesianTools::getSample(samples, coda = TRUE), method = 'gelman.plot')
  mean_estimates <- do.call(cbind, summary(samples_burned)[c('statistics', 'quantiles')])
  
  # include process error in error envelop (i.e. generate prediction interval)
  n_target <- 1000
  spec.length <- 2101
  param.samples <- do.call(rbind, samples_burned)
  param.samples <- param.samples[sample(nrow(param.samples), n_target), ]
  RT_pred <- array(data=NA,c(n_target,2101))
  print("*** Calculating error stats ***")
  for (r in seq_len(n_target)) {
    perturbed.prospect.ref <- rnorm(spec.length,prospect(param.samples[r,1:5], 5)[,1], param.samples[r,6])
    RT_pred[r,] <- perturbed.prospect.ref
  }
  
  # stats
  p.refl.stats$lower[i,] <- apply(RT_pred,2,quantile,probs=c(0.05), na.rm=T)
  p.refl.stats$upper[i,] <- apply(RT_pred,2,quantile,probs=c(0.95), na.rm=T)
  
  pdf(file = file.path(out.dir,paste0(unlist(refl_spec_info[i,title_var]),'_MCMC_trace_diag.pdf')), 
      width = 8, height = 6, onefile=T)
  par(mfrow=c(1,1), mar=c(2,2,2,2), oma=c(0.1,0.1,0.1,0.1)) # B, L, T, R
  plot(samples_burned)
  dev.off()
  
  # Generate modeled spectra
  input.params <- as.vector(unlist(mean_estimates[,1]))[1:prospect_ver]
  LRT <- prospect(param = input.params, version=prospect_ver)
  refl_spec_info2 <- refl_spec_info[i,c(1,2,3,5,6,7)]
  output_sample_info <- merge(refl_spec_info2, sample_info2, by="Sample_Barcode")
  output.LRT$Spec.Info[i,] <- as.character(output_sample_info)
  output.LRT$obs.Reflectance[i,] <- as.vector(unlist(sub_refl_data[i,]))
  output.LRT$mod.Reflectance[i,] <- LRT[,1]
  output.LRT$mod.Transmittance[i,] <- LRT[,2]
  
  # Extract results
  mod.params[i,] <- c(mean_estimates[row.names(mean_estimates)=="N",colnames(mean_estimates)=="Mean"],
                      mean_estimates[row.names(mean_estimates)=="N",colnames(mean_estimates)=="25%"],
                      mean_estimates[row.names(mean_estimates)=="N",colnames(mean_estimates)=="97.5%"],
                      mean_estimates[row.names(mean_estimates)=="Cab",colnames(mean_estimates)=="Mean"],
                      mean_estimates[row.names(mean_estimates)=="Cab",colnames(mean_estimates)=="25%"],
                      mean_estimates[row.names(mean_estimates)=="Cab",colnames(mean_estimates)=="97.5%"],
                      mean_estimates[row.names(mean_estimates)=="Car",colnames(mean_estimates)=="Mean"],
                      mean_estimates[row.names(mean_estimates)=="Car",colnames(mean_estimates)=="25%"],
                      mean_estimates[row.names(mean_estimates)=="Car",colnames(mean_estimates)=="97.5%"],
                      mean_estimates[row.names(mean_estimates)=="Cw",colnames(mean_estimates)=="Mean"],
                      mean_estimates[row.names(mean_estimates)=="Cw",colnames(mean_estimates)=="25%"],
                      mean_estimates[row.names(mean_estimates)=="Cw",colnames(mean_estimates)=="97.5%"],
                      mean_estimates[row.names(mean_estimates)=="Cm",colnames(mean_estimates)=="Mean"],
                      mean_estimates[row.names(mean_estimates)=="Cm",colnames(mean_estimates)=="25%"],
                      mean_estimates[row.names(mean_estimates)=="Cm",colnames(mean_estimates)=="97.5%"])
  
  setTxtProgressBar(pb, i)
  rm(samples,samples_burned,input.params,LRT,mean_estimates,
     param.samples,RT_pred,perturbed.prospect.ref,refl_spec_info2, output_sample_info)
  print(" ")
  print(" ")
  print("*** Starting new inversion ***")
  print(" ")
  print(" ")
  
  flush.console()
  
  
}) ## End of inversion loop
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Clean up
# names: N.mu, N.q25, N.q975, Cab.mu, Cab.q25, Cab.q75, Car.mu, Car.q25, Car.q75,
# Cw.mu, Cw.q25, Cw.q75, Cm.mu, Cm.q25, Cm.q75
mod.params <- as.data.frame(mod.params)
names(mod.params) <- c("N.mu", "N.q25", "N.q975", "Cab.mu", "Cab.q25", "Cab.q975", "Car.mu", "Car.q25",
                       "Car.q975","Cw.mu", "Cw.q25", "Cw.q975", "Cm.mu", "Cm.q25", "Cm.q975")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Plot comparison
#wl <- seq(400,2500,1)
#waves <- seq(500,2500,1)
pdf(file=file.path(out.dir,'PROSPECT5_Inversion_Diagnostics.pdf'),height=8,width=9)
par(mfrow=c(1,1), mar=c(4.3,4.3,1.0,4.3), oma=c(0.1,0.1,0.1,0.1)) # B L T R
for (i in seq_along(1:dim(sub_refl_data)[1] )) {
#for (i in seq_along(1:3 )) {
#for (i in seq_along(1:1 )) {
  plot(waves,output.LRT$obs.Reflectance[i,], type="l", col="black",xlab="Wavelength (nm)",ylab="Reflectance (0-1)",
       lwd=3,main=paste0(output.LRT$Spec.Info[i,4]," ", output.LRT$Spec.Info[i,7]) )
  lines(prospect_waves,prospect(param = mod.params[i,c(1,4,7,10,13)], version=prospect_ver)[,1],col="red",lwd=2)
  polygon(c(prospect_waves ,rev(prospect_waves)),c(p.refl.stats$upper[i,], rev(p.refl.stats$lower[i,])),
          col="grey70",border=NA)
  lines(waves,output.LRT$obs.Reflectance[i,],col="black",lwd=2)
  lines(prospect_waves,prospect(param = mod.params[i,c(1,4,7,10,13)], version=prospect_ver)[,1],col="red",lwd=2)
  legend("topright",legend=c("Observed","Modeled","Modeled 95% PI"),lty=1,col=c("black","red","grey70"),
         lwd=c(3,3,8),bty = "n")
  box(lwd=2.2)
  plot(waves,output.LRT$obs.Reflectance[i,], type="l", col="black",xlab="Wavelength (nm)",ylab="Reflectance (0-1)",
       lwd=3,ylim=c(0,1),main=paste0(output.LRT$Spec.Info[i,4]," ", output.LRT$Spec.Info[i,7]),cex.lab=1.7)
  lines(prospect_waves,prospect(param = mod.params[i,c(1,4,7,10,13)], version=5)[,1],col="red",lwd=3)
  lines(prospect_waves,1-prospect(param = mod.params[i,c(1,4,7,10,13)], version=5)[,2],col="grey70",lwd=3) 
  axis(4,labels = rev(seq(0.0, 1.0, 0.2)), at = seq(0.0, 1.0, 0.2))
  mtext("Transmittance (0-1)",side=4,line=3,cex=1.7)
  legend("right",legend=c("Meas. Reflectance","Mod. Reflectance","Mod. Transmittance"),
         lty=1,col=c("black","red","grey70"),
         lwd=c(3,3,8),bty = "n")
  box(lwd=2.2)
}
dev.off()
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Calculate absorption - use measured refl, modeled trans
full_spectrum_absorption <- data.frame(1 - output.LRT$mod.Transmittance - output.LRT$obs.Reflectance)
names(full_spectrum_absorption) <- paste0("Wave_",seq(Start.wave,End.wave,1))
plot(waves, full_spectrum_absorption[1,],type="l")
lines(waves, full_spectrum_absorption[2,])
lines(waves, full_spectrum_absorption[3,])

#abs.Start.wave <- 400  # start abs calc wavelength
#abs.End.wave <- 750    # end abs calc wavelength
keep_bands <- names(full_spectrum_absorption)[match(paste0("Wave_",seq(abs.Start.wave,abs.End.wave,1)),
                                                    names(full_spectrum_absorption))]
sub_spec_abs <- data.frame(droplevels(full_spectrum_absorption[,keep_bands]))
leaf_vis_absorption <- rowMeans(sub_spec_abs,na.rm = T)
leaf_vis_absorption_final <- data.frame(output.LRT$Spec.Info, leaf_vis_absorption)
names(leaf_vis_absorption_final) <- c(names_output_sample_info,"Leaf_VIS_Spectral_Absorption")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
## Output data

# Modeled Reflectance
output.refl <- data.frame(output.LRT$Spec.Info, output.LRT$mod.Reflectance)
names(output.refl) <- c(names_output_sample_info, paste0("Wave_",seq(Start.wave,End.wave,1)))
write.csv(output.refl,file=file.path(out.dir,'PROSPECT5_Modeled_Reflectance.csv'),
          row.names=FALSE)

output.resids <- data.frame(output.LRT$Spec.Info,output.LRT$mod.Reflectance-output.LRT$obs.Reflectance)
names(output.resids) <- c(names_output_sample_info, paste0("Wave_",seq(Start.wave,End.wave,1)))
write.csv(output.resids,file=file.path(out.dir,'PROSPECT5_Reflectance_Residuals.csv'),
          row.names=FALSE)

# Modeled Transmittance
output.trans <- data.frame(output.LRT$Spec.Info, output.LRT$mod.Transmittance)
names(output.trans) <- c(names_output_sample_info, paste0("Wave_",seq(Start.wave,End.wave,1)))
write.csv(output.trans,file=file.path(out.dir,'PROSPECT5_Modeled_Transmittance.csv'),
          row.names=FALSE)

# Modeled Absorption
output.abs <- data.frame(output.LRT$Spec.Info, full_spectrum_absorption)
names(output.abs) <- c(names_output_sample_info, paste0("Wave_",seq(Start.wave,End.wave,1)))
write.csv(output.abs,file=file.path(out.dir,'PROSPECT5_Modeled_Full_Spec_Absorption.csv'),
          row.names=FALSE)

write.csv(leaf_vis_absorption_final,file=file.path(out.dir,'PROSPECT5_Estimated_Leaf_VIS_Absorption.csv'),
          row.names=FALSE)

# Modeled params
output.params <- data.frame(output.LRT$Spec.Info,mod.params)
names(output.params) <-c(names_output_sample_info, names(mod.params))
write.csv(output.params,file=file.path(out.dir,'PROSPECT5_Modeled_Parameters.csv'),
          row.names=FALSE)
#--------------------------------------------------------------------------------------------------#
### EOF