

devtools::install_github("james-thorson/FishLife", ref="add-recruitment")

library(FishLife)
vignette("tutorial","FishLife")

# REPLACE THIS WITH LOADING Return.RData ON YOUR MACHINE
load( "D:/UW Hideaway (SyncBackFree)/AFSC/2018-10 -- steepness for pollock/Return.RData" )
# REPLACE THIS WITH LOADING Return.RData ON YOUR MACHINE

# Look at estimates ...
params = matrix(c("K", "M", "Winfinity", "Loo", "tmax", "tm", "Lm", "Temperature", "tm", "rho", "ln_margsd", "h"), ncol=2, byrow=TRUE)

# ... for pollock
#Cov_pred = Plot_taxa( Search_species(Genus="Theragra",Species="chalcogramma",add=FALSE)$match_taxonomy, params=params,
#  Cov_gjj=Return$Cov_gvv, Mean_gj=Return$beta_gv, ParentChild_gz=Return$ParentChild_gz, Y_ij=Return$Y_ij, mfrow=c(3,2) )

# ... for Pacific cod
Cov_pred = Plot_taxa( Search_species(Genus="Gadus",Species="predictive",add=FALSE)$match_taxonomy, params=params,
  Cov_gjj=Return$Cov_gvv, Mean_gj=Return$beta_gv, ParentChild_gz=Return$ParentChild_gz, Y_ij=Return$Y_ij, mfrow=c(3,2) )


# Mean steepness
Cov_pred[[1]]$Mean_pred['h']
sqrt( Cov_pred[[1]]$Cov_pred['h','h'] )

# Mean SigmaR
Cov_pred[[1]]$Mean_pred['ln_margsd']
sqrt( Cov_pred[[1]]$Cov_pred['ln_margsd','ln_margsd'] )

