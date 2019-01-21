# jonashaslbeck@gmail.com; November 2018

# Get iter from command line ...
iter <- commandArgs(trailingOnly=TRUE)
iter <- as.numeric(iter)

# -----------------------------------------------------------------------------
# ---------- Load Packages ----------------------------------------------------
# -----------------------------------------------------------------------------

# Parallelize
library(foreach)
library(doParallel)

# Data Generation
library(mlVAR)

source("aux_functions.R")


# -----------------------------------------------------------------------------
# ---------- Load Models ------------------------------------------------------
# -----------------------------------------------------------------------------

l_CM <- readRDS(file = "Models_60cells.RDS")
# l_CM <- readRDS(file = "/Volumes/Macintosh HD 2/Dropbox/VAR_vs_AR/3_code/VAR_vs_AR_code/Simulation_3_Variation_2Dim/l_CM_60cells.RDS")

# Fixed stuff
p <- 6
means <- rep(0, p) # mean vector
E <- matrix(0, p, p) # error covariance matrix
diag(E) <- 1

# n variation
n_seq <- 8:500

# -----------------------------------------------------------------------------
# ---------- Loop -------------------------------------------------------------
# -----------------------------------------------------------------------------

nClust <- 12 # We run the n-variations in parallel

# Fire up foreach
cl <- makeCluster(nClust)
registerDoParallel(cl)

# start clustered loops
outlist <- foreach(cell=1:60,
                   .packages=c("mlVAR"),
                   .verbose=FALSE,
                   .export=c("l_CM", "fSimulate", "fEstimate", "means", "E", "n_seq")) %dopar% {

                     # ----- Create Output List -----
                     
                     seed <- iter * 1000 + cell # unique seed from iter and cell
                     set.seed(seed)
                     
                     outlist <- fSimulate(Ar = l_CM[[cell]][[iter]],
                                          means = means,
                                          E = E,
                                          nIter = 100, # before 1
                                          verbatim = FALSE,
                                          n_seq = n_seq)
                     
                     # ----- Return -----
                     
                     return(outlist)
                     
                   } # end: foreach (n)


# Close down foreach
stopCluster(cl)


# --------- Save Results ---------

saveRDS(outlist, file = paste0('ARVAR_no3_Iter_', iter, '.RDS'))



