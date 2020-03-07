# jonashaslbeck@gmail.com; October 2018

mainDir <- "/Volumes/Macintosh HD 2/Dropbox/VAR_vs_AR/3_code/VAR_vs_AR_code/"

# -----------------------------------------------------------------------------------
# --------- Aux Functions -----------------------------------------------------------
# -----------------------------------------------------------------------------------

# Evaluation function

f_eval <- function(A) {
  
  # Stationary?
  eig <- eigen(A)
  eig_val <- abs(eig$values)
  stable <- ifelse(all(eig_val < 1), 1, 0)
  
  # Ratio
  R <- mean(abs(diag(A))) / mean(abs(A[upper.tri(A) | lower.tri(A)]))
  
  # Exp
  D <- mean(eig_val)
  
  # Exp
  eig_val_scale <- eig_val - min(eig_val)
  eig_val_scale <- eig_val_scale / max(eig_val_scale)
  realDim <- eig_val 
  realDim_sc <- var(eig_val_scale)
  
  # ADet
  Adet <- mean(abs(diag(A)))
  
  return(c(stable, R, D, realDim, realDim_sc, Adet))
  
} # eoF


# -----------------------------------------------------------------------------------
# --------- Mapping Out the Grid ----------------------------------------------------
# -----------------------------------------------------------------------------------

Laura_mm <- readRDS("files/MM_MindMaastricht.RDS")

nIter <- 10000
r_A <- list()
set.seed(5)
for(i in 1:nIter) r_A[[i]] <- apply(Laura_mm, 1:2, function(x) rnorm(1, x[1], x[2]))

# Evaluate
out <- do.call(rbind, lapply(r_A, f_eval))
out <- as.data.frame(out)
colnames(out) <- c("Stable", "R", "D")
out_stable <- out[out[, 1] == 1, 2:3]
saveRDS(out_stable, file = "files/out_stable.RDS")


# Extract Cells with at least one sample
n_cells <- (13-1)*(19-1)
m_cells <- matrix(NA, nrow = n_cells, ncol = 3)
x_axis <- seq(0, 9, length=18+1)
y_axis <- seq(0, .6, length=12+1)
m_cell_Bounds <- matrix(NA, n_cells, ncol=4)

counter <- 1
for(i in 1:(19-1)) {
  for(j in 1:(13-1)) {
    
    # at least one sampled?
    ind_x <- out_stable[, 1] >= x_axis[i] & out_stable[, 1] <= x_axis[i+1]
    ind_y <- out_stable[, 2] >= y_axis[j] & out_stable[, 2] <= y_axis[j+1]
    
    m_cells[counter, 3] <- any(ind_x & ind_y)
    
    # compute layout mean of cell for later plotting
    m_cells[counter, 1] <- (x_axis[i] + x_axis[i+1])/2
    m_cells[counter, 2] <- (y_axis[j] + y_axis[j+1])/2
    
    # save cell bounds
    m_cell_Bounds[counter, 1] <- x_axis[i]
    m_cell_Bounds[counter, 2] <- x_axis[i+1]
    m_cell_Bounds[counter, 3] <- y_axis[j]
    m_cell_Bounds[counter, 4] <- y_axis[j+1]
    
    counter <- counter + 1
  }
}

# Subset cells with at least one
m_cell_fill <- m_cells[m_cells[, 3] == 1, ]
saveRDS(m_cell_fill, file="files/cell_positions.RDS")

n_cells <- nrow(m_cell_fill)
n_cells

# Define Cells:
cell_def <- m_cell_Bounds[m_cells[, 3] == 1, ]


# -----------------------------------------------------------------------------------
# --------- Sampling the Grid -------------------------------------------------------
# -----------------------------------------------------------------------------------

# --- Meta Stuff ---

# Define Grid
l_CM <- vector("list", length = n_cells)
target_cell <- 100

# --- Repeat from here ---

set.seed(1)
for(d in 1:10000) { # d Batches
  
  # Simulate a couple
  nIter <- 5000
  
  r_A <- list()
  for(i in 1:nIter) r_A[[i]] <- apply(Laura_mm, 1:2, function(x) rnorm(1, x[1], x[2]))
  
  # Evaluate and sort in
  out_eval <- do.call(rbind, lapply(r_A, f_eval))
  
  for(i in 1:nIter) {
    
    if(out_eval[i, 1] == 1) {
      
      # Compute which rectangle it is (if any)
      horiz <- cell_def[, 1] <= out_eval[i, 2] & out_eval[i, 2] <= cell_def[, 2]
      vert <- cell_def[, 3] <= out_eval[i, 3] & out_eval[i, 3] <= cell_def[, 4]
      
      if(any(horiz & vert)) {
        
        ind_cell <- which(horiz & vert)
        length_ind_cell <- length(l_CM[[ind_cell]])
        if(length_ind_cell < 100) { # only ad dif not full yet
          l_CM[[ind_cell]][[length_ind_cell + 1]] <- r_A[[i]]
        } # end if: full?
        
      } # end if: in any rectangle?
      
    } # end if: stable?
    
  } # end for: sort in
  
  
  # --- Look at Stats ---
  
  sampled_frq <- do.call(c, lapply(l_CM, length))
  
  barplot(sampled_frq, ylim=c(0, 100))
  abline(h=100)
  
  print(d)
  if(all(sampled_frq == 100)) {
    print("FINISHED")
    break
  }
  
} # end for: 200 x 5000 batches

saveRDS(l_CM, file = "Simulation_3_Variation_2Dim/l_CM_60cells.RDS")


