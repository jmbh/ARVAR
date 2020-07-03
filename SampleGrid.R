# jonashaslbeck@gmail.com; May 2020

mainDir <- "Users/jonas/Dropbox/VAR_vs_AR/ARVAR/"

# -----------------------------------------------------------------------------------
# --------- Aux Functions -----------------------------------------------------------
# -----------------------------------------------------------------------------------

# Evaluation function

f_eval <- function(A) {
  
  # Stationary?
  eig <- eigen(A)
  eig_val <- abs(eig$values)
  stable <- ifelse(all(eig_val < 1), 1, 0)
  
  # Diagonal
  D <- mean(abs(diag(A)))
  
  # Off-diagonal
  OffD <- mean(abs(A[upper.tri(A) | lower.tri(A)]))
  
  return(c(stable, D, OffD))
  
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
colnames(out) <- c("Stable", "D", "OffD")
table(out$Stable)
out_stable <- out[out[, 1] == 1, 2:3]
saveRDS(out_stable, file = "files/out_stable.RDS")

n_box <- 15



# Extract Cells with at least one sample
n_cells <- (n_box-1)^2
m_cells <- matrix(NA, nrow = n_cells, ncol = 3)
m_cell_Bounds <- matrix(NA, n_cells, ncol=4)

counter <- 1
for(i in 1:(n_box-1)) {
  for(j in 1:(n_box-1)) {
    
    # at least one sampled?
    ind_x <- out_stable[, 1] >= x_seq[i] & out_stable[, 1] <= x_seq[i+1]
    ind_y <- out_stable[, 2] >= y_seq[j] & out_stable[, 2] <= y_seq[j+1]
    
    m_cells[counter, 3] <- any(ind_x & ind_y)
    
    # compute layout mean of cell for later plotting
    m_cells[counter, 1] <- (x_seq[i] + x_seq[i+1])/2
    m_cells[counter, 2] <- (y_seq[j] + y_seq[j+1])/2
    
    # save cell bounds
    m_cell_Bounds[counter, 1] <- x_seq[i]
    m_cell_Bounds[counter, 2] <- x_seq[i+1]
    m_cell_Bounds[counter, 3] <- y_seq[j]
    m_cell_Bounds[counter, 4] <- y_seq[j+1]
    
    counter <- counter + 1
  }
}

# Look at results
table(m_cells[, 3]) # about 1/3 have at least one 

# Subset cells with at least one
m_cell_fill <- m_cells[m_cells[, 3] == 1, ]
saveRDS(m_cell_fill, file="files/cell_positions.RDS")

n_cells <- nrow(m_cell_fill)
n_cells

# Define Cells:
cell_def <- m_cell_Bounds[m_cells[, 3] == 1, ]


# Flag nonempty cells

# Make a Grid
pdf("figures/App_SampledGrid_empty.pdf", height=6.2, width=6)

plot.new()
plot.window(xlim=c(0, .5), ylim=c(0, .2))
axis(1)
axis(2, las=2)
title(xlab = "Mean Absolute Diagonal")
title(ylab = "Mean Absolute Off-Diagonal")

points(out$D, out$OffD, pch = 20, col = adjustcolor('black', alpha = .5))

x_seq <- seq(0, .5, length=n_box)
y_seq <- seq(0, .2, length=n_box)

segments(x_seq, rep(0, n_box), x_seq, rep(.2, n_box))
segments(rep(0, n_box), y_seq, rep(.5, n_box), y_seq)

for(i in 1:74) text(m_cell_fill[i, 1], m_cell_fill[i,2], i, col="red")

dev.off()




# Plot the nonempty cells
# Make a Grid
pdf("figures/App_SampledGrid_empty_nonzero.pdf", height=6.2, width=6)

plot.new()
plot.window(xlim=c(0, .5), ylim=c(0, .2))
axis(1)
axis(2, las=2)
title(xlab = "Mean Absolute Diagonal")
title(ylab = "Mean Absolute Off-Diagonal")
# points(out$D, out$OffD)

x_seq <- seq(0, .5, length=n_box)
y_seq <- seq(0, .2, length=n_box)

x_length <- .5/n_box / 2
y_length <- .2/n_box / 2

rect(xleft = m_cell_fill[, 1] - x_length, 
     ybottom = m_cell_fill[, 2] - y_length, 
     xright = m_cell_fill[, 1] + x_length, 
     ytop = m_cell_fill[, 2] + y_length)

dev.off()


# -----------------------------------------------------------------------------------
# --------- Sampling the Grid -------------------------------------------------------
# -----------------------------------------------------------------------------------

# The goal here is to sample 100 models into each of the 74 non-empty cells above

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
        if(length_ind_cell < 100) { # only add if not full yet
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

saveRDS(l_CM, file = "files/Models_74cells.RDS")
