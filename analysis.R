source('aux_functions.R')

l_CM <- readRDS(file = "files/l_CM_60cells.RDS")
out_var <- readRDS('~/Desktop/ARVAR/out_VAR_500.RDS')

all_matrices <- vector('list', length = 6000)

ix <- 1
for (i in seq(60)) {
  for (j in seq(100)) {
    all_matrices[[ix]] <- l_CM[[i]][[j]]
    ix <- ix + 1
  }
}

models6000_100 <- list()

j <- 1
for (cell in seq(n_cells)) {
  for (iter in seq(100)) {
    models6000_100[[j]] <- out_var[[cell]][[iter]]
    j <- j + 1
  }
}

diff_first <- compute_diff_first(out_var)

ngaps_models <- sapply(models6000_100, function(m) {
  ee <- m$EE_all
  pe <- m$PE_all
  gaps <- sapply(seq(100), function(i) {
    ee_diff <- ee[, , i][, 1] - ee[, , i][, 2]
    pe_diff <- pe[, , i][, 1] - pe[, , i][, 2]
    get_gap_diff(pe_diff, ee_diff, n_seq = n_seq)
  })
  
  median(gaps)
})



X <- lapply(all_matrices, f_eval2)
X <- abs(do.call('rbind', X))

dat <- data.frame(gaps = ngaps_models, X)
dat <- na.omit(dat)

m <- lm(gaps ~ R + D, data = dat)
m1 <- lm(gaps ~ R * D, data = dat)

