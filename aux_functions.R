# jonashaslbeck@gmail.com; August 2018


# -----------------------------------------------------------------------------------
# --------- Aux Functions 1 -----------------------------------------------------------
# -----------------------------------------------------------------------------------

# --------- Simulate VAR ----------

simVAR <- function(pars, 
                   means = 0, 
                   lags = 1, Nt = 100, init, residuals = 0.1, 
                   burnin) 
{
  
  print(pars)
  
  if (is.matrix(pars)) 
    pars <- list(pars)
  if (any(sapply(pars, function(x) length(unique(dim(x))) > 
                 1))) 
    stop("non-square graph detected.")
  if (missing(burnin)) {
    burnin <- min(round(Nt/2), 100)
  }
  # Ni <- ncol(pars[[1]])
  
  residuals <- 1
  Ni <- 6
  
  
  if (length(means) == 1) 
    means <- rep(means, Ni)
  maxLag <- max(lags)
  if (length(residuals) == 1) {
    residuals <- diag(residuals, Ni)
  }
  else if (length(residuals) == Ni) {
    residuals <- diag(residuals)
  }
  if (!is.matrix(residuals) && ncol(residuals) != Ni && nrow(residuals) != 
      Ni) {
    stop("'residuals' is not a square matrix")
  }
  totTime <- Nt + burnin
  if (missing(init)) {
    init <- matrix(0, maxLag, Ni)
  }
  Res <- matrix(, totTime, Ni)
  Res[1:maxLag, ] <- init
  for (t in (maxLag + 1):(totTime)) {
    Res[t, ] <- means + rowSums(do.call(cbind, lapply(seq_along(lags), 
                                                      function(i) pars[[i]] %*% (Res[t - lags[i], ] - means)))) + 
      mvtnorm::rmvnorm(1, rep(0, Ni), residuals)
  }
  return(as.data.frame(Res[-(1:burnin), ]))
}



# --------- Evaluate characteristics of VAR matrix ----------

f_eval <- function(A) {
  
  # Stationary?
  eig <- eigen(A)
  eig_val <- abs(eig$values)
  stable <- ifelse(all(eig_val < 1), 1, 0)
  
  # Ratio
  R <- mean(abs(diag(A))) / mean(abs(A[upper.tri(A) | lower.tri(A)]))
  
  # Exp
  D <- mean(eig_val)
  
  return(c(stable, R, D))
  
} # eoF


f_eval2 <- function(A) {
  
  # Ratio
  eig_val <- eigen(A)$values
  R <- mean(abs(diag(A))) / mean(abs(A[upper.tri(A) | lower.tri(A)]))
  
  # Dim
  D <- mean(eig_val)
  
  # Var eigen
  VarEig <- var(eig_val)
  
  # Trace abs value
  AbsTr <- sum(abs(diag(A)))
  
  out <- c(R, D, VarEig, AbsTr)
  names(out) <- c("R", "D", "VarEig", "AbsTr")
  
  return(out)
  
} # eoF



# --------- Compute the Gap ----------

get_gap_diff <- function(pe_diff, ee_diff, n_seq = 8:200) {
  n_pe <- which(pe_diff > 0)[1]
  n_ee <- which(ee_diff > 0)[1]
  
  gap <- n_seq[n_ee] - n_seq[n_pe]
  gap
}

# --------- ? ----------

compute_average <- function(error_mat) {
  avg <- matrix(NA, ncol = 2, nrow = n_len)
  
  for (i in seq(2)) {
    for (j in seq(n_len)) {
      sum_j <- 0
      
      for (k in seq(n_rep)) {
        sum_j <- sum_j + error_mat[[k]][j, i]
      }
      
      avg[j, i] <- sum_j / n_rep
    }
  }
  
  avg
}


# --------- Computes average Gap for all Cells ----------

#' takes simulation output, returns list of 600.000 estimation error differences
get_all_ee_curves <- function(out_var, cells = seq(60), cols = 1:193) {
  nsim <- length(cells) * 100 * 100
  ee_diff <- matrix(NA, nrow = nsim, ncol = length(cols))
  
  index <- 1
  
  for (cell in cells) {
    for (model in seq(100)) {
      for (iter in seq(100)) {
        m <- out_var[[cell]][[model]]
        
        out <- m$EE_all[, , iter]
        ee_diff[index, ] <- out[cols, 1] - out[cols, 2]
        
        index <- index + 1
      }
    }
  }
  
  ee_diff
}


# Input: Simulation Results
# Output: average gap for each Cell
compute_diff_first <- function(out_var) {
  
  diff_first <- numeric(n_cells)
  
  for (cell in seq(n_cells)) {
    pe <- lapply(out_var[[cell]], function(x) x$PE[, 1:2])
    ee <- lapply(out_var[[cell]], function(x) x$EE[, 1:2])
    
    # compute the differences across all ns
    arvar_diff <- sapply(seq(n_rep), function(i) {
      pe_diff <- pe[[i]][, 1] - pe[[i]][, 2] # AR - VAR (PE)
      ee_diff <- ee[[i]][, 1] - ee[[i]][, 2] # AR - VAR (EE)
      
      gap <- get_gap_diff(pe_diff, ee_diff)
      gap
    })
    
    # compute the mean across differences
    diff_first[cell] <- round(mean(arvar_diff, na.rm = TRUE), 2)
  }
  
  diff_first
}


#' Weights the estimation error difference curves by resampling them according to their probability (based on cells)
get_cell_weighted_samples <- function(cell_weights, diffs, n_cells = 60, seed = 1) {
  nvar <- 10000
  draws <- round(cell_weights * nvar)
  
  # if the following is uncommented, renormalizes to sample only 100 per cell
  # x <- cell_weights * nvar
  # draws <- round((x / max(x)) * 100)
  samples <- c()
  
  set.seed(seed)
  for (i in seq(n_cells)) {
    d <- draws[i]
    
    if (d != 0) {
      lines <- diffs[[i]][sample(seq(100), size = d, replace = TRUE), ]
      samples <- rbind(samples, lines)
    }
    
  }
  
  samples
}


#' Weights the estimation error difference curves by resampling them according to their probability (based on the mixed model)
get_weighted_samples <- function(diffs_mat, model_weights, seed = 1) {
  set.seed(seed)
  nr_models <- nrow(diffs_mat)
  index <- sample(nr_models, replace = TRUE, size = nr_models, prob = model_weights)
  diffs_mat[index, ]
}


# -------- Compute Weight of each Model --------

get_model_weights <- function() {
  library('mvtnorm')
  l_CM <- readRDS(file = "files/Models_60cells.RDS")
  mm <- readRDS(file="files/MM_MindMaastricht.RDS")
  
  M <- as.numeric(mm[, , 1])
  S <- diag(as.numeric(mm[, , 2]))
  
  counter <- 1
  weights <- rep(NA, 60*100)
  
  for(d in 1:60) {
    for(i in 1:100) {
      
      x <- as.numeric(l_CM[[d]][[i]])
      
      weights[counter] <- dmvnorm(x, mean = M, sigma = S, log = FALSE)
      
      counter <- counter + 1
    }
  }
  
  weights / sum(weights)
}


# -----------------------------------------------------------------------------------
# --------- Aux Functions 2 -----------------------------------------------------------
# -----------------------------------------------------------------------------------

# --------- Shading for SEs ----------

library(scales)

shadeSE <- function(x, y_up, y_low, color, alpha) {
  
  col <- alpha(colour = color, alpha = alpha)
  n <- length(x)
  
  n <- length(x)
  p_x <- c(x, x[n:1])
  p_y <- c(y_up, y_low[n:1])
  polygon(x = p_x, y = p_y, col=col, border = NA)
  
} # eoF


# --------- Estimation Function ----------

fEstimate <- function(train, test, model = "AR", trueVAR) {
  
  # Aux variables
  p <- ncol(train)
  n_train <- nrow(train)
  n_test <- nrow(test)
  
  # Lag data
  train_dep <- as.matrix(train[-1, ])
  train_ind <- as.matrix(train[-n_train, ])
  test_dep <- as.matrix(test[-1, ])
  test_ind <- as.matrix(test[-n_test, ])
  
  # Storage
  m_predictions <- matrix(NA, nrow = n_test-1, ncol = p)
  m_VAR <- matrix(0, p, p)
  
  # Fit model
  if(model == "AR") {
    for(i in 1:p) {
      model_i <- lm(train_dep[, i] ~ train_ind[, i]) # fit
      m_VAR[i, i] <- model_i$coefficients[2] # save in parameter matrix
      m_predictions[, i] <- model_i$coefficients[1] + model_i$coefficients[2]*test_ind[, i] # save predictions
    }
  }
  
  if(model == "VAR") {
    for(i in 1:p) {
      model_i <- lm(train_dep[, i] ~ train_ind) # fit
      m_VAR[i, ] <- model_i$coefficients[-1] # save in parameter matrix
      m_predictions[, i] <- model_i$coefficients[1] + test_ind %*% matrix(model_i$coefficients[-1], nrow=p)  # save predictions
    }
  }
  
  if(model == "cAR") {
    for(i in 1:p) {
      model_i <- lm(train_dep[, i] ~ train_ind) # fit
      m_VAR[i, ] <- model_i$coefficients[-1] # save in parameter matrix
      m_VAR[i, -i] <- 0
      model_i$coefficients[-1][-i] <- 0
      m_predictions[, i] <- model_i$coefficients[1] + test_ind %*% matrix(model_i$coefficients[-1], nrow=p)  # save predictions
    }
  }
  
  # browser()
  
  ## Compute errors
  # Estimation Error
  MSE_EE <- mean((m_VAR - trueVAR)^2)
  MSE_EE_false <- mean((diag(m_VAR) - diag(trueVAR))^2)
  
  # Prediction Error: MSE
  v_MSE_PE <- rep(NA, p)
  for(i in 1:p) v_MSE_PE[i] <- mean((test_dep[, i] - m_predictions[, i])^2)
  MSE_PE <- mean(v_MSE_PE)
  
  # Prediction Error: R2
  v_R2_PE <- rep(NA, p)
  for(i in 1:p) v_R2_PE[i] <- 1 - var(m_predictions[, i]-test_dep[, i]) / var(test_dep[, i])
  
  
  # Output
  outlist <- list("estVAR" = m_VAR, 
                  "predsTest" = m_predictions, 
                  "trueTest" = test_dep, 
                  "MSE_PE" = MSE_PE,
                  "MSE_EE" = MSE_EE, 
                  "MSE_EE_false" = MSE_EE_false, 
                  "R2_PE" = v_R2_PE)
  
} # eoF

# # Testing
# train <- data[1:1500, ]
# test <- data[1501:2000, ]
# 

# out <- fEstimate(train = train, 
#                  test = test, 
#                  model = "AR", 
#                  trueVAR = Ar)


# --------- Simulation Function ----------


fSimulate <- function(Ar, 
                      means, 
                      E, 
                      nIter, 
                      verbatim = TRUE, 
                      n_seq) 
  
{
  
  l_data <- list()
  
  # Generate Testset
  testset <- simVAR(pars = Ar, 
                    means = means, 
                    lags = 1, 
                    Nt = 2000,
                    residuals = E)
  
  # Storage
  n_var <- length(n_seq)
  a_PE <- a_EE <- a_EEf <- array(NA, dim = c(n_var, 3, nIter))
  
  
  for(i in 1:nIter) {
    
    # Generate Data
    data <- simVAR(pars = Ar, 
                   means = means, 
                   lags = 1, 
                   Nt = n_seq[length(n_seq)],
                   residuals = E)

    l_data[[i]] <- data
    
    for(n in 1:n_var) {
      
      # ----- Fit AR model & make predictions -----
      out_AR <- fEstimate(train = data[1:n_seq[n], ], 
                          test = testset, 
                          model = "AR", 
                          trueVAR = Ar)
      
      a_PE[n, 1, i] <- out_AR$MSE_PE
      a_EE[n, 1, i] <- out_AR$MSE_EE
      a_EEf[n, 1, i] <- out_AR$MSE_EE_false
      
      
      # ----- Fit VAR model & make predictions -----
      out_VAR <- fEstimate(train = data[1:n_seq[n], ], 
                           test = testset, 
                           model = "VAR", 
                           trueVAR = Ar)
      
      a_PE[n, 2, i] <- out_VAR$MSE_PE
      a_EE[n, 2, i] <- out_VAR$MSE_EE
      
      # ----- Fit VAR / cAR model & make predictions -----
      out_cAR <- fEstimate(train = data[1:n_seq[n], ], 
                           test = testset, 
                           model = "cAR", 
                           trueVAR = Ar)
      
      a_PE[n, 3, i] <- out_cAR$MSE_PE
      a_EE[n, 3, i] <- out_cAR$MSE_EE
      
      
      if(verbatim) cat(paste0("Iteration = ", i, " n = ", n, "\n"))
    } # for: loop n
    
  } # for: iterations
  
  # --------- Aggregate over Iteration ----------
  
  # browser()
  
  # Means
  m_PE <- apply(a_PE, 1:2, median)
  m_EE <- apply(a_EE, 1:2, median)
  
  # Upper / Lower 25% quantile
  PE_qu <- apply(a_PE, 1:2, function(x) quantile(x, probs = c(.25, .75)))
  EE_qu <- apply(a_EE, 1:2, function(x) quantile(x, probs = c(.25, .75)))
  
  
  # --------- Return Results ----------
  
  outlist <- list("PE" = m_PE, 
                  "PE_qu" = PE_qu,
                  "EE" = m_EE, 
                  "EE_qu" = EE_qu, 
                  "PE_all" = a_PE, 
                  "EE_all" = a_EE, 
                  "EEf_all" = a_EEf)
  
  
} # end of Simulation Func


# removes unneeded stuff that costs too much memory to store
prune_simoutput <- function(ARVAR_iter, len = 300) {
  x <- ARVAR_iter
  
  for (i in seq(60)) {
    x[[i]]$EE_qu <- NULL
    x[[i]]$PE_qu <- NULL
    x[[i]]$EEf_all <- NULL
  }
  
  x
}

get_weighted_sd <- function(x, weights) {
  n <- length(x)
  
  weighted_mean <- sum(x * weights, na.rm = TRUE)
  num <- sum(weights * (x - weighted_mean)^2, na.rm = TRUE)
  denom <- (n - 1)/n * sum(weights, na.rm = TRUE)
  
  sqrt(num / denom)
}

plot_weights <- function(cell = 1) {
  lo <- (cell - 1) * 100
  hi <- cell * 100
  
  w <- model_weights[seq(lo, hi)]
  hist(w, breaks = 20)
  abline(v = mean(w), lwd = 2, col = 'green')
  abline(v = cell_weights[cell] / 100, lwd = 2)
}

# diff <- numeric(60)
# 
# for (d in seq(60)) {
#   lo <- (d - 1) * 100
#   hi <- d * 100
#   
#   w <- mean(model_weights[seq(lo, hi)])
#   cw <- cell_weights[d] / 100
#   
#   diff[d] <- w - cw
# }




#' Adapted from plotfunctions::gradientLegend
.gradientLegend <- function (valRange, color = "topo", nCol = 30, pos = 0.5, side = 4, 
                             length = 0.25, depth = 0.05, inside = TRUE, coords = FALSE, 
                             pos.num = NULL, n.seg = 3, border.col = "black", dec = NULL, 
                             fit.margin = TRUE, cex = 0.8) 
{
  loc <- c(0, 0, 0, 0)
  if (is.null(pos.num)) {
    if (side %in% c(1, 3)) {
      pos.num = 3
    }
    else {
      pos.num = side
    }
  }
  if (length(pos) == 1) {
    pos.other <- ifelse(side > 2, 1, 0)
    if (side %in% c(1, 3)) {
      switch <- ifelse(inside, 0, 1)
      switch <- ifelse(side > 2, 1 - switch, switch)
      loc <- getCoords(c(pos - 0.5 * length, pos.other - 
                           switch * depth, pos + 0.5 * length, pos.other + 
                           (1 - switch) * depth), side = c(side, 2, side, 
                                                           2))
    }
    else if (side %in% c(2, 4)) {
      switch <- ifelse(inside, 0, 1)
      switch <- ifelse(side > 2, 1 - switch, switch)
      loc <- getCoords(c(pos.other - switch * depth, pos - 
                           0.5 * length, pos.other + (1 - switch) * depth, 
                         pos + 0.5 * length), side = c(1, side, 1, side))
    }
  }
  else if (length(pos) == 4) {
    if (coords) {
      loc <- pos
    }
    else {
      loc <- getCoords(pos, side = c(1, 2, 1, 2))
    }
  }
  mycolors <- c()
  if (length(color) > 1) {
    mycolors <- color
  }
  else if (!is.null(nCol)) {
    if (color == "topo") {
      mycolors <- topo.colors(nCol)
    }
    else if (color == "heat") {
      mycolors <- heat.colors(nCol)
    }
    else if (color == "terrain") {
      mycolors <- terrain.colors(nCol)
    }
    else if (color == "rainbow") {
      mycolors <- rainbow(nCol)
    }
    else {
      warning("Color %s not recognized. A palette of topo.colors is used instead.")
      mycolors <- topo.colors(nCol)
    }
  }
  else {
    stop("No color palette provided.")
  }
  vals <- seq(min(valRange), max(valRange), length = length(mycolors))
  if (!is.null(dec)) {
    vals <- round(vals, dec[1])
  }
  im <- as.raster(mycolors[matrix(1:length(mycolors), ncol = 1)])
  ticks <- c()
  if (side%%2 == 1) {
    rasterImage(t(im), loc[1], loc[2], loc[3], loc[4], col = mycolors, 
                xpd = T)
    rect(loc[1], loc[2], loc[3], loc[4], border = border.col, 
         xpd = T)
    ticks <- seq(loc[1], loc[3], length = n.seg)
    segments(x0 = ticks, x1 = ticks, y0 = rep(loc[2], n.seg), 
             y1 = rep(loc[4], n.seg), col = border.col, xpd = TRUE)
  }
  else {
    rasterImage(rev(im), loc[1], loc[2], loc[3], loc[4], 
                col = mycolors, xpd = T)
    rect(loc[1], loc[2], loc[3], loc[4], border = border.col, 
         xpd = T)
    ticks <- seq(loc[2], loc[4], length = n.seg)
    segments(x0 = rep(loc[1], n.seg), x1 = rep(loc[3], n.seg), 
             y0 = ticks, y1 = ticks, col = border.col, xpd = TRUE)
  }
  determineDec <- function(x) {
    out = max(unlist(lapply(strsplit(x, split = "\\."), function(y) {
      return(ifelse(length(y) > 1, nchar(gsub("^([^0]*)([0]+)$", 
                                              "\\1", as.character(y[2]))), 0))
    })))
    return(out)
  }
  labels = sprintf("%f", seq(min(valRange), max(valRange), 
                             length = n.seg))
  if (is.null(dec)) {
    dec <- min(c(6, determineDec(labels)))
  }
  eval(parse(text = sprintf("labels = sprintf('%s', round(seq(min(valRange), max(valRange), length = n.seg), dec) )", 
                            paste("%.", dec, "f", sep = ""))))
  if (pos.num == 1) {
    if (fit.margin) {
      lab.height = max(strheight(labels)) * 0.8
      max.pos = getFigCoords()[3]
      if ((max.pos - loc[2]) < lab.height) {
        warning("Increase bottom margin, because labels for legend do not fit.")
      }
    }
    text(y = loc[2], x = ticks, labels = seq(min(valRange), 
                                             max(valRange), length = n.seg), col = border.col, 
         pos = 1, cex = 0.8, xpd = T)
  }
  else if (pos.num == 2) {
    if (fit.margin) {
      checkagain = TRUE
      while (checkagain == TRUE) {
        lab.width = (max(strwidth(labels)) + 0.5 * par()$cxy[1]) * 
          0.8
        min.pos = getFigCoords()[1]
        if ((loc[1] - min.pos) < lab.width) {
          if (!is.null(dec)) {
            dec = max(c(0, dec - 1))
            if (dec == 0) {
              warning("Decimal set to 0 (dec=0), but the labels still don't fit in the margin. You may want to add the color legend to another side, or increase the margin of the plot.")
              checkagain = FALSE
            }
          }
          else {
            tmp = max(unlist(lapply(strsplit(labels, 
                                             split = "\\."), function(x) {
                                               return(ifelse(length(x) > 1, nchar(x[2]), 
                                                             0))
                                             })))
            dec = max(c(0, tmp - 1))
            if (dec == 0) {
              warning("Decimal set to 0 (dec=0), but the labels still don't fit in the margin. You may want to add the color legend to another side, or increase the margin of the plot.")
              checkagain = FALSE
            }
          }
          eval(parse(text = sprintf("labels = sprintf('%s', round(seq(min(valRange), max(valRange), length = n.seg), dec) )", 
                                    paste("%.", dec, "f", sep = ""))))
        }
        else {
          checkagain = FALSE
        }
      }
    }
    text(y = ticks, x = loc[1], labels = labels, pos = 2, 
         cex = 0.8, col = border.col, xpd = T)
  }
  else if (pos.num == 3) {
    if (fit.margin) {
      lab.height = max(strheight(labels)) * 0.8
      max.pos = getFigCoords()[4]
      if ((max.pos - loc[4]) < lab.height) {
        warning("Increase top margin, because labels for legend do not fit.")
      }
    }
    text(y = loc[4], x = ticks, labels = seq(min(valRange), 
                                             max(valRange), length = n.seg), col = border.col, 
         pos = 3, cex = 0.8, xpd = T)
  }
  else if (pos.num == 4) {
    if (fit.margin) {
      checkagain = TRUE
      while (checkagain == TRUE) {
        lab.width = (max(strwidth(labels)) + 0.5 * par()$cxy[1]) * 
          0.8
        max.pos = getFigCoords()[2]
        if ((max.pos - loc[3]) < lab.width) {
          if (!is.null(dec)) {
            dec = max(c(0, dec - 1))
            if (dec == 0) {
              warning("Decimal set to 0 (dec=0), but the labels still don't fit in the margin. You may want to add the color legend to another side, or increase the margin of the plot.")
              checkagain = FALSE
            }
          }
          else {
            tmp = max(unlist(lapply(strsplit(labels, 
                                             split = "\\."), function(x) {
                                               return(ifelse(length(x) > 1, nchar(x[2]), 
                                                             0))
                                             })))
            dec = max(c(0, tmp - 1))
            if (dec == 0) {
              warning("Decimal set to 0 (dec=0), but the labels still don't fit in the margin. You may want to add the color legend to another side, or increase the margin of the plot.")
              checkagain = FALSE
            }
          }
          eval(parse(text = sprintf("labels = sprintf('%s', round(seq(min(valRange), max(valRange), length = n.seg), dec) )", 
                                    paste("%.", dec, "f", sep = ""))))
        }
        else {
          checkagain = FALSE
        }
      }
    }
    text(y = ticks, x = loc[3] - .2, labels = labels, pos = 4, 
         col = border.col, cex = cex, xpd = T)
  }
}
