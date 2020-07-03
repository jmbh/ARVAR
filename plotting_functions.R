# dablander.fabian@gmail.com June 2020

library(shape)
library(weights)
library(latex2exp)
library(RColorBrewer)
library(plotfunctions)
source('aux_functions.R')


# -------------------------------------------------------------
# ---------- Figure 1: plot_fig123 for Panels (a), (b), and (c)
# -------------------------------------------------------------
plot_fig123 <- function(
  diffs_mat, model_weights = NULL, figtype = 1, n_cells = 74, n_seq = 25:200,
  ylim = c(-.052, .052), main = '', letter = '', cells_coloured = c(48, 125, 412),
  ylab = TeX('$$EE_{AR} - EE_{VAR}$$'), xlab = "Number of observations n", seed = 1
) {
  
  if (figtype == 1) {
    # resample 740.000 according to their probabilities
    diffs_mat <- get_weighted_samples(diffs_mat, model_weights, seed = seed)
  }
  
  skip <- seq(n_seq[1] - 8)
  
  plot.new()
  n_len <- length(n_seq)
  xlim <- c(n_seq[1], n_seq[length(n_seq)])
  
  plot.window(
    xlim = xlim, ylim = ylim, yaxs = 'i'
  )
  x_labels <- scales::pretty_breaks(n = 10)(xlim)
  x_labels <- round(seq(25, 200, 25))
  title(main = main)
  title(xlab = xlab, cex.lab = 1.2)
  title(ylab = ylab, line = 1.8, cex.lab = 1.2)
  
  plot_line <- function(x, y, lo, hi, col, col_shade = 'grey',...) {
    lines(x, y, panel.first = polygon(
      c(n_seq, rev(n_seq)),
      c(lo, rev(hi)),
      col = grDevices::adjustcolor(col_shade, alpha = .4), border = NA
    ), col = col, ...
    )
  }
  
  if (figtype == 1) {
    diffs_mean <- apply(diffs_mat, 2, median)[-skip][seq(n_len)]
    diffs_sd <- apply(diffs_mat, 2, sd)[-skip][seq(n_len)]
    
    plot_line(n_seq, diffs_mean, diffs_mean - diffs_sd, diffs_mean + diffs_sd, col = 'black')
  }
  
  if (!is.null(model_weights)) {
    alpha <- (model_weights) / max(10 * model_weights)
  } else {
    alpha <- rep(.9, n_cells)
  }
  
  colours <- rep('grey', 6000)
  colours[cells_coloured] <- brewer.pal(3, 'Set1') # red, blue, green
  cells <- which(colours != 'grey')
  
  
  if (figtype == 2) {
    
    # for (cell in seq(length(model_weights))) {
    for (cell in seq(nrow(diffs_mat))) {
      ee_curve <- diffs_mat[cell, ][-skip][seq(n_len)]
      col <- colours[cell]
      colour <- ifelse(col == 'grey', adjustcolor(col, alpha = alpha[cell]), col)
      plot_line(n_seq, ee_curve, ee_curve, ee_curve, col = colour)
    }
    
    for (cell in cells) {
      ee_curve <- diffs_mat[cell, ][-skip][seq(n_len)]
      
      colour <- colours[cell]
      plot_line(n_seq, ee_curve, ee_curve, ee_curve, col = colour)
    }
  }
  
  if (figtype == 3) {
    for (cell in cells) {
      cell_models <- diffs_mat[[cell]]
      mean_cell <- apply(cell_models, 2, median)[-skip][seq(n_len)]
      sd_cell <- apply(cell_models, 2, sd)[-skip][seq(n_len)]
      
      col <- colours[cell]
      plot_line(n_seq, mean_cell, mean_cell - sd_cell, mean_cell + sd_cell, col = col, col_shade = col)
    }
  }
  
  # plot axes the last
  # axis(1, c(xlim[1], 50, 100, 150, 200))
  axis(1, x_labels, cex.axis = 1.2)
  axis(2, at = c(-.05, 0, .05), labels = c('', '0', ''), las = 2, cex.axis = 1.2)
  
  lines(c(0, 200), c(0, 0), lty = 2, col = 'grey')
  text(25, 0.048, letter, col = 'black', cex = 1.5, adj = 0)
  
  # draw VAR is better with arrows
  Arrows(x0 = 199, y0 = .017, x1 = 199, y1 = .003, arr.type = 'triangle', arr.length = .2, arr.width = .2)
  text(199, .018, "VAR better", col = 'black', cex = 1, adj = 0, srt = 90)
  Arrows(x0 = 199, y0 = .039, x1 = 199, y1 = .049, arr.type = 'triangle', arr.length = .2, arr.width = .2)
  
  # draw VAR is better with arrows
  Arrows(x0 = 199, y0 = -.015, x1 = 199, y1 = -.003, arr.type = 'triangle', arr.length = .2, arr.width = .2)
  text(199, -.034, "AR better", col = 'black', cex = 1, adj = 0, srt = 90)
  Arrows(x0 = 199, y0 = -.035, x1 = 199, y1 = -.049, arr.type = 'triangle', arr.length = .2, arr.width = .2)
}


# -----------------------------------------------------------------------------
# ---------- Figure 2: plot_gaps and plot_histogram (estimation error) --------
# -----------------------------------------------------------------------------

#' Plot the n at which the estimation error crosses (median across 100 models) for each cell
#' as a function of R and D
plot_gaps <- function(
  m_gap, main, m_cell_fill, n_cells = 74,
  cols = brewer.pal(11, "PiYG"), cex = 1, alpha = 1
) {
  
  # Set up Plot
  plot.new()
  plot.window(xlim=c(0, 0.50), ylim=c(0, 0.2))
  axis(1, at = seq(0, 0.50, .1))
  axis(2, seq(0, .2, 0.05), las=2, cex.axis = 1)
  title(main = main)
  title(ylab = "Mean Absolute Off-Diagonal", cex.lab = 1)
  title(xlab = "Mean Absolute Diagonal", cex.lab = 1)
  
  # Plot the cells
  m_gap_cols <- m_gap + abs(min(m_gap))
  ind_color <- round((m_gap_cols / max(m_gap_cols)) * 9) + 1
  col_gap <- cols[ind_color]
  
  hs_y <- 0.01428571/2
  hs_x <- 0.03571429/2
  
  for (k in seq(n_cells)) {
    rect(xleft = m_cell_fill[k, 1] - hs_x, 
         ybottom = m_cell_fill[k, 2] - hs_y,
         xright = m_cell_fill[k, 1] + hs_x,
         ytop = m_cell_fill[k, 2] + hs_y,
         col = grDevices::adjustcolor(col_gap[k], alpha = alpha)
    )
  }
  
  # Plot the Gap
  text(m_cell_fill[, 1], m_cell_fill[, 2], round(m_gap), cex = cex)
}


#' Plot the sampling distribution of the ns at which the estimation error difference crosses
#' (i.e., VAR outperforms AR); uses all 600.000 models
plot_histogram <- function(nswitch, med, xlim) {
  
  hist(
    nswitch, breaks = 50, main = '',
    xlab = TeX('$$n_{e}$$'), col = 'grey76', ylab = '', cex = 0.9,
    xlim = xlim, xaxt = 'n', las = 2, yaxt = 'n', ylim = c(0, 100000)
  )
  # x_labels <- scales::pretty_breaks(n = 8)(xlim)
  x_labels <- c(8, seq(40, 500, 40))
  x_labels <- c(8, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500)
  axis(1, x_labels, las = 1)
  axis(2, labels = FALSE)
  title(ylab = 'Frequency', line = 1.8)
  
  lines(c(med, med), c(0, 85000), lty = 2, lwd = 1.2)
}


# ---------------------------------------------
# ---------- Figure 3: plot_model_curves
# ---------------------------------------------

#' Given a model, plot the estimation and prediction error curves for AR and VAR
plot_model_curves <- function(models6000_100, model, ylim, lwd, legend.cex = .5,
                              legend, n_seq_ss, block_cols = c("darkgreen", "tomato"),
                              alpha_val = .2, n_seq = 8:500) {
  
  # ID subsequence
  ind_ss <- n_seq %in% n_seq_ss
  
  # Set up Plot
  plot.new()
  plot.window(xlim = range(n_seq[ind_ss]), ylim=ylim)
  
  rescale_errors <- function(ar, var) {
    ar_err <- matrix(NA, nrow = nrow(ar), ncol = ncol(ar))
    var_err <- matrix(NA, nrow = nrow(ar), ncol = ncol(ar))
    
    for (i in seq(nrow(ar))) {
      mmin <- min(c(ar[i, ], var[i, ]))
      mmax <- max(c(ar[i, ], var[i, ]))
      
      ar_err[i, ] <- (ar[i, ] - mmin) / mmax
      var_err[i, ] <- (var[i, ] - mmin) / mmax
    }
    
    list('ar' = ar_err, 'var' = var_err)
  }
  
  m <- models6000_100[[model]]
  
  ee_ar <- t(m$EE_all[, 1, ])
  ee_var <- t(m$EE_all[, 2, ])
  pe_ar <- t(m$PE_all[, 1, ])
  pe_var <- t(m$PE_all[, 2, ])
  
  ee <- rescale_errors(ee_ar, ee_var)
  ee$ar <- list(
    'mean' = apply(ee$ar, 2, mean),
    'sd' = apply(ee$ar, 2, sd),
    'lo' = apply(ee$ar, 2, quantile, .25),
    'hi' = apply(ee$ar, 2, quantile, .75)
  )
  
  ee$var <- list(
    'mean' = apply(ee$var, 2, mean),
    'sd' = apply(ee$var, 2, sd),
    'lo' = apply(ee$var, 2, quantile, .25),
    'hi' = apply(ee$var, 2, quantile, .75)
  )
  
  pe <- rescale_errors(pe_ar, pe_var)
  pe$ar <- list(
    'mean' = apply(pe$ar, 2, mean),
    'sd' = apply(pe$ar, 2, sd),
    'lo' = apply(pe$ar, 2, quantile, .25),
    'hi' = apply(pe$ar, 2, quantile, .75)
  )
  
  pe$var <- list(
    'mean' = apply(pe$var, 2, mean),
    'sd' = apply(pe$var, 2, sd),
    'lo' = apply(pe$var, 2, quantile, .25),
    'hi' = apply(pe$var, 2, quantile, .75)
  )
  
  # PE AR
  lines(
    n_seq[ind_ss], pe$ar$mean[ind_ss], type="l", lty = 3, lwd=lwd,
    panel.first =  polygon(
      c(n_seq[ind_ss], rev(n_seq[ind_ss])),
      c(pe$ar$lo[ind_ss], rev(pe$ar$hi[ind_ss])),
      col = grDevices::adjustcolor('gainsboro', alpha = .4), border = NA
    )
  )
  
  # PE VAR
  lines(
    n_seq[ind_ss], pe$var$mean[ind_ss], type="l", lty = 3, lwd=lwd, col = 'red',
    panel.first =  polygon(
      c(n_seq[ind_ss], rev(n_seq[ind_ss])),
      c(pe$var$lo[ind_ss], rev(pe$var$hi[ind_ss])),
      col = grDevices::adjustcolor('gainsboro', alpha = .4), border = NA
    )
  )
  
  # EE AR
  lines(
    n_seq[ind_ss], ee$ar$mean[ind_ss], type="l", lty = 1, lwd=lwd,
    panel.first = polygon(
      c(n_seq[ind_ss], rev(n_seq[ind_ss])),
      c(ee$ar$lo[ind_ss], rev(ee$ar$hi[ind_ss])),
      col = grDevices::adjustcolor('gainsboro', alpha = .4), border = NA
    )
  )
  
  # EE VAR
  lines(
    n_seq[ind_ss], ee$var$mean[ind_ss], type="l", lty = 1, lwd=lwd, col = 'red',
    panel.first = polygon(
      c(n_seq[ind_ss], rev(n_seq[ind_ss])),
      c(ee$var$lo[ind_ss], rev(ee$var$hi[ind_ss])),
      col = grDevices::adjustcolor('gainsboro', alpha = .4), border = NA
    )
  )
  
  n_var <- 193
  # Plot Box
  ## Plot Gap
  # Intersection point PE
  signs_y <- pe$ar$mean - pe$var$mean
  int_n <- sign(signs_y[-1]) == sign(signs_y[-n_var]) 
  n_PE <- n_seq[which(int_n == FALSE)+1] # first n variation after the switch
  
  # Intersection point EE
  signs_y <- ee$ar$mean - ee$var$mean
  int_n <- sign(signs_y[-1]) == sign(signs_y[-n_var]) 
  n_EE <- n_seq[which(int_n == FALSE)+1] # first n variation after the switch
  
  n_ax_cmb <- c(n_PE, n_EE)
  
  ifelse(n_EE > n_PE, col <- block_cols[1], col <- block_cols[2])
  col <- scales::alpha(col, alpha = alpha_val)
  
  rect(
    xleft = min(n_ax_cmb), ybottom = 0, xright = max(n_ax_cmb),
    ytop = .25, col = col, border = FALSE
  )
  
  if (legend) {
    legend("topright", 
           legend = c("Estimation Error AR", "Estimation Error VAR",
                      "Prediction Error AR", "Prediction Error VAR"), 
           lty=c(1,1,3,3), 
           col=c("black", "red", "black", "red"), 
           bty = "n", 
           lwd = rep(lwd, 4), 
           cex = legend.cex)
  }
  
  x_labels <- c(min(n_seq[ind_ss]), median(n_seq[ind_ss]), max(n_seq[ind_ss]))
  y_at <- seq(ylim[1], ylim[2], length.out = 6)
  axis(1, x_labels, las=1)
  axis(2, at = y_at, labels = c('0', '', '', '', '', ''), las = 1)
  
  title(xlab="Number of observations n")
  title(ylab="Scaled MSE", line = 1.8)
  
}


# ---------------------------------------------
# ---------- Figure 3: plot_ngaps_hist, 
# ---------------------------------------------

#' Histogram for the expected gap
plot_ngaps_hist <- function(ngaps_models, weight = NULL) {
  xlim <- c(min(ngaps_models) - .5, max(ngaps_models) + 1)
  xlim <- c(-30, 50)
  
  if (!is.null(weight)) {
    wtd.hist(
      ngaps_models, breaks = 100, main = '',
      xlab = TeX('$$Expected n_{gap}$$'), col = 'grey76', ylab = '',
      xlim = xlim, xaxt = 'n', las = 2, weight = weight, ylim = c(0, .3)
    )
    
    title(ylab = 'Weighted Frequency', line = 3)
    
  } else {
    hist(
      ngaps_models, breaks = 50, main = '',
      xlab = TeX('$$Expected n_{gap}$$'), col = 'grey76', ylab = '',
      xlim = xlim, xaxt = 'n', las = 2, weight = weight
    )
    
    title(ylab = 'Frequency', line = 3)
  }
  
  x_labels <- seq(xlim[1], xlim[2], 7)
  x_labels <- seq(-30, 50, 10)
  axis(1, x_labels, las = 1)
  axis(2, labels = FALSE)
  
  med <- median(ngaps_models)
  lines(c(med, med), c(0, 1500), lty = 2, lwd = 1.2)
}


# -----------------------------------------------------------------------------
# ---------- Figure 5: Two plots about the 1SE Rule ---------------------------
# -----------------------------------------------------------------------------

plot_1SE <- function(EE_comp_n_NA, model_weights, mul = 100) {
  # ee_comp_mean[is.na(ee_comp_mean)] <- 0 # there is no disagreement!
  # ee_comp_sd[is.na(ee_comp_sd)] <- 0 # there is no disagreement!
  
  ee_comp_mean <- apply(EE_comp_n_NA, 2, function(col) {
    sel <- which(!is.na(col))
    # if (length(sel) == 0) { return(NA) }
    sum(col[sel] * model_weights[sel] / sum(model_weights[sel]))
  })[seq(200)]
  
  ee_comp_sd <- apply(EE_comp_n_NA, 2, function(col) {
    sel <- which(!is.na(col))
    # if (length(sel) == 0) { return(NA) }
    get_weighted_sd(col[sel], model_weights[sel] / sum(model_weights[sel]))
  })[seq(200)]
  
  par(mfrow = c(1, 2))
  
  plot.new()
  plot.window(
    xlim = c(8, 200), ylim = c(-.005, .02) * mul, yaxs = 'i'
  )
  
  lo <- (ee_comp_mean - ee_comp_sd) * mul
  hi <- (ee_comp_mean + ee_comp_sd) * mul
  lines(ee_comp_mean * mul)
  lines(lo, col = 'skyblue')
  lines(hi, col = 'skyblue')
  
  x_labels <- c(8, 50, 100, 150, 200)
  title(xlab = 'Number of observations n', cex.lab = 1)
  title(ylab = TeX('$$EE_{comp} x 100$$'), line = 2.5, cex.lab = 1)
  
  axis(1, x_labels, cex.axis = 1.2)
  axis(2, at = c(-.005, 0, .01, .02, .03)*100, las = 2, cex.axis = 1)
  
  lines(c(0, 500), c(0, 0), lty = 2, col = 'grey')
  text(180, 1.9, '(a)', col = 'black', cex = 1.2, adj = 0)
  
  ee_prop_n <- apply(EE_comp_n_NA, 2, function(col) mean(col > 0, na.rm = TRUE))[seq(193)]
  
  plot.new()
  plot.window(
    xlim = c(8, 200), ylim = c(0, 1), yaxs = 'i'
  )
  
  lines(ee_prop_n)
  
  x_labels <- c(8, 50, 100, 150, 200)
  title(xlab = 'Number of observations n', cex.lab = 1)
  title(ylab = TeX('$$EE_{comp} > 0$$'), line = 2.5, cex.lab = 1)
  
  axis(1, x_labels, cex.axis = 1.2)
  axis(2, seq(0, 1, .2), las = 2, cex.axis = 1)
  
  lines(c(0, 500), c(.5, .5), lty = 2, col = 'grey')
  text(180, .97, '(b)', col = 'black', cex = 1.2, adj = 0)
}


#' Adapted from plotfunctions::gradientLegend (used in Figure 2)
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
    text(y = ticks, x = loc[3] - 0.005, labels = labels, pos = 4, 
         col = border.col, cex = cex, xpd = T)
  }
}
