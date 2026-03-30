# source("copilot_fn.R")

###################################################################
### my fn
###################################################################

###################################################################
plot_gam1 <- function(ft1,stpdf=NULL) {
  
  if(!is.null(stpdf)) pdf(file=stpdf, width=12, height=9)

  plot(ft1, pages=1, shade=TRUE)
  
  q1       <- predict(ft1) 
  iy2001   <- which(trunc(data01$time)==2001 & data01$isobs==1)
  i0x      <- which.max(q1[ iy2001] )
  i0n      <- which.min(q1[ iy2001] )
  doyn     <- data01$doy[ iy2001][i0n]
  sdoyn    <- data01$sdoy[iy2001][i0n]
  doyx     <- data01$doy[ iy2001][i0x]
  sdoyx    <- data01$sdoy[iy2001][i0x]
  nd       <- data01
  nd$sdoy  <- sdoyn
  q1.n     <- predict(ft1,  newdata=nd )
  nd$sdoy  <- sdoyx
  q1.x     <- predict(ft1,  newdata=nd )
  # remove interannual variability
  nd       <- data01
  nd$fYear <- "2003"
  q1.2003  <- predict(ft1, newdata=nd) 
  
  up.1()
  plot(ft1$y, pch=20, cex=.3, main=paste("x ~ ",ft1$formula[3]), cex.main=0.9, xlab="index", ylab="Temperature")
  lines(q1,   col=2, lwd=2)
  lines(q1.n, col=4, lwd=2)
  lines(q1.x, col=3, lwd=2)
  ix <- which(data01$doy==doyx)
  lines(ix, q1.2003[ix], col=6, lwd=2)
  im <- which(data01$doy==doyn)
  lines(im, q1.2003[im], col=6, lwd=2)
  legend("topleft", legend=c("OBS/GCM","median(year,doy)", "Winter min", "Winter min @2003 IAV", "Summer max", "Summer max @2003 IAV"), 
                          col=c(1,2,4,6,3,6), lwd=c(NA,2,2,2,2,2), pch=c(20,NA,NA,NA,NA,NA), bty="n", cex=1)
  grid()
  if(is.null(stpdf)) readline("continue?")
  
  up.1()
  ix <- 16000:18000-185
  plot(ft1$y[ix], pch=20, cex=.3, main=paste("x ~ ",ft1$formula[3]), cex.main=0.9, xlab="index", ylab="Temperature")
  lines(q1[ix],   col=2, lwd=2)
  lines(q1.n[ix], col=4, lwd=2)
  lines(q1.x[ix], col=3, lwd=2)
  iy <- which(data01$doy[ix]==doyx)
  lines(iy, q1.2003[ix][iy], col=6, lwd=2); points(iy, q1.2003[ix][iy], col=6, pch=20, cex=2.0)
  iz <- which(data01$doy[ix]==doyn)
  lines(iz, q1.2003[ix][iz], col=6, lwd=2); points(iz, q1.2003[ix][iz], col=6, pch=20, cex=2.0)
  legend("topleft", legend=c("OBS/GCM","median(year,doy)", "Winter min", "Winter min @2003 IAV", "Summer max", "Summer max @2003 IAV"), 
                          col=c(1,2,4,6,3,6), lwd=c(NA,2,2,2,2,2), pch=c(20,NA,NA,20,NA,20), bty="n", cex=1)
  grid()
  if(is.null(stpdf)) readline("continue?")
  
  # residuals
  up.1()
  plot(ft1$y-q1, pch=20, cex=.3, main=paste("Residuals"), cex.main=0.9, xlab="index", ylab="Temperature (C)")
  grid()
  if(!is.null(stpdf)) dev.off()
} 

### legacy
plot_regam <- plot_gam1


plot_gam2  <- function(ft1,stpdf=NULL) {
  
  if(!is.null(stpdf)) pdf(file=stpdf, width=12, height=9)

  plot(ft1, pages=1, shade=TRUE)
  
  q1       <- predict(ft1) 
  iy2001   <- which(trunc(data01$time)==2001 & data01$isobs==1)
  i0x      <- which.max(q1[ iy2001] )
  i0n      <- which.min(q1[ iy2001] )
  doyn     <- data01$doy[ iy2001][i0n]
  sdoyn    <- data01$sdoy[iy2001][i0n]
  doyx     <- data01$doy[ iy2001][i0x]
  sdoyx    <- data01$sdoy[iy2001][i0x]
  nd       <- data01
  nd$sdoy  <- sdoyn
  q1.n     <- predict(ft1,  newdata=nd )
  nd$sdoy  <- sdoyx
  q1.x     <- predict(ft1,  newdata=nd )
  # remove interannual variability
  nd       <- data01
  nd$fYear <- "2003"
  q1.2003  <- predict(ft1, newdata=nd) 
  
  up.1()
  plot(ft1$y, pch=20, cex=.3, main=paste("x ~ ",ft1$formula[3]), cex.main=0.9, xlab="index", ylab="Temperature")
  lines(q1,   col=2, lwd=2)
  legend("topleft", legend=c("OBS/GCM","prediction"), col=c(1,2), lwd=c(NA,2), pch=c(20,NA), bty="n", cex=1)
  grid()
  if(is.null(stpdf)) readline("continue?")
  
  up.1()
  ix <- 16000:18000-185
  plot(ft1$y[ix], pch=20, cex=.3, main=paste("x ~ ",ft1$formula[3]), cex.main=0.9, xlab="index", ylab="Temperature")
  lines(q1[ix],   col=2, lwd=2)
  legend("topleft", legend=c("OBS/GCM","prediction"), col=c(1,2), lwd=c(NA,2), pch=c(20,NA), bty="n", cex=1)
  grid()
  if(is.null(stpdf)) readline("continue?")
  
  # residuals
  up.1()
  plot(ft1$y-q1, pch=20, cex=.3, main=paste("Residuals"), cex.main=0.9, xlab="index", ylab="Temperature (C)")
  grid()
  if(!is.null(stpdf)) dev.off()
} 


###################################################################
plot_regam_RE <- function(ft1, do.select=NULL, stpdf=NULL, idx=1:4, do.years=1980:2025) {
    
    library(gratia)
    library("patchwork")

    if(!is.null(stpdf)) pdf(file=stpdf, width=12, height=9)

    sm1       <- smooth_estimates(ft1)
    names_ft1 <- unique(sm1[[1]])
    p1 <- draw(smooth_estimates(ft1, select=names_ft1[idx]))
    print(p1)
    if(is.null(stpdf)) readline("continue?")

    # sm.re.o <- smooth_estimates(ft1, select="s(sdoy,fYear):ftypeobs", data=data01[which(data01$isobs==1),])
    # draw(sm.re.o)
    # sm.re.m <- smooth_estimates(ft1, select="s(sdoy,fYear):ftypemod", data=data01[which(data01$isobs==0),])
    # draw(sm.re.m)

    ncol    <- 6
    nrow    <- 4
    newpage <- TRUE
    for( i in do.years ) {
      iy2  <- which(data01$year == i & data01$isobs==1)
      smo <- smooth_estimates(ft1, select=do.select, data=data01[iy2,])
      if(newpage) {
        pp <- draw(smo) 
        newpage <- FALSE
      } else {
        pp <- pp + draw(smo)
      }
      iy2  <- which(data01$year == i & data01$isobs==0)
      smm <- smooth_estimates(ft1, select=do.select, data=data01[iy2,])
      pp <- pp + draw(smm)
      # print(p)
      # p1 + p2 +plot_layout(ncol = 2)
      # cat(data01$year[iy2[1]],length(pp), newpage, "\n")
      if(length(pp)==(ncol*nrow)) {
        print(pp + plot_layout(ncol=ncol, nrow=nrow))
        newpage <- TRUE
        readline("Continue?")
      }
    }
    print(pp + plot_layout(ncol=ncol, nrow=nrow))

    if(is.null(stpdf)) readline("continue?")
    if(!is.null(stpdf)) dev.off()
} 




###################################################################
# copilot fn
###################################################################
plot_gam_terms <- function(fit, terms=NULL, data=NULL, sort_by="sdoy", add=FALSE, col="blue", lwd=2, LINES=TRUE, SE=TRUE,...) {

  # Use model data unless custom data provided
  if (is.null(data)) data <- fit$model

  # Predict all term contributions
  p <- predict(fit, data, type = "terms", se.fit = SE)

  # List available terms
  if(SE) term_names <- colnames(p$fit) else term_names <- colnames(p)

  if (is.null(terms)) {
    # message("Available terms:\n", paste(term_names, collapse = "\n"))
    # stop("Specify one or more terms = c('term1','term2',...).")
    terms <- term_names
  }

  # Check supplied term names
  bad <- setdiff(terms, term_names)
  if (length(bad) > 0) {
    stop("These terms were not found:\n", paste(bad, collapse = "\n"),
         "\n\nAvailable terms:\n", paste(term_names, collapse = "\n"))
  }

  # Extract selected contributions
  if(SE) {
    effect <- rowSums(p$fit[, terms, drop = FALSE])
    se     <- sqrt(rowSums(p$se.fit[, terms, drop = FALSE]^2))
  } else {
    effect <- rowSums(p[, terms, drop = FALSE])
    se <- NULL
  } 

  # Sorting axis
  if (!sort_by %in% colnames(data)) {
    # stop("sort_by must be a column in the data.")
    ord <- seq_len(nrow(data))
    xx  <- ord
  } else {
    ord <- order(data[[sort_by]])
    xx  <- data[[sort_by]][ord]
    }

  yy <- effect[ord]
  yy_up <- yy + 2 * se[ord]
  yy_dn <- yy - 2 * se[ord]

  # Plot
  if (!add) {
    if(LINES) plot(xx, yy, type = "l", lwd = lwd, col=col,
                        xlab = sort_by, ylab = "Contribution to linear predictor", ...) else 
          plot(xx, yy, col=col, pch=20, cex=0.2, xlab = sort_by, ylab = "Contribution to linear predictor", ...)
  } else {
    if(LINES) lines(xx, yy, lwd = lwd, col=col) else points(xx, yy, col=col, pch=20, cex=0.2)
  }

  # Add optional CI ribbon
  if(SE) {
    lines(xx, yy_up, col=col, lty = 2)
    lines(xx, yy_dn, col=col, lty = 2)
  }

  invisible(list(x = xx, fit = yy, upper = yy_up, lower = yy_dn))
}


###################################################################
plot_random_effect <- function(fit, fyear, ftype, n = 200) {

  # Sequence for scaled day-of-year
  newdat <- data.frame(
    sdoy  = seq(0, 1, length.out = n),
    fYear = factor(fyear, levels = levels(fit$model$fYear)),
    ftype = factor(ftype, levels = levels(fit$model$ftype)),
    gmst  = mean(fit$model$gmst, na.rm = TRUE)  # dummy value — removed later
  )

  # Predict ONLY the random-effect term
  # type="terms" returns a matrix of individual smooth contributions
  pred_terms <- predict(fit, newdat, type = "terms", se.fit = TRUE)

  # Find the correct column: mgcv names it like "s(sdoy,fYear,ftype)"
  term_name <- grep("sdoy.*fYear.*ftype", colnames(pred_terms$fit), value = TRUE)

  if (length(term_name) != 1) {
    stop("Could not uniquely identify the random-effect smooth term.")
  }

  reffect  <- pred_terms$fit[, term_name]
  rese     <- pred_terms$se.fit[, term_name]

  # Plot with confidence interval
  plot(
    newdat$sdoy, reffect, type = "l", lwd = 2,
    xlab = "Scaled day of year (sdoy)",
    ylab = "Random-effect smooth contribution",
    main = paste("Random-effect smooth for", ftype, ", Year =", fyear)
  )

  lines(newdat$sdoy, reffect + 2 * rese, lty = 2)
  lines(newdat$sdoy, reffect - 2 * rese, lty = 2)
}


###################################################################
plot_random_effect_time <- function(ft1, data01, ftype='obs', n = 200) {

  # Sequence for scaled day-of-year
  newdat <- data.frame(
    sdoy  = seq(0, 1, length.out = n),
    fYear = factor(fyear, levels = levels(ft1$model$fYear)),
    ftype = factor(ftype, levels = levels(ft1$model$ftype)),
    gmst  = mean(ft1$model$gmst, na.rm = TRUE)  # dummy value — removed later
  )

  # Predict ONLY the random-effect term
  # type="terms" returns a matrix of individual smooth contributions
  pred_terms <- predict(ft1, data01, type = "terms")

  # Find the correct column: mgcv names it like "s(sdoy,fYear,ftype)"
  term_name <- grep("sdoy.*fYear.*ftype", colnames(pred_terms))

  if (length(term_name) != 1) {
    stop("Could not uniquely identify the random-effect smooth term.")
  }

  reffect  <- pred_terms[, term_name]

  # Plot with confidence interval
  plot(
    data01$time, reffect, type = "l", lwd = 2,
    xlab = "Date",
    ylab = "Random-effect smooth contribution",
    main = paste("Random-effect smooth")
  )

  lines(newdat$sdoy, reffect + 2 * rese, lty = 2)
  lines(newdat$sdoy, reffect - 2 * rese, lty = 2)
}

# plot(pred_terms[, 1], pch=20, cex=0.3)

###################################################################
# Helper: identify the column in predict(type="terms") that corresponds to the sz random-effect term
.get_sz_term_name <- function(fit) {
  # Terms matrix for a tiny newdata to read column names safely
  tmp <- predict(fit, fit$model[1, , drop = FALSE], type = "terms")
  nm  <- colnames(tmp)
  # Look for the term that contains all three: sdoy, fYear, ftype
  cand <- grep("sdoy.*fYear.*ftype", nm, value = TRUE)
  cand <- grep("sdoy.*fYear", nm, value = TRUE)
  if (length(cand) != 1L) {
    stop("Could not uniquely identify the sz random-effect term. Found: ",
         paste(cand, collapse = ", "))
  }
  cand
}

###################################################################
# Base function: get random-effect contributions for (years × sdoy) at a fixed ftype
get_random_effect_curves <- function(fit, years, ftype, n = 200) {
  # Ensure factors match model levels
  ftype <- factor(ftype, levels = levels(fit$model$ftype))
  years <- factor(years, levels = levels(fit$model$fYear))

  # Build prediction grid for all requested years
  newdat <- expand.grid(
    sdoy  = seq(0, 1, length.out = n),
    fYear = years,
    ftype = ftype,
    KEEP.OUT.ATTRS = FALSE
  )

  # gmst is required but we only need the sz-term; set to a benign value
  # (It won't affect the extracted term for sz.)
  newdat$gmst <- mean(fit$model$gmst, na.rm = TRUE)

  # Predict all terms and then extract the sz term
  sz_name <- .get_sz_term_name(fit)
  p <- predict(fit, newdat, type = "terms", se.fit = TRUE)

  out <- data.frame(
    sdoy   = newdat$sdoy,
    fYear  = newdat$fYear,
    ftype  = newdat$ftype,
    effect = p$fit[, sz_name],
    se     = p$se.fit[, sz_name]
  )
  rownames(out) <- NULL
  out
}

###################################################################
# ---------- 1) BASE R OVERLAY PLOT ----------
plot_random_effect_overlay_base <- function(fit, years, ftype, n = 200, col_pal = NULL) {
  df <- get_random_effect_curves(fit, years, ftype, n = n)

  # Colors
  years_chr <- as.character(unique(df$fYear))
  if (is.null(col_pal)) {
    col_pal <- grDevices::hcl.colors(length(years_chr), "Dark3", rev = FALSE)
  } else if (length(col_pal) < length(years_chr)) {
    stop("col_pal must have at least as many colors as years requested.")
  }

  # Set up empty plot
  rng_y <- range(df$effect + 2*df$se, df$effect - 2*df$se, na.rm = TRUE)
  plot(NA, xlim = c(0,1), ylim = rng_y,
       xlab = "Scaled day of year (sdoy)",
       ylab = "Random-effect smooth contribution",
       main = paste0("Year-specific seasonal anomalies (", ftype, ")"))

  # Draw each year
  for (i in seq_along(years_chr)) {
    yi  <- years_chr[i]
    dyy <- subset(df, fYear == yi)
    lines(dyy$sdoy, dyy$effect, col=col_pal[i], lwd = 2)
    # Optional CI band as dashed lines (comment out if cluttered)
    # lines(dyy$sdoy, dyy$effect + 2*dyy$se, col=col_pal[i], lty = 3)
    # lines(dyy$sdoy, dyy$effect - 2*dyy$se, col=col_pal[i], lty = 3)
  }

  legend("topleft", legend = years_chr, col=col_pal, lwd = 2, bty = "n", ncol = 1,
         title = "Year")
}

# Example (Base R):
# plot_random_effect_overlay_base(fit1, years = c("1990","1998","2005","2010"), ftype = "obs")


###################################################################
# ---------- 2) GGPLOT2 OVERLAY PLOT ----------
# (Nice for publication; supports CI ribbon toggle)
plot_random_effect_overlay_gg <- function(fit, years, ftype, n = 200,
                                          show_ci = FALSE, alpha_ci = 0.15) {
  # Require ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.")
  }

  df <- get_random_effect_curves(fit, years, ftype, n = n)
  df$lower <- df$effect - 2*df$se
  df$upper <- df$effect + 2*df$se

  gg <- ggplot2::ggplot(df, ggplot2::aes(sdoy, effect, colour = fYear)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
    (if (show_ci) ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = fYear),
                                       colour = NA, alpha = alpha_ci, show.legend = FALSE) else NULL) +
    ggplot2::geom_line(size = 0.9) +
    ggplot2::scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    ggplot2::labs(
      x = "Scaled day of year (sdoy)",
      y = "Random-effect smooth contribution",
      colour = "Year",
      title = paste0("Year-specific seasonal anomalies (", ftype, ")")
    ) +
    ggplot2::theme_minimal(base_size = 12)

  return(gg)
}

# Examples (ggplot2):
# plot_random_effect_overlay_gg(fit1, years = c("1990","1998","2005","2010"), ftype = "obs")
# plot_random_effect_overlay_gg(fit1, years = levels(fit1$model$fYear)[1:10], ftype = "mod", show_ci = TRUE)





###################################################################
########## sim ####################################################
###################################################################

###################################################################
mk_newdata_for_date <- function(future_date, gmst_value, ftype_level,
                                origin_date, tz = "UTC") {
  # Ensure Date/POSIX handling is consistent with your original pre-processing
  future_date <- as.Date(future_date, tz)
  sdoy <- as.integer(strftime(future_date, format = "%j", tz = tz))
  
  # stime in the SAME units/origin used in model fit
  stime <- as.numeric(future_date - as.Date(origin_date, tz))
  
  data.frame(
    sdoy  = sdoy,
    gmst  = gmst_value,
    ftype = factor(ftype_level, levels = levels_used_in_fit),  # supply levels from your fit
    stime = stime
  )
}

###################################################################
simulate_gam <- function(model, newdata, nsim = 1000,
                         component = c("full", "s(stime)", "smooth"),
                         ftype_level = NULL,
                         response = TRUE,
                         include_residual = FALSE,
                         seed = NULL) {
  # component:
  #   "full"     -> all terms
  #   "s(stime)" -> only the temporal GP smooth columns (optionally for a specific ftype level)
  #   "smooth"   -> all smooth terms (drops parametric/fixed effects)
  

    # What it does - from copilot:
    #     - Uses predict(..., type="lpmatrix") and vcov() to simulate from the posterior of the GAM coefficients—this naturally includes your s(stime, bs="gp", by=ftype) term.
    #     - component="full" gives everything.
    #     - component="s(stime)" isolates only the GP temporal random effect (and can target a specific ftype level with ftype_level="obs" or "model" etc.).
    #     - include_residual=TRUE adds Gaussian residual noise for a predictive distribution.


  component <- match.arg(component)
  if (!is.null(seed)) set.seed(seed)
  
  # 1) Design matrix for newdata
  X <- predict(model, newdata, type = "lpmatrix")
  cn <- colnames(X)
  
  # 2) Posterior covariance of coefficients & mean
  #    unconditional=TRUE includes smoothing parameter uncertainty (recommended)
  V <- vcov(model, unconditional = TRUE) # theo: model$Vp
  beta_hat <- coef(model)
  
  # 3) Restrict to chosen component
  if (component != "full") {
    keep <- rep(TRUE, length(cn))
    
    if (component == "s(stime)") {
      # Keep only columns that belong to s(stime), optionally for an ftype level
      keep <- grepl("s\\(stime\\)", cn)
      if (!is.null(ftype_level)) {
        # mgcv typically encodes by-factor in the column name as ":ftype<level>"
        keep <- keep & grepl(paste0(":ftype", ftype_level), cn)
      }
      if (!any(keep)) {
        stop("Couldn't find columns for s(stime) ",
             if (!is.null(ftype_level)) paste0("with ftype level '", ftype_level, "'") else "",
             ". Inspect colnames(predict(..., type='lpmatrix')).")
      }
    } else if (component == "smooth") {
      # Drop parametric (non-smooth) terms; smooth columns usually contain "s(" or "ti("
      keep <- grepl("s\\(|ti\\(", cn)
    }
    
    # Zero-out everything else
    X[, !keep] <- 0
  }
  
  # 4) Draw coefficients from posterior
  # NOTE: MASS::mvrnorm is the standard choice
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Please install.packages('MASS')")
  }
  beta_sim <- MASS::mvrnorm(nsim, mu = beta_hat, Sigma = V)
  # theo: b_sims <- rmvn(n.sims, coef(model), model$Vp )

  # 5) Simulated linear predictors
  eta_sim <- beta_sim %*% t(X)   # nsim x nrow(newdata)
  # theo: SC_sims     <- tcrossprod( X[,DoYindex],  b_sims[,DoYindex] )

  # 6) Transform to response scale if requested
  if (response) {
    eta_sim <- model$family$linkinv(eta_sim)
  }
  
  # 7) Optionally add residual noise (for predictive distribution)
  if (include_residual) {
    fam <- family(model)$family
    if (grepl("Gaussian", fam, ignore.case = TRUE)) {
      sigma2 <- summary(model)$scale
      eps <- matrix(stats::rnorm(length(eta_sim), sd = sqrt(sigma2)),
                    nrow = nrow(eta_sim), ncol = ncol(eta_sim))
      eta_sim <- eta_sim + eps
    } else {
      warning("include_residual=TRUE currently implemented for Gaussian family only.")
    }
  }
  
  eta_sim
}

