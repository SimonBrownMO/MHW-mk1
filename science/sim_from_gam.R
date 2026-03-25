# source("sim_from_gam.R")

library(gratia)



###################################################################
# Example assumptions:
# - 'date' was converted to 'stime' as days since origin_date (numeric)
# - 'sdoy' = day-of-year 1..366 used with bs="cc"
# - 'gmst' is provided externally for the future date (scenario/projection)
# - 'ftype' has levels e.g., c("obs", "model")

# Important: Replace levels_used_in_fit with levels(data01$ftype) from the data you used 
# to fit the model, and pass the same origin_date/transforms you used when you created 
# stime for training.
# If you centered/scaled stime or gmst before fitting, apply the same transforms to 
# the future values.


###################################################################
## End‑to‑end example: simulate many samples for observations at a future date


fmla     <- x ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts') +ti(sdoy,gmst,bs=c("cc","ts")) +s(stime,bs="gp",by=ftype, m=1, k=140*8, sp=c(40,10))
load(file="ft1_regamJ2m1.RData",verbose=TRUE)
ft1 <- ft1.regamJ2m1
[1] "s(sdoy):ftypeobs"  "s(sdoy):ftypemod"  "s(gmst)"          
[4] "ti(sdoy,gmst)"     "s(stime):ftypeobs" "s(stime):ftypemod"

## Get factor levels as used in the fit
levels_used_in_fit <- levels(data01$ftype)

iobs      <- which(data01$isobs == 1)
r.stime.o <- range(data01$stime[iobs])

## Build newdata for OBSERVATIONS branch
new_obs_1d  <- data.frame(sdoy=0.6,                 gmst=0.1, ftype=factor("obs",levels=levels_used_in_fit), stime=-0.1)
new_obs_1y  <- data.frame(sdoy=seq(0,1,length=365), gmst=0.1, ftype=factor("obs",levels=levels_used_in_fit), stime=-0.1)
new_obs_1yb <- data.frame(sdoy=seq(0,1,length=365), gmst=0.1, ftype=factor("obs",levels=levels_used_in_fit), stime=0.25)
new_obs_1t  <- data.frame(sdoy=0.6,                 gmst=0.1, ftype=factor("obs",levels=levels_used_in_fit), stime=seq(r.stime.o[1],r.stime.o[2],length=100) )

### (A) Full predictive distribution INCLUDING residual noise

### 1d
sim_obs_1d_full      <- simulate_gam(model=ft1, newdata=new_obs_1d, component="full",     include_residual=TRUE)
sim_obs_1d_stim      <- simulate_gam(model=ft1, newdata=new_obs_1d, component="s(stime)", include_residual=TRUE)
sim_obs_1d_full_nres <- simulate_gam(model=ft1, newdata=new_obs_1d, component="full",     include_residual=FALSE)
sim_obs_1d_stim_nres <- simulate_gam(model=ft1, newdata=new_obs_1d, component="s(stime)", include_residual=FALSE)
up.2()
hist(sim_obs_1d_full,      breaks=100, col=rgb(0,0,1,0.3), main="FULL Simulated obs for sdoy=0.6 stime=-0.1", xlab="x")
hist(sim_obs_1d_stim,      breaks=100, col=rgb(0,1,1,0.3), main="s(stime) Simulated obs for sdoy=0.6 stime=-0.1", xlab="x")
# hist(sim_obs_1d_full_nres, breaks=100, col=rgb(0,0,1,0.3), main="FULL Simulated obs (no residual) for sdoy=0.6 stime=-0.1", xlab="x")
# hist(sim_obs_1d_stim_nres, breaks=100, col=rgb(0,1,1,0.3), main="s(stime) Simulated obs (no residual) for sdoy=0.6 stime=-0.1", xlab="x")

readline("stop1")

### 1y
sim_obs_1y_full  <- simulate_gam(model=ft1, newdata=new_obs_1y,  component="full",     include_residual=TRUE)
sim_obs_1yb_full <- simulate_gam(model=ft1, newdata=new_obs_1yb, component="full",     include_residual=TRUE)
sim_obs_1y_stim  <- simulate_gam(model=ft1, newdata=new_obs_1y,  component="s(stime)", include_residual=TRUE)
sim_obs_1yb_stim <- simulate_gam(model=ft1, newdata=new_obs_1yb, component="s(stime)", include_residual=TRUE)
up.2()
plot(1:365, sample(sim_obs_1y_full, 365), type="n", main="FULL Simulated obs for 1 year stime=-0.1", xlab="Day of Year", ylab="x", ylim=range(sim_obs_1y_full))
for(i in 1:365) points(rep(i,1000),     sim_obs_1y_full[,i],  col=rgb(0,0,1,0.05), pch=3, cex=0.3)
for(i in 1:365) points(rep(i,1000)+0.5, sim_obs_1yb_full[,i], col=rgb(1,0,0,0.05), pch=3, cex=0.3)
lines(1:365, apply(sim_obs_1y_full,  2, mean), col=4, lwd=2)
lines(1:365, apply(sim_obs_1yb_full, 2, mean), col=rgb(0.4,0,0,1), lwd=2)
grid()

plot(1:365, sample(sim_obs_1y_stim, 365), type="n", main="s(stime) Simulated obs for 1 year stime=-0.1", xlab="Day of Year", ylab="x", ylim=range(sim_obs_1y_stim))
for(i in 1:365) points(rep(i,1000),      sim_obs_1y_stim[,i], col=rgb(.4,0,0,0.1), pch=3, cex=0.3)
for(i in 1:365) points(rep(i,1000)+0.5, sim_obs_1yb_stim[,i], col=rgb(0,0,0.8,0.1), pch=3, cex=0.3)
grid()

readline("stop2")

### 1t
sim_obs_1t_full  <- simulate_gam(model=ft1, newdata=new_obs_1t,  component="full",     include_residual=TRUE)
sim_obs_1t_stim  <- simulate_gam(model=ft1, newdata=new_obs_1t,  component="s(stime)", include_residual=TRUE)

plot(new_obs_1t$stime, sample(sim_obs_1t_full, 100), type="n", main="Simulated obs for sdoy=0.6 stime=obs_range_100", xlab="stime", ylab="x", ylim=range(sim_obs_1t_full))
for(i in 1:100) points(rep(new_obs_1t$stime[i],1000), sim_obs_1t_full[,i], col=rgb(0,0,1,0.1), pch=20, cex=0.3)
for(i in 1:100) points(new_obs_1t$stime[i], mean(sim_obs_1t_full[,i]), col=1, pch=20)

# NEED TO TRY CHANGES WITH GMST TOO




## Summaries
mean_pred   <- mean(sim_full_obs)
pi90        <- quantile(sim_full_obs, probs = c(0.05, 0.95))
c(mean = mean_pred, `5%` = pi90[1], `95%` = pi90[2])

## (B) Only the GP temporal random effect for the 'obs' branch
sim_gp_obs <- simulate_gam(
  model = m,
  newdata = new_obs,
  nsim = 5000,
  component = "s(stime)",
  ftype_level = "obs",
  response = TRUE,        # link-inv; for Gaussian with identity link, same as link scale
  include_residual = FALSE,
  seed = 42
)

## Compare distributions
hist(sim_full_obs, breaks = 40, col = rgb(0,0,1,0.3),
     main = "Full vs GP-only simulations (obs, future date)",
     xlab = "x (response scale)")
hist(sim_gp_obs,  breaks = 40, col = rgb(1,0,0,0.3), add = TRUE)
legend("topright", legend = c("Full (with residual)", "GP s(stime) only"),
       fill = c(rgb(0,0,1,0.3), rgb(1,0,0,0.3)), bty = "n")



### Notes & pitfalls (worth a skim)

    # Outside-range behaviour (future dates):
    # The gp smoother is stationary; as you predict further from the data, the temporal GP effect 
    # tends to shrink toward its prior mean and the posterior variance can increase. The simulation 
    # code above respects that automatically because it draws from the fitted coefficients with 
    # their covariance.

    # By‑factor GP (by=ftype):
    # The function detects columns belonging to s(stime) and, if ftype_level is provided, restricts 
    # to that branch (e.g., "obs" vs "model"). If your column labels differ, run:

    #     > colnames(predict(m, new_obs, type = "lpmatrix"))

    # to confirm the exact pattern (they’re typically like s(stime):ftypeobs.1, etc.).

    # Fixed sp values:
    # You’ve provided sp=c(40,10) within the GP smooth to strongly penalise it. That’s fine; just 
    # know that unconditional=TRUE in vcov() won’t add smoothing‑parameter uncertainty for fixed 
    # sp. You can set unconditional=FALSE if you’ve fully fixed all sps, but leaving it TRUE is 
    # harmless.

    # Scaling/centering:
    # If you centered/scaled gmst or stime when fitting, do the same transform for future 
    # predictions—otherwise your GP/ts terms won’t behave as intended.

    # Residual noise:
    # Set include_residual=TRUE if you want the predictive distribution (process + observation noise). 
    # Leave it FALSE to simulate only the latent process defined by the smooths and parametric terms.

    # Memory/performance:
    # For large nsim, multiplying beta_sim %*% t(X) is fast for a single date (since X is a row), but 
    # if you simulate many future dates at once, consider chunking or vectorising over dates.
    