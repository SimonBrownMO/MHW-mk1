
# source("doFitInitValue.R")

list2env(HDconfig$initV , envir = .GlobalEnv)

load(st_events,verb=TRUE)

gpd.th.u  <- initV$th.u

# ### ensure the doy basis fn is consistent with hotseason if present
# if(events01$info$USEHOTSEASON & events01$hotseas$IsHotSeason) {
#     fmla.IVgpd <- initV$fmla.hseas.IVgpd
# } else {
#     fmla.IVgpd <- initV$fmla.ally.IVgpd
# }

# allIVfits <- list()
# for(i in seq_along(fmla.IVgpd)) {
#     allIVfits[[i]]      <- fn_jointFitInitVal(events01, gpd.th.u, fmla.IVgpd[[i]], IsJoint=TRUE)
#     allIVfits[[i]]$name <- names(fmla.IVgpd)[i]
#     cat("Fitted IV",i,names(fmla.IVgpd)[i],cr)
# }

# ### filter out bad fits
# isFitGood         <- NULL
# for(i in seq_along(fmla.IVgpd)) {
#     x1 <- summary(allIVfits[[i]])
#     pval <- c( x1[[1]]$logscale$'Pr(>|t|)', x1[[2]]$logscale$'Pr(>|t|)', x1[[2]]$shape$'Pr(>|t|)', x1[[1]]$shape$'Pr(>|t|)')
#     if(all(pval < 1e-6)) {isFitGood[i] <- FALSE; cat(i,names(fmla.IVgpd)[i],pval,cr)} else isFitGood[i] <- TRUE
# }
# allIVfits  <-  allIVfits[isFitGood]
# fmla.IVgpd <- fmla.IVgpd[isFitGood]

# ### save statistics of each fit
# st.InitVal.stats <- paste(dirname(HDconfig$initV$st_InitValue),'diag',sub('.RData','.stats',basename(HDconfig$files$st_InitValue)),sep='/')
# for(i in seq_along(allIVfits)) {
#     if(i==1) fn_diagInitVal(allIVfits[[i]], st.GpdInitDiag=st.InitVal.stats, SAVE=TRUE, NEWFILE=TRUE) else
#                 fn_diagInitVal(allIVfits[[i]], st.GpdInitDiag=st.InitVal.stats, SAVE=TRUE, NEWFILE=FALSE)
# }

# ### choose best fit
# allIVfits.aic        <- lapply(allIVfits, AIC)
# names(allIVfits.aic) <- names(fmla.IVgpd)
# cat("AIC",cr)
# for(i in seq_along(fmla.IVgpd)) cat(names(fmla.IVgpd)[i],allIVfits.aic[[i]],cr)

# for(i in seq_along(fmla.IVgpd)) {
#     cat(cr,'################',cr,names(fmla.IVgpd)[i],cr)
#     print(summary(allIVfits[[i]]))
# }

# # I dont trust the AIC as it does not pick out the failed fits
# cat(cr,cr,"##############################",cr,"Best IV GPD fit",cr)
# ix1 <- which.min(allIVfits.aic)
# cat(names(fmla.IVgpd)[ix1],allIVfits.aic[[ix1]],cr)
# chosen.IVgpd.name <- names(fmla.IVgpd)[ix1]
# chosen.IVgpd      <- allIVfits[[ix1]]

#################################################################################################
### now simple setup
chosen.IVgpd<- fn_jointFitInitVal(events01, gpd.th.u, initV$fmla, IsJoint=TRUE)

### save
save(file=HDconfig$initV$st_InitValue, chosen.IVgpd)

# Plot chosen fit
st.pdf.IV     <- paste(dirname(HDconfig$initV$st_InitValue),'plots',sub('.RData','.pdf',basename(HDconfig$initV$st_InitValue)),sep='/')
if(DOINITVALUEPLOT) fn_plotInitVal(chosen.IVgpd, chosen.IVgpd$fit.data, xplot=1, retp=NULL, st.pdf=st.pdf.IV)

cat("FitInitValue all done",cr)
