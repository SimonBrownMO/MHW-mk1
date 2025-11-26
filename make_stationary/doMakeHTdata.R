# source("doMakeHTdata.R")

sysdate   <- Sys.time()
cat("###############################################################################",cr)
stsysdate <- print(sysdate)
cat("Start doMakeStationary",cr)
cat("###############################################################################",cr,cr)


### Calculate annualcycle and climate change signal with QGAM
# fit QGAM

if(MSconfig$DOQGAM) {

    cat(cr,"Fitting QGAM",cr)

    if(length(data01$time)>MSconfig$nmax_qgam) {
        set.seed(0)
        i0 <- sample(seq_along(data01$time), MSconfig$nmax_qgam, replace=FALSE)
    } else i0 <- seq_along(data01$time)

    ### need to add here byYear covariates
    y.all <- trunc(data01$time)
    y.u   <- unique(y.all)
    cov.y <- array(0, dim=c(length(data01$time), length(y.u)) ) 
    fac.y <- NULL
    for(i in seq_along(y.u)) {
        iy <- which(y.all==y.u[i])
        cov.y[iy,i] <- 1
        fac.y       <- c(fac.y, rep(i, length(iy)) )
    }

readline("stop1")
    # fn_fit_qgam(MSconfig$fmla.MSqgam, data01[i0,], st.qgam=MSconfig$files$st_qgam, ptiles=MSconfig$do.ptiles)

    ### add years as factors to data01
    data01$year  <- factor(fac.y,labels=y.u) # weired factor thing labels have to be the ranked order of the data
    data01$ryear <- fac.y/max(fac.y)
    data01$x0 <- data01$x
    ### plan 1
        # fit just for gmst and remove
        # fit for year to get year to year variability and remove
        # fit for sdoy to get annual cycle +gmst interaction
    
    # gmst
    fmla.gmst       <- list(x0 ~ s(gmst, bs='tp', k=4), ~ 1 )
    fit.gmst        <- mqgam(fmla.gmst, data=data01, qu=c(0.5,0.9))
    fit.gmst$fmla   <- fmla.gmst
    fit.gmst$ptiles <- c(0.5,0.9)
    q50       <- qdo(fit.gmst, 0.5, predict)
    data01$x1 <- data01$x0 - q50
        # plot(data01$time, data01$x0, pch=20, cex=.3)
        # points(data01$time, data01$x2+mean(data01$x), pch=20, cex=.3, col=4)

    # year as factor, using 90th quantile to represent summer annual variability
    fmla.year       <- list(x1 ~ year, ~ 1 ) 
    fit.year        <- mqgam(fmla.year, data=data01, qu=c(0.5,0.9))
    fit.year$fmla   <- fmla.year
    fit.year$ptiles <- c(0.5,0.9)
    q90y      <- qdo(fit.year, 0.9, predict)

    data01$x <- data01$x1 - q90y    # BUT this makes winter worse interannual variability
        # plot(data01$time, data01$x0, pch=20, cex=.3)
        # points(data01$time, q90y+mean(data01$x), pch=20, cex=.3, col=2)
        # plot(data01$time, data01$x,  pch=20, cex=.3, col=2)
        # abline(h=0,col=4)

    readline("stop1b") 

    

    # ? just do summer?
    # annual cycle
    fmla.htdata <- list(x ~ s(sdoy, bs='cc',k=ms.k$doy) +ti(sdoy, gmst, bs=c('cc','tp')), ~ s(sdoy)) 
    fit.htdata  <- mqgam(fmla.htdata, data=data01, qu=do.ptiles)
    fit.htdata$fmla   <- fmla.htdata
    fit.htdata$ptiles <- do.ptiles
        # q90d      <- qdo(fit.htdata, 0.9, predict)
        # plot(data01$time, data01$x, pch=20, cex=.3)
        # points(data01$time, q90d, pch=20, cex=.3, col=2)
        # abline(h=0,col=1)


    save(file=MSconfig$files$st_htdata_qgam, fit.htdata, fit.year, fit.gmst, do.ptiles)
        
    cat(cr)
    cat(cr,"Generating q2p",cr)
    ptiles <- do.ptiles
    fn_gen_q2p(data01, st.q2p=MSconfig$files$st_htdata_q2p, fit.qgam=fit.htdata, st.qgam=MSconfig$files$st_htdata_qgam)
    cat(cr,"Generating p2q",cr)
    load(MSconfig$files$st_htdata_q2p)
    fn_gen_p2q(data01, qgam.q2p.fn, ms.thgpd.u, st.p2q=MSconfig$files$st_htdata_p2q, st.q2p=MSconfig$files$st_htdata_q2p)

    gc()

} # if doing qgam+q2p+p2q
# readline("stop2")


HERE HERE

cat(cr,"doMakeStationary: (re)loading pre-calculated QGAM",cr,MSconfig$files$st_qgam,cr)
load(MSconfig$files$st_qgam, verb=TRUE)
load(MSconfig$files$st_q2p,  verb=TRUE)

cat(cr,"Applying q2p",cr)
data01$uqgam <- fn_apply_q2p(data01, qgam.q2p.fn, st.q2p=NULL, SANITY=FALSE)

# calc MSgpd
cat(cr,"Fitting EVGAM",cr)
fn_fit_MSgpd(data01, ms.thgpd.u, fmla.MSgpd, chosen.MSgpd.name=chosen.MSgpd.name, fit.qgam=NULL, st.qgam=MSconfig$files$st_qgam, st.msgpd=MSconfig$files$st_msgpd)
load(MSconfig$files$st_msgpd, verb=TRUE) # if doing only apply MSgpd

# apply MSgpd
cat(cr,"Applying EVGAM",cr)
u1 <- fn_apply_MSgpd(data01, ms.thgpd.u, fit.qgam=NULL, st.qgam=st_qgam, st.msgpd=st_msgpd)

data01$u     <- u1

### general testing plots
    stin   <- MSconfig$files$st_qgam
    stdiag <- paste(dirname(stin),'plots',sub('MSqgam.RData','general.pdf',basename(stin)),sep='/')
    fn_diag_general(savefile=stdiag, width=10, height=7, DOPAUSE=FALSE)


# st.u0 <- paste(MSconfig$files$MSSAVEDIR,sub('.RData','_u0.RData',basename(MSconfig$files$st_base)),sep='' )
save(file=st_msdata01, data01, data01.std.param)
cat(cr,"Saved msdata",cr,st_msdata01,cr)


if(DODIAGQGAM) {
cat(cr,"DODIAGQGAM",cr)
    stin <- MSconfig$files$st_qgam
    stdiag <- paste(dirname(stin),'plots',sub('.RData','.pdf',basename(stin)),sep='/')
    fn_diag_qgam(fit.qgam, savefile=stdiag, width=10, height=7, DOPAUSE=FALSE)
}

if(DODIAGMSGPD) {
cat(cr,"DODIAGMSGPD",cr)
    stin <- MSconfig$files$st_msgpd
    stdiag <- paste(dirname(stin),'plots',sub('.RData','.pdf',basename(stin)),sep='/')
    fn_diag_evgam(savefile=stdiag, st.msgpd=stin, width=10, height=7, DOPAUSE=FALSE)
}

if(DODIAGQ2P) {
cat(cr,"DODIAGQ2P",cr)
    stin <- MSconfig$files$st_q2p
    stdiag <- paste(dirname(stin),'plots',sub('.RData','.pdf',basename(stin)),sep='/')
    fn_diag_q2p(savefile=stdiag, width=10, height=7, DOPAUSE=FALSE)
}

tidy()

cat(cr,cr)
cat("###############################################################################",cr)
sysdate   <- Sys.time()
stsysdate <- print(sysdate)
cat("All Done doMakeStationary",cr)
cat("###############################################################################",cr,cr)






#
