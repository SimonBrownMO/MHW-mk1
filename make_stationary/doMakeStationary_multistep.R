# source("doMakeStationary_multistep.R")

### 2026.03.26 
###    - have settled on a multistep approach
###    - step 1: remove annual cycle and climate change term from the mean
###         - option 1: Two step
###                     A: x      ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=5, m=2) +ti(sdoy,gmst,bs=c("cc","ts")) 
###                     B: residA ~ s(stime, bs='ts', k=2*ny.om, by=ftype)
###         - option 2: One step with by-year random effect
###                     AB: x     ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs='ts', k=5, m=2) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYearOM, bs="sz")
###    - step 2: QGAM as normal
###    - step 3: EVGAM as normal

### general setup
# data01[,ftype   := factor(isobs, levels=c(1,0), labels=c("obs","mod"))]
# data01[,fYearOM := factor(c(paste("o",data01$year[iob],sep=""), paste("m",data01$year[-iob],sep="")))]
data01$ftype   <- factor(data01$isobs, levels=c(1,0), labels=c("obs","mod"))
iob            <- which(data01$ftype=="obs")
d01.years      <- trunc(data01$time)
data01$fYearOM <- factor(c(paste("o",d01.years[iob],sep=""), paste("m",d01.years[-iob],sep="")))
ny.om          <- length(unique(data01$fYearOM[iob])) + length(unique(data01$fYearOM[-iob])) # total number of obs years + model years

##################################################################################
###    - step 1: remove annual cycle and climate change term from the mean
###         - option 1: Two step
###                     A: x      ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs="ts", k=5, m=2) +ti(sdoy,gmst,bs=c("cc","ts")) 
###                     B: residA ~ s(stime, bs='ts', k=2*ny.om, by=ftype)
###         - option 2: One step with by-year random effect
###                     AB: x     ~ ftype +s(sdoy,bs="cc",k=28,by=ftype) +s(gmst,bs="ts", k=5, m=2) +ti(sdoy,gmst,bs=c("cc","ts"))   +s(sdoy, fYearOM, bs="sz")
##################################################################################
if(DOSTEP1) {
    
    # if(STEP1 == "TWOSTEP") {
        ## fit 1
        ft1A.gam       <- gam(fmla1A, data=data01, method="GCV.Cp", select=TRUE) 

        ## fit 2
        data01$resid1A    <- resid(ft1A.gam)
        ft1B.gam          <- gam(fmla1B, data=data01, method="GCV.Cp", select=TRUE) 
        data01$resid1B    <- resid(ft1B.gam)

        cat("################################",cr)
        cat("Summary ft1A.gam",cr)
        print(summary(ft1A.gam))
        cat("################################",cr,cr)
        cat("################################",cr)
        cat("Summary ft1B.gam",cr)
        print(summary(ft1B.gam))
        cat("################################",cr,cr)

        plot_gam1(ft1A.gam, stpdf=paste(MSSAVEDIR,"step1/plots/GAM_step1_ft1A.pdf",sep=''))
        plot_gam2(ft1B.gam, stpdf=paste(MSSAVEDIR,"step1/plots/GAM_step1_ft1B.pdf",sep=''))

        # save(ft1A.gam, ft1B.gam, data01, file=paste(MSSAVEDIR,"step1/ft1_step1.RData",sep=''))

    # } else 
    # if(STEP1 == "ONESTEP") {

        ft0       <- gam(fmla0, data=data01, method="GCV.Cp", select=TRUE, fit=FALSE)
        sp.sz     <- rep(1, length(which(grepl("sdoy,fYearOM", names(ft0$sp)))))                  
        rm(ft0)
        ft1AB.gam         <- gam(fmla1AB, data=data01, method="GCV.Cp", select=TRUE) 
        data01$residStep1 <- resid(ft1AB.gam)

        cat("################################",cr)
        cat("Summary ft1AB.gam",cr)
        print(summary(ft1AB.gam))
        cat("################################",cr,cr)
                  
        plot_gam1(ft1AB.gam, stpdf=paste(MSSAVEDIR,"step1/plots/GAM_step1_ft1AB.pdf",sep=''))

        save(ft1A.gam, ft1B.gam, ft1AB.gam, data01, file=paste(MSSAVEDIR,"step1/ft1_step1.RData",sep=''))

    # }

} else {
    load(paste(MSSAVEDIR,"step1/ft1_step1.RData",sep=''),verb=TRUE)
}


##################################################################################
###    - step 2: QGAM as normal
##################################################################################
### subsample if too large as qgam fitting can be slow
SUBOPTION <- 2
if(length(data01$x)>nmax_qgam) {
    cat(cr,"Subsampling data for QGAM fitting to",nmax_qgam,"points",cr)
    ## what the policy should be is up for debate
    ## to characterise biases in climate models we need to maximise sampling in the overlap pierod
    ## to characterise climate chage terms we need to maximise sampling in the future period
    ## data i the early observed period is less useful for either
    ## OPTION 1 - changing this to 50% in observed period of which 50:50 is obs model, and 50% future
    ## OPTION 2 - 10% obs pre-model, 50% obs model overlap, 40% future model
    if(SUBOPTION==1) {
        # obs period
        r.o.time <- range(data01$time[which(data01$isobs==1)])
        io       <- which(data01$isobs==1 & data01$time>=r.o.time[1] & data01$time<=r.o.time[2])
        im       <- which(data01$isobs==0 & data01$time>=r.o.time[1] & data01$time<=r.o.time[2])
        i1.o     <- sample(io, (length(io)/length(c(io,im))) * nmax_qgam/2, replace=FALSE)
        i1.m     <- sample(im, (length(im)/length(c(io,im))) * nmax_qgam/2, replace=FALSE)
        # future period
        im       <- which(data01$isobs==0 & data01$time>r.o.time[2])
        i2.m     <- sample(im, nmax_qgam/2, replace=FALSE)
        i0qgam   <- c(i1.o,i1.m,i2.m)
    } else if(SUBOPTION==2) {
        # obs pre-model period
        r.o.time <- range(data01$time[which(data01$isobs==1)])
        r.m.time <- range(data01$time[which(data01$isobs==0)])
        io       <- which(data01$isobs==1 & data01$time>=r.o.time[1] & data01$time<r.m.time[1])
        nsamp    <- nmax_qgam/10
        if(length(io)>nsamp) i1.o <- sample(io, nsamp, replace=FALSE) else i1.o <- io
        # obs model overlap period
        io       <- which(data01$isobs==1 & data01$time>=r.m.time[1] & data01$time<=r.o.time[2])
        im       <- which(data01$isobs==0 & data01$time>=r.m.time[1] & data01$time<=r.o.time[2])
        nsamp    <- nmax_qgam/4  # 25% + 25% = 50%
        if(length(io)>nsamp) i2.o <- sample(io, nsamp, replace=FALSE) else i2.o <- io
        if(length(im)>nsamp) i2.m <- sample(im, nsamp, replace=FALSE) else i2.m <- im
        # future period
        im       <- which(data01$isobs==0 & data01$time>r.o.time[2])
        nsamp    <- nmax_qgam - sum(c(length(i1.o), length(i2.o), length(i2.m)))
        if(length(im)>nsamp) i3.m <- sample(im, nsamp, replace=FALSE) else  i3.m <- im
        # combine
        i0qgam   <- c(i1.o,i2.o,i2.m,i3.m)
    }
} else i0qgam <- 1:length(data01$x)

### fit QGAM
fit.qgam <- fn_fit_qgam(fmla.qgam, data01[i0qgam,], ptiles=do.ptiles)
    invisible(qdo(fit.qgam, ms.thgpd.u, plot, page=1, shade=TRUE, seWithMean=TRUE, main="90th percentile"))

fit.qgam$i0qgam <- i0qgam
if (!dir.exists(dirname(st_qgam))) dir.create(dirname(st_qgam), recursive = TRUE)
save(file=st_qgam, fit.qgam)

### generate q2p functions
cat(cr,"Generating q2p",cr)
fn_gen_q2p(data01, MSconfig$files$st_q2p, fit.qgam=fit.qgam)

### generate p2q function
cat(cr,"Generating p2q",cr)
fn_gen_p2q(data01, NULL, ms.thgpd.u, MSconfig$files$st_p2q, st.q2p=MSconfig$files$st_q2p)

### apply qgam
cat(cr,"doMakeStationary: (re)loading pre-calculated QGAM",cr,MSconfig$files$st_qgam,cr)
load(MSconfig$files$st_qgam, verb=TRUE)
load(MSconfig$files$st_q2p,  verb=TRUE)

cat(cr,"Applying q2p",cr)
data01$uqgam <- fn_apply_q2p(data01$residStep1, qgam.q2p.fn, st.q2p=NULL, SANITY=FALSE)


gc()


##################################################################################
###    - step 3: EVGAM as normal
##################################################################################
fn_fit_MSgpd(data01, ms.thgpd.u, fmla.MSgpd, chosen.MSgpd.name=chosen.MSgpd.name, fit.qgam=NULL, st.qgam=MSconfig$files$st_qgam, st.msgpd=MSconfig$files$st_msgpd)
    load(MSconfig$files$st_msgpd, verb=TRUE) # if doing only apply MSgpd



# apply MSgpd
    cat(cr,"Applying EVGAM",cr)
    u1       <- fn_apply_MSgpd(data01, st.msgpd=MSconfig$files$st_msgpd)
    data01$u <- u1

    # st.u0 <- paste(MSconfig$files$MSSAVEDIR,sub('.RData','_u0.RData',basename(MSconfig$files$st_base)),sep='' )
    save(file=st_msdata01, data01)
    cat(cr,"Saved msdata",cr,st_msdata01,cr)
readline("stop)")

    ### general testing plots
        stin   <- MSconfig$files$st_qgam
        stdiag <- paste(dirname(stin),'plots',sub('MSqgam.RData','general.pdf',basename(stin)),sep='/')
        fn_diag_general_twostep(savefile=stdiag, width=10, height=7, DOPAUSE=FALSE)


    if(DODIAGQGAM) {
    cat(cr,"DODIAGQGAM",cr)
        stin <- MSconfig$files$st_qgam
        stdiag <- paste(dirname(stin),'plots',sub('.RData','.pdf',basename(stin)),sep='/')
        fn_diag_qgam_twostep(fit.qgam, savefile=stdiag, width=10, height=7, DOPAUSE=FALSE)
    }

    if(DODIAGMSGPD) {
    cat(cr,"DODIAGMSGPD",cr)
        stin <- MSconfig$files$st_msgpd
        stdiag <- paste(dirname(stin),'plots',sub('.RData','.pdf',basename(stin)),sep='/')
        fn_diag_evgam_twostep(savefile=stdiag, st.msgpd=stin, width=10, height=7, DOPAUSE=FALSE)
    }

    if(DODIAGQ2P) {
    cat(cr,"DODIAGQ2P",cr)
        stin <- MSconfig$files$st_q2p
        stdiag <- paste(dirname(stin),'plots',sub('.RData','.pdf',basename(stin)),sep='/')
        fn_diag_q2p_twostep(savefile=stdiag, width=10, height=7, DOPAUSE=FALSE)
    }

    tidy()
