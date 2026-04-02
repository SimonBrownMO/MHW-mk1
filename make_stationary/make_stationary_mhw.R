# source("make_stationary_mhw.R")

# sbatch --time=360 --mem=30000 -o $PWD/main_MakeStationary.R.rout --wrap="Rscript $PWD/main_MakeStationary.R"

sysdate0   <- Sys.time()
cat("###############################################################################",cr)
stsysdate0 <- print(sysdate0)
cat("starting make_stationary_mhw.R",cr)
cat("###############################################################################",cr,cr)
cr     <- '\n'
st.pwd <- system("pwd", intern=TRUE)

datestamp  <- "2026-03-26" 
st_version <-  'v5' # "v_Hobday" # 'v3' # "v1" # 
    # v_Hobday: match Hobday as best we can, event threshold 0.90, climC term linear with time
    # v1: as v_Hobday but allowing linear trend with time but fixed annual cycle
    # v2: as for v1 but allowing linear trend with GMST but fixed annual cycle  x ~ gmst  +s(sdoy, bs='cc',k=ms.k$doy) 
    # v3: standard LST model for reference - not currently advocating it
    # v4: seasonal make stationayr
    # v5: multistep approach - remove annual cycle and climate change term from the mean, then QGAM, then EVGAM
do.region  <- "UKV"  # one of:
                    # "Global_.90S_to_90N."  
                    # "Global_.60S_to_60N."  
                    # "North_Atlantic_.CR." 
                    # "North_Atlantic_.MO."  
                    # "UKV"                  
                    # "NWS"                 
                    # "North_Sea"            
                    # "N_Hem_.20N_to_60N."   
                    # "Tropics_.20S_to_20N."
                    # "S_Hem_.60S_to_20S." 

iregion <- 12  # model equivalent region in digest_mass_files.R - magic number

source(paste(st.pwd,"/setup_MakeStationary.R",sep=''))

### pre-proc data ###################################################
    
    # source(paste(st.pwd,"/pre-proc-data.R",sep=''))
    # load(file=MSconfig$files$st_preproc, verb=TRUE)

  ### Observations
    load(MSconfig$files$st_infile_o,verbose=TRUE) # l.mhw
    o.time <- l.mhw$date
    # if(max(o.x, na.rm=TRUE)>100) {o.x <- o.x - deg0C; cat('Converting OBS to DegC',cr)}
    o.info        <- list()
    o.info$years  <- as.integer(format(o.time,"%Y"))
    o.info$month  <- as.integer(format(o.time,"%m"))
    o.info$date   <- as.integer(format(o.time,"%d"))
    o.info$uyears        <- sort(unique(o.info$years))
    o.info$nyears        <- length(o.info$uyears)
    o.info$ref_year      <- o.info$uyears[1]
    o.info$doy           <- as.integer(strftime(o.time, format = "%j"))
    o.info$maxdoy        <- 'crash' # max(o.info$doy,na.rm=TRUE)
    o.info$sdoy          <- o.info$doy/(unlist(lapply(o.info$years,days_in_year))+1)  # +1 so last day of year sdoy<1
    o.info$udoy          <- sort(unique(o.info$doy))
    o.info$imidyear      <- which(o.info$doy == 182)
    o.time.std           <- as.double(o.time)
    o.info$dtx.std.param <- list(mean=mean(o.time.std), max=max(o.time.std), min=min(o.time.std))
    o.info$time.std      <- (o.time.std - mean(o.time.std))/(max(o.time.std) - min(o.time.std))
    o.info$regions       <- names(l.mhw)[which(names(l.mhw)!="date")]

    # read observed GMST
    source(st_obs_gmst)
    # years_out          <- o.info$years + o.info$doy/max(o.info$doy)
    years_out          <- o.info$years + o.info$sdoy
    gmst.o             <- get_obs_global_annual_mean_temperature(update_download=FALSE, years_out=years_out)
    igmst              <- which(trunc(gmst.o$date) %in% gmst_ref_period)
    gmst.o$global.temp <- gmst.o$global.temp - mean(gmst.o$global.temp[igmst],na.rm=TRUE)

    # # make required data structures and save pre-proc file
    # l.o     <- list(time=o.time, data=l.mhw[which(names(l.mhw)!="date")], gmst=gmst.o$global.temp, info=o.info)
    # om_data <- list(obs=l.o, mod=data01.m)

    ### prepare data01 for MakeStationary ###############################################
    o.t2      <- o.info$years
    for(iy in 1:o.info$nyears) {
        yr             <- o.info$uyears[iy]
        iiy            <- which(o.info$years == o.info$uyears[iy])
        o.t2[iiy]      <- o.info$years[iiy] + (o.info$doy[iiy] - 0.5)/days_in_year(yr)
        # cat(iy,yr,days_in_year(yr),cr)
    }
    # o.time.lim       <- range(trunc(c(o.t2)))
    # o.time.std.param <- list(mean=mean(o.time.lim), max=max(o.time.lim), min=min(o.time.lim))
    # o.time           <- c(o.t2)
    # o.time.std       <- (o.time - o.time.std.param$mean)/(o.time.std.param$max - o.time.std.param$min)

    data01.o <- data.frame(time=o.info$years+o.info$sdoy, x=c(l.mhw[[do.region]]), gmst=gmst.o$global.temp, sdoy=o.info$sdoy, doy=o.info$doy, isobs=1)

    ### diagnostic checks for pre-proc data
    # if(DODIAGPRE) {
    #     par(mfrow=c(2,3))
    #     plot(om.time,                   om_data$obs$data[[do.region]], pch=46, main=paste("Obs data:", do.region))
    #     plot(om_data$obs$info$time.std, om_data$obs$data[[do.region]], pch=46, main=paste("Obs data:", do.region))
    #     plot(om_data$obs$info$doy,      om_data$obs$data[[do.region]], pch=46, main=paste("Obs data:", do.region))
    #     plot(om_data$obs$info$sdoy,     om_data$obs$data[[do.region]], pch=46, main=paste("Obs data:", do.region))
        
    #     plot(om.time, om_data$obs$gmst, pch=20, cex=.3, main="Obs GMST")
    #     par(mfrow=c(1,1))
    # }

    # ### make data01 data frames
    # data01.regions       <- data.frame(x=c(om_data$obs$data), time=om.time, stime=om.time.std)
    # data01.regions$doy   <- c(om_data$obs$info$doy )
    # data01.regions$sdoy  <- c(om_data$obs$info$sdoy)
    # om.doy.std.param     <- list(o_max=max(om_data$obs$info$doy,na.rm=TRUE) )
    # data01.regions$gmst  <- c(om_data$obs$gmst)
    # data01.regions$isobs <- 1
    # data01.std.param     <- list(time=om.time.std.param, doy=om.doy.std.param)
    # data01.regions$uqgam <- NA

  ### End Observations


  ### CPM data
    # currently only works with one reion
    load(MSconfig$files$st_infile_m,verbose=TRUE)

    # read model GMST
    load(st_mod_gmst, verb=TRUE)  # cpm_gmst$m001$gmst

    m.sdoy <- m.sst$doy / (max(m.sst$doy) +1)  # can do this as model has no leap year # +1 so last day of year sdoy<1) 
    # CPM data only goes to 2080-11-30 12:00:00 so need to crop m.sst to match
    i1       <- which(m.sst$date         %in% cpm_gmst$m001$time) 
    i2       <- which(cpm_gmst$m001$time %in% m.sst$date[i1])
    data01.m <- data.frame(time=m.sst$time[i1]+m.sdoy[i1], x=m.sst$sst[iregion,i1], gmst=cpm_gmst$m001$gmst[i2], sdoy=m.sdoy[i1], doy=m.sst$doy[i1] )
    data01.m <- data.table(data01.m)
    data01.m[,isobs := 0]
  ### end CPM data

    save(file=                 MSconfig$files$st_preproc, data01.o, data01.m)
    cat("Pre-proc saved to :", MSconfig$files$st_preproc, cr)

### END pre-proc data ###############################################

data01 <- rbind(data01.o, data01.m)
data01 <- data.table(data01)
# change gmst reference daate to be the last observation year
iob               <- which(data01$isobs==1)
data01$gmst[iob]  <- data01$gmst[iob] - tail(data01$gmst[iob],1)
iref              <- which.min(abs(data01$time[-iob]-tail(data01$time[iob],1)))
data01$gmst[-iob] <- data01$gmst[-iob] - data01$gmst[-iob][iref]

data01$stime      <- (data01$time - 2000)/100

# plot(data01$time, data01$gmst, pch=20, cex=.3, col=data01$isobs+1)
plot(data01$time, data01$x, pch=20, cex=.3, col=data01$isobs+1)
# readline("Stop")

### call doMakeStationary.R ###############################################
source(paste(st.pwd,"/doMakeStationary_multistep.R",sep=''))
# source(paste(st.pwd,"/doMakeHTdata.R",sep=''))


cat("###############################################################################",cr)
sysdate    <- Sys.time()
stsysdate  <- print(sysdate)

cat("Time taken",cr)
print(sysdate-sysdate0)
cat("All Done main_MakeStationary.R",cr)
cat("###############################################################################",cr,cr)

#
