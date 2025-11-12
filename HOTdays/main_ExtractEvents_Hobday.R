# source("main_ExtractEvents_Hobday.R")

### main_ExtractEvents_Hobday.R extracts events from anomalies  on the observed scale
### how the anomalies are calculated is specified below

st.pwd    <- system("pwd", intern=TRUE)

### setup variables #####################################
source(paste(st.pwd,"/../setup_all.R",sep=''))
event.th.u    <- 0.90    # definition of when a heatwave event begins/ends
event.length  <- 80      # maximum allowed length of heatwave.  For wrt-peak this is fwd and bwd so events can be up to 2x longer than this (days)
USEHOTSEASON  <- TRUE
FINDHOTSEASON <- TRUE    # TRUE== determine hot season based on annual cycle greater than thHotSeason
thHotSeason   <- 5.0     # annual range degC must exceed this if a hot season is to be defined
                        # for regions with a strong annual cycle it is a good idea to exclude weather types that are not associated with hot days in an absolute sense
facHotSeas    <- 0.25   # fraction of the year that defines the hot season. 
                        # E.g. 0.5 = half the year, 0.25 = quarter of the year eg ~JJA
eventsDOY     <- NULL   # Specify here desired DOY for events to be extracted if find FINDHOTSEASON is FALSE
DOEVENTSPLOT  <- TRUE

clim.year     <- 1996  # middle year of climatology period

events                <- list()
events$DOEVENTSPLOT   <- DOEVENTSPLOT
events$event.th.u     <- event.th.u
events$event.length   <- event.length
events$USEHOTSEASON   <- USEHOTSEASON
events$FINDHOTSEASON  <- FINDHOTSEASON
events$thHotSeason    <- thHotSeason

### currently missing from MSconfig
do.region <- "NWS"  

### v_Hobday
# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v_Hobday/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
# st_sub     <- "_Hobday"
# DOHOBWTIME <- FALSE

### v1  _hobwt    trend in annual cycle
MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v1/MSref/ostia_cdr_nrt_regions.MSref.2025-11-07.RData"
st_sub     <- NULL
DOHOBWTIME <- TRUE

### v2
# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v2/MSref/ostia_cdr_nrt_regions.MSref.2025-10-29.RData"

### v3
# MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/NWS/v3/MSref/ostia_cdr_nrt_regions.MSref.2025-11-04.RData"


cat(cr,"MSref.file",cr)
cat(MSref.file,cr,cr)
load(MSref.file, verb=TRUE)
list2env(MSconfig , envir = .GlobalEnv)
list2env(MSconfig$files , envir = .GlobalEnv)

st_events <- sub('.RData','_EventsThMMM.RData',basename(st_base))
st_events <- paste(MSSAVEDIR,'Events/',sub('MMM',event.th.u,st_events),sep='')
st_events <- sub('.RData',paste(st_sub,'.RData',sep=''), st_events)

readline("Continue?")

# if(st_version=='v1') {
#     load(st_msdata01, verb=TRUE)
#     data01$doy <- data01$doy + 0.5
#     st_msdata01 <- sub('/MSdata01/','/v1/MSdata01/',st_msdata01)
#     st_qgam     <- sub('/Qgam/','/v1/Qgam/',         st_qgam)
#     MSSAVEDIR   <- sub('RRR',do.region, sub('VVV',  st_version,    MSSAVEDIR0))
# }

### setup HotDays variables #####################################
    # source(paste(st.pwd,"/setup_HotDays.R",sep=''))
    # cat("HotDays main: Completed setup_HotDays.R",cr,"#########################",cr,cr)

# if(st_version=='v1') {
#     st_msdata01 <- sub('/MSdata01/','/v1/MSdata01/',st_msdata01)
#     MSSAVEDIR   <- sub('RRR',do.region, sub('VVV',st_version,    MSSAVEDIR0))
# }


    load(st_msdata01, verb=TRUE)
    cat("Loaded msdata01 from ", st_msdata01, cr,cr)

    
### extract events #####################################
    if(!DOHOBWTIME)  { ################################################################################
        source(paste(st.pwd,"/doExtractEvents_Hobday.R",sep=''))

        ### apply Hobday catogries ##############################
        # use the th90-th50 range to define categories

        load(st_events,   verb=TRUE)
        th.range <- th.hob$th90 - th.hob$th50
            # plot(data01$time, th.range, pch=20, cex=.3, main='Hobday th90-th50 range')

        # go through each event and assign categories to each day of the event
        all.events <- events01$events01_ally$obs
        hob_cat    <- all.events$ch.st.o
        hob_cat[]  <- NA
        for(k in seq_along(all.events$ch.sev)) {
            j1          <- all.events$ch.i[,k]
            hob_cat[,k] <- ceiling(all.events$ch.st.o[,k]/th.range[j1])
        }
        events01.hob$events01_ally$obs$hob_cat <- hob_cat
        save(file=st_events, events01.hob, th.qgam, th.hob, data01_an_hob)
        cat("Saved events to :", st_events, cr, cr )

        ###  plots  #####################################################
        ###  wrt time 
        load(st_events,   verb=TRUE)
        c21      <- distinct21colours()
        # all years observed values
        isum <- seq_along(data01$doy)
        r1   <- range(data01$x[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
        plot(data01$time[isum], data01$x[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Observed temperature RRR') )    
        lines(data01$time[isum], th.hob$th90[isum], col='red', lwd=1)
        for(k in seq_along(events01$obs$ch.sev)) {
            ch.t <- events01$obs$ch.time[[k]]
            points(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1])
            lines(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
        }
        grid()
        readline("Continue?")

        # 1994 2003 observed values
        isum <- which(data01$doy >0 & data01$time>=1995 & data01$time<2004)
        r1   <- range(data01$x[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
        plot(data01$time[isum], data01$x[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Observed temperature RRR') )    
        lines(data01$time[isum], th.hob$th90[isum], col='red', lwd=1)
        for(k in seq_along(events01$obs$ch.sev)) {
            ch.t <- events01$obs$ch.time[[k]]
            points(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1])
            lines(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
        }
        grid()
        readline("Continue?")

        # all year anomalies
        isum <- which(data01$doy >0 )
        r1   <- range(data01$an_hob[finite(events01$events01_ally$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
        plot(data01$time[isum], data01$an_hob[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Observed temperature RRR') )    
        for(k in seq_along(events01$events01_ally$obs$ch.sev)) {
            ch.t <- events01$events01_ally$obs$ch.time[[k]]
            points(ch.t, data01$an_hob[finite(events01$events01_ally$obs$ch.i[,k])], col=c21[k %% 21 +1])
            lines(ch.t, data01$an_hob[finite(events01$events01_ally$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
        }
        grid()
        readline("Continue?")

        # 1994 2003 anomalies
        isum <- which(data01$doy >0 & data01$time>=1995 & data01$time<2004)
        r1   <- range(data01$an_hob[finite(events01$events01_ally$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
        plot(data01$time[isum], data01$an_hob[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Observed temperature RRR') )    
        for(k in seq_along(events01$events01_ally$obs$ch.sev)) {
            ch.t <- events01$events01_ally$obs$ch.time[[k]]
            points(ch.t, data01$an_hob[finite(events01$events01_ally$obs$ch.i[,k])], col=c21[k %% 21 +1])
            lines(ch.t, data01$an_hob[finite(events01$events01_ally$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
        }
        grid()

        
    } else { ################################################################################
        source(paste(st.pwd,"/doExtractEvents_Hobday_wtime.R",sep=''))

        ### apply Hobday catogries ##############################
        # use the th90-th50 range to define categories

        load(st_events,   verb=TRUE)
        th.range <- th.hobwt$th90 - th.hobwt$th50
            # plot(data01$time, th.range, pch=20, cex=.3, main='Hobday th90-th50 range')

        # go through each event and assign categories to each day of the event
        all.events <- events01$events01_ally$obs
        hobwt_cat    <- all.events$ch.st.o
        hobwt_cat[]  <- NA
        for(k in seq_along(all.events$ch.sev)) {
            j1          <- all.events$ch.i[,k]
            hobwt_cat[,k] <- ceiling(all.events$ch.st.o[,k]/th.range[j1])
        }
        events01.hobwt$events01_ally$obs$hobwt_cat <- hobwt_cat
        save(file=st_events, events01.hobwt, th.hobwt, data01_an_hobwt)
        cat("Saved events to :", st_events, cr, cr )

        ###  plots  #####################################################
        ###  wrt time 
        load(st_events,   verb=TRUE)
        c21      <- distinct21colours()
        # all years observed values
        isum <- seq_along(data01$doy)
        r1   <- range(data01$x[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
        plot(data01$time[isum], data01$x[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Observed temperature RRR') )    
        lines(data01$time[isum], th.hobwt$th90[isum], col='red', lwd=1)
        for(k in seq_along(events01$obs$ch.sev)) {
            ch.t <- events01$obs$ch.time[[k]]
            points(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1])
            lines(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
        }
        grid()
        readline("Continue?")

        # 1994 2003 observed values
        isum <- which(data01$doy >0 & data01$time>=1995 & data01$time<2004)
        r1   <- range(data01$x[finite(events01$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
        plot(data01$time[isum], data01$x[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Observed temperature RRR') )    
        lines(data01$time[isum], th.hobwt$th90[isum], col='red', lwd=1)
        for(k in seq_along(events01$obs$ch.sev)) {
            ch.t <- events01$obs$ch.time[[k]]
            points(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1])
            lines(ch.t, data01$x[finite(events01$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
        }
        grid()
        readline("Continue?")

        # all year anomalies
        isum <- which(data01$doy >0 )
        r1   <- range(data01$an_hobwt[finite(events01$events01_ally$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
        plot(data01$time[isum], data01$an_hobwt[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Observed temperature RRR') )    
        for(k in seq_along(events01$events01_ally$obs$ch.sev)) {
            ch.t <- events01$events01_ally$obs$ch.time[[k]]
            points(ch.t, data01$an_hobwt[finite(events01$events01_ally$obs$ch.i[,k])], col=c21[k %% 21 +1])
            lines(ch.t, data01$an_hobwt[finite(events01$events01_ally$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
        }
        grid()
        readline("Continue?")

        # 1994 2003 anomalies
        isum <- which(data01$doy >0 & data01$time>=1995 & data01$time<2004)
        r1   <- range(data01$an_hobwt[finite(events01$events01_ally$obs$ch.i)], na.rm=TRUE )   # range(data01$x[isum] )
        plot(data01$time[isum], data01$an_hobwt[isum], pch=20, cex=.3, ylim=r1, main=sub('RRR', do.region, 'Observed temperature RRR') )    
        for(k in seq_along(events01$events01_ally$obs$ch.sev)) {
            ch.t <- events01$events01_ally$obs$ch.time[[k]]
            points(ch.t, data01$an_hobwt[finite(events01$events01_ally$obs$ch.i[,k])], col=c21[k %% 21 +1])
            lines(ch.t, data01$an_hobwt[finite(events01$events01_ally$obs$ch.i[,k])], col=c21[k %% 21 +1],lwd=2)
        }
        grid()

    }

    # iobs <- which(data01$isobs==1)
    # # if(DOEVENTSPLOT) fn_plotEvents(events01$obs, data01[iobs,], nplot=3, savefile=sub('Events/ukgd','Events/plots/ukgd',sub('.RData','.pdf',st_events)) )
    # if(DOEVENTSPLOT) fn_plotEvents(events01$obs, data01[iobs,], nplot=3, savefile=NULL )
    # cat("HotDays main: Completed doExtractEvents.R",cr,"#########################",cr,cr)







