# source("explore_HT.R")


### Use this file to determine the form of Initialising an event model
st.pwd <- system("pwd", intern=TRUE)
# source(paste(st.pwd,"/../../setup_all.R",sep=''))
source("setup_dev.R")


    MSref.file <- "/home/users/simon.brown/extremes/heatwaves/mhw/DATA/UKV/v5/MSref/ostia_cdr_nrt_regions.MSref.2026-03-26.RData"
    load(MSref.file, verb=TRUE)

    list2env(MSconfig ,       envir = .GlobalEnv)
    list2env(MSconfig$files , envir = .GlobalEnv)

    load(st_msdata01, verb=TRUE)

    source(paste(st.pwd,"/../../make_stationary/setup_MakeStationary.R",sep=''))

    st_events  <- list.files(glue(dirname(MSref.file),"/../Events"), pattern='_EventsTh', full.names=TRUE)
    st_events  <- st_events[grep('.RData', st_events)]
    load(st_events, verb=TRUE)
    list2env(events01$info , envir = .GlobalEnv)

    source(paste(st.pwd,"/../setup_HotDays.R",sep=''))

readline("Stop1")

### doFitHt.R ##################################### FAST

ht.th.u   <- fitht$ht.th.u
cr.th.u   <- fitht$cr.th.u
scale.doy <- list(obs=(366+1), mod=(360+1)) # bit of a bodge as assumes all obs are from leap years
stcovs    <- c('Index', 'CC', 'DOY', 'Age')

htJ <- fn_fitJointHt(events01, ht.th.u, cr.th.u, scale.doy, IsJoint=TRUE)

# save(file=HDconfig$fitht$st_Ht, htJ, ht.th.u, cr.th.u, scale.doy)

## not needed
# SjbPlotHtParametersCov(htJ$htN,new=TRUE)

## not needed
#  lag01 <- fn_lag_events01(events01, scale.doy)
#     list2env(lag01$om,env=parent.frame())
# SjbPlotHtModelPoints(htJ$htN, om12[om.1E1,])
st.pdf.ht <- paste(dirname(HDconfig$fitht$st_Ht),'plots',sub('.RData','.pdf',basename(HDconfig$fitht$st_Ht)),sep='/')

if(HDconfig$fitht$DOHTPLOT) fn_htJointPlot(htJ, stplot=st.pdf.ht)

















#