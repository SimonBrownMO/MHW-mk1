# source("doFitInitEvent.R")

list2env(HDconfig$initI , envir = .GlobalEnv)

# load(st_events, verb=TRUE)
# > str1(events01)
    # List of 5
    #  $ obs          :List of 16
    #  $ mod          :List of 16
    #  $ events01_ally:List of 2
    #  $ info         :List of 6
    #  $ hotseas      :List of 3
    # > str1(events01$obs)
        # List of 16
        #  $ ch.pk       :List of 6
        #  $ ch.st.l     : num [1:36, 1:114] 1.44 1.73 2.37 2.88 3.73 ...
        #  $ ch.st.o     : num [1:36, 1:114] 16.3 16.4 16.6 16.7 16.9 ...
        #  $ ch.doy      : int [1:36, 1:114] 226 227 228 229 230 231 232 233 234 235 ...
        #  $ ch.age      : num [1:36, 1:114] 0 1 2 3 4 5 6 7 8 9 ...
        #  $ ch.i        : int [1:36, 1:114] 5705 5706 5707 5708 5709 5710 5711 5712 5713 5714 ...
        #  $ ch.time     :List of 114
        #  $ ch.sev      : num [1:114] 36.3 39.3 62.2 53.7 49.3 ...
        #  $ ch.mean.sev : num [1:114] 2.79 2.45 2.39 2.83 2.24 ...
        #  $ ch.pkval    : num [1:114] 7.27 6.78 6.07 5.98 5.92 ...
        #  $ ch.pkval.day: num [1:114] 8 9 14 11 11 4 14 9 13 8 ...
        #  $ ch.duration : num [1:114] 13 16 26 19 22 12 25 19 27 13 ...
        #  $ ch.gmst     : num [1:114] -0.743 -0.4405 -0.382 -0.0668 -0.9448 ...
        #  $ ch.th.u     : num 0.9
        #  $ hotseas     :List of 3
        #  $ time        : num [1:16723] 1980 1980 1980 1980 1980 ...

### fit init event model to obs & model jointly
InitEventModel <- fn_fitInitEvent(initI$fmla, events01$events01_ally, data01, initI, DOHOTSEASON=DOHOTSEASON, IsJoint=TRUE, SAVEINITEVENTDIAG=initI$SAVEINITEVENTDIAG)
save(file=HDconfig$files$st_InitEvent, InitEventModel, initI)

load(st_InitEvent, verb=TRUE)
st.pdf.initev <- paste(dirname(HDconfig$initI$st_InitEvent),'plots',sub('.RData','.pdf',basename(HDconfig$initI$st_InitEvent)),sep='/')
if(DOINITEVENTPLOT) fn_plotInitEvent(InitEventModel$Model, main0='Init ~', st.pdf=st.pdf.initev)

#
