# source("doFitInitEvent.R")

list2env(HDconfig$initI , envir = .GlobalEnv)

load(st_events, verb=TRUE)
# > str1(events)
# List of 4
#  $ obs    :List of 13
#  $ mod    :List of 13
#  $ hotseas:List of 3
#  $ info   :List of 5
#        > str1(events$obs)
#        List of 13
#         $ ch.pk       :List of 6
#         $ ch.st.l     : num [1:19, 1:749] 0.446 2.349 8.846 2.515 0.79 ...
#         $ ch.st.o     : num [1:19, 1:749] 15.4 18.3 23.1 18.1 15.5 ...
#         $ ch.doy      : int [1:19, 1:749] 302 303 304 305 306 NA NA NA NA NA ...
#         $ ch.age      : num [1:19, 1:749] 0 1 2 3 0 NA NA NA NA NA ...
#         $ ch.i        : int [1:19, 1:749] 20026 20027 20028 20029 20030 NA NA NA NA NA ...
#         $ ch.time     :List of 749
#         $ ch.sev      : num [1:749] 7.35 16.64 18.27 14.63 5.65 ...
#         $ ch.mean.sev : num [1:749] 2.45 3.33 2.03 3.66 5.65 ...
#         $ ch.pkval    : num [1:749] 8.85 8.75 8.42 7.97 7.77 ...
#         $ ch.pkval.day: num [1:749] 2 3 7 4 1 2 2 3 5 3 ...
#         $ ch.duration : num [1:749] 3 5 9 4 1 4 2 4 17 5 ...
#         $ ch.th.u     : num 0.94

### fit init event model to obs & model singly and jointly
# covariate         <- initI$ievent.covaraite
# SAVEINITEVENTDIAG <- initI$SAVEINITEVENTDIAG
HS_InitEvent <- fn_fitInitEvent(events01, data01, initI, IsJoint=TRUE, SAVEINITEVENTDIAG=SAVEINITEVENTDIAG)
InitModel    <- HS_InitEvent$InitModel
initI        <- HS_InitEvent$initI
save(file=HDconfig$files$st_InitEvent, InitModel, initI)

load(st_InitEvent, verb=TRUE)
st.pdf.initev <- paste(dirname(HDconfig$initI$st_InitEvent),'plots',sub('.RData','.pdf',basename(HDconfig$initI$st_InitEvent)),sep='/')
if(DOINITEVENTPLOT) fn_plotInitEvent(InitModel$joint$BDT, ylim=c(0,0.1), main0='Init ~', st.pdf=st.pdf.initev)

#
