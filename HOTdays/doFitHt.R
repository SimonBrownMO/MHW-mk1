# source("doFitHt.R")

# source("../libs/fn_JointHotDays.R")

# load(st_msdata01, verb=T)
# load(st_events,verb=TRUE)

ht.th.u   <- fitht$ht.th.u
cr.th.u   <- fitht$cr.th.u
stcovs    <- c('Index', 'CC', 'DOY', 'Age')

htJ <- fn_fitJointHt(events01, ht.th.u, cr.th.u, IsJoint=TRUE)

save(file=HDconfig$fitht$st_Ht, htJ, ht.th.u, cr.th.u)

## not needed
# SjbPlotHtParametersCov(htJ$htN,new=TRUE)

## not needed
#  lag01 <- fn_lag_events01(events01)
#     list2env(lag01$om,env=parent.frame())
# SjbPlotHtModelPoints(htJ$htN, om12[om.1E1,])
st.pdf.ht <- paste(dirname(HDconfig$fitht$st_Ht),'plots',sub('.RData','.pdf',basename(HDconfig$fitht$st_Ht)),sep='/')

if(HDconfig$fitht$DOHTPLOT) fn_htJointPlot(htJ, stplot=st.pdf.ht)

cat("HT all done",cr)
