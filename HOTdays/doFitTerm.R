# source("doFitTerm.R")
# library(Rcpp) #, lib="/home/h03/hadsx/extremes/R/packages")
# library(evgam) #, lib="/home/h03/hadsx/extremes/R/packages")
# library(mgcv)
# library(qgam)
# library(ncdf4)
# library(PCICt)
# source("../libs/fn_JointMakeStationary.R")
# source("../libs/fn_JointHotDays.R")
# source("../libs/lib_HotDay.R")
# source("/home/users/simon.brown/code/R/libs/Rutils/sjb_colours.R")
# source("/home/users/simon.brown/extremes/conditional_extremes/HeffTawn_Hierarchical/BV_HT_Hier.R")

load(st_msdata01, verb=TRUE)
load(st_events,verb=TRUE)

scale.doy    <- list(obs=data01.std.param$doy$o_max, mod=data01.std.param$doy$m_max)
# stcovs       <- c('class', 'CC', 'DOY', 'Age')

### do GAM
# fmla1N       <- list(day1=termGAM.fmla$hseas$day1$TICDiCD, dayN=termGAM.fmla$hseas$dayN$TICDAiCD) # fit chosen
fmla1N       <- termGAM.fmla$hseas # fit all
st.diag.term <- paste(dirname(HDconfig$term$st_Term),'diag',sub('.RData','.GAMdiag.txt',basename(HDconfig$term$st_Term)),sep='/')
termGAM      <- fn_fitJointTermGAM(events01, fmla1N, scale.doy, stdiag=st.diag.term)

### do GLM
# fmla1N       <- list(day1=termGLM.fmla$hseas$day1$TICD, dayN=termGLM.fmla$hseas$dayN$TICDA) # fit chosen
fmla1N       <- termGLM.fmla$hseas # fit all
st.diag.term <- paste(dirname(HDconfig$term$st_Term),'diag',sub('.RData','.GLMdiag.txt',basename(HDconfig$term$st_Term)),sep='/')
termGLM      <- fn_fitJointTermGLM(events01, fmla1N, scale.doy, stdiag=st.diag.term)

# save
save(file=HDconfig$term$st_Term, termGAM, termGLM)

# plot
st.pdf.term <- paste(dirname(HDconfig$term$st_Term),'plots',sub('.RData','.GAMdiag.pdf',basename(HDconfig$term$st_Term)),sep='/')
if(HDconfig$term$DOTERMPLOT) fn_termJointGAMPlot2(termGAM, st.pdf=st.pdf.term)

st.pdf.term <- paste(dirname(HDconfig$term$st_Term),'plots',sub('.RData','.GLMdiag.pdf',basename(HDconfig$term$st_Term)),sep='/')
if(HDconfig$term$DOTERMPLOT) fn_termJointGLMPlot2(termGLM, st.pdf=st.pdf.term)

cat("doFitTerm all done",cr)
