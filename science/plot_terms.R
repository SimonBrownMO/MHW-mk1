# source("plot_terms.R")

library(mgcv)
library(data.table)
library(gratia)
source("copilot_fn.R")

# load( "./RData/ft1_regamJ2m1c.RData", verb=TRUE)
# pdf("./plots/terms_regamJ2m1c.pdf", width=12, height=10)
# ft1 <- ft1.regamJ2m1c
load( "./RData/ft1_regamJ4o.RData", verb=TRUE)
pdf("./plots/terms_regamJ4o.pdf", width=12, height=10)
ft1 <-         ft1.regamJ4o

# ftype
# s(sdoy):ftypeobs
# s(sdoy):ftypemod
# s(gmst)
# ti(sdoy,gmst)
# s(sdoy,fYear)
st.terms <- c("ftype",smooths(ft1))
iob      <- which(data01$isobs==1)

plot_gam_terms(ft1, terms=st.terms[2:5], main='All terms',data=data01[-iob,], sort_by="time", col="red",  SE=FALSE, xlim=range(data01$time))
plot_gam_terms(ft1, terms=st.terms[2:5], main='All terms',data=data01[iob,],  sort_by="time", col="blue", SE=FALSE, add=TRUE)
grid()

plot_gam_terms(ft1, terms="s(gmst)", data=data01[iob,],  sort_by="time", add=FALSE, col="blue", ylim=c(-5,6), xlim=range(data01$time) )
plot_gam_terms(ft1, terms="s(gmst)", data=data01[-iob,], sort_by="time", add=TRUE,  col="red")
grid()

plot_gam_terms(ft1, terms="ti(sdoy,gmst)",    data=data01[iob,],  sort_by="time", add=FALSE, col="blue", ylim=c(-5,6), xlim=range(data01$time) )
plot_gam_terms(ft1, terms="ti(sdoy,gmst)",    data=data01[-iob,], sort_by="time", add=TRUE,  col="red")
plot_gam_terms(ft1, terms="s(sdoy):ftypemod", data=data01[-iob,], sort_by="time", add=TRUE, LINES=FALSE, SE=FALSE,  col=rgb(0.7,0,0,0.2))
grid()

plot_gam_terms(ft1, terms="ti(sdoy,gmst)", data=data01[iob,],  sort_by="sdoy", add=FALSE, LINES=FALSE, SE=FALSE, col="blue", ylim=c(-5,6), xlim=range(data01$sdoy) )
plot_gam_terms(ft1, terms="ti(sdoy,gmst)", data=data01[-iob,], sort_by="sdoy", add=TRUE, LINES=FALSE,  SE=FALSE, col="red")
grid()

plot_gam_terms(ft1, terms="s(sdoy):ftypeobs", data=data01[iob,],  sort_by="sdoy", add=FALSE, col="blue", ylim=c(-5,6), xlim=range(data01$sdoy) )
plot_gam_terms(ft1, terms="s(sdoy):ftypemod", data=data01[-iob,], sort_by="sdoy", add=TRUE,  col="red")
i1 <- which(data01$isobs==1 & data01$year %in% c(1981,2024,2079))
i2 <- which(data01$isobs==0 & data01$year %in% c(1981,2024,2079))
plot_gam_terms(ft1, terms="ti(sdoy,gmst)", data=data01[i1,], sort_by="sdoy", LINES=FALSE, SE=FALSE, add=TRUE, col=rgb(0,0.6,0.6,0.2))
plot_gam_terms(ft1, terms="ti(sdoy,gmst)", data=data01[i2,], sort_by="sdoy", LINES=FALSE, SE=FALSE, add=TRUE, col=rgb(0.7,0,0,0.2) )
grid()


plot_gam_terms(ft1, terms=c("s(gmst)","ti(sdoy,gmst)"), data=data01[-iob,], sort_by="time", add=FALSE, col="red", xlim=range(data01$time) )
plot_gam_terms(ft1, terms=c("s(gmst)","ti(sdoy,gmst)"), data=data01[iob,],  sort_by="time", add=TRUE,  col="blue")
grid()


plot_gam_terms(ft1, terms=st.terms[-1], data=data01[-iob,], sort_by="time", col="red", xlim=range(data01$time), SE=FALSE)
plot_gam_terms(ft1, terms=st.terms[-1], data=data01[iob,],  sort_by="time", col="blue", SE=FALSE, add=TRUE)
grid()

plot_gam_terms(ft1, terms="s(sdoy,fYear,ftype)", data=data01[-iob,], sort_by="time", col="red", xlim=range(data01$time), SE=FALSE)
plot_gam_terms(ft1, terms="s(sdoy,fYear,ftype)", data=data01[iob,],  sort_by="time", col="blue", SE=FALSE, add=TRUE)
grid()

plot_gam_terms(ft1, terms="s(sdoy,fYear,ftype)", data=data01, sort_by="index", col=1,  SE=FALSE, LINES=F)
grid()

plot_gam_terms(ft1, terms=NULL, data=data01, sort_by="index", col=1,  SE=FALSE)

# plot_gam_terms(ft1, terms=st.terms[7], data=data01[-iob,], sort_by="time", col="red", xlim=range(data01$time), SE=FALSE)
# plot_gam_terms(ft1, terms=st.terms[6], data=data01[iob,],  sort_by="time", col="blue", SE=FALSE, add=TRUE)
# grid()

dev.off()
