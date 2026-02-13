# source("theo_annual_random_effect.R")

# setwd("/home/users/simon.brown/extremes/heatwaves/mhw/science")

library(mgcv)
library(data.table)

# pdf("theo_annual_random_effect.pdf", width=14, height=10)

# load the data
load("data01_NWS.RData")
data01 <- data.table(data01)
data01

### Fit a Gaussian gam where year is a random effect
# introduce years
data01[,year := trunc(time,0)]
# factor year
data01[,fYear := factor(year)]

model <- gam(x ~ s(doy,bs="cc") # seasonal cycle
               + s(doy,fYear,bs="sz",id=1,xt=list(bs="cc")),
               knots=list(doy=c(0,366)),data=data01
)
# which GAM coefficients relate to doy?
DoYindex <- grep("doy",names(coef(model)))
# model matrix
X <- predict(model,type="lpmatrix")  # [1:366, 1:414] [knots, coefs]
# GAM coefficients
b <- coef(model)                     # [414]    
# estimated seasonal cycles
SC <- X[,DoYindex] %*%  b[DoYindex]  # [1:16723, 1]
# add the intercept to put it on the scale of the data
SC <- SC + b[1]
# Check the seasonal cycle fit
plot(x~time,data=data01,pch=20, main="s(doy,bs=cc)+s(doy,fYear,bs=sz)")
lines(data01$time,SC,col="red",lwd=2)
readline("Continue?2")
plot(data01$time, data01$x-SC,pch=20, cex=.3, main="x - s(doy,bs=cc)+s(doy,fYear,bs=sz)")
grid()
readline("Continue?2b")

## Now see about integrating out the year-specific seasonal cycle
## by simulating from the posterior of the coefficients.
# simulate from posterior
n.sims <- 1000
b_sims <- rmvn(n.sims, coef(model), model$Vp )
# coefficients index relating to just the doy-year interaction
DoYindex <- grep("doy",names(coef(model)))
# grid of years-doy 
newd <- data.table( expand.grid(1:366,levels(data01$fYear)) )
names(newd) <- c("doy","fYear")
X           <- predict(model,newd,type="lpmatrix")
SC_sims     <- tcrossprod( X[,DoYindex],  b_sims[,DoYindex] )
# now put these simulations into a table with DoY samples
# irrespective of year
SC_sims <- data.table(SC_sims)
# include the "observed" DoYs from the data
SC_sims[,doy := newd$doy]
# list of 366 matrices
SC_list <- split(SC_sims[, -"doy"], SC_sims[["doy"]], drop = TRUE)
# list of 366 vactors
SC_list <- lapply(SC_list, function(x){unlist(x)})
# put in a table
SC_table <- do.call(rbind,SC_list)

# Now let's visualise the seasonal cycle with 
# the integrated out doy
plot(1:366,apply(SC_table,1,mean),type="l",ylim=c(-6,6), main="SC_" )
lines(1:366,apply(SC_table,1,quantile,probs=0.025),lty=2)
lines(1:366,apply(SC_table,1,quantile,probs=0.975),lty=2)
# readline("Continue?3")
# compare with the "main effect" of doy
# plot(model,select=1,ylim=c(-6,6))
plot(model)
readline("Continue?4")


# Now predict 2026
# The trick is to predict all terms from the gam
# and then exclude the doy ones. The year needs to
# be a year already in the data so predict() does
# not complain about a factor level that is not
# in the data
newD <- data.table(doy=1:366,fYear="2025")
X    <- predict(model,newD,type="lpmatrix")
pred <- tcrossprod( X[,-DoYindex] ,b_sims[,-DoYindex])
# Now add n.sims of the doy simulations
pred <- pred + SC_table[,sample(1:ncol(SC_table),n.sims)]

pred0 <- tcrossprod( X[,-DoYindex] ,b_sims[,-DoYindex])
pred1 <- SC_table[,sample(1:ncol(SC_table),n.sims)]
# Now add n.sims of the sdoy simulations
pred <- pred0 + pred1
plot(x~time,data=data01,pch=20,xlim=c(1980,2026))
lines(data01$time,SC,col="red")
lines((1:366 )/366 + 2026,apply(pred,1,mean),col="blue",lwd=2)
lines((1:366 )/366 + 2026,apply(pred,1,quantile,probs=0.025),col="grey50",lty=2)
lines((1:366 )/366 + 2026,apply(pred,1,quantile,probs=0.975),col="grey50",lty=2)
readline("Continue?5")

### simon
# i1 <- which(data01$doy==230)
# hist(SC_table[230,]+b[1],20,freq=FALSE)
# hist(SC[i1],20,col=rgb(1,0,0,0.5),freq=FALSE,add=TRUE)
# # mean(SC_table[230,]+b[1]) # 16.17457
# # mean(SC[i1])              # 16.17444
# readline("Continue?6")

# look wrt doy
# fitted seasonal cycle
plot(x~doy,data=data01,pch=20, cex=.3)
for(y in unique(data01$year)) {
    iy <- which(data01$year == y)
    lines(data01$doy[iy],SC[iy],col=rgb(0,0,1,0.4))
}
lines(1:366,apply(SC_table,1,mean)+b[1],                      lwd=3,col="red2")
lines(1:366,apply(SC_table,1,quantile,probs=0.025)+b[1],lty=2,lwd=3,col="red2")
lines(1:366,apply(SC_table,1,quantile,probs=0.975)+b[1],lty=2,lwd=3,col="red2")
readline("Continue?7")

# predicted spread
plot(1:366,apply(pred,1,mean),type="n",ylim=range(pred),xlab="Day of Year",ylab="Seasonal Cycle")
for(i in seq_along(pred[1,])) lines(1:366, pred[,i],col=rgb(0,0,1,0.1))
lines(1:366,apply(pred,1,mean),                lwd=3,col="red2")
lines(1:366,apply(pred,1,quantile,probs=0.025),lwd=3,col="red2",lty=2)
lines(1:366,apply(pred,1,quantile,probs=0.975),lwd=3,col="red2",lty=2)
readline("Continue?8")


# discontinuities
# iy2 <- which(data01$year %in% 2000:2001)
# y   <- seq_along(iy2)
# plot(y, data01[iy2,x],pch=46, xlim=c(350,400), ylim=c(7,11))
# for(y1 in 1980:2000){
#     iy2 <- which(data01$year %in% y1:(y1+1))
#     y   <- seq_along(iy2)
#     lines(y,SC[iy2],col="red")
# }
# readline("Continue?9")

### focus on 1997-1998
iy2 <- which(data01$year %in% 1997:1999)
y   <- seq_along(iy2)
plot(y, data01[iy2,x],pch=20, cex=.5)
lines(y,SC[iy2],col="red")
grid()

# dev.off()