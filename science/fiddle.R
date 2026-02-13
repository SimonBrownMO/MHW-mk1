# source("fiddle.R")

# setwd("/home/users/simon.brown/extremes/heatwaves/mhw/science")

library(mgcv)
library(data.table)


# load the data
load("data01_NWS.RData")
data01 <- data.table(data01)
data01



### Fit a Gaussian gam where year is a random effect
# introduce years
data01[,year := trunc(time,0)]
# factor year
data01[,fYear := factor(year)]

model1 <- gam(x ~ s(sdoy,bs="cc") 
               + s(sdoy,fYear,bs="sz",id=1,xt=list(bs="cc")),
               knots=list(sdoy=c(0,1)),data=data01
)
model2 <- gam(x ~ s(sdoy,bs="cc") 
               + s(sdoy,fYear,bs="sz",id=1),
               knots=list(sdoy=c(0,1)),data=data01
)
model3 <- gam(x ~ s(sdoy,bs="cc", k=12) 
               + s(sdoy,fYear,bs="sz",id=1),
               knots=list(sdoy=c(0,1)),data=data01
)
model4 <- gam(x ~ s(sdoy,bs="cc", k=12) 
                + s(sdoy,fYear,bs="sz",id=1),
                knots=list(sdoy=seq(0,1,length=12)),data=data01
)
model5 <- gam(x ~ s(sdoy,bs="cc", k=4) 
                + s(sdoy,fYear,bs="sz",id=1),
                knots=list(sdoy=seq(0,1,length=4)),data=data01
)
model6 <- gam(x ~ s(sdoy,bs="cc", k=4) 
                + s(sdoy,fYear,bs="sz",id=1),
                ,data=data01
)
model <- model1
# which GAM coefficients relate to doy?
DoYindex <- grep("sdoy",names(coef(model)))
# model matrix
X <- predict(model,type="lpmatrix")  # [1:366, 1:414] [knots, coefs]
# GAM coefficients
b <- coef(model)                     # [414]    
# estimated seasonal cycles
SC <- X[,DoYindex] %*%  b[DoYindex]  # [1:16723, 1]
# add the intercept to put it on the scale of the data
SC <- SC + b[1]

iy2 <- which(data01$year %in% 1989:1992)
y   <- seq_along(iy2)

# plot(x~time,data=data01,pch=20)
# lines(y,SC[iy2],col=2,lwd=4, lty=2)
# readline("Continue?1")

# Check the seasonal cycle fit
plot(x~time,data=data01,pch=20,xlim=c(1980,2026))
lines(data01$time,SC,col="red",lwd=2)
readline("Continue?2")

plot(y, data01[iy2,x],pch=20, cex=.5)
lines(y,SC[iy2],col="red")
grid()

readline("Continue?3")
### end fit #########################################

## Now see about integrating out the year-specific seasonal cycle
## by simulating from the posterior of the coefficients.
# simulate from posterior
n.sims <- 1000
b_sims <- rmvn(n.sims, coef(model), model$Vp )
# coefficients index relating to just the sdoy-year interaction
DoYindex <- grep("sdoy",names(coef(model)))
# grid of years-sdoy 
newd <- data.table( expand.grid(1:366,levels(data01$fYear)) )
names(newd) <- c("sdoy","fYear")
X           <- predict(model,newd,type="lpmatrix")
SC_sims     <- tcrossprod( X[,DoYindex],  b_sims[,DoYindex] )
# now put these simulations into a table with DoY samples
# irrespective of year
SC_sims <- data.table(SC_sims)
# include the "observed" DoYs from the data
SC_sims[,sdoy := newd$sdoy]
# list of 366 matrices
SC_list <- split(SC_sims[, -"sdoy"], SC_sims[["sdoy"]], drop = TRUE)
# list of 366 vactors
SC_list <- lapply(SC_list, function(x){unlist(x)})
# put in a table
SC_table <- do.call(rbind,SC_list)

# Now let's visualise the seasonal cycle with 
# the integrated out sdoy
 plot(1:366,apply(SC_table,1,mean),type="l",ylim=c(-6,6) )
lines(1:366,apply(SC_table,1,quantile,probs=0.025),lty=2)
lines(1:366,apply(SC_table,1,quantile,probs=0.975),lty=2)
readline("Continue?")
# compare with the "main effect" of sdoy
plot(model,select=1,ylim=c(-6,6))
readline("Continue?")


# Now predict 2026
plot(x~time,data=data01,pch=20,xlim=c(1980,2026))
lines(data01$time,SC,col="red")
# The trick is to predict all terms from the gam
# and then exclude the sdoy ones. The year needs to
# be a year already in the data so predict() does
# not complain about a factor level that is not
# in the data
newD <- data.table(sdoy=1:366,fYear="2025")
X <- predict(model,newD,type="lpmatrix")
pred <- tcrossprod( X[,-DoYindex] ,b_sims[,-DoYindex])
# Now add n.sims of the sdoy simulations
pred <- pred + SC_table[,sample(1:ncol(SC_table),n.sims)]
# add to the plot
xval <- (1:366 )/366 + 2026
lines(xval,apply(pred,1,mean),col="blue")
lines(xval,apply(pred,1,quantile,probs=0.025),col="blue")
lines(xval,apply(pred,1,quantile,probs=0.975),col="blue")

### simon
i1 <- which(data01$doy==230)
hist(SC_table[230,]+b[1],20,freq=FALSE)
hist(SC[i1],20,col=rgb(1,0,0,0.5),freq=FALSE,add=T)
# mean(SC_table[230,]+b[1]) # 16.17457
# mean(SC[i1])              # 16.17444

# look wrt doy
# fitted seasonal cycle
c1 <- rainbow(46)
plot(x~sdoy,data=data01,pch=46)
for(y in unique(data01$year)) {
    iy <- which(data01$year == y)
    # lines(data01$sdoy[iy],SC[iy],col=rgb(0,0,1,0.4))
    lines(data01$sdoy[iy],SC[iy],col=c1[y-1979])
}
legend("topleft",legend=1980:2025,col=c1,lty=1,cex=0.6)


# predicted spread
plot(1:366,apply(pred,1,mean),type="n",ylim=range(pred),xlab="Day of Year",ylab="Seasonal Cycle")
for(i in seq_along(pred[1,])) lines(1:366, pred[,i],col=rgb(0,0,1,0.1))
lines(1:366,apply(pred,1,mean),lwd=2,col="red4")



iy2 <- which(data01$year %in% 2000:2001)
y   <- seq_along(iy2)
plot(y, data01[iy2,x],pch=46, xlim=c(350,400), ylim=c(7,11))
for(y1 in 1980:2000){
    iy2 <- which(data01$year %in% y1:(y1+1))
    y   <- seq_along(iy2)
    lines(y,SC[iy2],col="red")
}

### focus on 1997-1998
iy2 <- which(data01$year %in% 1997:1998)
y   <- seq_along(iy2)
plot(y, data01[iy2,x],pch=20, cex=.5)
lines(y,SC[iy2],col="red")
grid()



c1 <- rainbow(7)
plot(x~sdoy,data=data01,pch=46)
for(y in 1982:1986) {
    iy <- which(data01$year == y)
    # lines(data01$sdoy[iy],SC[iy],col=rgb(0,0,1,0.4))
    lines(data01$sdoy[iy],SC[iy],col=c1[y-1981],lwd=4)
}
legend("topleft",legend=1982:1986,col=c1,lty=1,cex=1,lwd=4)