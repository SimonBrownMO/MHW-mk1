# plots looking at time differences that show SSTs are just AR processes?
plot(data01$time[-1], diff(data01$x0), pch=20, cex=.3, main='running median')
grid()
readline("Stop3 c")


plot(data01$sdoy[-1], diff(data01$x0), pch=20, cex=.3, main='running median')
grid()
readline("Stop3 c")


rm.x0 <- runmed(data01$x0, k=91)
rm.x1 <- runmed(data01$x1, k=91)

plot(data01$time, rm.x0, pch=20, cex=.3, main='running median')
grid()
readline("Stop3 c")

plot(data01$time, rm.x1, pch=20, cex=.3, main='running median')
grid()
readline("Stop3 c")


plot(data01$time, runmed(data01$x0, k=141), pch=20, cex=.3, main='running median'); grid()

rn.x0 <- zapsmall(convolve(data01$x0, rep(1/91,91), type="o"))


n1 <- 120
i1 <- seq((n1-1)/2+1, length=nrow(data01), by=1)
plot(data01$time, zapsmall(convolve(data01$x0, rep(1/n1,n1), type="o"))[i1], pch=20, cex=.3, main='running mean'); grid()