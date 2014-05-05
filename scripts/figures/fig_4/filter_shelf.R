svg("../../../output/figures/fig_4/filter_shelf.svg")

dist_m1 = c(0.00003726653,0.0001866447,0.008364835,0.0354626,0.1203860,0.3865920,1.110900,2.856550,6.572853,13.53353,24.93522,100,100,100,100,100,100,100,100,100,24.93522,13.53353,6.572853,2.856550,1.110900,0.3865920,0.1203860,0.0354626,0.008364835,0.0001866447,0.00003726653)

plot(c(6,36),c(0,355),type="n",xlab="Position in sequence",ylab="Start/stop position count")

dist_m1c = dist_m1*(355/100)
dist_m2c = dist_m1*(85/100)
dist_m3c = dist_m1*(100/100)


lines(c(16,16),c(0,355),lwd=2,col="darkgreen",lty=2)
lines(c(20,20),c(0,85),col="purple",lwd=2,lty=2)
lines(c(28,28),c(0,100),col="blue",lwd=2,lty=2)


polygon(seq(1,length(dist_m1c)),dist_m1c,col="#00FF0033",lty=0)
lines(dist_m1c,lwd=2,col="darkgreen")

polygon(seq(1,length(dist_m1c))+4,dist_m2c,col="#DD00DD33",lty=0)
lines(seq(1,length(dist_m1c))+4,dist_m2c,lwd=2,col="purple")

polygon(seq(1,length(dist_m1c))+12,dist_m3c,col="#0000FF33",lty=0)
lines(seq(1,length(dist_m1c))+12,dist_m3c,lwd=2,col="blue")

dev.off()
