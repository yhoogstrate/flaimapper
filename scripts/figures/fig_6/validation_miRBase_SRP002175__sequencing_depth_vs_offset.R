t = read.table("validation_miRBase_SRP002175__sequencing_depth_vs_offset.txt",stringsAsFactors=F,header=T)

t1 = subset(t,X5p_error <= 10)
t2 = subset(t,X3p_error >= -10)

svg("../../../output/figures/fig_6/validation_miRBase_SRP002175__sequencing_depth_vs_offset.svg",width=7,height=3.5)

plot(t1[,2],log(t1[,1]),col="red",cex=0.2,pch=19,xlim=c(-10,30),xaxt = "n",ylab="log(corresponding reads)",xlab="offset with miRBase")
points(t2[,4]+20,log(t2[,3]),col="blue",cex=0.2,pch=19,xlim=c(-10,10))
abline(v=10,col="gray")
axis(1, at=-10:30, labels=c(-10:9,"(-)10",-9:10))
text(-5,14,"5'")
text(15,14,"3'")

dev.off()
