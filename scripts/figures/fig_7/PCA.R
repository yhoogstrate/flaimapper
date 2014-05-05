## Used plot:

column_names = c("pre-miR","SNORD","SNORA","tRNA","SCARNA","MISC")

S2175   = c(947, 202, 56, 1122, 26, 107)
S6788_a = c(680, 140, 24, 957, 11, 755)
S6788_b = c(686, 181, 34, 1043, 24, 805)
S6788_c = c(560, 30, 15, 208, 9, 590)
S6788_d = c(517, 42, 17, 204, 6, 550)
S6788_e = c(738, 294, 112, 753, 47, 902)
S6788_f = c(649, 281, 87, 520, 32, 806)

# For proportional plot:
#S2175   = S2175   / sum(S2175)   * 100
#S6788_a = S6788_a / sum(S6788_a) * 100
#S6788_b = S6788_b / sum(S6788_b) * 100
#S6788_c = S6788_c / sum(S6788_c) * 100
#S6788_d = S6788_d / sum(S6788_d) * 100
#S6788_e = S6788_e / sum(S6788_e) * 100
#S6788_f = S6788_f / sum(S6788_f) * 100

S2175   = log(S2175)
S6788_a = log(S6788_a)
S6788_b = log(S6788_b)
S6788_c = log(S6788_c)
S6788_d = log(S6788_d)
S6788_e = log(S6788_e)
S6788_f = log(S6788_f)

table = t(data.frame(S2175,S6788_a,S6788_b,S6788_c,S6788_d,S6788_e,S6788_f))
colnames(table) = column_names
rownames(table) = c("","Untreated","RRP40","AGO","AGO & RRP40","XRN","Nucleus")

pdf("/dev/null")
pca = prcomp(table)
biplot(pca,xlim=c(-0.6,0.6),cex=c(1, 1/2))
dev.off()

svg("../../../output/figures/fig_7/PCA.svg")

x = pca$x[,1]
y = pca$x[,2]

plot(c(min(x)-1,max(x)+1),c(min(y),max(y)+0.05),type="n",xlab="PC 1",ylab="PC 2")
points(x,y,cex=0.85,pch=19,col=c("red","blue","blue","blue","blue","blue","blue"))
text(x,y,rownames(table),pos=2,col="darkgray",cex=2)

legend(0.5,-0.75,c("SRP002175","SRP006788"),col=c("red","blue"),pch=19,cex=1)


dev.off()
