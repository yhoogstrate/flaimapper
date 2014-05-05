t <- read.table("alignment_plot.txt",header=T)
t[,4] = (t[,2]+t[,3])/2
t_sorted <- t[order(-t[,4],-t[,2], t[,3]), ]

y <- 1

x_lim = c(1,93)
y_lim = c(1,sum(t_sorted[,1])+14)

svg("../../../output/figures/fig_1/alignment_plot__SNORD74_SRP006788_SRR207111_HeLa18-30.svg")

plot(x_lim,y_lim,type="n")

for(i in 1:nrow(t_sorted)) {
  element = t_sorted[i,]
  print(element)
  for(j in 1:element[,1]) {
    lines(c(element[,2],element[,3]),c(y,y),lwd=0.2,col="#0000FF")
    y = y + 1
  }
}

lines(c(1,93),c(y+11,y+11),lwd=3,col="#000000")

dev.off()
