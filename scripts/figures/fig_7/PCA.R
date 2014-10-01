#!/usr/bin/env R
args <- commandArgs(trailingOnly = TRUE)

output_dir = args[1]

precursors_names = c("pre-miR","SNORD","SNORA","tRNA","SCARNA","MISC")
precursors_total = c(1386,264,106,451,23,39)
datasets_names   = c("SRP002175* (Pigment)", "SRP006788 (HeLa)", "SRP028959 (HeLa)", "SRP034013* (B cells)", "SRP041082* (Prostate)")

samples_SRP002175_merged = c("SRP002175* (Pigment)")
samples_SRP006788 = c("SRR207111 HeLa (Total RNA)","SRR207112 HeLa (Total RNA)","SRR207113 HeLa (AGO pd)","SRR207114 HeLa (AGO pd)","SRR207115 HeLa (XRN)","SRR207116 HeLa (Nuclear)")
samples_SRP028959 = c("SRR954957 HeLa (Total RNA)", "SRR954958 HeLa (Nuclear)", "SRR954959 HeLa (Cytoplasmic)")
samples_SRP034013_merged = c("SRP034013* (B cells)")
samples_SRP041082_merged = c("SRP041082* (Prostate)")

datasets_length = c( length(samples_SRP002175_merged),
                     length(samples_SRP006788),
                     length(samples_SRP028959),
                     length(samples_SRP034013_merged),
                     length(samples_SRP041082_merged))

samples_names = c(
                  samples_SRP002175_merged,
                  samples_SRP006788,
                  samples_SRP028959,
                  samples_SRP034013_merged,
                  samples_SRP041082_merged)
samples_color = c(
                  rep("red",datasets_length[1]),
                  rep("blue",datasets_length[2]),
                  rep("green",datasets_length[3]),
                  rep("gray",datasets_length[4]),
                  rep("purple",datasets_length[5]))

samples_cex    = c(
                  rep(1.4,datasets_length[1]),
                  rep(1.4,datasets_length[2]),
                  rep(1.4,datasets_length[3]),
                  rep(1.4,datasets_length[4]),
                  rep(1.4,datasets_length[5]))


samples_pch    = c(
                  rep(16,datasets_length[1]),
                  rep(15,datasets_length[2]),
                  rep(18,datasets_length[3]),
                  rep(17,datasets_length[4]),
                  rep(6,datasets_length[5]))


fragments = data.frame(
                       c( 947,  202,   56, 1210,   26,  107), #SRP002175
                       
                       #c( 500,  104,   23,  937,   11,  567), #SRP002175_SRR038852
                       #c( 346,   64,    7,  554,    1,  360), #SRP002175_SRR038853
                       #c( 519,   93,   14,  810,    7,  547), #SRP002175_SRR038854
                       #c( 374,   69,    9,  687,    6,  394), #SRP002175_SRR038855
                       #c( 416,   60,    8,  767,    8,  436), #SRP002175_SRR038856
                       #c( 546,  119,   33, 1094,   14,  586), #SRP002175_SRR038857
                       #c( 404,   61,    7,  700,    6,  427), #SRP002175_SRR038858
                       #c( 572,  115,   28,  939,   14,  607), #SRP002175_SRR038859
                       #c( 452,  138,   34,  927,   13,  488), #SRP002175_SRR038860
                       #c( 515,  131,   26, 1033,   13,  574), #SRP002175_SRR038861
                       #c( 446,   80,   11,  750,   10,  466), #SRP002175_SRR038862
                       
                       c( 680,  140,   24,  998,   11,  755), #SRP006788_SRR207111
                       c( 686,  181,   34, 1104,   24,  806), #SRP006788_SRR207112
                       c( 560,   30,   15,  208,    9,  590), #SRP006788_SRR207113
                       c( 517,   42,   17,  209,    6,  550), #SRP006788_SRR207114
                       c( 738,  294,  112,  760,   47,  903), #SRP006788_SRR207115
                       c( 649,  280,   87,  519,   32,  806), #SRP006788_SRR207116
                       
                       c( 203,   15,    0,  229,    2,  210), #SRP028959_SRR954957
                       c( 169,   33,    0,  126,    4,  176), #SRP028959_SRR954958
                       c( 197,    2,    0,  274,    1,  204), #SRP028959_SRR954959
                       
                       c(2230,  670,  384, 1478,  105,  391), #SRP034013
                       
                       #c(1644,  521,  202, 1317,   59, 1956), #SRP034013_SRR1049397
                       #c(1629,  614,  303, 1389,   92, 1969), #SRP034013_SRR1049398
                       #c(1595,  549,  259, 1360,   72, 1939), #SRP034013_SRR1049399
                       
                       c(1454,  552,   201, 1222,   99, 304 ) #SRP041082
                       
                       #c(1205,  438,   92, 1175,   68, 1452), #SRP041082_SRR1232072
                       #c(1249,  507,  184, 1166,   82, 1525)  #SRP041082_SRR1232073
                       
                      )

fragments = t(fragments)
colnames(fragments) = precursors_names
rownames(fragments) = samples_names

fragments_samples_sum = apply(fragments, 1, sum)

fragments_normalized_1 = fragments / fragments_samples_sum

pca = prcomp(fragments_normalized_1)

x = pca$x[,1]
y = pca$x[,2]

print(paste(output_dir,"/PCA_merged.svg",sep=""))
svg(paste(output_dir,"/PCA_merged.svg",sep=""))

dy = abs(max(y) - min(y))
dx = abs(max(x) - min(x))
Sys.sleep(0.1);
plot(c(min(x)*1.35,max(x)*1.85),c(min(y)*1.05,max(y)*1.05),type="n",xlab="PC 1",ylab="PC 2")
legend(min(x)-(0.20*dx),max(y) - (0.65*dy),c(datasets_names,"*merged"),col=c(unique(samples_color),'black'),pch=c(unique(samples_pch),-1),cex=c(samples_cex,samples_cex[1])*0.8)
points(x,y,col=samples_color,cex=samples_cex*0.8,pch=samples_pch)
text(x+(0.05*dx),y+(0.025*dy),samples_names,col=samples_color,cex=samples_cex*0.8)

dev.off()

#Variance(s)
n = (pca$sdev)^2 / sum(pca$sdev^2)
m = cumsum((pca$sdev)^2) / sum(pca$sdev^2)
paste(round(100*m, 2), "%", sep="")
