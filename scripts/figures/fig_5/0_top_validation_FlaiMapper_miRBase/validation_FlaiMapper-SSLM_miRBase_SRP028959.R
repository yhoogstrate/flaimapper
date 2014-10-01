#!/usr/bin/env R
args <- commandArgs(trailingOnly = TRUE)

data_dir = args[1]
output_dir = args[2]
output_dir_extra = args[3]

experiments = c('SRR954957', 'SRR954958', 'SRR954959')

for(experiment in experiments) {
  f1 = paste(data_dir,"/validation_FlaiMapper_miRBase_SRP028959_",experiment,"_sensitivity.txt",sep="")
  f2 = paste(data_dir,"/validation_FlaiMapper_miRBase_SRP028959_",experiment,"_offset.txt",sep="")
  
  f3 = paste(output_dir_extra,"/validation_FlaiMapper_miRBase_SRP028959_",experiment,"_sensitivity_barplot.svg",sep="")
  f4 = paste(output_dir_extra,"/validation_FlaiMapper_miRBase_SRP028959_",experiment,"_sensitivity_piechart.svg",sep="")
  
  f5 = paste(output_dir_extra,"/validation_FlaiMapper_miRBase_SRP028959_",experiment,"_offset_with_miRBase_start-positions.svg",sep="")
  f6 = paste(output_dir_extra,"/validation_FlaiMapper_miRBase_SRP028959_",experiment,"_offset_with_miRBase_stop-positions.svg",sep="")
  
  table_overlap <- read.table(f1,header=T,row.names=1)
  table_error <- read.table(f2,header=T,row.names=1)
  
  # Overlap plot
  svg(f3)
  column <- 3
  x_axes = barplot(table_overlap[,column],ylim=c(0,650),names.arg=c("predicted","not predicted\nsupporting reads","not predicted\nno supporting reads"),ylab="Annotated miRNAs in miRBase 20")
  percentages <- round(table_overlap[,column]/sum(table_overlap[,column])*100,1)
  texting <- paste(table_overlap[,column],rep("\n(",3),percentages,rep("%)",3),sep="")
  text(x_axes,table_overlap[,column]+36,texting)
  dev.off()
  
  
  svg(f4)
  column <- 3
  pie(table_overlap[,column],labels=paste(c("predicted","not predicted","no reads"),texting,sep=": "),col=c("darkgreen","red","gray"),main="miRBase predicted by FlaiMapper\nsample SRP028959_4-1x_Ion_Torrent_PGM_HeLa")
  dev.off()
  
  
  # Error plot 5p
  svg(f5)
  column <- 5
  percentages <- round(table_error[,column]/sum(table_error[,column])*100,1)
  x_axes = barplot(percentages,ylim=c(0,110),names.arg=row.names(table_error),
  main="5'-error between miRBase miRNAs and FlaiMapper",ylab="% of occurance",xlab="Difference in 5' position between predicted fragment and annotated miRNA")
  texting <- paste(table_error[,column],rep("\n",21),percentages,rep("%",12),sep="")
  text(x_axes,percentages+8,texting)
  dev.off()
  
  
  # Error plot 3p
  svg(f6)
  column <- 6
  percentages <- round(table_error[,column]/sum(table_error[,column])*100,1)
  x_axes = barplot(percentages,ylim=c(0,110),names.arg=row.names(table_error)
  ,main="3'-error between miRBase miRNAs and FlaiMapper",ylab="% of occurance",xlab="Difference in 3' position between predicted fragment and annotated miRNA")
  texting <- paste(table_error[,column],rep("\n",21),percentages,rep("%",12),sep="")
  text(x_axes,percentages+8,texting)
  dev.off()
  
}
