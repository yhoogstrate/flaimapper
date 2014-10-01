#!/usr/bin/env R
args <- commandArgs(trailingOnly = TRUE)

data_dir = args[1]
output_dir = args[2]
output_dir_extra = args[3]

table_sensitivity <- read.table(paste(data_dir,"/validation_FlaiMapper-SSLM_miRBase_SRP002175_sensitivity.txt",sep=""),header=T,row.names=1)
table_offset <- read.table(paste(data_dir,"/validation_FlaiMapper-SSLM_miRBase_SRP002175_offset.txt",sep=""),header=T,row.names=1)

# Overlap plot
svg(paste(output_dir_extra,"validation_FlaiMapper_miRBase_SRP002175_sensitivity_barplot.svg",sep="/"))
column <- 3
x_axes = barplot(table_sensitivity[,column],ylim=c(0,950),names.arg=c("predicted","not predicted\nsupporting reads","not predicted\nno supporting reads"),ylab="Annotated miRNAs in miRBase 20")
percentages <- round(table_sensitivity[,column]/sum(table_sensitivity[,column])*100,1)
texting <- paste(table_sensitivity[,column],rep("\n(",3),percentages,rep("%)",3),sep="")
text(x_axes,table_sensitivity[,column]+55,texting)
dev.off()


svg(paste(output_dir_extra,"validation_FlaiMapper_miRBase_SRP002175_sensitivity_piechart.svg",sep="/"))
column <- 3
pie(table_sensitivity[,column],labels=paste(c("predicted","not predicted","no reads"),texting,sep=": "),col=c("darkgreen","red","gray"),main="miRBase predicted by FlaiMapper\nsample SRP002175")
dev.off()


# Error plot 5p
print(paste(output_dir,"0_top_left_validation_FlaiMapper_miRBase_SRP002175_offset_start-positions.svg",sep="/"))
svg(paste(output_dir,"0_top_left_validation_FlaiMapper_miRBase_SRP002175_offset_start-positions.svg",sep="/"),height=4.5)
column <- 5
percentages <- round(table_offset[,column]/sum(table_offset[,column])*100,1)
x_axes = barplot(percentages,ylim=c(0,110),names.arg=row.names(table_offset),
main="5'-error between miRBase miRNAs and FlaiMapper",ylab="% of occurance",xlab="Difference in 5' position between predicted fragment and annotated miRNA")
texting <- paste(table_offset[,column],rep("\n",21),percentages,rep("%",12),sep="")
text(x_axes,percentages+8,texting)
dev.off()


# Error plot 3p
print(paste(output_dir,"0_top_right_validation_FlaiMapper_miRBase_SRP002175_offset_stop-positions.svg",sep="/"))
svg(paste(output_dir,"0_top_right_validation_FlaiMapper_miRBase_SRP002175_offset_stop-positions.svg",sep="/"),height=4.5)
column <- 6
percentages <- round(table_offset[,column]/sum(table_offset[,column])*100,1)
x_axes = barplot(percentages,ylim=c(0,110),names.arg=row.names(table_offset),
main="3'-error between miRBase miRNAs and FlaiMapper",ylab="% of occurance",xlab="Difference in 3' position between predicted fragment and annotated miRNA")
texting <- paste(table_offset[,column],rep("\n",21),percentages,rep("%",12),sep="")
text(x_axes,percentages+8,texting)
dev.off()
