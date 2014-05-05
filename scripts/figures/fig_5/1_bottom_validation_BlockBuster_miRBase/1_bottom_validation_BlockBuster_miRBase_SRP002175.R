#for(_sd in seq(0.05,0.1,0.2,0.35,0.5,1,1.5,2)) {
for(x_sd in c(0.05)) {
  #  for(_dist in seq(1,50)) {
  for(x_dist in c(26)) {
    
    param = paste('min-dist_',x_dist,"_sd_",x_sd,sep='')
    
    table_sensitivity <- read.table("1_bottom_validation_BlockBuster_min-dist_26_sd_0.05_miRBase_SRP00217501_sensitivity.txt",header=T,row.names=1)
    table_offset <- read.table("1_bottom_validation_BlockBuster_min-dist_26_sd_0.05_miRBase_SRP00217501_offset.txt",header=T,row.names=1)
    
    # Overlap plot
    svg(paste("../../../../output/figures/fig_5_extra/validation_BlockBuster_",param,"_miRBase_SRP002175_sensitivity_barplot.svg",sep=""))
    column <- 3
    x_axes = barplot(table_sensitivity[,column],ylim=c(0,650),names.arg=c("predicted","not predicted\nsupporting reads","not predicted\nno supporting reads"),ylab="Annotated miRNAs in miRBase 20")
    percentages <- round(table_sensitivity[,column]/sum(table_sensitivity[,column])*100,1)
    texting <- paste(table_sensitivity[,column],rep("\n(",3),percentages,rep("%)",3),sep="")
    text(x_axes,table_sensitivity[,column]+36,texting)
    dev.off()
    
    
    svg(paste("../../../../output/figures/fig_5_extra/validation_BlockBuster_",param,"_miRBase_SRP002175_sensitivity_piechart.svg",sep=""))
    column <- 3
    pie(table_sensitivity[,column],labels=paste(c("predicted","not predicted","no reads"),texting,sep=": "),col=c("darkgreen","red","gray"),main="miRBase predicted by BlockBuster\nsample SRP006788")
    dev.off()
    
    
    # Error plot 5p
    svg(paste("../../../../output/figures/fig_5/1_bottom_left_validation_BlockBuster_",param,"_miRBase_SRP00217501_offset_start_positions.svg",sep=""),height=4.5)
    column <- 5
    percentages <- round(table_offset[,column]/sum(table_offset[,column])*100,1)
    x_axes = barplot(percentages,ylim=c(0,110),names.arg=row.names(table_offset),
    main="5'-error between miRBase miRNAs and BlockBuster",ylab="% of occurance",xlab="Difference in 5' position between predicted fragment and annotated miRNA")
    texting <- paste(table_offset[,column],rep("\n",21),percentages,rep("%",12),sep="")
    text(x_axes,percentages+8,texting)
    dev.off()
    
    
    
    svg(paste("../../../../output/figures/fig_5/1_bottom_right_validation_BlockBuster_",param,"_miRBase_SRP00217501_offset_stop_positions.svg",sep=""),height=4.5)
    column <- 6
    percentages <- round(table_offset[,column]/sum(table_offset[,column])*100,1)
    x_axes = barplot(percentages,ylim=c(0,110),names.arg=row.names(table_offset)
    ,main="3'-error between miRBase miRNAs and BlockBuster",ylab="% of occurance",xlab="Difference in 3' position between predicted fragment and annotated miRNA")
    texting <- paste(table_offset[,column],rep("\n",21),percentages,rep("%",12),sep="")
    text(x_axes,percentages+8,texting)
    dev.off()
  }
}

