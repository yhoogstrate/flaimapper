library(scatterplot3d)

p5 = read.table("validation_BlockBuster_miRBase__optimal_settings__root_square_error_plateau_start_positions.txt",header=T,row.names=1)
p3 = read.table("validation_BlockBuster_miRBase__optimal_settings__root_square_error_plateau_stop_positions.txt",header=T,row.names=1)

colnames(p5) = as.numeric(gsub("X", "", colnames(p5)))
colnames(p3) = as.numeric(gsub("X", "", colnames(p3)))


persp(seq(10, 300, 5), seq(10, 300, 5), z, phi = 45, theta = 45,
      xlab = "X Coordinate (feet)", ylab = "Y Coordinate (feet)",
      main = "Surface elevation data"
)


persp(sort(as.numeric(rownames(p5))), sort(as.numeric(colnames(p5)))  *100, round( p5), phi = 45, theta = 45,
      xlab = "X Coordinate (feet)", ylab = "Y Coordinate (feet)",
      main = "Surface elevation data"
)

x = rep(seq(1,50),8)
y = rep(seq(1,8),50)

scatterplot3d(x,y,as.list(as.matrix(p5)))
scatterplot3d(x,y,as.list(as.matrix(p3)))
