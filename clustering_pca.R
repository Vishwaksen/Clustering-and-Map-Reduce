# R Script for K-means, Heirarchical, DBSCAN PCA Visualization

# Set the current working directory here

setwd("C:/Users/Vishwaksen/workspace/Clustering")

dat = read.table("PCA_Dbscan_iyer.csv", header = TRUE)

pc <- princomp(dat)
plot(pc,type="lines")
biplot(pc)

prin_comp <- prcomp(dat, scale. = T)
prin_comp <- prcomp(dat)
biplot(prin_comp)

plot3d(prin_comp)