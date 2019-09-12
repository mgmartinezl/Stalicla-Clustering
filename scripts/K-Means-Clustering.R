# Clustering ASD patients using K-means over logistic PCA

# Author: Gabriela Martinez
# airamgabriela17@gmail.com

# The input dataset for this exercise is the binary matrix generated from the PBPM protocol. 


# ------ Define user parameters ------ #

input_path = "~//data//raw//binary-matrix-2019-08-23_PTV.csv"
output_folder = "~//reports//KMeans//"

patients_clust = "PTV-kMeans-patients.csv"
km_clust = "PTV-kMeans-patients.png"
elbow_clust = "elbow-optimal-clusters.png"
silhouette_clust = "silhouette-optimal-clusters.png"
silhouette_val = "silhouette-validation.png"

num_principal_components = 4
num_clusters = 4

library(glue)
patients_clust = glue(output_folder, patients_clust)
km_clust = glue(output_folder, km_clust)
elbow_clust = glue(output_folder, elbow_clust)
silhouette_clust = glue(output_folder, silhouette_clust)
silhouette_val = glue(output_folder, silhouette_val)


#----- Data reading and processing ------#

patients = read.csv(input_path, header=TRUE, sep=",", skip = 1, row.names = 1)

# Drop features with only 0s or 1s
patients <- patients[, colSums(patients != 0) > 0]


#----- Logistic PCA ------#

# Principal component analysis (PCA) for binary data, known as logistic PCA, 
# has become a popular alternative to dimensionality reduction of binary data.
# https://arxiv.org/pdf/1510.06112v1.pdf

library(logisticPCA)

# Exponential Family PCA
logsvd_model = logisticSVD(patients, k = num_principal_components)

# From this model, it is possible to extract the score matrix A and loading matrix B
# We will use the scores to apply partition-based algorithms
PCA <- cbind(patients, logsvd_model$A)
PCA <- select(PCA, '1' ,'2', '3', '4')

# Info about the model: 3 components explain 65% of the variance, 4 (76%) and 6 of them explain 96%.
# logsvd_model
# summary(logsvd_model)

# Two-dimensional components plot
# plot(logsvd_model, type = "scores") # 50% of variance explained

# Three-dimensional components plot
# library(plotly)
# logsvdpc1 <- logsvd_model$A[,1]
# logsvdpc2 <- logsvd_model$A[,2]
# logsvdpc3 <- logsvd_model$A[,3]
# logsvdPCs <- cbind(patients, logsvdpc1, logsvdpc2, logsvdpc3)
# logsvd.plot <- plot_ly(logsvdPCs, x=~logsvdpc1, y=~logsvdpc2, z=~logsvdpc3)
# logsvd.plot


#----- Optimal number of clusters and algorithms -----#

library(clValid)
intern <- clValid(as.matrix(PCA), 
                  nClust = 2:9, 
                  metric = "correlation",
                  clMethods = c("hierarchical","kmeans","pam","clara","agnes","fanny"), 
                  validation = "internal")

# Hierarchical clustering with k=2 is suggested as the best option 
# Do not trust this blindly!
summary(intern)


#----- K-Means -----#

# Silhouette method to compute optimal clusters >> 9 clusters
# Do not trust this blindly!
library(factoextra)
sil_clust <- fviz_nbclust(PCA, kmeans, method = "silhouette") + labs(subtitle = "Silhouette method")

# Save plot
library(cowplot)
save_plot(silhouette_clust, sil_clust)

# Elbow method to compute optimal clusters >> 2 clusters
# Do not trust this blindly!
elbow <- fviz_nbclust(PCA, kmeans, method = "wss") + geom_vline(xintercept = 2, linetype = 2) + labs(subtitle = "Elbow method")

# Save plot
library(cowplot)
save_plot(elbow_clust, elbow)

# K-means 
km <- kmeans(PCA, centers = num_clusters, nstart = 30)
km_plot <- fviz_cluster(km, data = PCA, ellipse.type = "convex")

# Save plot
library(cowplot)
save_plot(km_clust, km_plot)

# Internal validation for K-means with Silhouete method >> good classification according to Silhouette
# Do not trust this blindly!
silkm <- silhouette(km$cluster, dist(PCA))
summary(silkm)
sil_plot <- fviz_silhouette(silkm)

# Save plot
library(cowplot)
save_plot(silhouette_val, sil_plot)

# Export clusters
patients.clusters <- cbind(patients, cluster = km$cluster)
write.csv(patients.clusters, file = patients_clust, row.names = TRUE)
