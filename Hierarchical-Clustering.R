# Clustering ASD patients and pathways with hierarchical agglomerative methods

# Author: Gabriela Martinez
# airamgabriela17@gmail.com

# The input dataset for this exercise is the binary matrix generated from the PBPM protocol. 

# This script only runs agglomerative hierarchical techniques for clustering.

# ------ Define user parameters ------ #

input_path = "C:\\Users\\gabim\\Documents\\Master\\Stalicla\\Work\\Repos\\Clustering-PTV\\Part_II\\Data\\binary-matrix-2019-08-23_13???08???35.csv"
output_folder = "C:\\Users\\gabim\\Documents\\Master\\Stalicla\\Work\\Repos\\Clustering-PTV\\Part_III\\Results\\Hierarchical\\"

patients_clust = "PTV-agglomhc-patients.csv"
pathways_clust = "PTV-agglomhc-pathways.csv"
patients_dengrogram = "PTV-patients-dendrogram.png"
pathways_dengrogram = "PTV-pathways-dendrogram.png"
patients_dengrogram_info = "PTV-patients-dendrogram-labels-height.csv"
pathways_dengrogram_info = "PTV-pathways-dendrogram-labels-height.csv"

min_clusters_patients = 2
max_clusters_patients = 5

min_clusters_pathways = 2
max_clusters_pathways = 5

library(glue)
patients_clust = glue(output_folder, patients_clust)
pathways_clust = glue(output_folder, pathways_clust)
patients_dengrogram = glue(output_folder, patients_dengrogram)
pathways_dengrogram = glue(output_folder, pathways_dengrogram)
patients_dengrogram_info = glue(output_folder, patients_dengrogram_info)
pathways_dengrogram_info = glue(output_folder, pathways_dengrogram_info)


#----- Data reading and processing ------#

patients = read.csv(input_path, header=TRUE, sep=",", skip = 1, row.names = 1)

# Drop features with only 0s or 1s
patients <- df[, colSums(df != 0) > 0]


#----- Similarity matrix for patients -----#

# Compute similarity matrix for patients, using Jaccard's distance for binary dichotomic variables
patients.simmilarity <- dist(patients, method="binary")


#----- Hierarchical clustering for patients -----#

# Compute clusters
library(FactoClass)
patients.clustering <- ward.cluster(patients.simmilarity, plots=True, h.clust = 1)

# Compute and export dendrogram
library(ggdendro)
patients.clustering.dendro <- ggdendrogram(patients.clustering)

library(cowplot)
save_plot(patients_dengrogram, patients.clustering.dendro)

# Extract clusters >> min_clusters_patients to max_clusters_patients
patients.clusters <- cutree(patients.clustering, k=min_clusters_patients:max_clusters_patients)
patients.clusters <- data.frame(patients.clusters)

# Extract the labels and the height of patients in the dendrogram
patients.dendro.labels <- patients.clustering$label[patients.clustering$order]
patients.dendro.height <- patients.clustering$height[patients.clustering$order]
patients.dendro.summary <- data.frame(cbind(patients.dendro.labels, patients.dendro.height))

# Append clusters to patients IDs
library(data.table)
setnames(patients.clusters,paste0(names(patients.clusters),"_clusters"))

# Count the occurrences for any number of clusters
# Replace the number according to the number of clusters >> X2 for 2 clusters
# library(dplyr)
# count(patients.clusters, patients.clusters$X2)
# count(patients.clusters, patients.clusters$X3)

# Export results to a csv
write.csv(patients.clusters, file = patients_clust, row.names = TRUE)
write.csv(patients.dendro.summary, file = patients_dengrogram_info, row.names = FALSE)


#----- Similarity matrix for pathways -----#

# Compute phi coefficient of correlation for pathways
pathways.correlation <- cor(patients, use = "pairwise.complete.obs")

# Compute similarity matrix for pathways
pathways.simmilariy = 1 - pathways.correlation
pathways.distance = as.dist(pathways.simmilariy)


#----- Hierarchical clustering for pathways -----#

# Compute clusters
library(FactoClass)
pathways.clustering <- ward.cluster(pathways.distance, plots=True, h.clust = 1)

# Compute and export dendrogram
library(ggdendro)
pathways.clustering.dendro <- ggdendrogram(pathways.clustering)

library(cowplot)
save_plot(pathways_dengrogram, pathways.clustering.dendro)

# Extract clusters >> min_clusters_pathways to max_clusters_pathways
pathways.clusters <- cutree(pathways.clustering, k=min_clusters_pathways:max_clusters_pathways)
pathways.clusters <- data.frame(pathways.clusters)

# Extract the labels and the height of pathways in the dendrogram
pathways.dendro.labels <- pathways.clustering$label[pathways.clustering$order]
pathways.dendro.height <- pathways.clustering$height[pathways.clustering$order]
pathways.dendro.summary <- data.frame(cbind(pathways.dendro.labels, pathways.dendro.height))

# Append clusters to pathways names
library(data.table)
setnames(pathways.clusters,paste0(names(pathways.clusters),"_clusters"))

# Export clusters and dendrogram information to a csv
write.csv(pathways.clusters, file = pathways_clust, row.names = TRUE)
write.csv(pathways.dendro.summary, file = pathways_dengrogram_info , row.names = FALSE)
