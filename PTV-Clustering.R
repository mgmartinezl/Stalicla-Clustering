# Clustering PTV patients

# Author: Gabriela Martinez
# airamgabriela17@gmail.com

# The input dataset for this exercise is the binary matrix generated
# from the PBPM protocol. In this particular run, only PTV patients
# were considered.

# ------ Define parameters ------ #

input_path = "C:\\Users\\gabim\\Documents\\Master\\Stalicla\\Work\\Repos\\Clustering-PTV\\Part_II\\Data\\binary-matrix-2019-08-23_13???08???35.csv"
output_folder = "C:\\Users\\gabim\\Documents\\Master\\Stalicla\\Work\\Repos\\Clustering-PTV\\Part_III\\Results\\"

patients_file = "PTV-agglomhc-patients.csv"
pathways_file = "PTV-agglomhc-pathways.csv"

min_clusters_patients = 2
max_clusters_patients = 10

min_clusters_pathways = 2
max_clusters_pathways = 10

library(glue)
patients_filename = glue(output_folder, patients_file)
pathways_filename = glue(output_folder, pathways_file)

#----- Data reading and processing ------#

df = read.csv(input_path, header=TRUE, sep=",", skip = 1)

# Drop features with only 0s or 1s
df <- df[, colSums(df != 0) > 0]

#----- Similarity matrix for patients -----#

# Compute similarity matrix for patients, using Jaccard's distance for binary dichotomic variables
patients <- df[,2:ncol(df)]
patients.simmilarity <- dist(patients, method="binary")
#patients.simmilarity.matrix <- as.matrix(patients.simmilarity)

#----- Hierarchical clustering for patients -----#

# Compute clusters
library(FactoClass)
patients.clustering <- ward.cluster(patients.simmilarity, plots=True, h.clust = 1)

# Visualize dendrogram
library(ggdendro)
ggdendrogram(patients.clustering)

# Extract clusters >> min_clusters_patients to max_clusters_patients
patients.clusters <- cutree(patients.clustering, k=min_clusters_patients:max_clusters_patients)

# Append clusters to patient IDs
x <- data.frame(patients.clusters)
library(data.table)
setnames(x,paste0(names(x),"_clusters"))
patients.results <- cbind(child_id = df$child_id, x)

# Count the occurrences for any number of clusters
# Replace the number according to the number of clusters >> X2 for 2 clusters
# library(dplyr)
# count(results, results$X2)
# count(results, results$X3)

# Export results to a csv
write.csv(patients.results, file = patients_filename, row.names=FALSE)


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

# Visualize dendrogram
library(ggdendro)
ggdendrogram(pathways.clustering)

# Extract clusters >> min_clusters_pathways to max_clusters_pathways
pathways.clusters <- cutree(pathways.clustering, k=min_clusters_pathways:max_clusters_pathways)

# Append clusters to pathways names
pathways.results <- data.frame(pathways.clusters)
library(data.table)
setnames(pathways.results,paste0(names(pathways.results),"_clusters"))

# Export results to a csv
write.csv(pathways.results, file = pathways_filename, row.names = TRUE)
