# Clustering ASD patients using Latent Class Analysis (LCA)

# Author: Gabriela Martinez
# airamgabriela17@gmail.com

# The input dataset for this exercise is the binary matrix generated from the PBPM protocol. 


# ------ Define user parameters ------ #

input_path = "~//data//raw//binary-matrix-2019-08-23_PTV.csv"
output_folder = "~//reports//LCA//"
patients_clust = "PTV-LCA.csv"

# LCA requires a function specifying the column names of the dataset
# Change this manually to include or delete more features
f <- cbind(pathway_1,pathway_10,pathway_11,pathway_12,pathway_13,pathway_14,pathway_15,pathway_17,pathway_18,pathway_19,pathway_2,pathway_21,pathway_22,pathway_23,pathway_24,pathway_25,pathway_26,pathway_27,pathway_28,pathway_29,pathway_3,pathway_30,pathway_31,pathway_32,pathway_33,pathway_34,pathway_35,pathway_36,pathway_37,pathway_38,pathway_39,pathway_4,pathway_40,pathway_41,pathway_42,pathway_43,pathway_44,pathway_45,pathway_46,pathway_47,pathway_48,pathway_49,pathway_5,pathway_50,pathway_51,pathway_52,pathway_53,pathway_54,pathway_55,pathway_56,pathway_57,pathway_58,pathway_59,pathway_6,pathway_60,pathway_61,pathway_62,pathway_63,pathway_64,pathway_65,pathway_66,pathway_67,pathway_68,pathway_69,pathway_7,pathway_70,pathway_71,pathway_72,pathway_73,pathway_74,pathway_75,pathway_76,pathway_77,pathway_78,pathway_79,pathway_8,pathway_80,pathway_9)~ 1

library(glue)
patients_clust = glue(output_folder, patients_clust)


#----- Data reading and processing ------#

patients = read.csv(input_path, header=TRUE, sep=",", skip = 1, row.names = 1)

# Drop features with only 0s or 1s
patients <- patients[, colSums(patients != 0) > 0]

# For LCA to work, the dataset cannot have 0s
# We will replace all 0s by 2s
patients[patients == 0] <- 2 


#----- Latent Class Analysis ------#

# Latent class analysis is a technique used to classify observations 
# based on patterns of categorical responses.

library(poLCA)

# LCA models with 2, 3, 4, 5 and 6 clusters
M2 <- poLCA(f, patients, nclass=2, nrep=10, graphs=TRUE) #BIC(2): 7487
M3 <- poLCA(f, patients, nclass=3, nrep=10, graphs=TRUE) #BIC(3): 7285
#M4 <- poLCA(f, patients, nclass=4, nrep=10, graphs=TRUE) #Not possible
#M5 <- poLCA(f, patients, nclass=5, nrep=10, graphs=TRUE) #Not possible
#M6 <- poLCA(f, patients, nclass=6, nrep=10, graphs=TRUE) #Not possible


#----- Latent Class Analysis Validation ------#

# These clustering technique is unsupervised and model-based
# Rather than looking at simmilarities amongst data points,
# it looks for common underlying probabilistic distributions
# to form clusters.

# Then, the number of clusters is evaluated with estimators
# In this case, we wil consider the BIC, which, the lower, the better
# This particular case fits best with k=3 classes
# See BIC's above


#----- Export results ------#

# Append clusters to original data
LCA_clusters <- cbind(patients, C3 = M3$predclass, C2 = M2$predclass)
write.csv(LCA_clusters, file = patients_clust, row.names = TRUE)
