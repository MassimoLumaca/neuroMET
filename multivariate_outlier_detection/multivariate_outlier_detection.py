import scipy.io
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

# This script first normalizes the data using z-score normalization, then applies
# PCA. It calculates the Euclidean distance of each point from the origin in the 
# space of the first two principal components. Outliers are identified as points
# whose distance is more than two standard deviations away from the mean distance.
# The script then plots the scores on the first two principal components and highlights
# the outliers.

# Load the MATLAB .mat file
file_path = 'data/behavioural/table.covariates.mat'
mat = scipy.io.loadmat(file_path)

# Extract the data including: MET percentage scores (melodic, rhythmic, Total) and Gold-msi scores
table_covariate = pd.DataFrame(mat['table_covariate'])

# Selecting specific columns (4 = MelodicPercentageScore; 6 = RhythmicPercentageScore;
# 8 = TotalMETpercentageScore; 9 = F1; 10 = F2; 11 = F3; 12 = F4; 13 = F5; 14 = GS)
table_covariate_reduced = table_covariate.iloc[:, [4, 6, 8, 9, 10, 11, 12, 13, 14]]

# Normalize the data
scaler = StandardScaler()
X = scaler.fit_transform(table_covariate_reduced)

# Perform PCA
pca = PCA(n_components=2)
score = pca.fit_transform(X)

# Calculate the distance of each point from the origin
distances = np.sqrt(np.sum(score**2, axis=1))

# Identify outliers
n_std = 2
outliers = distances > np.mean(distances) + n_std * np.std(distances)

# Plot the scores and highlight the outliers
plt.figure()
plt.scatter(score[:, 0], score[:, 1], label='Normal Points')
plt.scatter(score[outliers, 0], score[outliers, 1], color='r', label='Outliers')
plt.xlabel('First Principal Component')
plt.ylabel('Second Principal Component')
plt.title('PCA-Based Outlier Detection')
plt.legend()
plt.show()
