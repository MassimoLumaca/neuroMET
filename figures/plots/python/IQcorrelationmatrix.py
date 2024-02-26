import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests
from itertools import product

# Load the data
data = pd.read_csv('MET_IQ_scores.csv')

# Compute the correlation matrix
correlation_matrix = data.corr()

# Extract the lower triangle of the correlation matrix to avoid duplicate calculations
lower_triangle_indices = np.tril_indices_from(correlation_matrix, k=-1)

# Generate column pairs for the lower triangle
column_names = correlation_matrix.columns
columns_lower_triangle = [(column_names[i], column_names[j]) for i, j in zip(*lower_triangle_indices)]

# Calculate p-values for all pairwise correlations in the lower triangle
p_values = [pearsonr(data[col1], data[col2])[1] for col1, col2 in columns_lower_triangle]

# Apply FDR correction
_, p_values_corrected, _, _ = multipletests(p_values, method='fdr_bh')

# Dictionary to store corrected p-values for each pair of variables
p_value_dict = {(col1, col2): p_val for (col1, col2), p_val in zip(columns_lower_triangle, p_values_corrected)}

def annotate_plot_with_p_value_and_fdr(x, y, **kws):
    """Function to annotate the regression plots with R^2 value, p-value, and asterisk for FDR"""
    r, _ = pearsonr(x, y)
    ax = plt.gca()
    col1, col2 = ax.get_xlabel(), ax.get_ylabel()
    corrected_p = p_value_dict.get((col2, col1), 1)  # Reverse order because of how pairplot is plotted
    asterisk = '*' if corrected_p < 0.05 else ''
    ax.text(0.05, 0.9, 'R = {:.2f}'.format(r), transform=ax.transAxes, fontsize=12)
    ax.text(0.05, 0.8, 'pFDR = {:.3f}{}'.format(corrected_p, asterisk), transform=ax.transAxes, fontsize=12)

# Plot settings
plot_kwargs = dict(marker='o', scatter_kws={"color": "black"}, line_kws={'color': 'red'}, ci=95)
diag_kwargs = dict(fill=True, color="grey")  # Updated to use `fill` instead of `shade`

# Create a pairplot with regression lines and grey KDE plots
g = sns.pairplot(data, kind='reg', corner=True, diag_kind="kde", 
                 markers="o", plot_kws=plot_kwargs, diag_kws=diag_kwargs)

# Annotate each scatterplot with the R value, p-value, and asterisk for FDR
g.map_lower(annotate_plot_with_p_value_and_fdr)

plt.savefig('MET_IQ_scores.pdf', format='pdf')
plt.show()
