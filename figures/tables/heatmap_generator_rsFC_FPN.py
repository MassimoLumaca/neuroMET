import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
data_wmi = pd.read_csv('PercentageTotalMETScore2rsFCfrontoparietalnetworkTable.csv')

# Define the label mappings
x_label_mapping = {
    "E_global": "Global Efficiency",
    "Local_Efficiency": "Local Efficiency",
    "Clustering_Coefficient": "Clustering Coefficient",
    "Betweenness_Centrality": "Betweenness Centrality"
}

# Labels

y_label_mapping = {
    "lG_front_middle": "l MFG",
    "lG_front_sup": "l SupFG",
    "lG_pariet_inf.Angular": "l AngG",
    "lG_pariet_inf.Supramar": "l SuMarG",
    "lG_parietal_sup": "l SupPL",
    "lS_front_middle": "l MFS",
    "lS_front_sup":"l SupFs",
    "lS_intrapariet_and_P_trans": "l IntPS",
    "rG_front_middle": "r MFG",
    "rG_front_sup": "r SupFG",
    "rG_pariet_inf.Angular": "r AngG",
    "rG_pariet_inf.Supramar": "r SuMarG",
    "rG_parietal_sup": "r SupPL",
    "rS_front_middle": "r MFS",
    "rS_front_sup":"r SupFs",
    "rS_intrapariet_and_P_trans": "r IntPS"
}


# Set NodeNames as the index and rename columns/rows based on label mappings
data_wmi.set_index('NodeNames', inplace=True)
data_wmi.rename(columns=x_label_mapping, index=y_label_mapping, inplace=True)

# Reorder the columns
column_order = ["Global Efficiency", "Local Efficiency", "Clustering Coefficient", "Betweenness Centrality"]
data_wmi = data_wmi[column_order]

# Create the heatmap with a greenish color scheme
# For fronto-parietal networks use cmap="YlGnBu"; for visual networks use cmap="YlGn"
# Create the heatmap with a greenish color scheme and set color bar range
plt.figure(figsize=(12, 10))
sns.heatmap(data_wmi, cmap="YlGnBu", annot=True, fmt=".2f", linewidths=.5, vmin=-3.5, vmax=3.5)
plt.title("Metrics Heatmap for Nodes")
plt.savefig('FrontoParietal_percentageTOTMET_rs-FC_Heatmap.pdf', format='pdf', bbox_inches='tight')  # Save the heatmap as a PDF figure
plt.show()


