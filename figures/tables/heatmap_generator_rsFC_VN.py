import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
data_wmi = pd.read_csv('PercentageTotalMETScore2rsFCvisualnetworkTable.csv')

# Define the label mappings
x_label_mapping = {
    "E_global": "Global Efficiency",
    "Local_Efficiency": "Local Efficiency",
    "Clustering_Coefficient": "Clustering Coefficient",
    "Betweenness_Centrality": "Betweenness Centrality"
}

# Labels

 y_label_mapping = {
     "lG_and_S_occipital_inf": "l InfOcG/S",
     "lG_occipital_middle": "l MOcG",
     "lG_occipital_sup": "l SupOcG",
     "lG_oc.temp_lat.fusifor": "l FuG",
     "lG_oc.temp_med.Lingual": "l Cos/Lins",
     "lPole_occipital": "l OcPo",
     "lS_occipital_ant":"l AOcS",
     "lS_parieto_occipital": "l PoCs",
     "rG_and_S_occipital_inf": "r InfOcG/S",
     "rG_occipital_middle": "r MOcG",
     "rG_occipital_sup": "r SupOcG",
     "rG_oc.temp_lat.fusifor": "r FuG",
     "rG_oc.temp_med.Lingual": "r Cos/Lins",
     "rPole_occipital": "r OcPo",
     "rS_occipital_ant":"r AOcS",
     "rS_parieto_occipital": "r PoCs"
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
plt.savefig('Visual_percentageTOTMET_rs-FC_Heatmap.pdf', format='pdf', bbox_inches='tight')  # Save the heatmap as a PDF figure
plt.show()


