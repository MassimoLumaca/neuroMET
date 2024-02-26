library(ggplot2)
library(dplyr)
setwd("/Users/au540169/Desktop/MINDLAB2016_MR-SensCognFromNeural ")

graph_metrics_rois <- read.csv("graph_metrics_rois.SC.csv")
behavioural_covariates <- read.csv("behavioural_covariates.SC.csv")

# Assuming you want to merge on row names
merged_data <- merge(graph_metrics_rois, behavioural_covariates, by = "row.names")

# Select the target graph theory metrics (e.g., Betweenness_Centrality) and node (e.g., roi_10) in the following part of the script
# List of nodes: 
# Left hemisphere
#roi_1 =  middle frontal gyrus; roi_2 =  superior frontal gyrus; roi_3 =  angular gyrus
# roi_4 =  supramarginal gyrus; roi_5 =  superior parietal gyrus; roi_6 =  middle frontal sulcus;
# roi_7 =  superior frontal sulcus; roi_8 =  intraparietal sulcus

# Right hemisphere
#roi_9 = middle frontal gyrus; roi_10 = superior frontal gyrus; roi_11 =  angular gyrus
# roi_12 =  supramarginal gyrus; roi_13 =  superior parietal gyrus; roi_14 =  middle frontal sulcus;
# roi_15 =  superior frontal sulcus; roi_16 = intraparietal sulcus

correlation_test <- cor.test(merged_data$Betweenness_Centrality_roi_10, merged_data$PercentageOfTotalMETScore)

# Extract the correlation coefficient (R)
R_value <- correlation_test$estimate

print(paste("Correlation (R):", R_value))

# Create the plot
p <- ggplot(data = merged_data, aes(x = Betweenness_Centrality_roi_10, y = PercentageOfTotalMETScore)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", level = 0.95) +
  labs(x = "Eglobal rMFG", y = "Percentage Total MET Score", title = "Scatterplot with 95% Confidence Interval") +
  theme_minimal() +
  scale_x_continuous(limits = c(min(merged_data$Betweenness_Centrality_roi_10) * 0.95, max(merged_data$Betweenness_Centrality_roi_10) * 1.05)) +
  theme(aspect.ratio = 1)

# Print the plot
print(p)









