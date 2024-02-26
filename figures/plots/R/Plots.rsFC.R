library(ggplot2)
library(dplyr)

setwd("/Users/au540169/Desktop/MINDLAB2016_MR-SensCognFromNeural")

graph_metrics_rois <- read.csv("graph_metrics_rois.rsFC.csv")
behavioural_covariates <- read.csv("behavioural_covariates.rsFC.csv")
# Assuming you want to merge on row names
merged_data <- merge(graph_metrics_rois, behavioural_covariates, by = "row.names")

# Create the plot
p <- ggplot(data = merged_data, aes(x = E_global_roi_9, y = PercentageOfTotalMETScore)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", level = 0.95) +
  labs(x = "Eglobal rMFG", y = "Percentage Total MET Score", title = "Scatterplot with 95% Confidence Interval") +
  theme_minimal() +
  scale_x_continuous(limits = c(min(merged_data$Clustering_Coefficient_roi_10) * 0.95, max(merged_data$Clustering_Coefficient_roi_10) * 1.05)) +
  theme(aspect.ratio = 1)

# Print the plot
print(p)
