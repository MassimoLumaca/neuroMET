# Load the required libraries
if (!requireNamespace("lavaan", quietly = TRUE)) {
  install.packages("lavaan")
}
if (!requireNamespace("semPlot", quietly = TRUE)) {
  install.packages("semPlot")
}
library(lavaan)
library(semPlot)

##### Specify main variables and parameters #####
WAIS_name <- "WMI"
NeuralMetricName <- "E_global_roi_9"

# Load Data
graph_metrics_rois <- read.csv("graph_metrics_rois_iq.csv")
NeuralMetric <- graph_metrics_rois[[NeuralMetricName]]
behavioural_covariates <- read.csv("behavioural_covariates_iq.csv")

# Define variables
Musicality <- behavioural_covariates$PercentageOfTotalMETScore
WAIS_Subscale <- behavioural_covariates[[WAIS_name]]
Age <- behavioural_covariates$Age
Sex <- behavioural_covariates$Gender
MusicalTraining <- behavioural_covariates$F3
mydata <- data.frame(NeuralMetric, WAIS_Subscale, Musicality, Age, Sex, MusicalTraining)

# SEM model definition
model <- '
  WAIS_Subscale ~ c*NeuralMetric 
  Musicality ~ a*WAIS_Subscale + b*NeuralMetric 
  indirect := c*a
  total := c*a + b
'

# Fit the SEM model
fit <- sem(model, data = mydata, se = "bootstrap", bootstrap = 1000)

# Capture output
parameter_output <- capture.output(print(parameterEstimates(fit, boot.ci.type="bca.simple")))
param_estimates <- parameterEstimates(fit, standardized=TRUE)

# Create a new data frame with path labels and standardized coefficients
path_coefficients_df <- data.frame(
  Path = paste(param_estimates$lhs, param_estimates$op, param_estimates$rhs),
  Standardized_Coefficient = param_estimates$std.all
)

# Print the data frame to view the labeled path coefficients
print(path_coefficients_df)

# Optionally, save the labeled path coefficients to a CSV file
csv_filename <- paste0("labeled_path_coefficients_", WAIS_name, ".csv")
write.csv(path_coefficients_df, csv_filename, row.names = FALSE)


# Define the PDF filenames based on the WAIS subscale
filename1 <- paste0("mediationanalysis_mediator_", WAIS_name, "_NoCorrection1.pdf")
filename2 <- paste0("mediationanalysis_mediator_", WAIS_name, "_NoCorrection2.pdf")
final_filename <- paste0("mediationanalysis_mediator_", WAIS_name, "NoCorrection.pdf")

# Start the PDF device for text
pdf(filename1, width = 10, height = 11)

# Create a blank plot for textual output
plot(1, type="n", axes=FALSE, xlab="", ylab="", xlim=0:1, ylim=0:1)
text(0.05, 0.9, paste(parameter_output, collapse="\n"), adj=c(0,1), cex=0.65, family="mono")

# Close the PDF device
dev.off()

# Save the SEM paths plot to a separate PDF
pdf(filename2, width = 8, height = 11)
print(semPaths(fit, whatLabels = "est.std", layout = "tree", style="lisrel", 
               curveAdjacent = TRUE, edge.label.cex = 0.7))
dev.off()

# Combine the two PDFs (requires the 'pdftk' command line tool)
system(paste("pdftk", filename1, filename2, "cat output", final_filename))

# Delete the intermediate PDFs
file.remove(filename1, filename2)

