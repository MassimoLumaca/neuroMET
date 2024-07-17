import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

# Get the absolute path of the script file
script_path = os.path.abspath(__file__)

# Extract the directory path
script_dir = os.path.dirname(script_path)

# Change the current working directory to the script directory
os.chdir(script_dir)

# Load the data from Excel instead of CSV
file_path = 'filteredGoldmsi-behavioural_covariates.241.xlsx'  # Updated file path
data = pd.read_excel(file_path)

# Define mappings for various behaviors
mappings = {
    'msi32': {1: '0', 2: '1', 3: '2', 4: '3', 5: '4-5', 6: '6-9', 7: '10 or more'},  # Duration of practice
    'msi33': {1: '0', 2: '0.5', 3: '1', 4: '1.5', 5: '2', 6: '3-4', 7: '5 or more'},  # Hours of Daily Practice
    'msi35': {1: '0', 2: '0.5', 3: '1', 4: '2', 5: '3', 6: '4-6', 7: '7 or more'},  # Music theory
    'msi36': {1: '0', 2: '0.5', 3: '1', 4: '2', 5: '3-5', 6: '6-9', 7: '10 or more'},  # Years of formal training
    'msi37': {1: '0', 2: '1', 3: '2', 4: '3', 5: '4', 6: '5', 7: '6 or more'}  # Instruments Played
}

# Specify the behavior to analyze (e.g., 'msi36')
behavior = 'msi36'  # Change this variable to switch behaviors
# Apply the mapping to the specified behavior column
data[f'{behavior}_mapped'] = data[behavior].replace(mappings[behavior])

# Calculate the percentage distribution for the mapped data
distribution = data[f'{behavior}_mapped'].value_counts(normalize=True).sort_index()

# Create a DataFrame for the visualization
distribution_table = pd.DataFrame({
    'Category': distribution.index,
    'Percentage': (distribution.values * 100).round(2)
})

# Sort the DataFrame, moving 'or more' categories to the end
more_categories = distribution_table['Category'].str.contains('or more')
distribution_table = pd.concat([
    distribution_table[~more_categories],
    distribution_table[more_categories]
]).reset_index(drop=True)

# Create a pie chart with the ordered categories
fig, ax = plt.subplots(figsize=(10, 6))  # Adjust the figure size here
wedges, texts = ax.pie(distribution_table['Percentage'], startangle=90,
                       colors=sns.color_palette("Blues", n_colors=len(distribution_table)),
                       textprops={'fontsize': 8})

# Create the legend with percentage values
labels = [f'{label} : {pct:.1f}%' for label, pct in zip(distribution_table['Category'], distribution_table['Percentage'])]
ax.legend(wedges, labels, title="Category and percentage", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

ax.set_ylabel('')  # Remove the y-label
ax.set_title(f'Pie Chart of {behavior}')

# Define the output PDF file path using the script directory
output_pdf_path = os.path.join(script_dir, f'Pie_Chart_{behavior}.pdf')

# Save the figure to PDF in the specified path, making sure to include the whole figure
fig.savefig(output_pdf_path, format='pdf', bbox_inches='tight')