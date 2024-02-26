import numpy as np
import nibabel as nib
from nilearn import plotting, datasets
from scipy.io import loadmat
import matplotlib.pyplot as plt

def load_adjacency_matrix(filename, variable_name):
    """
    Load adjacency matrix from .mat file
    """
    mat_data = loadmat(filename)
    adjacency_matrix = mat_data[variable_name]
    return adjacency_matrix

def amplify_differences(matrix, power=2):
    """
    Apply a power transformation to amplify differences.
    """
    amplified_matrix = np.power(matrix, power)
    return amplified_matrix

def get_atlas_data():
    """
    Load Destrieux atlas data
    """
    destrieux = datasets.fetch_atlas_destrieux_2009()
    atlas_filename = destrieux['maps']
    coordinates = plotting.find_parcellation_cut_coords(labels_img=atlas_filename)
    return coordinates, atlas_filename

def get_node_attributes(matrix):
    """
    Create list of node sizes, colors and alpha values
    """
    large_node_value = 20  # Increased node size
    small_node_value = 0  # Nodes with no edges will not be visible
    large_node_color = (1, 0, 0, 1)
    small_node_color = (0, 0, 0, 0)
    
    # Compute node sizes
    node_sizes = np.where(np.sum(matrix, axis=1) > 0, large_node_value, small_node_value)

    # Compute node colors. Initialize with small_node_color and then modify based on condition.
    node_colors = np.full((matrix.shape[0], 4), small_node_color)
    node_colors[np.where(np.max(matrix, axis=1) == 1)] = large_node_color
    
    return node_sizes, node_colors

def plot_connectome(matrix, coordinates, node_sizes, node_colors):
    """
    Plot 2D and 3D representation of the connections and nodes
    """
    # Only consider nodes that have at least one edge
    nodes_to_include = np.where(node_sizes > 0)[0]
    matrix = matrix[nodes_to_include][:, nodes_to_include]
    coordinates = coordinates[nodes_to_include]
    node_sizes = node_sizes[nodes_to_include]
    node_colors = node_colors[nodes_to_include]

    # Plotting the 2D connectome
    display = plotting.plot_connectome(matrix, coordinates,
                                       node_size=node_sizes,
                                       node_color=node_colors,
                                       display_mode='y',
                                       colorbar=True,
                                       edge_cmap='Blues',
                                       edge_threshold='0%')  # edges are color-coded based on their strength
    plt.show()

    # Plotting the 3D connectome
    view = plotting.view_connectome(matrix, coordinates, node_size=node_sizes, edge_cmap='Blues', edge_threshold='0%', colorbar=True)
    view.open_in_browser()

if __name__ == "__main__":
    adjacency_matrix = load_adjacency_matrix('AvAdjacencyMatprop.mat', 'Csym')
    
    # Amplify differences
    amplified_matrix = amplify_differences(adjacency_matrix, power=2)
    
    coordinates, atlas_filename = get_atlas_data()
    node_sizes, node_colors = get_node_attributes(amplified_matrix)  # Use amplified matrix
    plot_connectome(amplified_matrix, coordinates, node_sizes, node_colors)

def plot_matrix_heatmap(matrix, indices):
    """
    Plot the connectome matrix as a heatmap.
    """
    # Subset matrix based on desired indices
    subset_matrix = matrix[np.ix_(indices, indices)]
    
    # Plot heatmap
    plt.imshow(subset_matrix, cmap='Blues', interpolation='nearest')
    plt.colorbar()
    plt.title('Connectome Matrix Heatmap')
    plt.show()

if __name__ == "__main__":
    adjacency_matrix = load_adjacency_matrix('AvAdjacencyMatprop.rsFC.FPN.mat', 'Csym')
    
    # Amplify differences
    amplified_matrix = amplify_differences(adjacency_matrix, power=2)
    
    coordinates, atlas_filename = get_atlas_data()
    node_sizes, node_colors = get_node_attributes(amplified_matrix)  # Use amplified matrix
    
    # Indices of interest from the Destrieux labeling
    desired_indices = [14, 52, 15, 53, 24, 25, 26, 55, 88, 126, 89, 127, 98, 99, 100, 129]  # Python uses 0-based indexing
    
    # Plot the matrix heatmap
    plot_matrix_heatmap(amplified_matrix, desired_indices)
    
    # Plot the connectome
    plot_connectome(amplified_matrix, coordinates, node_sizes, node_colors)