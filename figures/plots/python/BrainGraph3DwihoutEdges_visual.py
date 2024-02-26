import numpy as np
import nibabel as nib
from nilearn import plotting, datasets
from scipy.io import loadmat
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from nilearn import plotting, image
from nilearn.input_data import NiftiLabelsMasker

def load_adjacency_matrix(filename, variable_name):
    """
    Load adjacency matrix from .mat file
    """
    mat_data = loadmat(filename)
    adjacency_matrix = mat_data[variable_name]
    return adjacency_matrix

def normalize_matrix(matrix):
    """
    Normalize a matrix to have values between 0 and 1
    """
    return matrix / np.max(matrix)

def get_atlas_data():
    """
    Load Destrieux atlas data
    """
    destrieux = datasets.fetch_atlas_destrieux_2009()
    atlas_filename = destrieux['maps']
    masker = NiftiLabelsMasker(labels_img=atlas_filename, standardize=True)
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

def plot_connectome(matrix, coordinates, node_sizes, node_colors, show_edges=False):
    """
    Plot 2D and 3D representation of the connections and nodes
    """
    # Only consider nodes that have at least one edge
    nodes_to_include = np.where(node_sizes > 0)[0]
    matrix = matrix[nodes_to_include][:, nodes_to_include]
    coordinates = coordinates[nodes_to_include]
    node_sizes = node_sizes[nodes_to_include]
    node_colors = node_colors[nodes_to_include]

    display = plotting.plot_connectome(matrix, coordinates,
                                       title='Connectome (axial view)',
                                       node_size=node_sizes,
                                       node_color=node_colors,
                                       display_mode='x',
                                       colorbar=True,
                                       edge_cmap='Blues',
                                       edge_threshold=1,
                                       edge_kwargs={'linewidth': 3})
    plt.show()

    # If show_edges is False, set the matrix to all zeros to not display any edges
    if not show_edges:
        matrix = np.zeros_like(matrix)
        min_edge_value = 0  # or any default value you'd like to set
    else:
        min_edge_value = np.min(matrix[matrix > 0])

    custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom_viridis', plt.cm.viridis(np.linspace(0, 1, 256)), N=256)

    view = plotting.view_connectome(matrix, coordinates, node_size=node_sizes, edge_cmap=custom_cmap, edge_threshold=0)
    view.open_in_browser()

if __name__ == "__main__":
    adjacency_matrix = load_adjacency_matrix('symmetric_matrix_visual.mat','A')
    adjacency_matrix_norm = normalize_matrix(adjacency_matrix)
    coordinates, atlas_filename = get_atlas_data()
    node_sizes, node_colors = get_node_attributes(adjacency_matrix_norm)
    plot_connectome(adjacency_matrix_norm, coordinates, node_sizes, node_colors)
