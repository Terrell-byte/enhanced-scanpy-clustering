import numpy as np
import scanpy as sc
import anndata as ad
import scanpy_clustering.clustering as cl
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

print('Running')
# Load example data
adata = sc.read_h5ad('data\\symsim_observed_counts_5000genes_5000cells_complex.h5ad')[:750, :30]  # Example dataset from Scanpy

# Load example dataset (replace with your AnnData object)
#adata = sc.datasets.pbmc3k()
X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X  # Convert sparse if needed

'''
# Set k (min_samples for DBSCAN, usually 4 or 5)
k = 5

# Compute k-nearest neighbors using KDTree
tree = KDTree(X)
distances, _ = tree.query(X, k=k+1)  # k+1 because first neighbor is itself

# Get the k-th nearest neighbor distance for each point
k_distances = distances[:, -1]  # Last column gives k-th neighbor distance

# Sort distances in ascending order
k_distances.sort()

# Plot k-distance graph
plt.figure(figsize=(8, 6))
plt.plot(np.arange(len(k_distances)), k_distances, marker='o', linestyle='-')
plt.xlabel("Data Point Index (sorted)")
plt.ylabel(f"Distance to {k}-th Nearest Neighbor")
plt.title(f"k-Distance Graph for k={k}")
plt.grid()
plt.show() '''

cl.cluster(adata, algorithm='DBSCANOPT', key_added='dbscan_labels', eps=25, min_samples=4, metric='euclidean')
# Check results
print(adata.obs)