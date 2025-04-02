import importlib
import numpy as np
import scanpy as sc
import anndata as ad
import scanpy_clustering.clustering as cl
import time
from sklearn.metrics.cluster import adjusted_rand_score

print('Running')
# Load example data
adata = sc.read_h5ad('data\\symsim_observed_counts_5000genes_5000cells_complex.h5ad')[:750, :30]  # Example dataset from Scanpy

# Start timer
tic = time.perf_counter()

# Run DBScan clustering with custom implementation
cl.cluster(adata, algorithm='DBScan_Base', key_added='dbscan_labels', eps=25, min_samples=4, metric='euclidean')

# Get time taken for custom implementation
toc = time.perf_counter()
print(f"Ran in: {toc - tic:0.4f} seconds")

# Check results
print(adata.obs)

# Check rand_score
from sklearn.metrics.cluster import adjusted_rand_score
adjusted_rand_score(adata.obs['group'], adata.obs['dbscan_labels'])
