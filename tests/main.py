import numpy as np
import scanpy as sc
import anndata as ad
from scanpy_clustering.algorithms.dbScan import db
from sklearn.cluster import DBSCAN

print('Running')
# Load example data
adata = sc.read_h5ad('data\symsim_observed_counts_5000genes_5000cells_complex.h5ad')  # Example dataset from Scanpy

# Apply DBScan
adata = db(adata, eps=0.5, min_samples=10)

# Check results
print(adata.obs['dbscan_labels'].value_counts())