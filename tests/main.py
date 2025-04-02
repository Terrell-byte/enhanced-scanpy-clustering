import numpy as np
import scanpy as sc
import anndata as ad
from scanpy_clustering import clustering as cl

print('Running')
# Load example data
adata = sc.read_h5ad('data\\symsim_observed_counts_5000genes_5000cells_complex.h5ad')[:750, :30]  # Example dataset from Scanpy

# Apply DBScan
#cl.enable_scanpy_integration()
#print(cl.list_algorithms())
#sc.tl.DBScan(adata, key_added='dbscan_labels', eps=0.5, min_samples=10, metric='euclidean')
#sc.tl.DBScan.cluster(adata, eps=0.5, min_samples=10, metric='euclidean')
cl.cluster(adata, algorithm='DBSCANOPT', key_added='dbscan_labels', eps=2500, min_samples=30, metric='euclidean')
# Check results
#print(adata.obs['dbscan_labels'].value_counts())