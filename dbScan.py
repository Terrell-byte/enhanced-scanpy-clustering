import numpy as np
import scanpy as sc
import anndata as ad
from sklearn.cluster import DBSCAN

def dbscan_ann_data(adata: ad.AnnData, eps=0.5, min_samples=5, metric='euclidean'):
    """
    Applies DBScan clustering to an AnnData object.

    Parameters:
    - adata: AnnData object
    - eps: The maximum distance between two samples for them to be considered as in the same neighborhood.
    - min_samples: The number of samples required to form a dense region.
    - metric: The distance metric to use (default: 'euclidean').

    Returns:
    - The modified AnnData object with cluster labels stored in `adata.obs['dbscan_labels']`
    """

    # Extract feature matrix
    X = adata.X
    if isinstance(X, np.ndarray):
        data = X
    else:
        data = X.toarray()  # Convert sparse matrix to dense

    # Apply DBScan
    clustering = DBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit(data)

    # Store results
    adata.obs['dbscan_labels'] = clustering.labels_.astype(str)  # Convert to string to avoid categorical issues

    return adata

print('Running')
# Load example data
adata = sc.datasets.pbmc3k()  # Example dataset from Scanpy

# Apply DBScan
adata = dbscan_ann_data(adata, eps=0.5, min_samples=10)

# Check results
print(adata.obs['dbscan_labels'].value_counts())
