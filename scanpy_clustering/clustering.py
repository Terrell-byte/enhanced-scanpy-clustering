"""
Core API for scanpy_clustering
"""
from typing import Optional
from anndata import AnnData

from .algorithms import get_algorithm


def cluster(
    adata: AnnData,
    algorithm: str = "default",
    key_added: str = 'cluster',
    n_clusters: Optional[int] = None,
    **kwargs
) -> None:
    """
    Perform clustering on the data.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    algorithm : str
        The algorithm to use for clustering.
    key_added : str, default: 'cluster'
        Key under which to add the cluster labels to adata.obs.
    n_clusters : int, optional
        Number of clusters (algorithm-specific).
    **kwargs
        Additional arguments to pass to the algorithm.
    """
    # Get the appropriate algorithm implementation
    algo_impl = get_algorithm(algorithm)
    
    # Run the algorithm
    algo_impl.cluster(adata, key_added=key_added, n_clusters=n_clusters, **kwargs)
    
    return None 