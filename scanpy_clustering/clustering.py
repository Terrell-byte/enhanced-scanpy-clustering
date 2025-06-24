"""
Core API for scanpy_clustering
"""
import scanpy as sc
import inspect
from anndata import AnnData
from scanpy_clustering.algorithms import _ALGORITHMS

def enable_scanpy_integration():
    """
    Monkey-patches Scanpy to include the new clustering algorithm in sc.tl.
    Users must explicitly call this function to enable it.
    """
    for name, algorithm in _ALGORITHMS.items():
        setattr(sc.tl, name, algorithm().cluster)

def list_algorithms() -> list:
    """Returns the clustering algorithms and their 'cluster' method parameters."""
    algorithms_info = []
    
    for name, algo_class in _ALGORITHMS.items():
        if hasattr(algo_class, 'cluster'):
            sig = inspect.signature(algo_class.cluster)
            params = {
                k: v.default
                for k, v in sig.parameters.items()
                if k not in ('self', 'adata')
            }
            algorithms_info.append({
                'algorithm': name,
                'parameters': params
            })
    
    return algorithms_info

def cluster(
    adata: AnnData,
    algorithm: str = "default",
    key_added: str = 'cluster',
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
    **kwargs
        Additional arguments to pass to the algorithm.
        """
    
    if algorithm not in _ALGORITHMS:
        raise ValueError(f"Unknown clustering method: {algorithm}. Available options: {list_algorithms()}")

    clustering_algo = _ALGORITHMS[algorithm]()
    clustering_algo.cluster(adata, key_added=key_added, **kwargs)