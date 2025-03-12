"""
Core API for scanpy_clustering
"""
from typing import Optional
from anndata import AnnData

from .algorithms import get_algorithm


def neighbors(
    adata: AnnData,
    algorithm: str = "default",
    n_neighbors: int = 15,
    **kwargs
) -> None:
    """
    Compute a neighborhood graph of observations.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    algorithm : str
        The algorithm to use for computing the neighborhood graph.
    n_neighbors : int, default: 15
        Number of neighbors to use.
    **kwargs
        Additional arguments to pass to the algorithm.
    """
    # Get the appropriate algorithm implementation
    algo_impl = get_algorithm(algorithm)
    
    # Run the algorithm (when implemented)
    algo_impl.compute_neighbors(adata, n_neighbors=n_neighbors, **kwargs)
    
    return None 