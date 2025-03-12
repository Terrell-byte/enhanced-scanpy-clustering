"""
Base class for clustering algorithms
"""
from abc import ABC, abstractmethod
from anndata import AnnData
from typing import Optional


class BaseAlgorithm(ABC):
    """
    Base class for clustering algorithms.
    
    All algorithm implementations must inherit from this class
    and implement the required methods.
    """
    
    @abstractmethod
    def cluster(
        self,
        adata: AnnData,
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
        key_added : str, default: 'cluster'
            Key under which to add the cluster labels to adata.obs.
        n_clusters : int, optional
            Number of clusters to find (algorithm-specific).
        **kwargs
            Additional arguments specific to the algorithm.
            
        Returns
        -------
        Updates `adata.obs[key_added]` with cluster assignments.
        """
        pass 