"""
Base class for clustering algorithms
"""
from abc import ABC, abstractmethod
from anndata import AnnData


class BaseAlgorithm(ABC):
    """
    Base class for clustering algorithms.
    
    All algorithm implementations must inherit from this class
    and implement the required methods.
    """
    
    @abstractmethod
    def compute_neighbors(
        self,
        adata: AnnData,
        n_neighbors: int = 15,
        **kwargs
    ) -> None:
        """
        Compute a neighborhood graph of observations.
        
        Parameters
        ----------
        adata : AnnData
            Annotated data matrix.
        n_neighbors : int, default: 15
            Number of neighbors to use.
        **kwargs
            Additional arguments specific to the algorithm.
            
        Returns
        -------
        Updates `adata` with neighbors data.
        """
        pass 