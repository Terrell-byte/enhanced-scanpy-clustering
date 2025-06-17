"""
Base class for clustering algorithms
"""
from abc import ABC, abstractmethod
from anndata import AnnData
from typing import Optional

import numpy as np


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
        **kwargs
    ) -> None:
        """
        Perform clustering on the data.
        
        Parameters
        ----------
        adata : AnnData
            Annotated data matrix.
        key_added : str, default: 'cluster' (Naming for the cluster labels)
            Key under which to add the cluster labels to adata.obs.
        **kwargs
            Additional arguments specific to the algorithm.
            
        Returns
        -------
        Updates `adata.obs[key_added]` with cluster assignments.
        """
        pass 

    @classmethod
    def register(cls):
        """
        Register the algorithm class with the global registry.
        """
        from scanpy_clustering.algorithms import register_algorithm
        register_algorithm(cls.__name__, cls)

    @staticmethod
    def _convert_anndata(adata: AnnData) -> np.ndarray:
        """
        Convert the AnnData object to a NumPy array suitable for clustering.

        Parameters
        ----------
        adata : AnnData
            Annotated data matrix.

        Returns
        -------
        np.ndarray
            The data matrix as a dense NumPy array.
        """
        x = adata.X
        if not isinstance(x, np.ndarray):
            x = x.toarray()  # Convert sparse matrix to dense if needed
        return x
