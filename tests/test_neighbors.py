"""
Tests for the neighbors functionality
"""
import pytest
import numpy as np
import scanpy as sc
from anndata import AnnData

from scanpy_clustering import neighbors


def test_neighbors_basic():
    """Test basic functionality of the neighbors function."""
    # Create a simple test dataset
    X = np.random.random((100, 20))
    adata = AnnData(X)
    
    # Pre-process
    sc.pp.pca(adata)
    
    # Run our neighbors function
    neighbors(adata, algorithm="algo1", n_neighbors=10)
    
    # Check that the neighbors were computed
    assert "neighbors" in adata.uns
    assert "connectivities" in adata.obsp
    assert "distances" in adata.obsp 