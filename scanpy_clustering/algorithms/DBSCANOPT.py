import numpy as np
from anndata import AnnData
from typing import List
from collections import deque
from scanpy_clustering.algorithms.BaseAlgorithm import BaseAlgorithm

class DBSCANOPT(BaseAlgorithm):
    """
    Custom DBSCAN clustering algorithm implementation.
    """

    def _euclidean_distance(self, p1: np.ndarray, p2: np.ndarray) -> float:
        """Compute the Euclidean distance between two points."""
        p1 = p1.toarray().flatten() if hasattr(p1, "toarray") else np.array(p1)
        p2 = p2.toarray().flatten() if hasattr(p2, "toarray") else np.array(p2)
        return np.sqrt(np.sum(((p1 - p2) ** 2)))

    def _region_query(self, X: np.ndarray, point_idx: int) -> List[int]:
        """Find all neighbors of a given point within epsilon distance."""
        neighbors = []
        for i in range(X.shape[0]):
            if self._euclidean_distance(X[point_idx], X[i]) < self.eps:
                neighbors.append(i)
        return neighbors

    def _expand_cluster(self, X: np.ndarray, labels: np.ndarray, point_idx: int, cluster_id: int, NeighborPts: List[int]) -> None:
        """Expand the cluster from the given core point."""
        queue = deque(NeighborPts)
        labels[point_idx] = cluster_id

        while queue:
            neighbor_idx = queue.popleft()
            if labels[neighbor_idx] == -1:  # Previously labeled as noise
                labels[neighbor_idx] = cluster_id
            if labels[neighbor_idx] != 0:  # Already processed
                continue

            labels[neighbor_idx] = cluster_id
            new_neighbors = self._region_query(X, neighbor_idx)
            if len(new_neighbors) >= self.min_samples:
                queue.extend(new_neighbors)

    def cluster(
        self,
        adata: AnnData,
        key_added: str = 'cluster',
        eps: float = 0.5,
        min_samples: int = 5
    ) -> None:
        """
        Perform DBSCAN clustering on the data.

        Parameters
        ----------
        adata : AnnData
            Annotated data matrix.
        key_added : str, default: 'cluster'
            Key under which to add the cluster labels to adata.obs.
        eps : float, default: 0.5
            The maximum distance between two samples for them to be considered as in the same neighborhood.
        min_samples : int, default: 5
            The number of samples required to form a dense region.

        Returns
        -------
        None
            Updates `adata.obs[key_added]` with cluster assignments.
        """
        if adata.X is None:
            raise ValueError("AnnData object does not contain data in `.X`.")

        X = adata.X

        self.eps = eps
        self.min_samples = min_samples

        labels = np.zeros(X.shape[0], dtype=int)  # 0 = unvisited
        cluster_id = 0
        for i in range(X.shape[0]):
            if labels[i] != 0:  # Already visited
                continue

            neighbors = self._region_query(X, i)
            if len(neighbors) < self.min_samples:
                labels[i] = -1  # Mark as noise
            else:
                cluster_id += 1
                self._expand_cluster(X, labels, i, cluster_id, neighbors)

        # Store results in AnnData
        adata.obs[key_added] = labels.astype(str)
