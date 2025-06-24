from sklearn.cluster import DBSCAN
from scanpy_clustering.algorithms.base import BaseAlgorithm

class DBScan_Base(BaseAlgorithm):
    def cluster(self, adata, key_added = 'cluster', **kwargs):
        """
        Perform DBSCAN clustering on the data.

        Parameters
        ----------
        adata : AnnData
            Annotated data matrix.
        key_added : str, default: 'cluster'
            Key under which to add the cluster labels to adata.obs.
        eps : float, default=0.5
            The maximum distance between two samples for one to be considered as in the neighborhood of the other.
        min_samples : int, default=5
            The number of samples (or total weight) in a neighborhood for a point to be considered as a core point.
        metric : str or callable, default='euclidean'
            The metric to use when calculating distance between instances in a feature array.
        metric_params : dict, default=None
            Additional keyword arguments for the metric function.
        algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, default='auto'
            The algorithm to be used by the NearestNeighbors module to compute pointwise distances.
        leaf_size : int, default=30
            Leaf size passed to BallTree or KDTree.
        p : float, default=2
            The power of the Minkowski metric to be used to calculate distance between points.
        n_jobs : int, default=None
            The number of parallel jobs to run for neighbors search.
        """

        # Extract feature matrix
        data = self._convert_anndata(adata)

        # Apply DBScan
        clustering = DBSCAN(**kwargs).fit(data)

        # Store results
        adata.obs[key_added] = clustering.labels_.astype(str)  # Convert to string to avoid categorical issues