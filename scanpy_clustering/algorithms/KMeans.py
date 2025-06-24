from sklearn.cluster import KMeans as km
from scanpy_clustering.algorithms.BaseAlgorithm import BaseAlgorithm

class KMeans(BaseAlgorithm):
    def cluster(self, 
                adata, 
                key_added = 'cluster', 
                **kwargs) -> None:
        """
        Perform KMeans clustering on the data.

        Parameters
        ----------
        adata : AnnData
            Annotated data matrix.
        key_added : str, default: 'cluster'
            Key under which to add the cluster labels to adata.obs.
        n_clusters : int, default=8
            The number of clusters to form as well as the number of centroids to generate.
        init : {'k-means++', 'random'} or ndarray, default='k-means++'
            Method for initialization.
        n_init : int, default=10
            Number of time the k-means algorithm will be run with different centroid seeds.
        max_iter : int, default=300
            Maximum number of iterations of the k-means algorithm for a single run.
        tol : float, default=1e-4
            Relative tolerance with regards to inertia to declare convergence.
        verbose : int, default=0
            Verbosity mode.
        random_state : int, RandomState instance or None, default=None
            Determines random number generation for centroid initialization.
        copy_x : bool, default=True
            When pre-computing distances it is more numerically accurate to center the data first.
        algorithm : {"auto", "full", "elkan"}, default="auto"
            K-means algorithm to use.

        Returns
        -------
        None
            The cluster labels are added to adata.obs[key_added].
        """
        # Extract feature matrix
        data = self._convert_anndata(adata)

        # Apply KMeans
        clustering = km(**kwargs).fit(data)
        # Store results
        adata.obs[key_added] = clustering.labels_.astype(str)  # Convert to string to avoid categorical issues