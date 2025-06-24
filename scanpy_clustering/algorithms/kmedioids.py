from sklearn_extra.cluster import KMedoids as kmedoids
from scanpy_clustering.algorithms.BaseAlgorithm import BaseAlgorithm

class KMedoids(BaseAlgorithm):
    def cluster(self, 
                adata, 
                key_added = 'cluster', 
                **kwargs) -> None:
        """
        Perform KMedoids clustering on the data.

        Parameters
        ----------
        adata : AnnData
            Annotated data matrix.
        key_added : str, default: 'cluster'
            Key under which to add the cluster labels to adata.obs.
        **kwargs
            Additional arguments to pass to sklearn_extra.cluster.KMedoids.

        Returns
        -------
        None
            The cluster labels are added to adata.obs[key_added].
        """
        # Extract feature matrix
        data = self._convert_anndata(adata)

        # Apply KMedoids
        clustering = kmedoids(**kwargs).fit(data)
        # Store results
        adata.obs[key_added] = clustering.labels_.astype(str)  # Convert to string to avoid categorical issues