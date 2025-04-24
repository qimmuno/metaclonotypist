import numpy as np
import pandas as pd
from statsmodels.stats.proportion import test_proportions_2indep
from statsmodels.stats.multitest import fdrcorrection
from pyrepseq import *
import scipy.sparse

def metaclonotypist(df, chain='beta',
                    max_edits=2, max_tcrdist=20,
                    clustering='cc', clustering_kwargs=dict(),
                    node_label='bioidentity'):
    """Integrated pipeline for metaclonotype identification.

    chain: str
        'alpha' or 'beta'
        used for identifying df columns to use in `pyrepseq.nearest_neighbor_tcrdist`
    max_edits: int
        maximal number of edits in CDR3 to consider in prefiltering step
        The default value of 2 should provide a reasonable starting point
        Set to 1 to speed up neighbor finding in large datasets.
        Caution: Values >2 might slow down neighbor finding a lot.
    max_tcrdist: int
        maximal TCRdist score used to prune the initial edit distance based neighbors
        Optimal choices might be between 0 (nearly exact matches only in large datasets)
        to ~100 (to find any clustering in small datasets)
    clustering: str
        choice of clustering algorithm
        see `pyrepseq.graph_clustering` for details
    node_label : str
        column name to use for labelling nodes in the resulting graph
        
    """
    neighbors = nearest_neighbor_tcrdist(df, max_edits=max_edits, max_tcrdist=max_tcrdist, chain=chain)
    return graph_clustering(neighbors, df[node_label], clustering=clustering, **clustering_kwargs)

def metaclonotypist_sceptr(df, chain='beta',
                    max_edits=2, max_sceptrdist=20,
                    clustering='cc', clustering_kwargs=dict(),
                    node_label='bioidentity'):
    """Integrated pipeline for metaclonotype identification.

    chain: str
        'alpha' or 'beta'
        used for identifying df columns to use in `pyrepseq.nearest_neighbor_tcrdist`
    max_edits: int
        maximal number of edits in CDR3 to consider in prefiltering step
        The default value of 2 should provide a reasonable starting point
        Set to 1 to speed up neighbor finding in large datasets.
        Caution: Values >2 might slow down neighbor finding a lot.
    max_sceptrdist: int
        maximal SCEPTR score used to prune the initial edit distance based neighbors
        Optimal choices might be between 0 (nearly exact matches only in large datasets)
        to ~1.5 (to find any clustering in small datasets)
    clustering: str
        choice of clustering algorithm
        see `pyrepseq.graph_clustering` for details
    node_label : str
        column name to use for labelling nodes in the resulting graph
        
    """
    neighbors = nearest_neighbor_sceptrdist(df, max_edits=max_edits, max_sceptrdist=max_sceptrdist, chain=chain)
    return graph_clustering(neighbors, df[node_label], clustering=clustering, **clustering_kwargs)

def flatten_hlas(hla_table):
    """Reshape a HLA table into a Boolean presence/absence format."""
    all_hlas = list(sorted(set(np.asarray(hla_table).flatten())-set([np.nan])))
    columns = []
    for hla in all_hlas:
        columns.append(hla_table.apply(lambda row: hla in set(row), axis=1))
    columns_arr = np.asarray(columns).T
    return pd.DataFrame(data=columns_arr, columns=all_hlas, index=columns[0].index)

def hla_association(clusters, hla_table, fdr_alpha=0.1, method='fisher',
                    sampleid='Sample.ID'):
    """Cluster HLA association test.

    clusters : pd.DataFrame
        needs columns 'cluster' and 'Sample.ID'
    hla_table: pd.DataFrame
        each row is a sample, each column presence/absence of a HLA allele
    fdr_alpha: float
        False Discovery Rate threshold for Benjamini/Hochberg procedure
    method: ['fisher' or 'agresti-caffo']
        test procedure to determine p-values
    sampleid: Sample ID column name
    """
    results = []
    hla_all_counts = hla_table.sum(axis=0)
    for cid, cluster in clusters.groupby('cluster'):
        cluster_individuals = cluster[sampleid].unique()
        counts = hla_table.loc[cluster_individuals].sum(axis=0)
        hla_other_counts = hla_all_counts - counts
        for hla in counts.index:
            nobs1, nobs2 = len(cluster_individuals), len(hla_table)-len(cluster_individuals)
            count1, count2 = counts.loc[hla], hla_other_counts.loc[hla]
            with np.errstate(all='ignore'):
                if method == 'fisher':
                    table = [[count1, nobs1-count1], [count2, nobs2-count2]]
                    pvalue = scipy.stats.fisher_exact(table, alternative='greater').pvalue
                    odds_ratio = scipy.stats.contingency.odds_ratio(
                                np.array(table, dtype=int), kind='sample').statistic
                elif method == 'agresti-caffo':
                    res = test_proportions_2indep(count1, nobs1,
                                                  count2, nobs2,
                                          alternative='larger',
                                          method='agresti-caffo', compare='diff')
                    pvalue = res.pvalue
                    pvalue = res.pvalue
                    odds_ratio = res.odds_ratio
                else:
                    raise NotImplementedError(f'{method} is not implemnted')
            result = cid, hla, count1, nobs1, count2, nobs2, pvalue, odds_ratio
            results.append(result)
    cluster_association = pd.DataFrame(results,
                                columns=['cluster', 'hla', 'count_allele', 'total_allele', 'count_other', 'total_other', 'pvalue', 'odds_ratio'])
    cluster_association['significant'] = fdrcorrection(cluster_association['pvalue'], alpha=fdr_alpha)[0]
    return cluster_association

def contingency_matrix(labels_true, labels_pred, eps=None, sparse=False, dtype=np.int64):
    """Build a contingency matrix describing the relationship between labels.

    Parameters
    ----------
    labels_true : array-like of shape (n_samples,)
        Ground truth class labels to be used as a reference.

    labels_pred : array-like of shape (n_samples,)
        Cluster labels to evaluate.

    eps : float, default=None
        If a float, that value is added to all values in the contingency
        matrix. This helps to stop NaN propagation.
        If ``None``, nothing is adjusted.

    sparse : bool, default=False
        If `True`, return a sparse CSR continency matrix. If `eps` is not
        `None` and `sparse` is `True` will raise ValueError.

    dtype : numeric type, default=np.int64
        Output dtype. Ignored if `eps` is not `None`.

    Returns
    -------
    contingency : {array-like, sparse}, shape=[n_classes_true, n_classes_pred]
        Matrix :math:`C` such that :math:`C_{i, j}` is the number of samples in
        true class :math:`i` and in predicted class :math:`j`. If
        ``eps is None``, the dtype of this array will be integer unless set
        otherwise with the ``dtype`` argument. If ``eps`` is given, the dtype
        will be float.
        Will be a ``csr_matrix`` if ``sparse=True``.

    Notes
    ------
    This code is adapted from a corresponding scikit-learn function (BSD-3 licensed)
    """

    if eps is not None and sparse:
        raise ValueError("Cannot set 'eps' when sparse=True")

    classes, class_idx = np.unique(labels_true, return_inverse=True)
    clusters, cluster_idx = np.unique(labels_pred, return_inverse=True)
    n_classes = classes.shape[0]
    n_clusters = clusters.shape[0]
    # Using coo_matrix to accelerate simple histogram calculation,
    # i.e. bins are consecutive integers
    # Currently, coo_matrix is faster than histogram2d for simple cases
    contingency = scipy.sparse.coo_matrix(
        (np.ones(class_idx.shape[0]), (class_idx, cluster_idx)),
        shape=(n_classes, n_clusters),
        dtype=dtype,
    )
    if sparse:
        contingency = contingency.tocsr()
        contingency.sum_duplicates()
    else:
        contingency = contingency.toarray()
        if eps is not None:
            # don't use += as contingency is integer
            contingency = contingency + eps
    return contingency

def entropy(labels):
    """Calculate the entropy for a labeling.

    Parameters
    ----------
    labels : array-like of shape (n_samples,), dtype=int
        The labels.

    Returns
    -------
    entropy : float
       The entropy for a labeling.

    Notes
    -----
    The logarithm used is the natural logarithm (base-e).
    """
    if len(labels) == 0:
        return 1.0

    pi = np.unique(labels, return_counts=True, equal_nan=False)[1].astype(np.float64)
    pi /= np.sum(pi)

    if pi.size == 1:
        return 0.0

    return -np.sum(pi * np.log(pi))

def conditional_entropies(labels_true, labels_pred):
    """Calculate conditional entropies of ground truth classes and predicted clusters.

    Parameters
    ----------
    labels_true : array-like of shape (n_samples,)
        Ground truth class labels to be used as a reference.

    labels_pred : array-like of shape (n_samples,)
        Cluster labels to evaluate.

    Returns
    -------
    entropy_joint, entropy_CK, entropy_KC, entropy_C, entropy_K
    """


    entropy_C = entropy(labels_true)
    entropy_K = entropy(labels_pred)

    contingency = contingency_matrix(labels_true, labels_pred, sparse=True)

    if isinstance(contingency, np.ndarray):
        # For an array
        nzx, nzy = np.nonzero(contingency)
        nz_val = contingency[nzx, nzy]
    else:
        # For a sparse matrix
        nzx, nzy, nz_val = scipy.sparse.find(contingency)

    nsum = np.sum(nz_val)
    entropy_joint = -np.sum(nz_val * np.log(nz_val))/nsum + np.log(nsum)

    entropy_CK = entropy_joint - entropy_K
    entropy_KC = entropy_joint - entropy_C
    return entropy_joint, entropy_CK, entropy_KC, entropy_C, entropy_K

def compression_score(labels_true, labels_predicted):
    entropy_joint, entropy_CK, entropy_KC, entropy_C, entropy_K = conditional_entropies(labels_true, labels_predicted)

    return 1-entropy_KC/(np.log(len(labels_true))-entropy_C)

