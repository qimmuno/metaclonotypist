from .main import *
import argparse
import os
import matplotlib.pyplot as plt
import seaborn as sns

from importlib.metadata import version


__version__ = version("metaclonotypist")

def script():
    """Entry point for the console script"""

    parser = argparse.ArgumentParser()
    
    parser.add_argument("--tcrpath",
                        required=True,
                        help="Path to input TCR data (CSV file)")
    
    parser.add_argument("--hlapath",
                        required=True,
                        help="Path to input HLA metadata (CSV file)")
    
    parser.add_argument('-o', '--output-dir',
                        required=True,
                        help='Path to the output directory')
    
    parser.add_argument("--chain", default='beta', help="chain to use (default: beta)",
                        choices=['alpha', 'beta'])
    
    parser.add_argument("--tcrdistmethod", default='tcrdist',
                        help="TCR distance method (default: tcrdist)",
                        choices=['tcrdist', 'sceptr'])
    
    parser.add_argument("--mincount", type=int, default=None,
                        help="Minimum count for clones (default: None, no filtering)")

    parser.add_argument("--maxtcrdist", type=float, default=15,
                        help="Maximum TCR distance (default: 15)")

    parser.add_argument("--clustering", default='leiden',
                        help="Clustering algorithm (default: leiden)",
                        choices=['leiden', 'multilevel'])
    
    parser.add_argument("--hlatest", default='fisher',
                        help="Statistical test method for HLA association (default: fisher)",
                        choices=['fisher', 'agresti-caffo'])
    
    parser.add_argument("--mindonors", type=int, default=4,
                        help="Minimum number of donors for HLA filtering (default: 4)")

    parser.add_argument("--maxedits", type=int, default=2,
                        help="Maximum edits for TCR distance (default: 2)")

    parser.add_argument("--version", action='version', version=f'%(prog)s {__version__}',
                        help="Show the version of Metaclonotypist")
    

    args = parser.parse_args()

    tcrpath = args.tcrpath
    hlapath = args.hlapath
    outputdir = args.output_dir
    os.makedirs(outputdir, exist_ok=True)
    chain = args.chain
    chain_suffix = 'A' if chain == 'alpha' else 'B'
    tcrdistmethod = args.tcrdistmethod
    mincount = args.mincount
    maxtcrdist = args.maxtcrdist
    clustering = args.clustering
    hlatest = args.hlatest

    maxedits = args.maxedits
    mindonors = args.mindonors

    if clustering == 'multilevel':
        clustering_kwargs = dict(resolution=10)
    elif clustering == 'leiden':
        clustering_kwargs = dict(resolution=0.1,
                             objective_function='CPM',
                             n_iterations=4)
    else:
        raise NotImplementedError(f'Unknown clustering method: {clustering}')


    print('Loading data')
    df = pd.read_csv(tcrpath)
    if not 'bioidentity' in df.columns:
        df['bioidentity'] = (df[f'TR{chain_suffix}V'] + '_'
                             + df[f'CDR3{chain_suffix}'] + '_'
                             + df[f'TR{chain_suffix}J'])
    hlas = flatten_hlas(pd.read_csv(hlapath, index_col=0))

    # filter clones < mincount
    if mincount is not None:
        df = df[df['clonal_count']>=mincount]
    # only keep samples found in both datasets
    df = df[df['Sample.ID'].isin(hlas.index)]
    hlas = hlas.loc[list(set(df['Sample.ID']))]
    # filter hlas < min_donors
    hlas = hlas[hlas.columns[hlas.sum(axis=0)>=mindonors]]

    print('Starting clustering')
    if tcrdistmethod == 'sceptr':
        clusters = metaclonotypist_sceptr(df, chain=chain,
                               max_sceptrdist=maxtcrdist, max_edits=maxedits,
                               clustering=clustering, clustering_kwargs=clustering_kwargs)
    elif tcrdistmethod == 'tcrdist':
        clusters = metaclonotypist(df, chain=chain,
                               max_tcrdist=maxtcrdist, max_edits=maxedits,
                               clustering=clustering, clustering_kwargs=clustering_kwargs)
    else:
        raise NotImplementedError(f'Unknown TCR distance method: {tcrdist}')
    clusters['Sample.ID'] = df.loc[clusters.index]['Sample.ID']


    print('Starting HLA association')

    # filter clusters with too few donors
    ndonors = clusters.groupby('cluster').apply(lambda cluster: len(cluster['Sample.ID'].unique()))
    clusters = clusters[clusters['cluster'].isin(ndonors[(ndonors >= mindonors)].index)]
       
    cluster_association = hla_association(clusters, hlas, method=hlatest)
    hla_metaclones = cluster_association[cluster_association['significant']]
    hla_metaclones = hla_metaclones.drop(columns=['significant'])

    params_str = f'{tcrdistmethod}_{chain}_td{maxtcrdist:.3f}_{hlatest}_{clustering}'
    params_str += f'_mincount{mincount}' if mincount is not None else ''
    params_str += f'_maxedits{maxedits}' if maxedits != 2 else ''
    #clusters.to_csv(f'{outputdir}/clustering_{params_str}.csv')
    #cluster_association.to_csv(f'{outputdir}/clusters_associations_{params_str}.csv', index=False)
    hla_metaclones.to_csv(f'{outputdir}/cluster_associations_{params_str}.csv', index=False)

    sig_clusters = clusters[(clusters['cluster'].isin(hla_metaclones['cluster']))].reset_index()
    sig_clusters = sig_clusters.merge(hla_metaclones, on='cluster')
    hla_match = [hlas.loc[row['Sample.ID']][row['hla']] for ind, row in sig_clusters.iterrows()]
    sig_clusters = sig_clusters.iloc[hla_match]
    sig_clusters.to_csv(f'{outputdir}/clustering_{params_str}.csv', index=False)


    # shuffle hlas
    hlas_shuffled = hlas.copy()
    hlas_shuffled.index = np.random.permutation(hlas_shuffled.index)
    cluster_association_shuffled = hla_association(clusters, hlas_shuffled, method=hlatest)
    #cluster_association_shuffled.to_csv(f'{outputdir}/clusters_association_shuffled_{params_str}.csv', index=False)
    hla_metaclones_shuffled = cluster_association_shuffled[cluster_association_shuffled['significant']]
    sig_clusters_shuffled = clusters[(clusters['cluster'].isin(hla_metaclones_shuffled['cluster']))].reset_index()
    hla_match_shuffled = [hlas.loc[row['Sample.ID']][row['hla']] for ind, row in sig_clusters_shuffled.iterrows()]
    sig_clusters_shuffled = sig_clusters_shuffled.iloc[hla_match_shuffled]

    data = [['tcrdistmethod', tcrdistmethod],
            ['chain', chain],
            ['maxtcrdist', maxtcrdist],
            ['mincount', mincount],
            ['hlatest', hlatest],
            ['clustering', clustering],
            ['nassociations', len(hla_metaclones)],
            ['nassociations_shuffled', len(hla_metaclones_shuffled)],
            ['nmetaclones', len(hla_metaclones['cluster'].unique())],
            ['nmetaclones_shuffled', len(hla_metaclones_shuffled['cluster'].unique())],
            ['clustered_fraction', len(clusters)/len(df)],
            ['sig_clonotype_fraction', len(sig_clusters)/len(df)],
            ['sig_clonotype_fraction_shuffled', len(sig_clusters_shuffled)/len(df)],
            ['sig_read_fraction', df.loc[sig_clusters['index']]['clonal_count'].sum()/df['clonal_count'].sum()],
            ['sig_read_fraction_shuffled', df.loc[sig_clusters_shuffled['index']]['clonal_count'].sum()/df['clonal_count'].sum()],
            ['id_fraction', len(sig_clusters['Sample.ID'].unique())/len(df['Sample.ID'].unique())],
            ['id_fraction_shuffled', len(sig_clusters_shuffled['Sample.ID'].unique())/len(df['Sample.ID'].unique())],
           ]
    index, values = list(zip(*data))
    pd.Series(index=index, data=values, name='results').to_csv(f'{outputdir}/stats_{params_str}.csv', index=True)

    print('Starting plotting')

    cluster_association.replace(np.inf, 200, inplace=True)
    cluster_association.replace(-np.inf, -200, inplace=True)
    cluster_association_shuffled.replace(np.inf, 200, inplace=True)
    cluster_association_shuffled.replace(-np.inf, -200, inplace=True)

    fig, axes = plt.subplots(ncols=2, sharex=True, sharey=True)
    sns.scatterplot(ax=axes[0], data=cluster_association,
                    x='odds_ratio',
                    y=-np.log10(cluster_association['pvalue']),
                    hue='significant',
                    s=5)
    sns.scatterplot(ax=axes[1], data=cluster_association_shuffled,
                    x='odds_ratio',
                    y=-np.log10(cluster_association_shuffled['pvalue']),
                    hue='significant',
                    s=5)
    axes[0].set_title('Data')
    axes[1].set_title('Shuffled HLA')
    for ax in axes:
        ax.set_xscale('log')
        ax.set_ylabel('p value')
        ax.set_xlabel('odds ratio')
    fig.tight_layout()
    fig.savefig(f'{outputdir}/volcano_plot_{params_str}.png')
