import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import sys
import os
sys.path.insert(
    0,
    os.path.abspath(
        os.path.join(
            os.path.dirname(__file__), '..'
        )
    )
)

class ClusterViz:
    def show_pmra_pmdec_plot(df, xlim=None, ylim=None):
        fig, ax = plt.subplots(figsize=(6,6))
        gr = ax.scatter(df['pmra'], df['pmdec'], s=10, c=df['label'])

        fig.colorbar(gr, ax=ax)
        ax = plt.gca()
        ax.invert_yaxis()
        if xlim is not None:
            plt.xlim(*xlim)
        if ylim is not None:
            plt.ylim(*ylim)

        plt.xlabel(r'$\mu_{\alpha*}$ (mas/yr)')
        plt.ylabel(r'$\mu_{\delta}$ (mas/yr)')

        plt.title("Proper motions of the clusters")

        plt.show()

    def show_label_hist(df):
        plt.figure(figsize=(6, 4))
        plt.hist(df['label'])

        plt.xlabel('Label of Cluster')
        plt.ylabel('Number of Sources')

        plt.title("Cluster label histogram")

        plt.show()

    def show_locations(init_df, clustered_df, clusterer="HDBSCAN"):
        rows = 5
        columns = 5
        grid = plt.GridSpec(rows, columns, wspace = .4, hspace = .4)
        plt.figure(figsize=(7, 7))
        ra_c = np.mean(clustered_df['ra'])
        ra_l = '${:.2f}\degree$'.format(ra_c)
        de_c = np.mean(clustered_df['dec'])
        de_l = '${:.2f}\degree$'.format(de_c)
        plt.subplot(grid[0, 0:-1])
        plt.title("Astrometric location of the cluster members")
        plt.hist(clustered_df['ra'], bins = 30, color = 'royalblue', alpha = .7)
        plt.axvline(ra_c,c='r', label=ra_l)
        plt.legend()
        plt.subplot(grid[1:rows+1, 0:-1])
        plt.plot(init_df['ra'], init_df['dec'], '.', mec='silver', mfc='darkgray', markersize=1., label="All sources")
        plt.plot(clustered_df['ra'], clustered_df['dec'], 'o', mfc='tab:orange', markersize=2., label=clusterer)
        plt.axis('equal')
        plt.xlabel(r'$\alpha$ (deg)')
        plt.ylabel(r'$\delta$ (deg)')
        plt.legend()
        plt.subplot(grid[1:rows+1, -1])
        plt.hist(clustered_df['dec'], bins = 30, orientation='horizontal', color = 'royalblue', alpha = .7)
        plt.axhline(de_c,c='r', label=de_l)
        plt.legend()
        
        plt.show()

    def show_proper_motions(init_df, clustered_df, xlim=None, ylim=None, clusterer="HDBSCAN"):
        rows = 5
        columns = 5
        grid = plt.GridSpec(rows, columns, wspace = .4, hspace = .4)
        plt.figure(figsize=(7, 7))
        pmra_c = np.mean(clustered_df['pmra'])
        pmra_l = '${:.2f}\degree$'.format(pmra_c)
        pmde_c = np.mean(clustered_df['pmdec'])
        pmde_l = '${:.2f}\degree$'.format(pmde_c)
        plt.subplot(grid[0, 0:-1])
        plt.title("Proper motion of the cluster members")
        plt.hist(clustered_df['pmra'], bins = 30, color = 'royalblue', alpha = .7)
        plt.axvline(pmra_c,c='r', label=pmra_l)
        plt.legend()
        plt.subplot(grid[1:rows+1, 0:-1])
        plt.plot(init_df['pmra'], init_df['pmdec'], '.', mec='silver', mfc='darkgray', markersize=1., label="All sources")
        plt.plot(clustered_df['pmra'], clustered_df['pmdec'], 'o', mfc='tab:orange', markersize=2., label=clusterer)
        plt.axis('equal')
        plt.xlabel(r'$\mu_{\alpha*}$ (mas/yr)')
        plt.ylabel(r'$\mu_{\delta}$ (mas/yr)')
        plt.legend()
        if xlim is not None:
            plt.xlim(*xlim)
        if ylim is not None:
            plt.ylim(*ylim)
        plt.subplot(grid[1:rows+1, -1])
        plt.hist(clustered_df['pmdec'], bins = 30, orientation='horizontal', color = 'royalblue', alpha = .7)
        plt.axhline(pmde_c,c='r', label=pmde_l)
        plt.legend()
        plt.show()

    def show_cmd(init_df, clustered_df, xlim=None, clusterer="HDBSCAN"):
        plt.plot(init_df['bp_rp'], init_df['phot_g_mean_mag'], '.', mec='silver', mfc='darkgray', markersize=2., label="All sources")
        plt.plot(clustered_df['bp_rp'], clustered_df['phot_g_mean_mag'], 'o', color='tab:orange', markersize=2., label=clusterer)

        plt.xlabel(r'$G_{BP}-G_{RP}$')
        plt.ylabel(r'$G$ (mag)')

        if xlim is not None:
            plt.xlim(*xlim)

        plt.title("Color magnitude diagram of the cluster")
        
        plt.gca().invert_yaxis()
        plt.legend()
        plt.show()

    def show_parallax_hist(init_df, clustered_df, xlim=None, clusterer="HDBSCAN", show_all=True):
        bins_all = np.arange(init_df['parallax'].min(), init_df['parallax'].max(), .01)
        bins_sam = np.arange(clustered_df['parallax'].min(), clustered_df['parallax'].max(), .01)

        plt.figure(figsize=(6, 4))
        if show_all:
            init_df.parallax.hist(bins=bins_all, color='gray', label="All Sources")
        clustered_df.parallax.hist(bins=bins_sam, color='orange', label=clusterer)
        plt.axvline(np.mean(clustered_df['parallax']),c='r', label=f'Mean: {np.mean(clustered_df["parallax"]):.2f}')

        plt.xlabel(r'$\omega$ (mas)')
        plt.ylabel('Number of Sources')

        plt.title("Parallax histogram of the cluster")

        if xlim is not None:
            plt.xlim(*xlim)

        plt.xticks()
        plt.yticks()

        plt.legend()
        plt.show()

    def get_cluster_info(clustered_df):
        cluster_info = {}
        print("Cluster info:\n")

        for metrics in ['ra', 'dec', 'pmra', 'pmdec', 'parallax']:
            cluster_info[f"{metrics}_c"] = np.mean(clustered_df[metrics])
            print(f'{metrics}_c\t\t\t= ', '{:.3f}'.format(cluster_info[f"{metrics}_c"]))
        
        cluster_info["dist_c"] = 1/cluster_info["parallax_c"] if cluster_info["parallax_c"] > 0 else -1
        print('dist_c\t\t\t= ', '{:.3f}'.format(cluster_info[f"dist_c"]))

        return cluster_info