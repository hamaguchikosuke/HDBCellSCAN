# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 16:15:49 2017

@author: hamag
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.datasets as data
# matplotlib inline
sns.set_context('poster')
sns.set_style('white')
sns.set_color_codes()
plot_kwds = {'alpha' : 0.5, 's' : 80, 'linewidths':0}

moons, _ = data.make_moons(n_samples=50, noise=0.05)
blobs, _ = data.make_blobs(n_samples=50, centers=[(-0.75,2.25), (1.0, 2.0)], cluster_std=0.25)
test_data = np.vstack([moons, blobs])
plt.scatter(test_data.T[0], test_data.T[1], color='b', **plot_kwds)

#%%

import sys
#sys.path.append('C:\\Users\\hamag\\.conda\\envs\\khtest\\Lib\\site-packages')
sys.path.append('C:\\Users\\hammer\\AppData\\Local\\conda\\conda\\envs\\CaImaging\\Lib\\site-packages')
import hdbscan

clusterer = hdbscan.HDBSCAN(min_cluster_size=5, gen_min_span_tree=True)
clusterer.fit(test_data)

#HDBSCAN(algorithm='best', alpha=1.0, approx_min_span_tree=True,
#    gen_min_span_tree=True, leaf_size=40, memory=Memory(cachedir=None),
#    metric='euclidean', min_cluster_size=5, min_samples=None, p=None)

#%%
clusterer.minimum_spanning_tree_.plot(edge_cmap='viridis',
                                      edge_alpha=0.6,
                                      node_size=80,
                                      edge_linewidth=2)

#%%
clusterer.single_linkage_tree_.plot(cmap='viridis', colorbar=True)

#%%
clusterer.condensed_tree_.plot()

#%%
palette = sns.color_palette()
cluster_colors = [sns.desaturate(palette[col], sat)
                  if col >= 0 else (0.5, 0.5, 0.5) for col, sat in
                  zip(clusterer.labels_, clusterer.probabilities_)]
plt.scatter(test_data.T[0], test_data.T[1], c=cluster_colors, **plot_kwds)

