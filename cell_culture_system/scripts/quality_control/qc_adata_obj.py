import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import bbknn

# Install packages for analysis and plotting
from scipy.io import mmread
from sklearn.decomposition import TruncatedSVD
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.sparse import csr_matrix
matplotlib.rcParams.update({'font.size': 22})

### Basic QC, Check the cells in 2D space

# Read in the count matrix that was output by `kb`.
mtx = mmread("/content/counts_unfiltered/cells_x_genes.mtx")

# Plot the cells in the 2D PCA projection
fig, ax = plt.subplots(figsize=(10, 7))

ax.scatter(X[:,0], X[:,1], alpha=0.5, c="green")

plt.axis('off')
plt.savefig("/content/pca_plot.pdf", bbox_inches='tight', pad_inches=0.1)


# density display for PCA plot
from scipy.interpolate import interpn

def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :
    """
    Scatter plot colored by 2d histogram
    """
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins)
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    sc = ax.scatter( x, y, c=z, **kwargs )
    return sc

fig, ax = plt.subplots(figsize=(7,7))

x = X[:,0]
y = X[:,1]

sc = density_scatter(x, y, ax=ax, cmap="Greens")

fig.colorbar(sc, ax=ax)
plt.axis('off')
plt.savefig("/content/density_scatter.pdf", bbox_inches='tight', pad_inches=0.1)

# Create sparse matrix representation of the count matrix
mtx = csr_matrix(mtx)


# Test for library saturation:
# Create a plot showing genes detected as a function of UMI counts.
fig, ax = plt.subplots(figsize=(10, 7))

ax.scatter(np.asarray(mtx.sum(axis=1))[:,0], np.asarray(np.sum(mtx>0, axis=1))[:,0], color="green", alpha=0.01)
ax.set_xlabel("UMI Counts")
ax.set_ylabel("Genes Detected")
ax.set_xscale('log')
ax.set_yscale('log', nonposy='clip')

ax.set_xlim((0.5, 4500))
ax.set_ylim((0.5,2000))

plt.show()
plt.savefig("/content/genes_vs_umi.pdf", bbox_inches='tight', pad_inches=0.1)

# Knee plot for library saturation
# Create the "knee plot"
knee = np.sort((np.array(mtx.sum(axis=1))).flatten())[::-1]
fig, ax = plt.subplots(figsize=(10, 7))

ax.loglog(knee, range(len(knee)),linewidth=5, color="g")

ax.set_xlabel("UMI Counts")
ax.set_ylabel("Set of Barcodes")

plt.grid(True, which="both")
plt.savefig("/content/knee_plot.pdf", bbox_inches='tight', pad_inches=0.1)

