# packages 
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

import matplotlib as mpl
from matplotlib.ticker import MultipleLocator,FormatStrFormatter,MaxNLocator
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
from matplotlib.transforms import ScaledTranslation

pd.set_option('display.max_colwidth', None)

# plot global parameters
mpl.rcParams['figure.figsize'] = [1, 0.8]

# display axis spines
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False

# font
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# font sizes
fontSize = 7
mpl.rcParams["axes.labelsize"]= fontSize
mpl.rcParams[ "xtick.labelsize"]= fontSize
mpl.rcParams["ytick.labelsize"]= fontSize
mpl.rcParams["legend.fontsize"]= fontSize
mpl.rcParams["font.size"]= fontSize

# line widths
mpl.rcParams["axes.linewidth"]= 0.5
mpl.rcParams["grid.linewidth"]= 0.5
mpl.rcParams["lines.linewidth"]= 1.
mpl.rcParams["lines.markersize"]= 3

# tick sizes
mpl.rcParams["xtick.major.size"]= 2
mpl.rcParams["ytick.major.size"]= 2
mpl.rcParams["xtick.major.width"]= 0.5
mpl.rcParams["ytick.major.width"]= 0.5
# mpl.rcParams["xtick.top"]= False
# mpl.rcParams["ytick.right"]= False

# save figure
mpl.rcParams["savefig.bbox"]= "tight"
mpl.rcParams["savefig.pad_inches"]= 0.05


# function here
def Gini(array):
    """Calculate the Gini coefficient of a numpy array."""
    array = array.to_numpy()
    array = array.flatten() #all values are treated equally, arrays must be 1d
    array = np.maximum(array, 0)
    array = np.add(array, 0.0000001, out=array, casting="unsafe") #values cannot be 0
    array = np.sort(array) # values must be sorted
    index = np.arange(1,array.shape[0]+1) #index per array element
    n = array.shape[0] # number of array elements
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array))) #Gini coefficient

def calculate_gini_index(values):
    """Calculate the Gini index of an array of values."""
    # Mean absolute difference
    mad = np.abs(np.subtract.outer(values, values)).mean()
    # Relative mean absolute difference
    rmad = mad / np.mean(values)
    # Gini coefficient
    g = 0.5 * rmad
    return g

# color
color = ["#F5E9CA", "#DDAAAB", "#B2CFEB", "#CDE6E0", "#C0B6D9"]

def plot_sig(xstart,xend, yend, sig):
    plt.hlines(yend, xstart, xend, color="black", linewidth = 0.5)
    plt.annotate(r'%s'%sig, xy=((xstart + xend)/2, yend+0.005), color="black", ha='center')