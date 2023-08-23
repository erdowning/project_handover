import numpy as np
import matplotlib.pyplot as plt

def hist(data, bins, weights, orientation='vertical', **kwargs):
    hist, bin_edges = np.histogram(data, bins=bins, weights=weights)
    centres = (bin_edges[:-1] + bin_edges[1:]) / 2
    if orientation=='vertical':
        plt.plot(centres, hist, **kwargs)
    elif orientation=='horizontal':
        plt.plot(hist, centres, **kwargs)
    return hist, centres