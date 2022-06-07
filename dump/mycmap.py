import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import pickle


####################### COLORMAP -- DIVERGENT
cdict1  = {'red':  ((0.0, 0.0, 0.0),
                   (0.35, 0.0, 0.0),
                   (0.48, 1.0, 1.0),
                   (0.5, 1.0, 1.0),
                   (0.52, 1.0, 1.0),
                   (0.65, 1.0, 1.0),
                   (0.8, 1.0, 1.0),
                   (1.0, 0.6, 0.6)),

         'green': ((0.0, 0.0, 0.0),
                   (0.2, 0.4, 0.4),
                   (0.35, 0.75, 0.75),
                   (0.48, 1.0, 1.0),
                   (0.5, 1.0, 1.0),
                   (0.52, 1.0, 1.0),
                   (0.65, 0.7, 0.7),
                   (0.8, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.6, 0.6),
                   (0.2, 0.8, 0.8),
                   (0.35, 0.6, 0.6),
                   (0.45, 1.0,1.0),
                   (0.5, 1.0, 1.0),
                   (0.55, 1.0, 1.0),
                   (0.65, 0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }

        
mycmaps = {
	'BuGnWYlRd': mcolors.LinearSegmentedColormap('BuGnWYlRd', cdict1),
}

