from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon



## HUJI Colour Palette ##
d_black = [35./255,31./255,32./255]			    #	Black				BlackC
d_olive	= [159.0/255.0,161.0/255.0,97.0/255.0] 	# 	Olive Green	5835C
d_blue	= [0,99.0/255.0,136.0/255.0] 			# 	Blue			634C
d_red	= [170.0/255.0,43.0/255.0,74.0/255.0] 	#   Red			201C
d_midblue = [145./255,184.0/255.0,189.0/255.0]	#   Mid Blue		5493C
d_pinck=[227.0/255.0, 37.0/255.0, 132.0/255.0]  #   Pinck
d_green=[84.0/255.0, 207.0/255.0, 146.0/255.0]  # Green
d_orange=[248.0/255.0, 69.0/255.0, 123.0/255.0] #Orange
dd_blue=[64.0/255.0, 210.0/255.0, 254.0/255.0]
d_green=[94.0/255.0, 184.0/255.0, 122.0/255.0]  #green

#Pastel colors
pastel_green=[201.0/255.0,255.0/255.0,178.0/255.0]
pastel_blue=[194.0/255.0,178.0/255.0,255.0/255.0]
pastel_yellow=[255.0/255.0,242.0/255.0,178.0/255.0]
pastel_light_blue=[178.0/255.0,225.0/255.0,255.0/255.0]
pastel_red=[255.0/255.0,178.0/255.0,185.0/255.0]
pastel_torqize=[178.0/255.0,255.0/255.0,218.0/255.0]
pastel_purpule=[255.0/255.0,178.0/255.0,249.0/255.0]





# update matplotlib fonts etc
plt.rc('font',**{'family':'Serif','serif':['Times New Roman']})
params={'axes.labelsize':16,'xtick.labelsize':16,'ytick.labelsize':16,'legend.fontsize': 16,'axes.linewidth':2}
plt.rcParams.update(params)


def gradient_fill(x, y, fill_color=None, ax=None, **kwargs):
    """
    Plot a line with a linear alpha gradient filled beneath it.

    Parameters
    ----------
    x, y : array-like
        The data values of the line.
    fill_color : a matplotlib color specifier (string, tuple) or None
        The color for the fill. If None, the color of the line will be used.
    ax : a matplotlib Axes instance
        The axes to plot on. If None, the current pyplot axes will be used.
    Additional arguments are passed on to matplotlib's ``plot`` function.

    Returns
    -------
    line : a Line2D instance
        The line plotted.
    im : an AxesImage instance
        The transparent gradient clipped to just the area beneath the curve.
    """
    if ax is None:
        ax = plt.gca()

    line, = ax.plot(x, y, **kwargs)
    if fill_color is None:
        fill_color = line.get_color()

    zorder = line.get_zorder()
    alpha = line.get_alpha()
    alpha = 1.0 if alpha is None else alpha

    z = np.empty((100, 1, 4), dtype=float)
    rgb = mcolors.colorConverter.to_rgb(fill_color)
    z[:,:,:3] = rgb
    z[:,:,-1] = np.linspace(0, alpha, 100)[:,None]

    xmin, xmax, ymin, ymax = x.min(), x.max(), y.min(), y.max()
    im = ax.imshow(z, aspect='auto', extent=[xmin, xmax, ymin, ymax],
                   origin='lower', zorder=zorder)

    xy = np.column_stack([x, y])
    xy = np.vstack([[xmin, ymin], xy, [xmax, ymin], [xmin, ymin]])
    clip_path = Polygon(xy, facecolor='none', edgecolor='none', closed=True)
    ax.add_patch(clip_path)
    im.set_clip_path(clip_path)

    ax.autoscale(True)
    return line, im