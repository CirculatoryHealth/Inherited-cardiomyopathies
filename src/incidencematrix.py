# A function to draw incidence matrices.
#
# It will make a n by m grid of lines  populating the intersections with dots.
# important intersection will be highlighted

# importing
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def draw_incidencematrix(
    data,
    fsize=(6,6),
    dot_colour=[('grey',0), ('black',1)],
    line_colour=['lightgrey', 'lightgrey'],
    dot_size=[4, 8],
    dot_transparency=[0.9, 1.0],
    lw=[0.3, 0.3],
    tick_labs=[4.5, 4.5],
    tick_len = [2, 2],
    tick_wid = [0.3, 0.3],
    margins = None,
    ax=None,
    break_limits = [-np.inf, np.inf]
):
    '''
    Creates a `categorical heatmap` = a visualization of an incidence matrix.
    
    Arguments
    ---------
    Data : pd.DataFrame
        A incidence matrix with index and column labels used as x-axis and
        y-axis tick labels. The matrix values will be plotted basd on the
        `dot_colour` breaks and colours, with a specified size and transparency.
    fsize : tuple
        A two element tuple, with the width and height in cm.
    dot_colour : `list` of `tuple`
        A list of arbitrary length, specifying the colour and upper bound
        the colour is applied to. Each tuple should have
        (<colour>, <upper bound>).
        
        The default: [('grey',0), ('black',1)], colours dots grey for value in
        (\infinity, 0], and colours dots black for values in (0, 1].
    dot_size : list
        A list of length equal to `dot_colour`. specifying the size of the dots.
    dot_transparency : list
        A list of length equal to `dot_colour`, specifying the alpha of the dots.
    line_colour : list
        A two element list specifying the horizontal and vertical line
        colours.
    lw : list
        A two element list specifying the horizontal and vertical line
        size.
    tick_labs : list
        A two element list specifying the label size of the x-, y-axis ticks.
    tick_len : list
        A two element list specifying the length of the x-, y-axis ticks.
    tick_wid : list
        A two element list specifying the width of the x-, y-axis ticks.
    margins : list
        A two element list specifying the margins of the x-, y-axis.
    ax : plt.axes
        An optional matplotlib axis -- will otherwise make one internally
    break_limits : list
        Currently used to specify the lower bound the first colour is applied
        to. Most likely you will never need to touch this.
        
    Returns
    -------
    Unpacks a matplotlib figure, axes, unless `ax` is supplied an plt.axis,
    in which case nothing is returned.
    '''
    
    # check inputs
    if not len(dot_colour) == len(dot_size):
        raise IndexError('The number of `dot_size` entries should equal `dot_colour`')
    if not len(dot_colour) == len(dot_transparency):
        raise IndexError('The number of `dot_transparency` entries should equal `dot_colour`')
    if not isinstance(data, pd.DataFrame):
        raise ValueError('`data` should be a pd.DataFrame')
    
    # do we need to make an axis
    if ax is None:
        cmtoinch = 0.393700787
        f, ax = plt.subplots(figsize=(fsize[0] * cmtoinch,
                                      fsize[1] * cmtoinch))
    else:
        f = None
    
    # the x and y coordinates
    M, N = data.shape
    x, y = np.meshgrid(np.arange(M), np.arange(N))
    
    ################
    # plot dots
    for it, value in enumerate(dot_colour):
        # unpack value
        col, cut = value
        # make breaks
        if it==0:
            sel = (data > break_limits[0]) & (data <= cut)
        # elif it==n:
        #     sel = (data >= cut) & (data < break_limits[1])
        #     print(sel)
        else:
            sel = (data > cut_old) & (data <= cut)
        
        # subset data
        xs = x[sel.to_numpy().T]
        ys = y[sel.to_numpy().T]
        
        # error out if there is no data.
        if xs.size == 0:
            raise ValueError('No data, for cut: ' + str(cut))
        
        # plot
        ax.scatter(xs.flat, ys.flat, c=col, edgecolor=(1, 1, 1, 0),
                   linewidths=0.0,
                   zorder=3, s=dot_size[it], alpha=dot_transparency[it]
                   )
        
        # store previous cut
        cut_old = cut
        # end loop
    ################
    
    # adding grid lines
    for xv in range(x.shape[1]):
        ax.axvline(x=xv, c=line_colour[0], linestyle='-', zorder=1,
                   linewidth=lw[0])
    
    for xy in range(x.shape[0]):
        ax.axhline(y=xy, c=line_colour[1], linestyle='-', zorder=1,
                   linewidth=lw[0])
    
    # ticks
    ax.set(xticks=np.arange(x.shape[1]), yticks=np.arange(x.shape[0]),
           xticklabels=data.index, yticklabels=data.columns)
    ax.tick_params(axis="x", labelsize=tick_labs[0], length=tick_len[0],
                   width=tick_wid[0], rotation=90)
    ax.tick_params(axis="y", labelsize=tick_labs[1], length=tick_len[1],
                   width=tick_wid[1])
    
    # trim margin
    if not margins is None:
        ax.margins(x=margins[0], y=margins[1])
    
    # return the figure and axes
    return f, ax
