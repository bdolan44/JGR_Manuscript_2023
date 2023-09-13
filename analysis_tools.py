# Brody Fuchs, CSU, July 2014
# brfuchs@atmos.colostate.edu

# functions to do analysis and plotting for CLEAR output
# and just a whole bunch of other damn things

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import matplotlib.ticker as plticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patheffects as PE
import datetime
import matplotlib.dates as md
from scipy import interpolate
from copy import deepcopy
import netCDF4

def convolve_smooth(xarr, yarr, window = 'boxcar', width = 5, mode = 'valid'):

    odd = np.fmod(width, 2) # checking to see if it's odd
    if not odd: width += 1

    cmode = deepcopy(mode)

    shorten =  int((width-1)/2)
    # this is filtering of y array stuff

    if window == 'boxcar':
        filt = np.ones((width,))/width

    if window == '121':
        filt = np.ones(width)*(1-(np.abs(np.arange(width)-(width-1)/2)/(width/2)))
        filt /= np.sum(filt) # this is to normalize


    if mode == 'extend':
        yarr = np.insert(yarr, 0, np.ones(shorten)*yarr[0])
        yarr = np.insert(yarr, -1, np.ones(shorten)*yarr[-1])
        cmode = 'valid'

    yout = np.convolve(yarr, filt, mode = cmode)
# this is for cutting out parts of the x array
    if mode == 'valid': 
        xout = xarr[shorten:-1*shorten]

    if mode == 'same' or mode == 'extend':
        xout = xarr

    



    return yout, xout

# put xout second in case you don't want it, can just throw away



def expand_axis(ax, which = 'x', factor = 1.05):
    if which == 'x':
        lims = ax.get_xlim()
        axis_range = lims[1]-lims[0]
        ax.set_xlim(lims[0]-axis_range*factor/2, lims[1]+axis_range*factor/2)	

    elif which == 'y':
        lims = ax.get_ylim()
        axis_range = lims[1]-lims[0]
        ax.set_ylim(lims[0]-axis_range*factor/2, lims[1]+axis_range*factor/2)	


def interpolate_data(x, y, new_x, method = 'spline'):
    if method == 'spline':
        tck = interpolate.splrep(x, y, s = 0)
        ynew = interpolate.splev(new_x, tck, der=0)

    elif method == 'interp1d':
        f = interpolate.interp1d(x, y)
        ynew = f(new_x)

    return ynew

def label_subplots(fig, xoff = 0.0, yoff = 0.02, nlabels = None,**kwargs):
    letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 
        'i', 'j', 'k', 'l', 'm','n','o']
    figaxes = fig.get_axes()
    if nlabels is None: nlabels = len(figaxes)

    for fa in range(nlabels):
        xbox = figaxes[fa].get_position()
        xmin, ymax = xbox.xmin, xbox.ymax
        # this is the position I want
        fig.text(xmin+xoff, ymax+yoff, '({})'.format(letters[fa]),**kwargs)
    

def set_ax_min(ax, value = 0, xory = 'x'):
    # set the minimum value of an axis, usually making the minimum 0, which is cumbersome
    if xory == 'x':
        axlim = ax.get_xlim()
        ax.set_xlim(value, axlim[1])     
    elif xory == 'y':
        axlim = ax.get_ylim()
        ax.set_ylim(value, axlim[1])     
    else:
        print ('xory needs to be x or y axis')
    return


def bp_custom(plot_handle):

    plt.setp(plot_handle['boxes'], color = 'black', linewidth = 1.5, facecolor='0.7')
    plt.setp(plot_handle['whiskers'], color = 'black', linewidth = 1.5)
    plt.setp(plot_handle['caps'], color = 'black', linewidth=1.5)
    plt.setp(plot_handle['medians'], color = 'darkorange', linewidth=1.5)


def histogram_average(values, hist_values, bins):
    # given a bunch of hist_values that put data into bins
    # give the average of values of all points within a bin
    dig = np.digitize(hist_values, bins)
    out = np.array( [np.average(values[dig == _]) for _ in range(bins.shape[0])] )
    return out


def axis_text(ax, xpos, ypos, axtext, **kwargs):
    ax_xlim, ax_ylim = ax.get_xlim(), ax.get_ylim()
    ax.text(ax_xlim[0]+(ax_xlim[1]-ax_xlim[0])*xpos, ax_ylim[0]+(ax_ylim[1]-ax_ylim[0])*ypos, axtext, **kwargs)

def edit_existing_colormap(cmap_string):
    # way to take an existing python colormap and edit it

    jetcmap = cm.get_cmap("jet", 10) #generate a jet map with 10 values 
    jet_vals = jetcmap(np.arange(10)) #extract those values as an array 
    jet_vals[0] = [0.1, 0, 0.1, 1] #change the first value 
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", jet_vals)
    return newcmap

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.itervalues():
        sp.set_visible(False)

def multiple_axes_one_legend(ax, handles, **kwargs):
    labs = [l.get_label() for l in handles]
    ax.legend(handles, labs, **kwargs)

def remove_top_and_right(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


def additional_yaxis(parent_ax, fig = plt.gcf(), side = 'right', xloc = None, color = 'black', min_val = None):
    # add another y axis to a plot, and change its color and position, returns axis so can actually plot to it
    new_ax = parent_ax.twinx()

    if xloc is not None:
        new_ax.spines[side].set_position(('axes',xloc))
        make_patch_spines_invisible(new_ax)
        new_ax.spines[side].set_visible(True)

    new_ax.yaxis.set_label_position(side)
    new_ax.yaxis.set_ticks_position(side)


    if color != 'black':
        for i in new_ax.get_yticklabels():
                i.set_color(color)

    return new_ax


def subplot_label(ax, letter, xpos = None, ypos = 1.03, **kwargs):
    axlim = ax.get_xlim()
    aylim = ax.get_ylim()
    if xpos is None: xpos = axlim[0]
    ax.text(xpos,aylim[1]*ypos,'(%s)'%(letter), horizontalalignment = 'left', fontsize = 10, **kwargs)


def custom_colorbar(handle, xpos = 0.95, xwidth = 0.02, ypos = 0.05, ywidth = 0.9, \
        figure = plt.gcf(), **kwargs):
    "A custom colorbar"
    
    cbar_ax = figure.add_axes([xpos, ypos, xwidth, ywidth])
    cbar = plt.colorbar(handle, cax = cbar_ax, **kwargs)
    cbar.set_alpha(1)
    cbar.draw_all()

    return cbar

@np.vectorize
def round_hour(dt_obj):
    # rounds the datetime object to the nearest hour and returns the new object
    if dt_obj.minute >= 30:
        new_dt = dt_obj + datetime.timedelta(hours = 1)
    else: new_dt = dt_obj
    out = datetime.datetime(new_dt.year, new_dt.month, new_dt.day, new_dt.hour, 0, 0)

    return out

@np.vectorize
def LST(dt, lon = 0):
	# this returns the local solar time when a time and longitude are input

    gamma = 2 * np.pi / 365 * (dt.timetuple().tm_yday - 1 + float(dt.hour - 12) / 24)
    eqtime = 229.18 * (0.000075 + 0.001868 * np.cos(gamma) - 0.032077 * np.sin(gamma) \
             - 0.014615 * np.cos(2 * gamma) - 0.040849 * np.sin(2 * gamma))
    decl = 0.006918 - 0.399912 * np.cos(gamma) + 0.070257 * np.sin(gamma) \
           - 0.006758 * np.cos(2 * gamma) + 0.000907 * np.sin(2 * gamma) \
           - 0.002697 * np.cos(3 * gamma) + 0.00148 * np.sin(3 * gamma)
    time_offset = eqtime + 4 * lon
    tst = dt.hour * 60 + dt.minute + dt.second / 60 + time_offset
    solar_time = datetime.datetime.combine(dt.date(), datetime.time(0)) + datetime.timedelta(minutes=tst)
    return solar_time

def day_of_year(dt_obj):
    tt = dt_obj.timetuple()
    float_tt = tt.tm_yday + dt_obj.hour/24. + dt_obj.minute/1440. + dt_obj.second/86400.
    # this actually includes the hours, minutes and seconds
    return float_tt

class TimeHistogram(object):

    def __init__(self, dates, times, bins = 30):

        self.date_array = np.array([ datetime.datetime( 2000+int(_[0:2]), int(_[3:5]), int(_[6:8]), \
        int(np.floor(__/3600.)), int(np.floor(np.mod(__, 3600)/60.0)) ) for _, __ in zip(dates, times) ])

        self.datenum = np.array([md.date2num(_) for _ in self.date_array])
        self.hist, self.edges = np.histogram(self.datenum, bins = bins)
        
        self.dtedges = [md.num2date(_) for _ in self.edges]


    def plot_hist(self, fig = None, ax = None, **kwargs):

        if fig is None:	fig = plt.figure()
        # make new figure object if one not provided
        if ax is None: ax = plt.gca()
        ax.plot(self.dtedges[:-1], self.hist, **kwargs)

        return fig



def filter_dict(dict, condition):

    "automatically filters dictionary for a given condition then returns filtered dict"
    out = {}

    for k in dict.keys():

        out[k] = dict[k][condition]

    return out

def dict_concatenate(old_dict, new_dict,mykeys):

    "loops thru keys of each dict and concatenates them"
#    assert old_dict.keys() == new_dict.keys()
    out = {}
    for key in mykeys:
        out[key] = np.concatenate((old_dict[key], new_dict[key]))

    return out

def dict_concat_auto(indict):
    # this will take a dictionary made of smaller dicts that all have matching keys
    # it'll then make one large dict of the concatenated smaller dicts
    # grab the keys from the first larger key
    print(list(indict.keys())[0])
    keys = indict[list(indict.keys())[0]].keys()
    print ('keys: {}'.format(keys))
    # now check to make sure all the keys are the same
#     for k in indict.keys():
#         assert indict[k].keys() == keys

    for k in indict.keys():
        try:
            print(keys)
            assert indict[k].keys() == keys
        except AssertionError as ae:
            ikeys = list(indict[k].keys())
        
            for test_key in ikeys:
                if test_key in keys:
                    pass
                else:
                    indict[k].pop(test_key, None)
            assert indict[k].keys() == keys

    out = indict[list(indict.keys())[0]]
    # now loop thru the keys and do some concatenation
    for _k in list(indict.keys())[1:]:
        out = dict_concatenate(out, indict[_k], keys)



    return out

def between_latlon(lat_arr, lon_arr, latmin, latmax, lonmin, lonmax):

    out = np.where( (lat_arr >= latmin) & (lat_arr < latmax) & (lon_arr >= lonmin) & \
            (lon_arr < lonmax) )
    return out

def between(arr1, arr2, min1, max1, min2, max2):

    out = np.where( (arr1 >= min1) & (arr1 < max1) & (arr2 >= min2) & (arr2 < max2) )
    return out

def between1(arr1, min1, max1):

    out = np.where( (arr1 >= min1) & (arr1 < max1) )
    return out

def surrounding(lats, lons, inlat_bounds, inlon_bounds, outlat_bounds, outlon_bounds):
    # want to return values that are within the outer bounds, but not the inner bounds
    # each of the bounds has to be a list of 2 values (min/max)
    #good = np.where( (lats >= outlat_bounds[0]) & (lats <= outlat_bounds[1]) & (lons >= outlon_bounds[0]) & \
    #		(lons <= outlon_bounds[1]) & (lats > inlat_bounds[1]) & (lats < inlat_bounds[0]) & \
    #		(lons > inlon_bounds[1]) & (lons < inlon_bounds[0]) )

    latgood = (lats >= outlat_bounds[0]) & (lats <= outlat_bounds[1])
    latbad =  (lats < inlat_bounds[1]) & (lats > inlat_bounds[0])

    lat_bool = np.bitwise_and( latgood, np.logical_not(latbad) )

    longood = (lons >= outlon_bounds[0]) & (lons <= outlon_bounds[1])
    lonbad = (lons < inlon_bounds[1]) & (lons > inlon_bounds[0])

    lon_bool = np.bitwise_and( longood, np.logical_not(lonbad) )
    
    good = np.bitwise_and(lat_bool, lon_bool)

    return good


def big_array(dict, cond = None):

    n_keys = len(dict.keys())
    out = np.zeros(0)
    if dict[dict.keys()[0]].ndim == 2:
        out = np.zeros((dict[dict.keys()[0]].shape[0],0))
    for i in range(n_keys):
        if dict.values()[i].any():
            if cond is not None:
                    out = np.concatenate((out, dict.values()[i][cond]))
            else:
                    out = np.concatenate((out, dict.values()[i]), axis = 1)

    return out


def big_ma_array(dict, cond = None): # for masked arrays

    n_keys = len(dict.keys())
    out = np.zeros(0)
    if dict[dict.keys()[0]].ndim == 2:
        out = np.zeros((dict[dict.keys()[0]].shape[0],0))
    for i in range(n_keys):
        if dict.values()[i].any():
            if cond is not None:
                    out = np.ma.concatenate((out, dict.values()[i][cond]))
            else:
                    out = np.ma.concatenate((out, dict.values()[i]), axis = 1)

    return out


def custom_hist (x, y, bins, min_samples = 5,r1=25,r2=75):
    hist, edges = np.histogram(x, bins = bins, density = True)
    hist = hist/np.sum(hist)
    y50 = np.zeros(edges.shape[0]-1)
    y75 = np.zeros(edges.shape[0]-1)
    y25 = np.zeros(edges.shape[0]-1)
    mad = np.zeros(edges.shape[0]-1)
    for h in range(edges[:-1].shape[0]):
        this_bin = np.where( (x >= edges[h]) & (x < edges[h+1]) )
        if this_bin[0].shape[0] >= min_samples:
            these_y = y[this_bin]
            y50[h] = np.median(these_y)
            y25[h] = np.percentile(these_y, r1)
            y75[h] = np.percentile(these_y, r2)
            mad[h] = np.median(np.abs(these_y-y50[h]))
        else:
            y50[h] = np.nan
            y25[h] = np.nan
            y75[h] = np.nan
            mad[h] = np.nan

    temp_mask = np.logical_not(np.isnan(y50))
    temp_mask = np.bitwise_and(temp_mask, y50 != 0)
    sl, int, r_val, p_val, std_err = stats.linregress(np.log(edges[:-1][temp_mask]), np.log(y50[temp_mask]))
    rho, p = stats.spearmanr(edges[:-1][temp_mask], y50[temp_mask])

    return hist, edges, y50, y25, y75, r_val, mad, rho, sl, int

def median_plot1d(xarr, yarr, xbins = np.arange(0,1,0.05), cmap = plt.cm.jet, cbar_label = '', \
        xlab = '', left_lab = '', right_lab = '', savename = '', min_pts = 5, \
        point_size = 250, e_width = 2, cdf = 0, regions = None, \
        colors = ['red','blue','DarkOrange','Green'], cdf_width = 2, xrange = None, \
        yrange = None, xscale = 'linear', yscale = 'symlog', yticks = None, verbose = 0, cond = None):

    "bin the x value and take the median y value corresponding to storms in each bin"

    if regions is not None:
        if len(regions) > 1: 
            X = big_array(xarr, cond = cond)
            Y = big_array(yarr, cond = cond)
        elif len(regions) == 1: 
            X = xarr[regions[0]]
            Y = yarr[regions[0]]
    else: # if not a dictionary, just a couple of arrays
        X = xarr
        Y = yarr


    # if more than 1 region, concatenate them together

    if verbose == 1:   
        print ('X: ', X )
        print ('Y: ', Y)

    hist, edges = np.histogram(X, bins = xbins, density = True)
    hist = 100.0*hist/np.sum(hist)
    out = np.zeros(hist.shape)
    out25 = np.zeros(hist.shape)
    out75 = np.zeros(hist.shape)
    for i in range(edges[:-1].shape[0]):
        this_bin = np.where( (X >= edges[i]) & (X < edges[i+1]) )
    #print this_bin
        if this_bin[0].shape[0] >= min_pts:
            these = Y[this_bin]
        #print 'these: ', these
            out[i] = np.median(these)
            out25[i] = np.percentile(these, 25)
            out75[i] = np.percentile(these, 75)

        else: out[i] = np.nan
#    bad = np.where(out25 == 0)
#    out25[bad] = 0.2

    fig = plt.figure()
    fig.subplots_adjust(bottom = 0.00, top = 0.93, left = 0.1, right = 0.9)
    sc = plt.scatter(edges[:-1], out, s = point_size, cmap = cmap, alpha = 0.9, c = hist, vmin = 0)
    err = plt.errorbar(edges[:-1], out, yerr = [out - out25, out75- out], ls = 'none', \
        color = 'k', elinewidth = e_width)

    good = np.logical_not(np.isnan(out))

    #print good

    sl, int, r_val, p_val, std_err = stats.linregress(edges[:-1][good], np.log(out[good]))
    rho, p = stats.spearmanr(edges[:-1][good], np.log(out[good]))

    cb = plt.colorbar(sc, orientation = 'horizontal', pad = 0.1)
    cb.set_label(cbar_label)
    cb.set_alpha(1)
    cb.draw_all()
    plt.axis('tight')
    plt.ylabel(left_lab)
    plt.xlabel(xlab)
    if yrange is not None:
    	plt.ylim(yrange[0], yrange[1])
    ax = plt.gca()
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if yticks is not None:
    	ax.set_yticks(yticks)
    	ax.set_yticklabels(yticks)
    ax.yaxis.set_major_formatter(plticker.FormatStrFormatter('%d'))

    plt.title(r'All Regions, R$^2$: %.2f; $\rho$: %.2f, N = %d'%(r_val**2, rho, X.shape[0]))

    ###### NOW DO RIGHT AXIS WITH CDF OF ENVIRONMENTAL VALUES

    if cdf == 1:
        ax2 = ax.twinx()
        ax2.set_ylim(0,1)
        ax2.set_xscale(xscale)
        if xrange is not None:
            ax2.set_xlim(xrange[0],xrange[1])
        for i, reg in enumerate(regions):

            x_temp, x_bins = xarr, xbins
            hist, edges = np.histogram(x_temp[reg], bins = x_bins, density = True)
            hist = np.cumsum(hist/np.sum(hist))
            ax2.plot(edges[:-1], hist, colors[i], label = reg, linewidth = cdf_width)

        ax2.legend(loc = 'upper left', prop={'size':10})
        ax2.set_ylabel(right_lab)

    else:
      ax2 = None

    if savename is not '':
        plt.savefig(savename)
        return None
    else:
        return fig, ax, ax2

def custom_hist2d(xarr, yarr, zarr, bins = [np.arange(0,1,0.05), np.arange(0,1,0.05)], min_pts = 5):


    X, Y, Z = xarr, yarr, zarr
    hist, xedges, yedges = np.histogram2d(X, Y, bins = bins)
    dx, dy = np.float(bins[0][1] - bins[0][0]), np.float(bins[1][1] - bins[1][0])
    hist = np.transpose(hist)
    f = np.zeros(hist.shape)
    for i in range(xedges[:-1].shape[0]):
        for j in range(yedges[:-1].shape[0]):
            this_bin = np.where( (X >= xedges[i]) & (X < xedges[i+1]) &
                            (Y >= yedges[j]) & (Y < yedges[j+1]) )
            if this_bin[0].shape[0] >= min_pts:
                these_f = Z[this_bin]
                f[j,i] = np.median(these_f)
            else: f[j,i] = np.nan


    f_m = np.ma.masked_array(f, np.isnan(f))
    return xedges, yedges, f_m, hist

def median_plot2d(xarr, yarr, zarr, bins = [np.arange(0,1,0.05), np.arange(0,1,0.05)], cmap = plt.cm.jet, \
        cbar_label = '', xlab = '', left_lab = '', right_lab = '', savename = '', min_pts = 5, \
            point_size = 250, e_width = 2, cdf = 0, regions = ['AL','DC','OK','CO'], \
            colors = ['red','blue','DarkOrange','Green'], cdf_width = 2, xrange = None, \
            yrange = None, text_color = 'white', verbose = 0):

    "bins both x and y variables, and plots a colormap of median z value for each 2d bin"

    ############## add in option for glowing text for numbers ################

    if len(regions) > 1:
        X = big_array(xarr)
        Y = big_array(yarr)
        Z = big_array(zarr)
    elif len(regions) == 1:
        X = xarr[regions[0]]
        Y = yarr[regions[0]]
        Z = zarr[regions[0]]

    # if more than 1 region, concatenate them together

    if verbose == 1:
        print ('X: ', X)
        print ('Y: ', Y)
        print ('Z: ', Z)

    hist, xedges, yedges = np.histogram2d(X, Y, bins = bins)
    dx, dy = np.float(bins[0][1] - bins[0][0]), np.float(bins[1][1] - bins[1][0])
    hist = np.transpose(hist)
    f = np.zeros(hist.shape)
    for i in range(xedges[:-1].shape[0]):
        for j in range(yedges[:-1].shape[0]):
            this_bin = np.where( (X >= xedges[i]) & (X < xedges[i+1]) &
                            (Y >= yedges[j]) & (Y < yedges[j+1]) )
            if this_bin[0].shape[0] >= min_pts:
                these_f = Z[this_bin]
                f[j,i] = np.median(these_f)
            else: f[j,i] = np.nan
    f_m = np.ma.masked_array(f, np.isnan(f))
    fig = plt.figure()
    pc = plt.pcolormesh(xedges, yedges, f_m, cmap = cmap)
    cb = plt.colorbar(pc)
    cb.set_label('Flash rate (/min)')
    # plotting white numbers on each box to show number of samples
    for i in range(xedges[:-1].shape[0]):
        for j in range(yedges[:-1].shape[0]):
            if hist[j,i] >= min_pts:
                plt.text(xedges[i]+dx/2., yedges[j]+dy/2., '%d'%hist[j,i], color = text_color,
                    horizontalalignment = 'center', fontsize = 9)
#    plt.xlabel('Normalized CAPE (J/kg*m)')
#    plt.ylabel('Cloud base height (m AGL)')
    plt.savefig(savename)


def stratify(x1, x2, x3, y, x1_vals = None, x2_vals = None):

#    x1_arg = np.argsort(x1)
#    x2_arg = np.argsort(x2)

#    x1_vals.insert(0,x1.min())
#    x2_vals.insert(0,x2.min())
#    x1_vals.append(x1.max())
#    x2_vals.append(x2.max())

    assert x1_vals is not None and x2_vals is not None

    s1 = x1.shape[0]
    s2 = x2.shape[0]
    s3 = x3.shape[0]

    print (x1_vals, x2_vals)

    x3_arrs = []
    y_arrs = []

    for i1 in range(len(x1_vals) - 1):
        for i2 in range(len(x2_vals) - 1):
            good = between(x1, x2, x1_vals[i1], x1_vals[i1+1], x2_vals[i2], x2_vals[i2+1])

            x3_arrs.append(x3[good])
            y_arrs.append(y[good])



    return x3_arrs, y_arrs



# BELOW IS WAY TO STRATIFY BY PERCENTILES, KEEPING IT AROUND JUST IN CASE

    # now have the stratifying points for x1, x2, go and grab x3
#    for i1 in range(len(x1_strat)-1):
#        for i2 in range(len(x2_strat)-1):
#            min_1 = np.round(x1_strat[i1]*s1/100.)
#            max_1 = np.round(x1_strat[i1+1]*s1/100.)-1
#
#            minval_1 = x1[x1_arg[min_1]]
#            maxval_1 = x1[x1_arg[max_1]]
#
#            min_2 = np.round(x2_strat[i2]*s2/100.)
#            max_2 = np.round(x2_strat[i2+1]*s2/100.)-1
#
#            minval_2 = x2[x2_arg[min_2]]
#            maxval_2 = x2[x2_arg[max_2]]
#
#            good = between(x1, x2, minval_1, maxval_1, minval_2, maxval_2)
#

def median_subplot(x, y, vals_1, vals_2, xbins, min_samples = 5, figname = '', \
        name1 = '', name2 = '', figsize = (6,6)):

    xs = len(x)
    v1s = len(vals_1)-1
    v2s = len(vals_2)-1

    fig, ax = plt.subplots(v1s, v2s, figsize = figsize)
    fig.subplots_adjust(top = 0.90, left = 0.05, right = 0.95, hspace = 0.4, bottom = 0.05)

    for i, a in enumerate(ax.flatten()):
        hist, edges, y50, y25, y75, r_val, mad, rho, sl, intercept = custom_hist(x[i], y[i], \
                xbins, min_samples = min_samples)

        a.scatter(edges[:-1], y50)
        a.axis('tight')
        a.set_title( '%s %.2f - %.2f\n%s %.2f - %.2f\nN = %d' % ( name1, vals_1[int(np.floor(i/v2s))], \
                vals_1[int(np.floor(i/v2s)+1)], name2, vals_2[np.mod(i,v2s)], vals_2[np.mod(i,v2s)+1], \
            x[i].shape[0] ), fontsize = 12 )

    if figname is not '':
        plt.savefig(figname)

def reflectivity_plot(x_arr, in_arr, condition, color = 'blue', lw = 3, label = ''):

    plt.plot(x_arr, np.array([np.average(in_arr[0,condition]), np.average(in_arr[1,condition]), \
    np.average(in_arr[2,condition]), np.average(in_arr[3,condition]), np.average(in_arr[4,condition]), \
    np.average(in_arr[5,condition])]), color, label = label, linewidth = lw)


#class Cells(region, ):

def mask_criteria(array, min_pct = 5., criteria = False):

    # loop thru array to see how many of the first 
    mask = np.ma.getmask(array)
    good = np.where(mask == criteria)

    out = np.ma.copy(array)

    for i in range(array.shape[0]):

        crit = np.where(good[0] == i)
        num_samples = crit[0].shape[0]
        pct = 100.0*float(num_samples)/(array.shape[1]) 

        #print i, num_samples, pct

    if pct < min_pct: # if not enough samples to count it, set the masks
        out[i,:] = np.ma.masked


    return out

#def glow(label_handle, color, width = 3): # make glowing text behind foreground text to make stand out

#    for tmp in label_handle:
#        tmp.set_path_effects([PE.withStroke(linewidth = width, foreground = color)])

def glow(ax, x, y, value, textcolor = 'black', width = 3, fontsize = 8, glow_color = 'white'):

    ax.text(x, y, value, color = textcolor, horizontalalignment = 'center', fontsize = 8, \
        path_effects = [PE.withStroke(linewidth = width, foreground = glow_color)])


def cdf(xarr, bins):

    hist, edges = np.histogram(xarr, bins = bins, density = True) # density makes it a float

    return np.cumsum(hist/np.sum(hist)), edges
#    hist = np.cumsum(hist/np.sum(hist))


def sloc (lat1,lon1,lat2,lon2):
    re = 6378.0
    c = np.pi/180.0
    lat1 = lat1*c
    lat2 = lat2*c
    lon1 = lon1*c
    lon2 = lon2*c

    d = re*np.arccos(np.sin(lat1)*np.sin(lat2)+np.cos(lat1)*np.cos(lat2)*np.cos(lon1-lon2))
    return np.round(d,decimals = 2)

def sph_ang (latref,lonref,latpoint,lonpoint):
    re = 6378.0
    c = np.pi/180.0
    latdiff = latpoint - latref
    ydiff = latdiff*(2*np.pi*re/360.0)
    londiff = lonpoint - lonref
    xdiff = londiff*(2*np.pi*re*np.cos(c*latref)/360.0)
   
    ang = np.arctan2(ydiff, xdiff)/c

    return np.round(ang,decimals = 2)

def copymask(oldarray, newarray):
    newarray = np.ma.masked_where(oldarray.mask, newarray)
    return newarray

def movingaverage(x,N, masked = False):
    mvav = np.convolve(x, np.ones((N,))/N)[(N-1):]
# transferring mask from original array to smoothed array 
    if masked: 
        #	mvav = np.ma.masked_where(x.mask, mvav)
        mvav = copymask(x, mvav)
        for pt in np.where(x.mask == True)[0]:
            if pt >= N: mvav.mask[pt-int(N):pt] = True
            else: mvav.mask[0:pt] = True    
        # also mask out points near the bad values of the raw 
        # array since they will be impacted
        # what about if close to the edge? might be wrapping around the edge

    return mvav

class BigArray(object):
    # take a dictionary with multiple keys and make it one big array with all data
    def __init__(self, dict):
        self.dict = dict
        self.array = self.big_array3d()
        self.keyarray = self.key_expand()


    def key_expand(self):
        out = np.zeros(self.array.shape[0], dtype = 'S8')
        ind = 0
        for key in self.dict.keys():
            a = self.dict[key].shape[0]
            out[ind:ind+a] = key
            ind += a

        return out

    def big_array3d(self):

        n_keys = len(self.dict.keys())
        if self.dict[self.dict.keys()[0]].ndim == 3:
            out = np.zeros((0,self.dict[self.dict.keys()[0]].shape[1], self.dict[self.dict.keys()[0]].shape[2]))
        elif self.dict[self.dict.keys()[0]].ndim == 2:
            out = np.zeros((0,self.dict[self.dict.keys()[0]].shape[1]))
        else: out = np.zeros(0)
        for i in range(n_keys):
            if self.dict.values()[i].any():
                out = np.concatenate((out, self.dict.values()[i]), axis = 0)

        return out

def bin_indices(vals, bins):
    if vals.ndim == 1:
        # this is just a reincarnation of digitize
        ovec = np.argsort(vals)
        ivec = np.searchsorted(vals, bins, sorter=ovec)
        out = []
        for b in range(bins.shape[0]-1):
            out.append(ovec[ivec[b] : ivec[b+1]])

    return out


def vert_cdf_barplot(vals, y, colors = None, ax = plt.gca()):

    for iy in range(vals.shape[0]):
        ax.barh(y[iy], vals[iy, 0], left = 0, height = 0.25, edgecolor = 'none', color = colors[0])
        for ix in range(1, vals.shape[1]):
                ax.barh(y[iy], vals[iy, ix], left = np.cumsum(vals, axis = 1)[iy, ix-1], \
                        color = colors[ix], height = 0.25, edgecolor = 'none')

    # no return value needed here I'm pretty sure, can add one later if needed


def write_dict_to_netcdf(filename, indict, dim_dict, attr_dict=None):

    # This will take in a dictionary and write it to a netcdf file
    # ***************************************************************
    # This will call write_nc_dict, which in turn calls write_nc_variable

    ncwid = netCDF4.Dataset(filename, 'w')

    # Then loop thru the dim_dict and make the appropriate dimensions

    for ddk in dim_dict.keys():
        ncwid.createDimension(ddk, dim_dict[ddk])

    #print dir(ncwid)
    #print ncwid.dimensions

    for ik in indict.keys():
        if 'type' in indict[ik].keys():
            var_type = indict[ik]['type']
        else:
            var_type = 'f'

        #print (ik, var_type, indict[ik].keys())

        long_name = indict[ik]['longname']
        units = indict[ik]['units']
        missing_val = indict[ik]['missing_val']
        write_nc_variable(ncwid, indict[ik]['data'], ik, indict[ik]['dims'], var_type,units=units,missing_val=missing_val,long_name=long_name)

    if attr_dict is not None:
        for ad in attr_dict.keys():
            ncwid.setncattr(ad, attr_dict[ad])


    ncwid.close()


def write_nc_variable(ncid, variable, name, dimension, vartype,units="",missing_val='-9999.',long_name = ''):
    #print 'writing %s variable to netcdf %s'%(name, variable)
    #if isinstance(variable, np.ndarray): print 'max value: %d'%variable.max()
    dummy = ncid.createVariable(name, vartype, dimension)
    print('np.sahpe dummy',np.shape(dummy))
    print('np.shape variable',np.shape(variable))
    print('name',name,dimension)
    dummy[:] = variable
    #print('attempting to write data')
    dummy.units = units
    dummy.missing_val = missing_val
    dummy.long_name = long_name

def delete_nc_dict(ncid, dict, vartype):
    for key in dict.keys():
        ncid.delncattr(key)

def into_dict(input_, variables): # this isnt working and may not work

    "takes multiple variables and puts them into a dictionary"

    out_dict = {}
    for v in variables:

        var_name = eval(v)
        out_dict['%s'%v] = np.array([input_[i].var_name for i in range(len(input))])


def write_nc_dict(dict, nc, dimension, vartype, prefix = '',units="",missing_val='-9999.',long_name = ''):

    "will loop thru dictionary keys to write nc variables using write_nc_variable function above"

    for key in dict.keys():

        write_nc_variable(nc, dict[key], '%s%s'%(prefix, key), dimension, vartype,units=units,missing_val=missing_val,long_name = long_name)







