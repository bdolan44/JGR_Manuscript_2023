# Originally written by Brenda Dolan
# Modified, starting in August 2018, by Brody Fuchs
# brfuchs@atmos.colostate.edu


import matplotlib
matplotlib.use('agg')
import h5py
import numpy as np
import matplotlib.pyplot as plt
import glob
import numpy as np
#import pycolors_new
import imp
import sys
import bottleneck as bn
#from gamma_dist import compute_gamma_dist
import oa_stats as OA 
import analysis_tools as AT

import dill as pickle
import os

import scipy.stats as stats
import scipy
import pandas
import dsd_helper_functions as dfunc
import warnings
import argparse
from matplotlib import colors
import Config
import steiner_houze_yuter_cs as shy
import RadarConfig


# Do I just wanna make this a script that's callable from the command line???
# that could be called by another driver script?

try:
    base_path = os.path.dirname(os.path.realpath(__file__))
except Exception as bpe:
    base_path = os.path.dirname(os.path.realpath('__file__'))
    print ('hit the exception: {}'.format(bpe))
    
print (base_path)

dbz_bins = np.arange(-10, 80, 5)
rr_bins = np.arange(0, 100, 2.5)
rr_log_bins = np.logspace(-3, 2, 15)


parser = argparse.ArgumentParser(description='namma-newtest-rain.h5')

parser.add_argument('--filename', action="store", dest="filename", type=str)
parser.add_argument('--exper', action="store", dest="exper", type=str)

pargs = parser.parse_args()

assert (pargs.filename is not None) and (len(glob.glob(pargs.filename)) == 1), 'Filename not found, check again'
assert (pargs.exper is not None), 'You need to specify --exper, so I know how to name the saved pickle file'
cfg = Config.Config('%s/simulation_settings.yaml'%(base_path))
ycfg = cfg.v

assert pargs.exper in ycfg['sims'].keys(), 'exper: {e} is not in the {y} config file, please check your spelling or make sure {e} is in {y}'.format(\
                e=pargs.exper, y='simulation_settings.yaml')

rcfg = RadarConfig.RadarConfig()


sim_cfg = ycfg['sims'][pargs.exper]
if 'd3'  in sim_cfg.keys():
    print('workign with 3D data!')
    d3 = True
else:
    d3 = False
filename = pargs.filename

if 'mu' in sim_cfg.keys():
    mu = sim_cfg['mu']
else:
    mu = 2

# Here we're reading in the h5 file using the h5py module
# will need to get familiar with this
rain_dsd_data = h5py.File(filename, 'r')
filtered_rain_dsd_dict = dfunc.format_subsample_dsd_data(rain_dsd_data, nu=mu, tskip=sim_cfg['tskip'], hskip=sim_cfg['hskip'],d3=d3)

out_filename = "{b}/pickles/{e}".format(b=base_path, e=sim_cfg['dsd_pickle'])

pickle.dump(filtered_rain_dsd_dict, open( out_filename, "wb" ) )
print ('saved {f}'.format(f=out_filename))


# these are the test parameters here for calculating dBZ
ttest = 30 # the time
ztest = 0 # the height index
mutest = 2 # the mu value to input

for iti in range(rain_dsd_data['rain_gam_lognw'].shape[0]): # loop thru the times

# let's calculate the reflectivity from the specified DSD
	dbz_test = dfunc.calc_dbz(10**rain_dsd_data['rain_gam_lognw'][iti, ztest, :, :], rain_dsd_data['rain_gam_d0'][iti, ztest, :, :], 
						np.zeros_like(rain_dsd_data['rain_gam_d0'][iti, ztest, :, :]) + mutest)
	dbz_test = np.ma.masked_invalid(dbz_test)



# now let's run convective stratiform partitioning on the data
# as if it were regular radar data

	Lon, Lat = np.meshgrid(rain_dsd_data['x_coords'][:], rain_dsd_data['y_coords'][:])

	yh_cs, yh_cc, yh_bkgnd = shy.conv_strat_latlon(dbz_test, Lat, Lon, CoreThresh=40.0, method='SYH', a=8, b=64, sm_rad=25)
	
	yh_cc_ma = np.ma.masked_where(yh_cc < 0, yh_cc)


	cs_arr = np.full(yh_cs.shape, np.nan)
	yh_conv = (yh_cs == 3) | (yh_cs == 4) | (yh_cs == 5)
	yh_mixed = (yh_cs == 1) | (yh_cs == 2)
	yh_strat = yh_cs == 0

	cs_arr[yh_conv] = 2
	cs_arr[yh_mixed] = 3
	cs_arr[yh_strat] = 1


# okay now let's make a figure
	fig, ax = plt.subplots(2, 2, figsize=(10, 10))

	axf = ax.flatten()

	dbzpc = axf[0].pcolormesh(rain_dsd_data['x_coords'][:], rain_dsd_data['y_coords'][:], dbz_test, vmin=rcfg.plot_params[rcfg.dz_name]['lims'][0], 
														vmax=rcfg.plot_params[rcfg.dz_name]['lims'][1], cmap=rcfg.plot_params[rcfg.dz_name]['cmap'])
	dbzcb = plt.colorbar(dbzpc, ax=axf[0], aspect=40)

	dbzspc = axf[1].pcolormesh(rain_dsd_data['x_coords'][:], rain_dsd_data['y_coords'][:], yh_bkgnd, vmin=rcfg.plot_params[rcfg.dz_name]['lims'][0], 
														vmax=rcfg.plot_params[rcfg.dz_name]['lims'][1], cmap=rcfg.plot_params[rcfg.dz_name]['cmap'])
	dbzscb = plt.colorbar(dbzspc, ax=axf[1], aspect=40)


	cspc = axf[2].pcolormesh(rain_dsd_data['x_coords'][:], rain_dsd_data['y_coords'][:], cs_arr, vmin=rcfg.plot_params['CS']['lims'][0], 
														vmax=rcfg.plot_params['CS']['lims'][1], cmap=rcfg.plot_params['CS']['cmap'])

	# cspc = axf[2].pcolormesh(rain_dsd_data['x_coords'][:], rain_dsd_data['y_coords'][:], yh_cs, vmin=rcfg.plot_params['CS']['lims'][0], 
	# 													vmax=rcfg.plot_params['CS']['lims'][1], cmap=plt.cm.jet)

	cscb = plt.colorbar(cspc, ax=axf[2], aspect=40)


	ccpc = axf[3].pcolormesh(rain_dsd_data['x_coords'][:], rain_dsd_data['y_coords'][:], yh_cc_ma, vmin=0, 
														vmax=6, cmap=plt.cm.jet)

	cccb = plt.colorbar(ccpc, ax=axf[3], aspect=40)


	plt.tight_layout()

	fig.suptitle('%s t=%03d'%(pargs.exper, iti))
	fig.subplots_adjust(top=0.93)

	plt.savefig('dummy_figs/%s_dbz_cs_%03d.png'%(pargs.exper, iti))

	plt.close(fig)


	# now add in some PDFs
	fig, ax = plt.subplots(1, 2, figsize=(11, 6))

	dbz_hist = np.histogram(dbz_test.compressed(), bins=dbz_bins)[0].astype(float)
	dbz_hist *= 100.0/dbz_hist.sum()

	rr_hist = np.histogram(rain_dsd_data['pcprr'][iti, :, :], bins=rr_log_bins)[0].astype(float)
	rr_hist *= 100.0/rr_hist.sum()



	ax[0].plot(dbz_bins[:-1], dbz_hist)
	ax[0].set_xlabel('dBZ')


	ax[1].plot(rr_log_bins[:-1], rr_hist)
	ax[1].set_xlabel('RR (mm/hr)')

	plt.tight_layout()
	
	fig.suptitle('%s t=%03d'%(pargs.exper, iti))
	fig.subplots_adjust(top=0.93)


	plt.savefig('dummy_figs/%s_dbz_rr_pdf_%03d.png'%(pargs.exper, iti))

	plt.close(fig)


