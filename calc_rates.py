####This is a program to read RAMS data, calculate the PC values given the 'RAMS world' input,
####Then plot the data and look at the contributing process rates.
### Brenda Dolan
### bdolan@atmos.colostate.edu
### June 2018


# Updates from Brody Fuchs, Aug 2018
# brfuchs@atmos.colostate.edu



import matplotlib
matplotlib.use('agg')

import xarray as xr
import pyart
pyart.config.get_field_colormap('CZ')
pycm = matplotlib.cm.get_cmap(pyart.config.get_field_colormap('CZ'))


import h5py
import numpy as np
import matplotlib.pyplot as plt

try:
    plt.style.use('presentation')
except Exception:
    pass

import glob
import numpy as np
import imp
import sys
#import bottleneck as bn
import oa_stats as OA 
import analysis_tools as AT

import Config

import dill as pickle

import scipy.stats as stats
import scipy
import pandas
import dsd_helper_functions as dfunc
from matplotlib import colors
import matplotlib as mpl
import pyart
import os
import argparse
import time
from matplotlib.colors import LogNorm


try:
    base_path = os.path.dirname(os.path.realpath(__file__))
except Exception as bpe:
    base_path = os.path.dirname(os.path.realpath('__file__'))
    print ('hit the exception: {}'.format(bpe))
    
print (base_path)


start = time.time()

parser = argparse.ArgumentParser(description='Put in a file to be processed')

parser.add_argument('--exper', action="store", dest="exper", type=str)

pargs = parser.parse_args()

cfg = Config.Config('%s/simulation_settings.yaml'%(base_path))
ycfg = cfg.v

assert pargs.exper in ycfg['sims'].keys(), 'exper: {e} is not in the {y} config file, please check your spelling or make sure {e} is in {y}'.format(\
				e=pargs.exper, y='simulation_settings.yaml')

sim_cfg = ycfg['sims'][pargs.exper]
var_cfg = ycfg['vars']
shortname = sim_cfg['shortname']
#pyart.config.get_field_colormap('CZ')
#pycm = matplotlib.cm.get_cmap(pyart.config.get_field_colormap('CZ'))
#print pycm


# These just look like a set of options for each different experiment
# **** These need to be in a centralized location and should just be one dictionary
# **** kinda like the RadarConfig.plot_params

 
w_bins = np.arange(-80, 85, 5)



plot_pca_flag = True
multi_panel = True

alldsd = pickle.load(open( "{b}/pickles/ramsworld_dsd.p".format(b=base_path), "rb" ) )

# so only one simulation is analyzed per run-thru of the script, and it's chosen right here

raindat = h5py.File('{b}/data/{f}'.format(b=base_path, f=sim_cfg['rain_file']), 'r')
microdat = h5py.File('{b}/data/{f}'.format(b=base_path, f=sim_cfg['micro_file']), 'r')
try:
    extradat = h5py.File('{b}/data/{f}'.format(b=base_path, f=sim_cfg['extra_file']), 'r')
    extra_flag = True
except Exception as ede:
    print ('Could not read in extradat file: {}'.format(ede))
    extra_flag = False

print ('after reading in files', time.time()-start)



if extra_flag:
        # do some analysis with the data in the extra file

    wdat = np.squeeze(extradat['w'][::sim_cfg['tskip'], :, ::sim_cfg['hskip'], ::sim_cfg['hskip']])
    
    # first, let's do a timeseries of a W PDF

    w_hist_ts = []

    for iti in range(wdat.shape[0]): # loop thru in time
        w_hist_tmp = np.histogram(wdat[iti].ravel(), bins=w_bins)[0].astype(float)
        w_hist_tmp *= 1.0/w_hist_tmp.sum()

        w_hist_ts.append(w_hist_tmp)

    w_hist_ts = np.array(w_hist_ts)




    # let's do PDF on the whole darn thing
    wdat_1d = wdat.ravel()
    w_hist = np.histogram(wdat_1d, bins=w_bins)[0].astype(float)
    w_hist *= 1.0/w_hist.sum()

    fig, ax = plt.subplots(2, 1, figsize=(8, 10))

    ax[0].plot(w_bins[:-1]+np.diff(w_bins)[0]/2.0, w_hist)

    ax[0].set_xlabel('W (m/s)')
    ax[0].set_ylabel('Frequency')

    ax[0].set_yscale('log')
    ax[0].grid()
    ax[0].set_title('{s} vertical motion PDF'.format(s=sim_cfg['shortname']))


    # now do a pcolormesh on the w pdf timeseries
    wts_pc = ax[1].pcolormesh(np.arange(w_hist_ts.shape[0]), w_bins[:-1], w_hist_ts.T, norm=LogNorm(vmin=1e-6, vmax=1), cmap=plt.cm.jet)

    wts_cb = plt.colorbar(wts_pc, ax=ax[1], orientation='horizontal', pad=0.1, fraction=0.1, aspect=40)

    ax[1].set_xlabel('Timestep')
    ax[1].set_ylabel('W (m/s)')


    plt.tight_layout()

    figname = '{b}/figures/{s}_w_pdf.png'.format(b=base_path, s=sim_cfg['shortname'])
    print ('saving %s'%figname)
    plt.savefig(figname, dpi=200)




#quit()

####Get the parameters back from the RAMS data. Here the data can be decimated in time or space
####Also determine which method to integrate the rates with int_type = 'sum' or int_type = 'max'
params2d = dfunc.get_rams_params_2d(raindat, microdat, tskip=sim_cfg['tskip'], hskip=sim_cfg['hskip'])

good_vars = var_cfg.keys()



pc1, pc2 = OA.get_eof_cast(params2d, alldsd, good_vars)
pc1 = np.reshape(pc1, np.shape(params2d['d00']))
pc2 = np.reshape(pc2, np.shape(params2d['d00']))

pc1 = np.ma.masked_where(params2d['d00']<=0, pc1)
pc2 = np.ma.masked_where(params2d['d00']<=0, pc2)

params2d['d00'] = np.ma.masked_where(params2d['d00'] <= 0, params2d['d00'])
params2d['nww'] = np.ma.masked_where(params2d['d00'] <= 0, params2d['nww'])

thresh = 1.5

ramsgroups = dfunc.get_groups(pc1, pc2, thresh, thresh, thresh, thresh)

groupvals = dfunc.make_groups(pc1, ramsgroups)


params2d['d00'] = np.ma.masked_where(params2d['d00'] <= 0, params2d['d00'])
params2d['nww'] = np.ma.masked_where(params2d['d00'] <= 0, params2d['nww'])
params2d['sigm'] = np.ma.masked_where(params2d['d00'] <= 0, params2d['sigm'])
params2d['rrr'] = np.ma.masked_where(params2d['d00'] <= 0, params2d['rrr'])



groupvals = np.ma.masked_where(params2d['d00']<= 0, groupvals)


##### Plot the recast data ################

########in D0-Nw space ######################

# ***** These plots don't need to be in here I'm pretty sure
# ***** I'd put these in the rams_world.ipynb

print ('plotting PCA results', time.time()-start)

if plot_pca_flag == True:

    fig, ax = plt.subplots(1, 1, figsize=(12,8))
    dfunc.plot_groups(params2d['d00'], params2d['nww'], ramsgroups, ax)
    plt.xlim(0, 4)
    plt.ylim(0, 8)
    plt.xlabel('D$_0$ (mm)')
    plt.ylabel('LogN$_w$')
    plt.title('{l} Global Groups'.format(l=sim_cfg['longname']))

    plt.tight_layout()
    figname = '{b}/figures/{s}_NwD0_globalgroups.png'.format(b=base_path, s=sim_cfg['shortname'])
    print( 'saving %s'%figname)
    plt.savefig(figname, dpi=200)

########in PC1PC2 space ######################
if plot_pca_flag == True:

    fig, ax = plt.subplots(1,1,figsize=(8,8))
    t = dfunc.get_2d_pchist(np.ravel(pc1),np.ravel(pc2),'{l}'.format(l=sim_cfg['longname']),ax)

    figname = '{b}/figures/{ex:}_PC1PC2_globalcast.png'.format(b=base_path, ex=sim_cfg['shortname'])

    print( 'saving %s'%figname)
    plt.savefig(figname, dpi=150)
    plt.clf()


#############################################################################
###  Set up an array of the microphysical processes so we can rank them
############################################################################

# BF: I'm actually not sure if this is needed at all...
print ('microphysics', time.time()-start)


microstack = params2d['agg'][np.newaxis,:]

microstack = np.vstack([microstack,params2d['c2r'][np.newaxis,:]])
microstack = np.vstack([microstack,params2d['i2r'][np.newaxis,:]])
microstack = np.vstack([microstack,params2d['melt'][np.newaxis,:]])
microstack = np.vstack([microstack,params2d['r2i'][np.newaxis,:]])
microstack = np.vstack([microstack,params2d['rime'][np.newaxis,:]])
microstack = np.vstack([microstack,params2d['vapi'][np.newaxis,:]])
microstack = np.vstack([microstack,params2d['vapl'][np.newaxis,:]])


hand = ['agg','c2r','i2r','melti','r2i','rime','vapi','vapl']

###Take the maximum value along the first dimension to give the highest rate in the colum
maxrate = np.argmax(microstack, axis=0)

microrank = np.zeros((8,len(np.ravel(params2d['d00']))))
microrank1 = np.zeros((len(np.ravel(params2d['d00']))))
microrank2 = np.zeros((len(np.ravel(params2d['d00']))))
microrank3 = np.zeros((len(np.ravel(params2d['d00']))))

dat = np.reshape(microstack,np.shape(microrank))


for i, a in enumerate(np.ravel(params2d['d00'])):
    dat1 = dat[:,i]
    rank = scipy.stats.rankdata(np.ravel(dat1),method='min')
    microrank[:,i] = rank
    try:
        microrank1[i] = np.argwhere(np.floor(rank)==1)
    except ValueError as ve:
        if len(np.argwhere(np.floor(rank)==1)) > 0:
            microrank1[i] = np.argwhere(np.floor(rank)==1)[0]

    try:
        microrank2[i] = np.argwhere(np.floor(rank)==2)
    except:
        if len(np.argwhere(np.floor(rank)==2)) > 0:
            microrank2[i] = np.argwhere(np.floor(rank)==2)[0]
    
    try:
        microrank3[i] = np.argwhere(np.floor(rank)==3)
    except:
        if len(np.argwhere(np.floor(rank)==3)) > 0:
            microrank3[i] = np.argwhere(np.floor(rank)==3)[0]

print ('microphysics ranks', time.time()-start)


shp3d = np.shape(params2d['d00'])

microrank = np.reshape(microrank,(8,shp3d[0],shp3d[1],shp3d[2]))
microrank1 = np.reshape(microrank1,(shp3d[0],shp3d[1],shp3d[2]))
microrank2 = np.reshape(microrank2,(shp3d[0],shp3d[1],shp3d[2]))
microrank3 = np.reshape(microrank3,(shp3d[0],shp3d[1],shp3d[2]))

microrank1 = np.ma.masked_where(params2d['d00']<=0,microrank1)
microrank2 = np.ma.masked_where(params2d['d00']<=0,microrank2)
microrank3 = np.ma.masked_where(params2d['d00']<=0,microrank3)

whbad = np.less_equal(params2d['d00'],0)

microrank1[whbad] =-1
microrank2[whbad] =-1
microrank3[whbad] =-1

rate_names = ['','agg','c2r','i2r','melti','r2i','rime','vapi','vapl']
rate_colors = [ 'white','pink','lightblue', 'limegreen', 'purple', 'blue',
              'red', 'gold','darkorange']

cmaprates = colors.ListedColormap(rate_colors)

# So these groups and colors and process rates and colors are gonna need to be in 
# a centralized location

hid_names = ['Amb','Grp1', 'Grp2', 'Grp3', 'Grp4',
              'Grp5', 'Grp6']
hid_colors = [ 'grey','Red', 'Green', 'GOld', 'Blue',
              'Orange', 'Purple']

cmaphid = colors.ListedColormap(hid_colors)

print ('plotting 4 panel groups', time.time()-start)



############################################################
if plot_pca_flag == True:

    fig, ax = plt.subplots(2,2,figsize=(12,12))
    axf = ax.flatten()


    dfunc.plot_groups(params2d['d00'], params2d['nww'], ramsgroups, axf[0])

    dfunc.contour_density(params2d['d00'], params2d['nww'], 'TOR2',axf[1])

    dfunc.plot_frequency(groupvals, axf[2], thresh)

    hist,eg = np.histogram(np.ravel(microrank1),bins=np.arange(-1,9,1))
    for i,h in enumerate(np.arange(0,8,1)):
        #print i
        dum = axf[3].bar(h, hist[i+1]/np.float(np.sum(hist[1:]))*100., width=0.5, color=rate_colors[i+1])
    axf[3].set_xlim(-0.5,7.5)
    axf[3].set_xticks(np.arange(0, 8, 1))
    axf[3].set_xticklabels(rate_names[1:],rotation=15)
    axf[3].set_title('RANK 1 Frequency')




    plt.tight_layout()
    plt.savefig('{b}/figures/{e}_PC{t}_4panel_groups.png'.format(b=base_path, t=thresh, e=shortname))

    plt.clf()

print ('done w/ 4 panel groups', time.time()-start)



########################################################################################################################
#####Now make a 6 panel figure with D0, Nw, Groups, and the highest 3 ranked rates in a column   ######################
########################################################################################################################


# Put this stuff in a centralized location as well!!

radarcbar = ['PeachPuff','Aqua','DodgerBlue','MediumBlue','Lime', \
    'LimeGreen','Green','Yellow','Orange','OrangeRed', \
    'Red', 'Crimson','Fuchsia','Indigo','DarkCyan','White']
temp_cmap = colors.ListedColormap(radarcbar)

rr_colors = ['#FFFFFF', '#ccd1d1', '#99a3a4', '#707b7c', '#f4ec3f', '#f3bd4a', '#f39c12', '#ef776a', '#C85046',
            '#9B2422', '#600000']
rr_cmap = mpl.colors.ListedColormap(rr_colors)
rr_bounds = np.array([0, 2, 3, 5, 10, 20, 30, 50, 75, 150, 225, 300])
rr_norm = mpl.colors.BoundaryNorm(rr_bounds, rr_cmap.N)


new_stern_cmap = dfunc.modify_cmap("gist_stern", {0: [0.3, 0.0, 0.0, 1.0], 1: [0.9, 0.0, 0.0, 1.0]}, 
            ncolors=8, new_name='newcmap')


groupvals = np.ma.masked_where(params2d['d00']<=0,groupvals)

if multi_panel == True:
    print ('Calculating 6panel rates....')
    #for t in range(20,21):

    for t in range(0,len(groupvals[:,0,0])):

        fig, ax = plt.subplots(3,4,figsize=(18,10))
        axf= ax.flatten()

    #t =20
        cb0 = axf[0].pcolor(params2d['d00'][t,:,:],vmin=0,vmax=3,cmap=new_stern_cmap)
        cb00 = plt.colorbar(cb0,ax=axf[0])
        axf[0].set_title('D$_0$ (mm)')

        cb1 = axf[1].pcolor(params2d['nww'][t,:,:],vmin=1,vmax=6,cmap=plt.cm.Spectral_r)
        cb11 = plt.colorbar(cb1,ax=axf[1])
        axf[1].set_title('Log(N$_w$)')


        cb = axf[2].pcolor(groupvals[t,:,:],vmin=0,vmax=7,cmap=cmaphid)
        cb22 = plt.colorbar(cb,ax=axf[2])
        cb22.ax.set_yticklabels(hid_names)
        axf[2].set_title('PCA Groups')

        ramssub = dfunc.get_groups(pc1[t,:,:],pc2[t,:,:],thresh,thresh,thresh,thresh)
        groupsub = dfunc.make_groups(pc1[t,:,:],ramssub)
        dfunc.plot_groups(params2d['d00'][t,:,:],params2d['nww'][t,:,:],ramssub,axf[3])
        axf[3].set_xlim(0,3)
        axf[3].set_ylim(0,7)
        axf[3].set_xlabel('D$_0$ (mm)')
        axf[3].set_ylabel('LogN$_w$')
    

        cb = axf[4].pcolor(params2d['dbz'][t,:,:],vmin=0,vmax=65,cmap=temp_cmap)
        cb23 = plt.colorbar(cb,ax=axf[4])
        #cb23.ax.set_yticklabels(hid_names)
        axf[4].set_title('Rain Reflectivity')

        cb = axf[5].pcolor(params2d['rrr'][t,:,:],cmap=rr_cmap,norm = rr_norm)
        cb24 = plt.colorbar(cb,ax=axf[5])
        axf[5].set_title('Rain Rate')
   
        cb = axf[6].pcolor(params2d['sigm'][t,:,:],vmin=0,vmax=3,
                       cmap=pycm)
        cb24 = plt.colorbar(cb,ax=axf[6])
        axf[6].set_title('Sigma_M')

    

        dfunc.plot_frequency(groupsub,axf[7],thresh)


    
        cb3 = axf[8].pcolor(microrank1[t,:,:],vmin=-1,vmax=8,cmap=cmaprates)
        #cb33 = plt.colorbar(cb3,ax=axf[3])
        #cb33.ax.set_yticklabels(rate_names)
        axf[8].set_title('ColRateSum 1')
        cb3= plt.colorbar(cb3,ax=axf[8])



        cb4 = axf[9].pcolor(microrank2[t,:,:],vmin=0,vmax=8,cmap=cmaprates)
        cb44 = plt.colorbar(cb4,ax=axf[9])

        #cb44 = plt.colorbar(cb4,ax=axf[4])
        #cb44.ax.set_yticklabels(rate_names)
        axf[9].set_title('ColRateSum 2')

        np.ma.masked_equal(microstack[:,t,:,:],-1)
        mssum =np.ma.sum(microstack[:,t,:,:],axis=2)
        print (np.shape(mssum))

        msum = np.ma.sum(mssum,axis=1)
        print (np.shape(msum))

        eg = np.arange(0,9,1)
        cb5 = axf[10].bar(eg[:-1],msum/np.ma.sum(msum)*100.,color=rate_colors[1:])
        axf[10].set_xlim(-0.5,7.5)
        axf[10].set_xticks(np.arange(0,8,1))
        axf[10].set_xticklabels(rate_names[1:],rotation=25)


        #cb44 = plt.colorbar(cb4,ax=axf[4])
        #cb44.ax.set_yticklabels(rate_names)
        axf[10].set_title('Total rates')

   
    
    #     cb5 = axf[10].pcolor(microrank3[t,:,:],vmin=0,vmax=8,cmap=cmaprates)
    #     cb55 = plt.colorbar(cb5,ax=axf[10])
    # #    cb55.ax.set_yticklabels(rate_names)
    #     axf[10].set_title('ColRateSum 3')
    #     cax = fig.add_axes([0.1, -0.01, 0.8, 0.02])
    #     fg = fig.colorbar(cb5,  cax=cax, orientation='horizontal',use_gridspec=True)

        #cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
    #     cax.set_xticks(np.arange(0, 11, 1))
    #     cax.set_xticklabels(rate_names,y=-0.01,horizontalalignment='center')  # vertically oriented colorbar
    #    fig.set_yticklabels(rate_names)

        hist,eg = np.histogram(microrank1[t,:,:],bins=np.arange(-1,9,1))
        for i,h in enumerate(np.arange(0,8,1)):
            print (i)
            dum = axf[11].bar(h,hist[i+1]/np.float(np.sum(hist[1:]))*100.,width=0.5,color=rate_colors[i+1])
        axf[11].set_xlim(-0.5,7.5)
        axf[11].set_xticks(np.arange(0,8,1))
        axf[11].set_xticklabels(rate_names[1:],rotation=25)
        axf[11].set_title('RANK 1 Frequency')
        #cb3= plt.colorbar(dum,ax=axf[5])



        plt.suptitle('time step: {t:02d}'.format(t=t),y=1.01,fontsize=25)

    #    plt.savefig('TOR2_Groups_surf_{t:02d}.png'.format(t=t),dpi=200,bbox_inches='tight')
        plt.tight_layout()
        plt.savefig('{b}/figures/{s}_Groups_surf_12panel_{t:02d}.png'.format(b=base_path, s=shortname, t=t), 
                            dpi=200, bbox_inches='tight')

        plt.clf()



########################################################################################################################



########################################################################################################################
#####Now plot histograms of rates for each column   ######################
fig, ax = plt.subplots(3,2,figsize=(12,8))
axf = ax.flatten()
gcolor={'Group1':'red','Group2':'green','Group3':'goldenrod','Group4':'blue','Group5':'orange','Group6':'purple'}
gcolor2={'Group1':'gray','Group2':'gray','Group3':'gray','Group4':'gray','Group5':'gray','Group6':'gray'}
gcolor3={'Group1':'black','Group2':'black','Group3':'black','Group4':'black','Group5':'black','Group6':'black'}


for i,g in enumerate(gcolor.keys()):
    whgr1 = ramsgroups[g]
#    print np.shape(whgr1)
    hist,eg = np.histogram(microrank1[whgr1],bins=np.arange(0,9,1))
    dum = axf[i].bar(eg[:-1]-0.2,hist/np.float(np.sum(hist))*100.,width=0.2,color=gcolor[g],label='Rank1')

    hist2,eg = np.histogram(microrank2[whgr1],bins=np.arange(0,9,1))
    dum = axf[i].bar(eg[:-1],hist2/np.float(np.sum(hist2))*100.,width=0.2,color=gcolor2[g],label='Rank2')

    hist3,eg = np.histogram(microrank3[whgr1],bins=np.arange(0,9,1))
    dum = axf[i].bar(eg[:-1]+0.2,hist3/np.float(np.sum(hist3))*100.,width=0.2,color=gcolor3[g],label='Rank3')

    axf[i].set_xticks(np.arange(0,9,1))
    axf[i].set_xticklabels(hand,rotation=20)
    axf[i].set_title('{g}'.format(g=g),fontsize=16)
    axf[i].legend(loc='best')
plt.tight_layout()
plt.suptitle('{s}Rank SUMS'.format(s=longname),fontsize = 30,y=1.02)
plt.savefig('{b}/figures/{s}_ranksums_process.png'.format(b=base_path, s=shortname), dpi=200, bbox_inches='tight')

########################################################################################################################
#####Look only at large D0 ######################
########################################################################################################################



fig, ax = plt.subplots(3,2,figsize=(12,8))
axf = ax.flatten()
# again, this needs to be somewhere else
gcolor={'Group1':'red','Group2':'green','Group3':'goldenrod','Group4':'blue','Group5':'orange','Group6':'purple'}
for i,g in enumerate(gcolor.keys()):
    whgr1 = ramsgroups[g]
    subdat = np.greater(params2d['d00'][whgr1],2.5,where=True)
    msubdat = microrank1[whgr1]
    msubdat2 = microrank2[whgr1]
    msubdat3 = microrank3[whgr1]
#    whgrf = np.where(subdat>2.5) 

    hist,eg = np.histogram(msubdat[subdat],bins=np.arange(0,9,1))
    dum = axf[i].bar(eg[:-1]-0.2,hist/np.float(np.sum(hist))*100.,width=0.2,color=gcolor[g],label='Rank1')

    hist2,eg = np.histogram(msubdat2[subdat],bins=np.arange(0,9,1))
    dum = axf[i].bar(eg[:-1],hist2/np.float(np.sum(hist2))*100.,width=0.2,color=gcolor2[g],label='Rank2')

    hist3,eg = np.histogram(msubdat3[subdat],bins=np.arange(0,9,1))
    dum = axf[i].bar(eg[:-1]+0.2,hist3/np.float(np.sum(hist3))*100.,width=0.2,color=gcolor3[g],label='Rank3')

    axf[i].set_xticks(np.arange(0,9,1))
    axf[i].set_xticklabels(hand,rotation=20)
    axf[i].set_title('{g}'.format(g=g),fontsize=16)
    axf[i].legend(loc='best')

#     except:
#         print 'no data'
plt.tight_layout()
plt.suptitle('{s} Rank SUMs D0>2.5'.format(s=shortname,fontsize = 30,y=1.02))
plt.savefig('{b}/figures/{s}_rank2sum_larged0_process.png'.format(b=base_path, s=shortname), dpi=200, bbox_inches='tight')

########################################################################################################################

########################################################################################################################
#####Look only at large Nw ######################
fig, ax = plt.subplots(3,2,figsize=(12,8))
axf = ax.flatten()
gcolor={'Group1':'red','Group2':'green','Group3':'goldenrod','Group4':'blue','Group5':'orange','Group6':'purple'}
for i,g in enumerate(gcolor.keys()):
    whgr1 = ramsgroups[g]
#    subdat = np.greater(params2d['d00'][whgr1],2.5,where=True)
    msubdat = microrank1[whgr1]
    msubdat2 = microrank2[whgr1]
    msubdat3 = microrank3[whgr1]
    whgrf = np.where(np.logical_and(params2d['nww'][whgr1]>5,params2d['d00'][whgr1]<1.0)) 

    hist,eg = np.histogram(msubdat[whgrf],bins=np.arange(0,9,1))
    dum = axf[i].bar(eg[:-1]-0.2,hist/np.float(np.sum(hist))*100.,width=0.2,color=gcolor[g],label='Rank1')

    hist2,eg = np.histogram(msubdat2[whgrf],bins=np.arange(0,9,1))
    dum = axf[i].bar(eg[:-1],hist2/np.float(np.sum(hist2))*100.,width=0.2,color=gcolor2[g],label='Rank2')

    hist3,eg = np.histogram(msubdat3[whgrf],bins=np.arange(0,9,1))
    dum = axf[i].bar(eg[:-1]+0.2,hist3/np.float(np.sum(hist3))*100.,width=0.2,color=gcolor3[g],label='Rank3')

    axf[i].set_xticks(np.arange(0,9,1))
    axf[i].set_xticklabels(hand,rotation=20)
    axf[i].set_title('{g}'.format(g=g),fontsize=16)
    axf[i].legend(loc='best')

#     except:
#         print 'no data'
plt.tight_layout()
plt.suptitle('{s} Rank SUMS nw>5,d0<1.0'.format(s=shortname),fontsize = 30,y=1.02)
plt.savefig('{b}/figures/{s}_ranksums_largenw_process.png'.format(b=base_path, s=shortname), dpi=200, bbox_inches='tight')
########################################################################################################################

########################################################################################################################
#####Look only at large Nw ######################
fig, ax = plt.subplots(3, 2, figsize=(12, 8))
axf = ax.flatten()
gcolor={'Group1':'red','Group2':'green','Group3':'goldenrod','Group4':'blue','Group5':'orange','Group6':'purple'}
for i,g in enumerate(gcolor.keys()):
    whgr1 = ramsgroups[g]
#    subdat = np.greater(params2d['d00'][whgr1],2.5,where=True)
    msubdat = microrank1[whgr1]
    msubdat2 = microrank2[whgr1]
    msubdat3 = microrank3[whgr1]
    whgr = np.logical_and(params2d['d00'][whgr1]>0.8,params2d['d00'][whgr1]<1.2) 
    whgr2 = np.logical_and(whgr,params2d['nww'][whgr1]>4.2)
    whgrf = np.where(whgr2)


    hist,eg = np.histogram(msubdat[whgrf],bins=np.arange(0,9,1))
    dum = axf[i].bar(eg[:-1]-0.2,hist/np.float(np.sum(hist))*100.,width=0.2,color=gcolor[g],label='Rank1')

    hist2,eg = np.histogram(msubdat2[whgrf],bins=np.arange(0,9,1))
    dum = axf[i].bar(eg[:-1],hist2/np.float(np.sum(hist2))*100.,width=0.2,color=gcolor2[g],label='Rank2')

    hist3,eg = np.histogram(msubdat3[whgrf],bins=np.arange(0,9,1))
    dum = axf[i].bar(eg[:-1]+0.2,hist3/np.float(np.sum(hist3))*100.,width=0.2,color=gcolor3[g],label='Rank3')

    axf[i].set_xticks(np.arange(0,9,1))
    axf[i].set_xticklabels(hand,rotation=20)
    axf[i].set_title('{g}'.format(g=g),fontsize=16)
    axf[i].legend(loc='best')

#     except:
#         print 'no data'
plt.tight_layout()
plt.suptitle('{s} Rank SUMs 0.8<D0<1.2, Nw>4.2'.format(s=shortname),fontsize = 30,y=1.02)
plt.savefig('{b}/figures/{s}_ranksum_horn_process.png'.format(b=base_path, s=shortname), dpi=200, bbox_inches='tight')

########################################################################################################################

fig, ax = plt.subplots(3,1,figsize=(14,10))
axf = ax.flatten()


print (np.shape(groupvals))

for i,v in enumerate(groupvals[:,0,0]):
#for i in [10,12]:

    p= plot_time_frequency(groupvals[i,:,:],axf[0],i,shortname)
 

    r = plot_time_frequency_rates(microstack[:,i,:,:],axf[1],i,shortname)
    rq =plot_time_frequency_rates(microrank1[i,:,:],axf[2],i,shortname,rank=1)

# axf[0].legend(p,hid_names)
# axf[1].legend(r,rate_names[1:])
# axf[2].legend(rq,rate_names[1:])

axf[0].set_xticks(np.arange(0,35,5))
axf[1].set_xticks(np.arange(0,35,5))
axf[2].set_xticks(np.arange(0,35,5))



axf[1].set_title('Total Rates')
axf[2].set_title('Rate Ranks1')



plt.tight_layout()

plt.savefig('{b}/figures/{s}_timeheight_freq.png'.format(b=base_path, s=shortname), dpi=200, bbox_inches='tight')
