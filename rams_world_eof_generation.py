# encoding: utf-8
######This is a program that reads in pickle files of RAMS DSD data in order to create
#####a global EOF.

#### B. Dolan
#### August 2018
#### Colorado STate University


## Modified by B. Fuchs, October 2018
## brfuchs@atmos.colostate.edu

import matplotlib
matplotlib.use('agg')

import h5py
import numpy as np
import matplotlib.pyplot as plt
try:
    plt.style.use('presentation')
except IOError as ImportError:
    pass

import glob
import numpy as np
#import pycolors_new
import imp
import sys
import bottleneck as bn
#from gamma_dist import compute_gamma_dist
import oa_stats as OA 
import analysis_tools as AT

import pickle as pickle
import io
import scipy.stats as stats
import scipy
import pandas
import dsd_helper_functions as dfunc
import os
import Config
from collections import OrderedDict
import analysis_tools as atools
#sys.setdefaultencoding('utf8')

def generate_random_color(redlims=[0, 255], greenlims=[0, 255], bluelims=[0, 255]):

    r = np.random.randint(*redlims)
    g = np.random.randint(*greenlims)
    b = np.random.randint(*bluelims)
    return '#%02x%02x%02x' % (r, g, b)




try:
    base_path = os.path.dirname(os.path.realpath(__file__))
except Exception as bpe:
    base_path = os.path.dirname(os.path.realpath('__file__'))
    print ('hit the exception: {}'.format(bpe))


cfg = Config.Config('simulation_settings_baseline3.yaml')
ycfg = cfg.v


# set up experiment dictionary to hold all the data
exper_dsd_data = OrderedDict()


expers = sorted(ycfg['sims'].keys())
n_expers = len(expers)


print ('experiments being included in the RAMS world PCA generation: {}'.format(expers))

# now load in the data
for ex in expers:
    print('ex', ex)
    with io.open('{bp}/pickles/{p}'.format(bp=base_path, p=ycfg['sims'][ex]['dsd_pickle']),'rb') as f:
#	heading_times, headings = np.load(ship_heading_files)
        print('pickle read:',f)
        exper_dsd_data[ex] =  pickle.load(f, encoding='latin1')


all_exper_data = atools.dict_concat_auto(exper_dsd_data)

all_exper_data = dfunc.remove_bad(all_exper_data)
all_exper_data = dfunc.remove_bad(all_exper_data,rainvar='lnt')


# First, let's just do a survey of the data included here

fig, ax = plt.subplots(1, 1, figsize=(8,8))

labels = exper_dsd_data.keys()
sizes = [len(_['lwcc']) for _ in exper_dsd_data.values()]
#'#%02x%02x%02x' % (0, 128, 64)

colors = [ycfg['sims'][_]['color'] if 'color' in ycfg['sims'][_].keys() else generate_random_color(redlims=[100,255], \
                        greenlims=[100,255], bluelims=[100,255]) for _ in exper_dsd_data.keys()]

#colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue']
#explode = (0.1, 0, 0, 0)  # explode 1st slice
 
# Plot
#plt.pie(sizes, explode=explode, labels=labels, colors=colors,
#        autopct='%1.1f%%', shadow=True, startangle=140)

wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%', 
                shadow=False, startangle=60, textprops={'fontsize': 8}, labeldistance=1.02, colors=colors)

for text in texts:
    text.set_color('black')
    text.set_size(8)

for w in wedges:
    w.set_linewidth(1.2)
    w.set_edgecolor('white')


ax.axis('equal')

plt.tight_layout()
fig.subplots_adjust(left=0.08, right=0.92, top=0.92, bottom=0.07)
fig.suptitle('Fraction of %d observations\nfrom %d RAMS simulations'%(len(all_exper_data['lwcc']), n_expers), fontsize=20)

plt.savefig('{bp}/figures/exper_fractions_{ne}.png'.format(bp=base_path, ne=n_expers), dpi=150)


##################################################################################################
##################################################################################################



# ### Go through and synthesizes the dictionaries from all projects into a single large one (data), but save the 
# ### sizes and locations of various datasets within the dataset.
# loca={}
# loca['TOR1']=0
# siz = {}
# siz['TOR1']=tor1sub['lrr'].shape[0]
# dat1,loca,siz = gather_exper(tor1sub,tor2sub,loca,siz,'TOR1','TOR2')
# dat2,loca,siz = gather_exper(dat1,tor4sub,loca,siz,'TOR2','TOR4')
# dat3,loca,siz = gather_exper(dat2,mc3esub,loca,siz,'TOR4','MC3E')
# dat4,loca,siz = gather_exper(dat3,nammabg,loca,siz,'MC3E','NAMMABG')
# data,loca,siz = gather_exper(dat4, nammaccn,loca,siz,'NAMMABG','NAMMACCN')

# for v in siz.keys():
#     print v, siz[v],loca[v]
    
# Instead of doing the chain calls of gather_exper, let's just use analysis_tools.concatenate_dict
# (or whatever it's called) to make them all into one!




good_vars = ['nww', 'dmm', 'sigm', 'llwc', 'lrr', 'lnt']
labels = [ycfg['vars'][_]['label'] for _ in good_vars]


# ### Now for the actual processing. call exper_eof from the helper functions, which calls the PCA from OA stats.
# ### Define the variables here.
# good_vars = ['nww', 'd00','sigm', 'llwc','lrr','lnt']#, 'dmaxx']#, 'sigdd']
# labels = ['logN$_w$', 'D$_m$','$\sigma$$_m$', 'logLWC', 'logRR','logNtt']#'D$_{max}$',]#, '$\sigma$']



alldsd = dfunc.exper_eof(all_exper_data, good_vars, labels)
### if the EOFS need to be flipped, do that here
#alldsd.flip_eof_and_pc(rank=1)
#alldsd.flip_eof_and_pc(rank=2)


#Save to pickle file for further use.

pickle.dump(alldsd, open( "{bp}/pickles/ramsworld3_dsd_dm.p".format(bp=base_path), "wb" ) )
print ('New RAMS world dsd pickled')


#### Now compare to the disdrometer world EOFs.

with io.open('{bp}/pickles/world_class.p'.format(bp=base_path),'rb') as f:
#	heading_times, headings = np.load(ship_heading_files)
    #print('pickle read:',f)
    world =  pickle.load(f, encoding='latin1')
# world = pickle.load(open('{bp}/pickles/world_class.p'.format(bp=base_path),'rb'))
    print ('loaded in the observation DSD pickle')



## This is where the EOFs are plotted.
pcoption = 'off'
exper = 'RAMS_WORLD'
fig, ax = dfunc.plot_eofs(alldsd, 'RAMS World', 'test', pcoption=pcoption)


if pcoption is not 'off':
    ax[0,0].plot(world.eof1, color='r')
    ax[1,0].plot(world.eof2, color='r')
    ax[2,0].plot(world.eof3, color='r')


    hist1,eg = np.histogram(world.pc1, bins=np.linspace(-5,5,50), density=True)
    hist1 *= 100./hist1.sum()
    ax[0,1].plot(eg[:-1]+np.diff(eg)[0]/2.,hist1, color='r', linewidth=1)


    hist2,eg = np.histogram(world.pc2, bins=np.linspace(-5,5,50), density=True)
    hist2 *= 100./hist2.sum()
    ax[1,1].plot(eg[:-1]+np.diff(eg)[0]/2.,hist2, color='r', linewidth=1)

    hist3,eg = np.histogram(world.pc3, bins=np.linspace(-5,5,50), density=True)
    hist3 *= 100./hist3.sum()
    ax[2,1].plot(eg[:-1]+np.diff(eg)[0]/2.,hist3, color='r', linewidth=1)

else:
    ax[0].plot(world.eof1, color='0.5', linewidth=2)
    ax[0].plot(0, 0, label='%d RAMS SIMS: %.1fM pts'%(len(exper_dsd_data.keys()), 1.0e-1*np.round(len(all_exper_data['d00'])/1.0e5)), 
                                                color='red', lw=4)
    ax[0].plot(0, 0, label='Global obs: %dK pts'%(np.round(world.size/1000.0)), color='0.5', lw=4)
    ax[0].legend(loc='best')

    ax[1].plot(world.eof2, color='0.5', linewidth=2)
    ax[2].plot(world.eof3, color='0.5', linewidth=2)

    for a in ax:
        a.set_ylim(-1.0, 1.0)


plt.savefig('{bp}/figures/{ex:}_eofs_{ne}.png'.format(bp=base_path, ex=exper, ne=n_expers), dpi=150)


### Make a plot of PC1 vs. PC2 using the dsd helper functions
fig, ax = plt.subplots(1, 1, figsize=(8,8))
t = dfunc.get_2d_pchist(alldsd.get_pc(rank=1),alldsd.get_pc(rank=2),'RAMS World',ax)
plt.savefig('{ex:}_PC1PC2.png'.format(ex=exper), dpi = 150)










