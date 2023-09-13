import matplotlib
import h5py
import numpy as np
import matplotlib.pyplot as plt
import glob
import numpy as np
import imp
import sys

import oa_stats as OA 
import analysis_tools as AT
import xarray as xr
import pickle

import scipy.stats as stats
import scipy
import pandas
from matplotlib import colors


import matplotlib as mpl
from scipy.special import gamma

def calc_gamma_opt(qr,nr,gnu): 
  #Mean mass diameter (meters), calcuated from m = cD^beta
    
    #print(np.shape(qr),np.shape(nr))
    ccfmas= 545.
    ppwmas = 3.
    dmean = (qr / (nr * ccfmas))**(1./ppwmas)
    #convert to Williams et al. units
    dmean = dmean * 1000. #(mm)
    #Hydrometeor density (kg/m3)
    dens = (6./3.14159)*ccfmas*(dmean**(ppwmas-3))
    #convert to Williams et al. units
    dens = dens / 1000. #(g/cm3)

    #Liquid/Ice water content temporary variable (g/m3)
    liwc  = qr * 1000.
    #Number mixing ratio temporary variable (#/m3)
    numc  = nr

    #****************************************************************************
    #Compute gamma distribution info
    #compute characteristic diameter (mm)
    dn = dmean / np.float(gnu)
    #print(dn)
    #compute gamma value
    
    gammlnx = gamma(gnu)
    gammaa = np.exp(gammlnx)
    #compute distribution bin increment (mm)
    dD = 0.02 * dn
    nbins = 1000
    db = np.zeros(nbins)
    fmg = np.zeros(nbins)
    #determine fraction per bin based on mean diameters and shape parameter
    dsizes = np.arange(0,nbins,dD)-0.5
    db = dD*dsizes
    fmg = dD*(db/dn)**(gnu-1.)/(dn*gammaa)*np.exp(-db/dn)
    #print(np.shape(fmg))
#     for q,ibin in enumerate(range(nbins-1)):
#         db[q] = dD * (float(ibin) - 0.5)
#         fmg[q] = dD * (db[q] / dn) ** (gnu - 1.) / (dn * gammaa) * np.exp(-db[q] / dn)
    #**************************************************************************
    #Work on bin values to get Dm
    SumMd = 0.0
    SumMdD =0.0
    
    Md = liwc*fmg
    #print(np.shape(Md))
    SumMd = np.sum(Md*dD)
    SumMdD = np.sum(Md*db*dD)
    
#     for q,ibin in enumerate(range(nbins-1)):
# 
#         Md      = liwc * fmg[q]      #Bin(mm) liquid/ice (g/m3)
#         D       = db[q]           #Bin(mm) diameter (mm)
#         SumMd   = SumMd + (Md*dD)
#         SumMdD  = SumMdD + (Md*D*dD)


    #Mass weighted mean diam (mm)
    Dm = SumMdD/SumMd
    Sum2nd = 0.0
    #Work on bin values of mass mixing ratio and compute total from distribution
#     for q,ibin in enumerate(range(nbins-1)):
# 
#         Md      = liwc * fmg[q]      #Bin(mm) liquid/ice (g/m3)
#         D       = db[q]           #Bin(mm) diameter (mm)
#         Sum2nd  = Sum2nd + (((D-Dm)**2)*Md*dD)    
# 
    #Standard deviation of mass spectrum (mm)
    Sum2nd  = np.sum(((db-Dm)**2)*Md*dD)    
    Sigma_m = np.sqrt(Sum2nd/SumMd)

    #Compute volumetic mean diameter (mm)
    #Use (gnu-1) since RAMS' gnu is different than Williams et al.
    D0 = Dm * (3.67 + (gnu-1.)) / (4.0 + (gnu-1.))

    #Compute normalized intercept parameter (1/mm * 1/m3) (Dolan et al. 2018)
    Nw = np.log10( 3.67**4. * 1000. / (3.1415 * dens) * (liwc/D0**4.) )

    return Nw, Sigma_m,Dm,liwc,D0


def calc_gamma(qr,nr,gnu): 
  #Mean mass diameter (meters), calcuated from m = cD^beta
    
    #print(np.shape(qr),np.shape(nr))
    ccfmas= 545.
    ppwmas = 3.
    dmean = (qr / (nr * ccfmas))**(1./ppwmas)
    #convert to Williams et al. units
    dmean = dmean * 1000. #(mm)
    #Hydrometeor density (kg/m3)
    dens = (6./3.14159)*ccfmas*(dmean**(ppwmas-3))
    #convert to Williams et al. units
    dens = dens / 1000. #(g/cm3)

    #Liquid/Ice water content temporary variable (g/m3)
    liwc  = qr * 1000.
    #Number mixing ratio temporary variable (#/m3)
    numc  = nr

    #****************************************************************************
    #Compute gamma distribution info
    #compute characteristic diameter (mm)
    dn = dmean / np.float(gnu)
    #print(dn)
    #compute gamma value
    gammlnx = gammln(gnu)
    gammaa = np.exp(gammlnx)
    #compute distribution bin increment (mm)
    dD = 0.02 * dn
    nbins = 1000
    db = np.zeros(nbins)
    fmg = np.zeros(nbins)
    #determine fraction per bin based on mean diameters and shape parameter
    dsizes = np.arange(0,nbins,dD)-0.5
    db = dD*dsizes
    fmg = dD*(db/dn)**(gnu-1.)/(dn*gammaa)*np.exp(-db/dn)

#     for q,ibin in enumerate(range(nbins-1)):
#         db[q] = dD * (float(ibin) - 0.5)
#         fmg[q] = dD * (db[q] / dn) ** (gnu - 1.) / (dn * gammaa) * np.exp(-db[q] / dn)
    #**************************************************************************
    #Work on bin values to get Dm
    SumMd = 0.0
    SumMdD =0.0
    
    Md = liwc*fmg
    SumMd = np.sum(Md*dD)
    SumMdD = np.sum(Md*db*dD)
    
#     for q,ibin in enumerate(range(nbins-1)):
# 
#         Md      = liwc * fmg[q]      #Bin(mm) liquid/ice (g/m3)
#         D       = db[q]           #Bin(mm) diameter (mm)
#         SumMd   = SumMd + (Md*dD)
#         SumMdD  = SumMdD + (Md*D*dD)


    #Mass weighted mean diam (mm)
    Dm = SumMdD/SumMd
    Sum2nd = 0.0
    #Work on bin values of mass mixing ratio and compute total from distribution
#     for q,ibin in enumerate(range(nbins-1)):
# 
#         Md      = liwc * fmg[q]      #Bin(mm) liquid/ice (g/m3)
#         D       = db[q]           #Bin(mm) diameter (mm)
#         Sum2nd  = Sum2nd + (((D-Dm)**2)*Md*dD)    
# 
    #Standard deviation of mass spectrum (mm)
    Sum2nd  = np.sum(((db-Dm)**2)*Md*dD)    
    Sigma_m = np.sqrt(Sum2nd/SumMd)

    #Compute volumetic mean diameter (mm)
    #Use (gnu-1) since RAMS' gnu is different than Williams et al.
    D0 = Dm * (3.67 + (gnu-1.)) / (4.0 + (gnu-1.))

    #Compute normalized intercept parameter (1/mm * 1/m3) (Dolan et al. 2018)
    Nw = np.log10( 3.67**4. * 1000. / (3.1415 * dens) * (liwc/D0**4.) )

    return Nw, Sigma_m,Dm,liwc,D0
def gammln(xx):
    
    
    cof= [76.18009173, -86.50532033, 24.01409822,-1.231739516, .120858003e-2, -.536382e-5]
    stp = 2.50662827465
    one = 1.0
    x=xx-1
    tmp=x+5.5
    tmp=(x+0.5)*np.log(tmp)-tmp
    ser=1
    for  j in range(5):
        x=x+one
        ser=ser+cof[j]/x
    
    gammln = tmp+np.log(stp*ser)

    return gammln




### The PCA crashes with the infinite values of RR created by taking the log of 0.0 points.

def remove_bad(dat, rainvar='lrr',nwvar ='nww'):

    whgood = np.isfinite(dat[rainvar])
    for k in dat.keys():
        #print('removing bad from ',k)
        dat[k] = dat[k][whgood]
    whgood = dat[nwvar]>0
    for k in dat.keys():
        dat[k] = dat[k][whgood]
    
    return dat


### This function will randomly subsample the data.
def random_sub(dat):
    ran = np.random.choice(np.arange(0,np.shape(dat['d00'])[0]),size=40000)
    for k in dat.keys():
        dat[k] = dat[k][ran]
    return dat

def format_subsample_dsd_data_diam(raindat,diamdat, tskip=1, hskip=1, nu=2,mu=0.0, d3=False):

    # BRF: not entirely sure what this function does right now,
    # but going to put comments in here so I can jot down notes
    # and whatnot, which will help me figure this out
    # Also, these functions should be in a different script. Will have to work on that
    # And need to change the name of this function to be more descriptive

    # I think this function takes in a dictionary of DSD params (at the surface)
    # and converts them into a nice format and also does some subsampling of the data

    # This should be done via some sort of dictionary so you don't have to write code for each variable

    # ***** MAYBE MAKE EACH OF THESE SIMULATIONS INTO AN OBJECT??? ******

    # this stuff here looks like a time and horizontal skip
    # np.squeeze will get rid of any dimensions that are size 1, which simplifies things

    # need to do this in a loop using a translation dictionary

    #zvals = raindat['z_coords'][:]
    #whrain =np.squeeze(np.where(zvals < 3000.))
    #print('shape of whrain',np.shape(whrain))

    nwd=np.squeeze(diamdat['rain_lognw_dm30'])[::tskip, ::hskip, ::hskip]#[0]#*1000.
    massd=np.squeeze(diamdat['rain_dm_d30'])[::tskip, ::hskip, ::hskip]#[0]#*1000.
    d0d = (3.67+nu)/((4+nu)*massd)
    #d0d = np.squeeze(diamdat['rain_d0_d30'][::tskip,whrain, ::hskip, ::hskip])
    print('npshape d0d',np.shape(d0d))
    d0=np.squeeze(diamdat['rain_d0_d43'])[::tskip, ::hskip, ::hskip]#[0]#*1000.
    nw=np.squeeze(diamdat['rain_lognw_d043'])[::tskip, ::hskip, ::hskip]#[0]#/1000.

    sigm=np.squeeze(diamdat['rain_sigma_mass'])[::tskip, ::hskip, ::hskip]#[0]#*1000.
    dm=np.squeeze(diamdat['rain_dm_d43'])[::tskip,  ::hskip, ::hskip]#[0]#*1000.
    d3 == True
    print('d3 is ',d3)
    if d3 == True:
        rr=np.squeeze(raindat['pcpvr'][::tskip, ::hskip, ::hskip])#[0]#*1000.
    else:
        rr=np.squeeze(raindat['pcprr'][::tskip,::hskip, ::hskip])#[0]#*1000.

    lwc=np.squeeze(raindat['rain_m3'])[::tskip,::hskip, ::hskip]#[0]#*1000.
    nt=np.squeeze(raindat['rain_concen_m3'])[::tskip,  ::hskip, ::hskip]#[0]#*1000.

    print('rr shape',np.shape(rr))
    # Here we're raveling and I'm not sure why? Also, why not pick out the 0th element
    # in the above lines of code? There might just be some extra dimension floating around?
    if d3 == True:
        print('saving the 3D version of data!')
        d0sfc = d0[:]
        dmsfc = dm[:]
        sigmsfc = sigm[:]
        rrsfc = rr[:]
        lwcsfc =lwc[:]
        nwdsfc = nwd[:]
        massdsfc = massd[:]
        d0dsfc = d0d[:]
        ntsfc = nt[:]
        nwsfc = nw[:]
        if mu == 0.0:
            muu = np.zeros_like(d0sfc)+nu
            
        else:
            muu = mu
        dat = {'d00': d0sfc,
                         'dmm': dmsfc,
                         'lwcc': lwcsfc,
                         'rrr': rrsfc,
                         'nww': nwsfc,
                         'sigm': sigmsfc,
                         'ntt': ntsfc,
                         'd0d':d0dsfc,
                         'nwd':nwdsfc,
                         'dmd':massdsfc,
                         'muu': muu}
    else:
        print(np.shape(d0))
        d0sfc = np.ravel(d0[:,:,:])
        dmsfc = np.ravel(dm[:,:,:])
        sigmsfc = np.ravel(sigm[:,:,:])
        nwdsfc = np.ravel(nwd[:,:,:])
        massdsfc = np.ravel(massd[:,:,:])
        d0dsfc =np.ravel(d0d[:,:,:])
        rrsfc = np.ravel(rr[:,:,:])
        lwcsfc = np.ravel(lwc[:,:,:])
        ntsfc = np.ravel(nt[:,:,:])
        nwsfc = np.ravel(nw[:,:,:])
    

        # okay, so specifying a mu/nu here, based on function input
        print(mu)
        if mu == 0.0:
            muu = np.zeros_like(d0sfc)+nu
        else:
            muu = mu


        # I'm guessing this is defining the good points here
        whgd1 = np.logical_and(rrsfc>0,ntsfc>0)
        whgd = np.squeeze(np.where(np.logical_and(whgd1,d0sfc>0)))
        awhgd = np.squeeze(np.argwhere(np.logical_and(whgd1,d0sfc>0)))
        shp = np.shape(d0[:,:,:])
        # then just go and grab points where the values are valid?
    #     dat = {'d00': d0sfc[whgd],
    #                  'dmm': dmsfc[whgd],
    #                  'lwcc': lwcsfc[whgd],
    #                  'rrr': rrsfc[whgd],
    #                  'nww': nwsfc[whgd],
    #                  'sigm': sigmsfc[whgd],
    #                  'ntt': ntsfc[whgd],
    #                  'muu': muu[whgd]}
        dat = {'d00': d0sfc[awhgd],
                     'dmm': dmsfc[awhgd],
                     'lwcc': lwcsfc[awhgd],
                     'rrr': rrsfc[awhgd],
                     'nww': nwsfc[awhgd],
                     'sigm': sigmsfc[awhgd],
                     'ntt': ntsfc[awhgd],
                     'muu': muu[awhgd],
                     'd0d':d0dsfc[awhgd],
                     'nwd':nwdsfc[awhgd],
                       'shp':shp,
                       'whgd':awhgd,
                     'dmd':massdsfc[awhgd]}




    ###Now let's make some of these in log-normal.
    #ramsdsddat['lnw'] = make_lognorm(ramsdsddat,'nww')

    # What is make_lognorm?
    # In dfunc, it's just taking the log10 of the values. So it's not actually making them lognormal

    dat['lnt'] = make_lognorm(dat, 'ntt')
    dat['lrr'] = make_lognorm(dat, 'rrr')
    dat['llwc'] = make_lognorm(dat, 'lwcc')

    return dat







def format_subsample_dsd_data(raindat, tskip=1, hskip=1, nu=2,d3=False,lev=0):

    # BRF: not entirely sure what this function does right now,
    # but going to put comments in here so I can jot down notes
    # and whatnot, which will help me figure this out
    # Also, these functions should be in a different script. Will have to work on that
    # And need to change the name of this function to be more descriptive

    # I think this function takes in a dictionary of DSD params (at the surface)
    # and converts them into a nice format and also does some subsampling of the data

    # This should be done via some sort of dictionary so you don't have to write code for each variable

    # ***** MAYBE MAKE EACH OF THESE SIMULATIONS INTO AN OBJECT??? ******

    # this stuff here looks like a time and horizontal skip
    # np.squeeze will get rid of any dimensions that are size 1, which simplifies things

    # need to do this in a loop using a translation dictionary

    zvals = raindat['z_coords'][:]
    whrain =np.squeeze(np.where(zvals < 3000.))
    print('shape of whrain',np.shape(whrain))

    d0=np.squeeze(raindat['rain_gam_d0'][::tskip,whrain, ::hskip, ::hskip])#[0]#*1000.
    nw=np.squeeze(raindat['rain_gam_lognw'][::tskip, whrain, ::hskip, ::hskip])#[0]#/1000.

    sigm=np.squeeze(raindat['rain_gam_sigma'][::tskip, whrain, ::hskip, ::hskip])#[0]#*1000.
    dm=np.squeeze(raindat['rain_gam_dm'][::tskip, whrain, ::hskip, ::hskip])#[0]#*1000.

    if d3 == True:
        try:
            rr=np.squeeze(raindat['pcpvr'][::tskip,whrain, ::hskip, ::hskip])#[0]#*1000.
        except:
            rr=np.squeeze(raindat['pcprr'][::tskip, ::hskip, ::hskip])#[0]#*1000.
        
    else:
        rr=np.squeeze(raindat['pcprr'][::tskip, ::hskip, ::hskip])#[0]#*1000.

    lwc=np.squeeze(raindat['rain_m3'][::tskip, whrain, ::hskip, ::hskip])#[0]#*1000.
    nt=np.squeeze(raindat['rain_concen_m3'][::tskip, whrain, ::hskip, ::hskip])#[0]#*1000.


    # Here we're raveling and I'm not sure why? Also, why not pick out the 0th element
    # in the above lines of code? There might just be some extra dimension floating around?
    if d3 == True:
        print('saving the 3D version of data!')
        d0sfc = d0[:]
        dmsfc = dm[:]
        sigmsfc = sigm[:]
        rrsfc = rr[:]
        lwcsfc =lwc[:]
        ntsfc = nt[:]
        nwsfc = nw[:]
        muu = np.zeros_like(d0sfc)+nu
        dat = {'d00': d0sfc,
                         'dmm': dmsfc,
                         'lwcc': lwcsfc,
                         'rrr': rrsfc,
                         'nww': nwsfc,
                         'sigm': sigmsfc,
                         'ntt': ntsfc,
                         'muu': muu}
    else:
        d0sfc = np.ravel(d0[:,lev,:,:])
        dmsfc = np.ravel(dm[:,lev,:,:])
        sigmsfc = np.ravel(sigm[:,lev,:,:])
        rrsfc = np.ravel(rr[:,:])
        lwcsfc = np.ravel(lwc[:,lev,:,:])
        ntsfc = np.ravel(nt[:,lev,:,:])
        nwsfc = np.ravel(nw[:,lev,:,:])
    
        shp = np.shape(d0[:,lev,:,:])
    # okay, so specifying a mu/nu here, based on function input
    muu = np.zeros_like(d0sfc)+nu

    # I'm guessing this is defining the good points here
    # okay, so specifying a mu/nu here, based on function input
    muu = np.zeros_like(d0sfc)+nu

    # I'm guessing this is defining the good points here
    whgd1 = np.logical_and(rrsfc>0,ntsfc>0)
    whgd = np.squeeze(np.where(np.logical_and(whgd1,d0sfc>0)))
    awhgd = np.squeeze(np.argwhere(np.logical_and(whgd1,d0sfc>0)))
    # then just go and grab points where the values are valid?
#     dat = {'d00': d0sfc[whgd],
#                  'dmm': dmsfc[whgd],
#                  'lwcc': lwcsfc[whgd],
#                  'rrr': rrsfc[whgd],
#                  'nww': nwsfc[whgd],
#                  'sigm': sigmsfc[whgd],
#                  'ntt': ntsfc[whgd],
#                  'muu': muu[whgd]}
    dat = {
            'd00': d0sfc[whgd],
                 'dmm': dmsfc[whgd],
                 'lwcc': lwcsfc[whgd],
                 'rrr': rrsfc[whgd],
                 'nww': nwsfc[whgd],
                 'sigm': sigmsfc[whgd],
                 'ntt': ntsfc[whgd],
                 'muu': muu[whgd],
                 'shp':shp,
                 'whgd':awhgd}

    print(np.shape(dat['d00']))
    ###Now let's make some of these in log-normal.
    #ramsdsddat['lnw'] = make_lognorm(ramsdsddat,'nww')

    # What is make_lognorm?
    # In dfunc, it's just taking the log10 of the values. So it's not actually making them lognormal

    dat['lnt'] = make_lognorm(dat, 'ntt')
    dat['lrr'] = make_lognorm(dat, 'rrr')
    dat['llwc'] = make_lognorm(dat, 'lwcc')

    return dat



def plot_time_frequency_rates(micro, ax, loc, exper, rank=0):
        # I'm guessing this function will go thru and calculate process rates and
        # plot them, but need to look back and confirm this


    totm=[]
    micro =np.ma.masked_equal(micro,-1)
    if rank == 0:

        for j in range(0,8):
            v= np.ma.sum(np.abs(micro[j,:,:]))
            totm.append(np.float(v))
    else:
        for j in range(0,8):
            try:
                v= len(np.squeeze(np.where(np.ma.ravel(micro)==j)))
                totm.append(np.float(v))
            except:
                totm.append(0.0)

    totmsum = np.float(np.ma.sum(np.abs(np.array(totm))))
#    print totm
    if totmsum>0:
        tfreq = np.zeros_like(np.array(totm))
        for j in range(0,8):
            tfreq[j] = np.float(totm[j])/totmsum*100.
        bot = 0
        r = []
        for i,t in enumerate(tfreq):
#            print bot, t
            rl = ax.bar(loc,t,0.9,bot,color=rate_colors[i+1])
            r.append(rl)
            bot = bot+t
    else:
        r = 0
#     ax.set_xticks(np.arange(0,7,1))
#     ax.set_xticklabels(hid_names,rotation=25)
#     ax.set_title('PC Thresh: {t} {e}'.format(t=thresh, e=exper))
    ax.set_xlabel('Time step')
    ax.set_ylabel('Relative Frequency')

    return r







def plot_time_frequency(groups, ax, loc, exper):

        # not sure what this function does yet
        # I think it goes thru and finds the fraction of values that correspond
        # to each group. But will have to confirm when I go thru the code itself
        # to see where/how it's called

    totm = []
    for j in range(0, 7):
        try:
            v = len(np.squeeze(np.where(np.ma.ravel(groups)==j)))
        except:
            v = 0
        totm.append(np.float(v))
    totmsum = np.float(np.sum(np.array(totm)))
#    print totm
    if totmsum>0:
        tfreq = np.zeros_like(np.array(totm))
        for j in range(0,7):
            tfreq[j] = np.float(totm[j])/totmsum*100.
        bot = 0
        p = []
        #p = np.zeros_like(totm)
        for i,t in enumerate(tfreq):
#            print bot, t
            pl = ax.bar(loc, t, 0.9, bot, color=hid_colors[i])
            p.append(pl)
            bot = bot+t
    else:
        p=0
#     ax.set_xticks(np.arange(0,7,1))
#     ax.set_xticklabels(hid_names,rotation=25)
    ax.set_title('PC Thresh: {t} {e}'.format(t=thresh, e=exper))
    ax.set_xlabel('Time step')
    ax.set_ylabel('Relative Frequency')
    return p






def modify_cmap(orig_cb_string, modify_dict, ncolors=8, new_name='newcmap'):
    
    oldcmap = plt.cm.get_cmap(orig_cb_string, ncolors) #generate a jet map with 10 values
    old_vals = oldcmap(np.arange(ncolors)) #extract those values as an array
    for k, v in modify_dict.iteritems():
        old_vals[k] = v
        new_cmap = mpl.colors.LinearSegmentedColormap.from_list(new_name, old_vals)
    return new_cmap

def modify_cmap(orig_cb_string, modify_dict, ncolors=8, new_name='newcmap'):
    
    oldcmap = plt.cm.get_cmap(orig_cb_string, ncolors) #generate a jet map with 10 values
    old_vals = oldcmap(np.arange(ncolors)) #extract those values as an array
    for k, v in modify_dict.items():
        old_vals[k] = v
        new_cmap = mpl.colors.LinearSegmentedColormap.from_list(new_name, old_vals)
    return new_cmap

def get_echo_top_height(refl,hts,thresh=0):
    
    maxvals = np.argmax(refl<thresh,axis=1)
    ets=hts[maxvals]
    return ets





def make_groups(pc1, ramsgroups):
    groupval = np.zeros_like(pc1)-1
    groupval[ramsgroups['Amb']]=0
    groupval[ramsgroups['Group1']] = 1
    groupval[ramsgroups['Group2']] = 2
    groupval[ramsgroups['Group3']] = 3
    groupval[ramsgroups['Group4']] = 4
    groupval[ramsgroups['Group5']] = 5
    groupval[ramsgroups['Group6']] = 6


    return groupval

def calc_dbz(nw,d0,mu):
    f_mu = (6./(3.67)**4.)*((3.67+mu)**(mu+4))/scipy.special.gamma(mu+4)
    gam = scipy.special.gamma(7+mu)
    end = d0**(7+mu)/(3.67+mu)**(7+mu)
    z = nw*f_mu*gam*end
    return 10.*np.log10(z)

def get_rams_params_2d(raindat,microdat,extra1=None,tskip=1,hskip=1, nu=2,xsize=400,ysize=500,tsize=31,zsize = 59,int_type = 'sum'):

    d0=np.squeeze(raindat['rain_gam_d0'][::tskip,:,::hskip,::hskip])#[0]#*1000.
    nw=np.squeeze(raindat['rain_gam_lognw'][::tskip,:,::hskip,::hskip])#[0]#/1000.

    sigm=np.squeeze(raindat['rain_gam_dm'][::tskip,:,::hskip,::hskip])#[0]#*1000.
    dm=np.squeeze(raindat['rain_gam_sigma'][::tskip,:,::hskip,::hskip])#[0]#*1000.

    rr=np.squeeze(raindat['pcprr'][::tskip,::hskip,::hskip])#[0]#*1000.
    lwc=np.squeeze(raindat['rain_m3'][::tskip,:,::hskip,::hskip])#[0]#*1000.
    nt=np.squeeze(raindat['rain_concen_m3'][::tskip,:,::hskip,::hskip])#[0]#*1000.

    d0sfc = (d0[:,0,:,:])
    print('d0 shape:',np.shape(d0sfc))
    dmsfc = (dm[:,0,:,:])
    sigmsfc = (sigm[:,0,:,:])
    rrsfc = (rr[:,:])
    lwcsfc = (lwc[:,0,:,:])
    ntsfc = (nt[:,0,:,:])
    nwsfc = (nw[:,0,:,:])
#    print np.shape(aggsum)
#    whgd = np.where(d0sfc>0)
#    whgdarg = np.argwhere(d0sfc>0)
    
    
    muu = np.zeros_like(d0sfc)+(nu-1)

    
    dbzsfc = calc_dbz(10.**nwsfc,d0sfc, muu)


#     xdat = raindat['x_coords']#[::hskip]
#     ydat = raindat['y_coords']#[::hskip]
    xsize = int(microdat['x_coords'][::hskip].shape[0])
    ysize = int(microdat['y_coords'][::hskip].shape[0])
    zsize = int(microdat['z_coords'].shape[0])
    tsize = int(microdat['t_coords'][::tskip].shape[0])
    
    xx = raindat['x_coords'][:]
    yy = raindat['y_coords'][:]


    
    x = np.arange(0,xsize,1)
    y = np.arange(0,ysize,1)
    z = np.arange(0,zsize,1)
    t = np.arange(0,tsize,1)
    agg = xr.DataArray(microdat['aggregatet'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])
    c2r = xr.DataArray(microdat['cld2raint'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])
    i2r = xr.DataArray(microdat['ice2raint'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])
    melti = xr.DataArray(microdat['melticet'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])
    r2i = xr.DataArray(microdat['rain2icet'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])
    rime = xr.DataArray(microdat['rimecldt'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])
    vapi = xr.DataArray(microdat['vapicet'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])
    vapl = xr.DataArray(microdat['vapliqt'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])

    if int_type == 'sum':
        aggsum = np.sum(agg.values,axis=1)
        c2rsum = np.sum(c2r.values,axis=1)
        i2rsum = np.sum(i2r.values,axis=1)
        meltisum = np.sum(melti.values,axis=1)
        r2isum = np.sum(r2i.values,axis=1)
        rimesum = np.sum(rime.values,axis=1)
        vapisum = np.sum(vapi.values,axis=1)
        vaplsum = np.sum(vapl.values,axis=1)
    else:
        aggsum = np.nanmax(agg.values,axis=1)
        c2rsum = np.nanmax(c2r.values,axis=1)
        i2rsum = np.nanmax(i2r.values,axis=1)
        meltisum = np.nanmax(melti.values,axis=1)
        r2isum = np.nanmax(r2i.values,axis=1)
        rimesum = np.nanmax(rime.values,axis=1)
        vapisum = np.nanmax(vapi.values,axis=1)
        vaplsum = np.nanmax(vapl.values,axis=1)

    dat = {'d00':d0sfc,
            'dbz':dbzsfc,
                 'dmm':dmsfc,
                 'dbz':dbzsfc,
                 'lwcc':lwcsfc,
                 'rrr':rrsfc,
                 'nww':nwsfc,
                 'sigm':sigmsfc,
                 'ntt':ntsfc,
                 'muu':muu,
                 'agg': aggsum,
                 'c2r': c2rsum,
                 'i2r': i2rsum,
                 'melt': meltisum,
                 'r2i': r2isum,
                 'rime': rimesum,
                 'vapi': vapisum,
                 'vapl': vaplsum}
#                 'x': xdat,
#                 'y':ydat}

    if extra1 is not None:

        try:
            lhf = xr.DataArray(extra1['latheatfrz'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])
            lhft = xr.DataArray(extra1['latheatfrzt'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])
            lhv = xr.DataArray(extra1['latheatvap'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])
            lhvt = xr.DataArray(extra1['latheatvapt'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])
        except:
            print("no latent heating information")


        refl = xr.DataArray(extra1['reflect_all'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])
        w = xr.DataArray(extra1['w'][::tskip,:,::hskip,::hskip], coords=[t,z,y,x], dims=['time', 'z','y','x'])


        hts=microdat['z_coords'][:]
        eth = get_echo_top_height(refl,hts)
        dbz30 = get_echo_top_height(refl,hts,thresh=30)
        dbz50 = get_echo_top_height(refl,hts,thresh=50)
        
        vta = xr.DataArray(extra1['vt_aggregatet'][::tskip,::hskip,::hskip], coords=[t,y,x], dims=['time','y','x'])
        vtc2r = xr.DataArray(extra1['vt_cld2raint'][::tskip,::hskip,::hskip], coords=[t,y,x], dims=['time','y','x'])
        vti2r = xr.DataArray(extra1['vt_ice2raint'][::tskip,::hskip,::hskip], coords=[t,y,x], dims=['time','y','x'])
        vtmi = xr.DataArray(extra1['vt_melticet'][::tskip,::hskip,::hskip], coords=[t,y,x], dims=['time','y','x'])
        vtr2i = xr.DataArray(extra1['vt_rain2icet'][::tskip,::hskip,::hskip], coords=[t,y,x], dims=['time','y','x'])
        vtrimec = xr.DataArray(extra1['vt_rimecldt'][::tskip,::hskip,::hskip], coords=[t,y,x], dims=['time','y','x'])
        vtvi = xr.DataArray(extra1['vt_vapicet'][::tskip,::hskip,::hskip], coords=[t,y,x], dims=['time','y','x'])
        vtvl = xr.DataArray(extra1['vt_vapliqt'][::tskip,::hskip,::hskip], coords=[t,y,x], dims=['time','y','x'])


        if int_type == 'sum':
            try:
                lhfsum = np.sum(lhf.values,axis=1)
                lhftsum = np.sum(lhft.values,axis=1)
                lhvsum = np.sum(lhv.values,axis=1)
                lhvtsum = np.sum(lhvt.values,axis=1)
            except:
                pass
            reflsum = np.sum(refl.values,axis=1)
            wsum = np.sum(w.values,axis=1)
#                 vtasum = np.sum(vta.values,axis=1)
#                 vtc2rsum = np.sum(vtc2r.values,axis=1)
#                 vti2rsum = np.sum(vti2r.values,axis=1)
#                 vtmisum = np.sum(vtmi.values,axis=1)
#                 vtr2isum = np.sum(vtr2i.values,axis=1)
#                 vtrimecsum = np.sum(vtrimec.values,axis=1)
#                 vtvisum = np.sum(vtvi.values,axis=1)
#                 vtvlsum = np.sum(vtvl.values,axis=1)

        else:
            try:
                lhfsum = np.nanmax(lhf.values,axis=1)
                lhftsum = np.nanmax(lhft.values,axis=1)
                lhvsum = np.nanmax(lhv.values,axis=1)
                lhvtsum = np.nanmax(lhvt.values,axis=1)
            except:
                pass
            reflsum = np.nanmax(refl.values,axis=1)
            wsum = np.nanmax(w.values,axis=1)
#                 vtasum = np.nanmax(vta.values,axis=1)
#                 vtc2rsum = np.nanmax(vtc2r.values,axis=1)
#                 vti2rsum = np.nanmax(vti2r.values,axis=1)
#                 vtmisum = np.nanmax(vtmi.values,axis=1)
#                 vtr2isum = np.nanmax(vtr2i.values,axis=1)
#                 vtrimecsum = np.nanmax(vtrimec.values,axis=1)
#                 vtvisum = np.nanmax(vtvi.values,axis=1)
#                 vtvlsum = np.nanmax(vtvl.values,axis=1)
        x, y = np.meshgrid(xx,yy)
        cs_arr = []
        for i in range(len(reflsum)):
            dbz_comp = reflsum[i,...]
            yh_cs, yh_cc, yh_bkgnd = shy.conv_strat_latlon(dbz_comp, y, x, CoreThresh=37.0, method='SYH', a=8, b=64, sm_rad=2,dist_thresh=5)

            yh_cc_ma = np.ma.masked_where(yh_cc < 0, yh_cc)


            cs_arr_s = np.full(yh_cs.shape, np.nan)
            yh_conv = (yh_cs == 3) | (yh_cs == 4) | (yh_cs == 5)
            yh_mixed = (yh_cs == 1) | (yh_cs == 2)
            yh_strat = yh_cs == 0

            cs_arr_s[yh_conv] = 2
            cs_arr_s[yh_mixed] = 2
            cs_arr_s[yh_strat] = 1

            mask = np.where(np.logical_or(np.isnan(dbz_comp),dbz_comp<=-10.))
            cs_arr_s[mask] = -9
    
            cs_arr.append(cs_arr_s)
            
        whc_check = np.where(np.logical_and(rrsfc >10., wsum >= 1.))
        whs_check = np.where(np.logical_and(rrsfc<=10., wsum < 1.))# 
        mask = np.where(np.logical_or(np.isnan(reflsum),reflsum<=-10.))
        rrw_arr = np.ones_like(rrsfc)
        rrw_arr[whc_check] =2
        rrw_arr[whs_check] =1

        rrw_arr[mask] =-9


        try:
            dat['lhf'] = lhfsum
            dat['lhft'] = lhftsum
            dat['lhv'] = lhvsum
            dat['lhvt'] = lhvtsum
        except:
            pass
        dat['refl'] = reflsum
        dat['w'] = wsum
        dat['vta'] = vta.values
        dat['vtc2r'] = vtc2r.values
        dat['vti2r'] = vti2r.values
        dat['vtmi'] = vtmi.values
        dat['vtr2i'] = vtr2i.values
        dat['vtrimec'] = vtrimec.values
        dat['vtvi'] = vtvi.values
        dat['vtvl'] = vtvl.values
        
        dat['eth'] = eth
        dat['dbz30'] = dbz30
        dat['dba50'] = dbz50
        
        dat['css'] = np.array(cs_arr)
        dat['csw'] = np.array(rrw_arr)


    ###Now let's make some of these in log-normal.
    #ramsdsddat['lnw'] = make_lognorm(ramsdsddat,'nww')
    dat['lnt'] = make_lognorm(dat,'ntt')
    dat['lrr'] = make_lognorm(dat,'rrr')
    dat['llwc'] = make_lognorm(dat,'lwcc')

#    dat = {}
    
    return dat

def contour_density(dat1, dat2, exper,ax):
    #histogram definition
    xyrange = [[0,6],[0,6]] # data range
    bins = [40,40] # number of bins
    thresh =2  #density threshold

    xdat = np.array(dat1)
    ydat = np.array(dat2)
    exper= exper
    d0=np.arange(0,6)
    d0l=np.arange(0,1.9,0.1)

#     print ('xdatshape', xdat.shape)
#     print ('ydatshape', ydat.shape)

    hh, locxy = np.histogramdd((np.ravel(xdat), np.ravel(ydat)), range=xyrange, bins=bins) #,normed=True)

    posx = np.digitize(xdat, locxy[0])
    posy = np.digitize(ydat, locxy[1])


    #select points within the histogram
    ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
    hhsub = hh[posx[ind] - 1, posy[ind] - 1] # values of the histogram where the points are
    xdat1 = xdat[ind][hhsub < thresh] # low density points
    ydat1 = ydat[ind][hhsub < thresh]
    hh[hh < thresh] = np.nan # fill the areas with low density by NaNs
    cb = ax.imshow(np.flipud(hh.T)/np.nansum(hh)*100.,cmap='CMRmap',extent=np.array(xyrange).flatten(),norm=colors.LogNorm(), interpolation='none', origin='upper')
         
    for label in (ax.get_yticklabels()):
        label.set_fontsize(18)
    for label in (ax.get_xticklabels()):
        label.set_fontsize(18)

    bringi_log_Nw = (-1.6*d0) + 6.3
    liz_log_Nw=0.0*(d0l)+3.85
    br,=ax.plot(d0,bringi_log_Nw,linestyle='--',color='silver',label='Bringi09',lw=5)
    lz,=ax.plot(d0l,liz_log_Nw,linestyle='-',color='silver',label='Thompson15',lw=5)


 
    ax.grid(True)
    ax.set_xlabel('D$_0$ (mm)',fontsize=20)
    ax.set_ylabel('logN$_w$',fontsize=20)
    ax.set_title('{e} '.format(e=exper),fontsize=22)
    ax.set_ylim(1,6)
    ax.set_xlim(0,4)
    #plt.plot(xdat1, ydat1, '.',color='darkblue')
#    plt.savefig('{e}_density_nwdm_norm.png'.format(e=exper),dpi=200)
    return cb


def plot_eofs(dsd,experiment,extra,pcoption='distribution'):
    fig, ax = dsd.plot_eofs(ranks = [1,2,3],pcoption=pcoption,colors='off')
    #ax[0,1].set_ylim(-5,5)
    fig.suptitle('Main EOFs of drop size distribution\nquantities during %s, %d%% variance explained'%
                     (experiment, dsd.get_variance_range(ranks = [1,2,3,4])), fontsize = 16)
    fig.subplots_adjust(top = 0.87)

    return fig, ax
    
hid_names = ['Amb','Grp1', 'Grp2', 'Grp3', 'Grp4',
              'Grp5', 'Grp6']
hid_colors = [ 'grey','Red', 'Green', 'GOld', 'Blue',
              'Orange', 'Purple']


def plot_frequency(groups,ax,thresh):
    totm=[]
    for j in range(0,7):
        try:
            v= len(np.squeeze(np.where(np.ma.ravel(groups)==j)))
        except TypeError as tp:
            print( tp)
            v = 0
        totm.append(float(v))
    #print totm
    totmsum = float(np.sum(np.array(totm)))
    #print totmsum
    if totmsum > 1:
        tfreq = np.zeros_like(np.array(totm))
        for j in range(0,7):
            tfreq[j] = float(totm[j])/totmsum*100.
     #   print tfreq,np.shape(tfreq)
        for i,t in enumerate(tfreq):
            ax.bar(i,t,color=hid_colors[i])
    ax.set_xticks(np.arange(0,7,1))
    ax.set_xticklabels(hid_names,rotation=25)
    ax.set_title('PC Thresh: {t}'.format(t=thresh))



    
def exper_eof(raindata,good_vars,labels):

    dsd = OA.Dataset(data = raindata, data_type = 'attributes')#, method = 'svd')
    dsd.into_matrix(varlist = good_vars)
    C = dsd.covariance_matrix()
    dsd.calculate_eofs()
    dsd.set_eof_labels(labels)
    dsd.auto_flip()
    return dsd

def make_lognorm(raindata,var):
    return np.log10(raindata[var])
    
import matplotlib.colors as colors


def get_2d_pchist(pc1,pc2,exper,ax, thresh=1.5):
    H,xedges,yedges = np.histogram2d(pc1,pc2,bins = [np.arange(-5,30,0.05),np.arange(-5,30,0.05)])
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
#    ax.axhline(y=1.5,color='yellow')
#    ax.axhline(y=-1.5,color='blue')
    #thresh = 1.5
    p = ax.axvspan(-10.0,-thresh, facecolor='green', alpha=0.3)
    p = ax.fill_between([-thresh, thresh], -thresh,-10,facecolor='blue', alpha=0.3)    
    p = ax.fill_between([-thresh,thresh], thresh,10,facecolor='gold', alpha=0.3)
    p = ax.fill_between([thresh, 10], thresh,10,facecolor='darkorange', alpha=0.3)
    p = ax.fill_between([thresh, 10], -thresh,thresh,facecolor='red', alpha=0.3)
    p = ax.fill_between([thresh, 10], -10,-thresh, facecolor='purple', alpha=0.3)
    p = ax.fill_between([-thresh, thresh], -1.0,thresh, facecolor='grey', alpha=0.3)

    p = ax.axhline(thresh, color='r', linestyle='dashed',lw=2)
    p = ax.axhline(-thresh, color='r', linestyle='dashed',lw=2)
    p = ax.axvline(thresh, color='r', linestyle='dashed',lw=2)
    p = ax.axvline(-thresh, color='r', linestyle='dashed',lw=2)




    cb =ax.imshow(Hmasked.T/float(np.sum(Hmasked))*100., interpolation='nearest', origin='lower',
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],cmap='CMRmap',zorder=10,norm=colors.LogNorm())

    cb1 = plt.colorbar(cb,ax = ax,shrink=0.4)
    cb1.set_label('Frequency (%)')

    p = ax.axvline(0, color='k', linestyle='dashed',lw=1,zorder=15)
    p = ax.axhline(0, color='k', linestyle='dashed',lw=1,zorder=15)


    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(15)

    ax.set_xticks(np.arange(-5,10+1, 1.))
    ax.set_yticks(np.arange(-5,10+1, 1.))

    ax.set_xlim(-4,8)
    ax.set_ylim(-5,5)
    ax.set_title(exper,fontsize=16)
    ax.set_xlabel('PC1 (Std. Anomaly)',fontsize=16)
    ax.set_ylabel('PC2 (Std. Anomaly)',fontsize=16)
    return cb
    
def get_groups(pc1, pc2, thresh1, thresh2, thresh3, thresh4):
        from collections import OrderedDict
        condc = pc1 > thresh1
        conds = pc1 < (-1.)*thresh2
        conda = pc1 > -999.

        condi = pc2 < (-1.)*thresh4
        condw = pc2 > thresh3

        condcc = np.logical_and(condc,condw)
        condci = np.logical_and(condc,condi)
    
        condc11 = np.logical_and(condc,~condcc)
        condc12 = np.logical_and(condc,~condci)
        condcon = np.logical_and(condc11,condc12)
        whmask = np.ma.getmask(pc1)

        condshal0 = np.logical_and(~whmask,~condcc)
        condshal = np.logical_and(condw,condshal0)

        condstrat1 = np.logical_and(~whmask,~condw)
        condstrat0 = np.logical_and(condstrat1,~condcc)
        condstrat = np.logical_and(conds,condstrat0)
    
        condice0 = np.logical_and(~whmask,~condci)
        condice = np.logical_and(condi,condice0)

        amb0 = np.logical_and(~whmask,~condc)
        amb1 = np.logical_and(amb0,~conds)
        amb2 = np.logical_and(~condi,~condw)
        amb = np.logical_and(amb1,amb2)

        amb0s = np.logical_and(~whmask,~condci)
        amb1s = np.logical_and(amb0s,~condcc)
        amb2s = np.logical_and(amb1s,~condshal)
        amb3s = np.logical_and(amb2s,~condcon)
        amb4s = np.logical_and(amb3s,~condci)
        amb5s = np.logical_and(amb4s,~conds)
        ambs = np.logical_and(amb5s,~condice)

        whnpc2c =np.where(condci)
        whppc2c =np.where(condcc)

        whsppc1 =np.where(condcon)
        whsnpc1 =np.where(condstrat)

        whppc1 =np.where(condc)
        whnpc1 =np.where(conds)

        whppc2 =np.where(condw)
        whnpc2 =np.where(condi)
    
        whsppc2 = np.where(condshal)
        whsnpc2 = np.where(condice)
    
        whamb = np.where(amb)
        whambs = np.where(ambs)
        whall = np.where(conda)
        
        grp1=whsppc1
        grp2=whsnpc1
        grp3=whsppc2
        grp4=whsnpc2
        grp5=whppc2c
        grp6=whnpc2c
        grpa = whambs
        groups =  OrderedDict()
        groups['Amb']=grpa
        groups['Group1']=grp1
        groups['Group2']=grp2
        groups['Group3']=grp3
        groups['Group4']=grp4
        groups['Group5']=grp5
        groups['Group6']=grp6
        return groups



def plot_hist(dat1,dat2,bins,ax,title1,title2,lw=4,color1='r',color2='g',alpha=0.9):
    hist1, edg = np.histogram(dat1,bins)
    hist2, edg = np.histogram(dat2,bins)
    
    ax.plot(edg[:-1],hist1*100./hist1.sum(),label=title1,lw=lw,color=color1,alpha=alpha)
    ax.plot(edg[:-1],hist2*100./hist2.sum(),label=title2,lw=lw,color=color2,alpha=alpha)
    return ax
    
def plot_groups(var1,var2,groups,ax,line_on=0):
    cols={'Group1':'r',
         'Group2':'g',
         'Group3':'gold',
         'Group4':'blue',
          'Group5':'orange',
         'Group6':'purple',
         'Amb':'gray'}
    ax.scatter(np.array(var1),np.array(var2),color='grey',alpha=0.2,marker='.')    
    for g in groups.keys():
        ax.scatter(np.array(var1)[groups[g]],np.array(var2)[groups[g]],color=cols[g],alpha=0.2,marker='.',label=g)
    if line_on == 1:
        d0=np.arange(0,5,0.2)
        d0l=np.arange(0,2.2,0.2)


        nwdum=np.arange(0,5,0.2)
        bringi_log_Nw = (-1.6*d0) + 6.3
        liz_log_Nw=0.0*(d0l)+3.85

        br,=ax.plot(d0,bringi_log_Nw,linestyle='--',color='gray',alpha=0.5,linewidth=4,label='Bringi 2009')



    return ax