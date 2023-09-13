##Written by Brenda Dolan
##Read in WRF data, interpolate to a grid, save the gamma parameters

##July 2020

import wrf
from netCDF4 import Dataset
import numpy as np
import glob
import pickle
import os

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
    #print('returning data!')
    return Nw, Sigma_m,Dm,liwc,D0

def make_lognorm(raindata,var):
    return np.log10(raindata[var])
   

##########################################################################################
##########################################################################################
##########################################################################################
# font = {'family' : 'normal',
#         'size'   : 16}
# 
# matplotlib.rc('font', **font)




exper = 'Thompson_mu1'


namma_names  = sorted(glob.glob(f'doe.wrf.supercell.{exper}/wrfout_d01_*nc'))
namma_norm = list()


mur = 1

nw =[]
sigm = []
dm = []
lwc = []
d0 = []
rr = []
for n in namma_names:
    print(f'working on {n}')
    tm = os.path.basename(n)
    #print(tm)
    time = tm[22:27]
    print(time)
    ncfile= Dataset(n)
    qr = wrf.getvar(ncfile, "QRAIN")
    nr = wrf.getvar(ncfile,"QNRAIN")
    rrt = wrf.getvar(ncfile,"RAINC").values
    u=wrf.getvar(ncfile,"U")
    v=wrf.getvar(ncfile,"V")
    w=wrf.getvar(ncfile,"W")
    
    #Need to interpret some of the variables to the regular grid heights
    zinterpH = np.arange(0,15.25,0.25)
    whz = np.where(zinterpH<3.0)
    qrain_cart = wrf.vinterp(ncfile,qr,'ght_msl',zinterpH).values
    nr_cart = wrf.vinterp(ncfile,nr,'ght_msl',zinterpH).values
#     ud_cart =wrf.destagger(u,'bottom_top')
#     u_cart = wrf.vinterp(ncfile,u,'ght_msl',zinterpH)
#     v_cart = wrf.vinterp(ncfile,v,'ght_msl',zinterpH)
#     w_cart = wrf.vinterp(ncfile,w,'ght_msl',zinterpH)
    interp_dat = {'qr':qrain_cart,'nr':nr_cart,'mu':mur,'rrt':rrt,'z':zinterpH}#'u':u_cart,'v':v_cart,'w':w_cart}
    
    base_path ='WRF'
    out_filename = f"{base_path}/pickles/{exper}_interpdat_wrf_{time}.p"

    pickle.dump(interp_dat, open( out_filename, "wb" ) )
    print ('saved {f}'.format(f=out_filename))




