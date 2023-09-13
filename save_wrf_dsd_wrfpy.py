##Written by Brenda Dolan
##Read in WRF data, interpolate to a grid, save the gamma parameters

##July 2020

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
    #print('gnu',gnu)
    dn = dmean / np.float(gnu)
    #print('dn',dn)
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
    #print('md',Md,'SumMd',SumMd)
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
    #print(D0)

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




exper = 'Morrison_mu0'

###Mu value
mur = 0
###Nu value
nur = mur+1



namma_names  = sorted(glob.glob(f'{exper}/INTERP/{exper}*interpdat*.p'))
namma_norm = list()
#print(namma_names)

nw =[]
sigm = []
dm = []
lwc = []
d0 = []
rr = []

hskip = 2

for n in namma_names:
    print(f'working on {n}')
    tm = os.path.basename(n)
    #print(tm)
    time = tm[-7:-2]
    #print(time)
    print(time)
    data= pickle.load(open(f"{n}","rb"))
    zinterpH = data['z']

    whz = np.where(zinterpH<3.0)
    qrain_cart = np.squeeze(data['qr'][whz,::hskip,::hskip])
    nr_cart = np.squeeze(data['nr'][whz,::hskip,::hskip])
    rrt = np.squeeze(data['rrt'][::hskip,::hskip])

    nwt = np.squeeze(np.zeros_like(qrain_cart))
    sigmt = np.squeeze(np.zeros_like(qrain_cart))
    dmt = np.squeeze(np.zeros_like(qrain_cart))
    lwdt = np.squeeze(np.zeros_like(qrain_cart))
    d0t = np.squeeze(np.zeros_like(qrain_cart))
    
    #print(np.shape(nwt),whz[0][-1])
    #print(np.nanmax(qrain_cart))
    if np.nanmax(qrain_cart)>0.00001:
        for t,hg in enumerate((nwt[:,0,0])):
            #print('working on time ',t)
            for i,xx in enumerate((qrain_cart[0,:,0])):
                for j,yy in enumerate((qrain_cart[0,0,:])):
                        #print(qr[t,i,j])
                        if qrain_cart[t,i,j]>0.00001:
                            #print('coing to calc qr!')
                            #print('TEST',nur)
                            nwt[t,i,j],sigmt[t,i,j],dmt[t,i,j],lwdt[t,i,j],d0t[t,i,j] = calc_gamma_opt(qrain_cart[t,i,j],nr_cart[t,i,j],nur)
            #print(np.max(d0t))
    muu = np.zeros_like(d0t)+mur
    if np.nanmax(d0t) == 0.0:
        print('Uh-oh, nanmax(d0t) == 0')
    dat = {'d00': np.array(d0t),
                 'dmm': np.array(dmt),
                 'lwcc': np.array(lwdt),
                 'rrr': np.array(rrt),
                 'nww': np.array(nwt),
                 'sigm': np.array(sigmt),
                 #'ntt': ntsfc[whgd],
                 'muu': np.array(muu)}

    dat['lrr'] = make_lognorm(dat, 'rrr')
    dat['llwc'] = make_lognorm(dat, 'lwcc')    

    base_path ='WRF'

    out_filenamesfc = f"{base_path}/pickles/{exper}_dsd_wrf_{time}_skip{hskip}.p"

    pickle.dump(dat, open( out_filenamesfc, "wb" ) )
    print ('saved {f}'.format(f=out_filenamesfc))





