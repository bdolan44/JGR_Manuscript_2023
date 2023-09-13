## Written by Brenda Dolan 
## December 2016
## bdolan@atmos.colostate.edu
1
## Data needs to already be processed and saved in common format (dictionary with different
## Variables, and saved in pickle files.
## Really should convert to project basis vectors on any dataset
##
##
#Here goes nothing...setting up a class for disdrometer data

import numpy as np 
import oa_stats as OA 
import matplotlib.pyplot as plt 
import scipy.stats as stats
import pickle as pickle
#from mpl_toolkits.basemap import Basemap
from scipy.io.idl import readsav
import analysis_tools as AT
from matplotlib.colors import LogNorm
import plot_oa as POA
from datetime import datetime, timedelta
from collections import OrderedDict




class Disdrometer(object):
    
    def __init__(self,distype='2dvd',color='black',region=0,exper='global',extra='test',loc=0,
                good_vars=[''],sav=1,thresh1=1.5,thresh2=1.5,thresh3=1.5,thresh4=1.5,pickdir='scratch/PCA_NEW/SIGM_TEST/'):
#        self.pccolors = {'conv': 'red','strat':'green','warm':'orange'}
        labels = {'nww':'N$_w$', 
                  'd00':'D$_0$',
                  'lnww':'logN$_w$',
                  'lnw':'logN$_w$',
                  'llwc':'logLWC',
                  'lrr':'logR',
                  'lnt':'logN$_t$',
                  'sigm':'$\sigma$$_m$',
                  'dmm':'D$_m$',
                  'dbzss':'dbz',
                  'muu':'$\mu$', 
                  'lwcc':'LWC', 
                  'rrr':'R',
                  'ntt':'N$_t$'}#'D$_{max}$',]#, '$\sigma$']

        eoflabels =OrderedDict()
        for k in good_vars:
            eoflabels[k]=labels[k]
        self.eoflabels=eoflabels

        labelsunits = {'nww':'N$_w$', 
                  'lnww':'logN$_w$',
                  'd00':'D$_0$ (mm)',
                  'dmm':'D$_m$ (mm)',
                  'dbzss':'dbz',
                  'muu':'$\mu$', 
                  'lwcc':'LWC (g m$^{-3}$)', 
                  'rrr':'R (mm hr$^{-1}$)',
                  'ntt':'N$_t$ (drops)'}#'D$_{max}$',]#, '$\sigma$']

#        world_raindata=pickle.load(open( "{s}world_raindata_dm.p".format(s=pickdir), "rb" ) )
#        world_dsd=pickle.load(open( "{s}world_dsd_dm.p".format(s=pickdir), "rb" ) )

        world_raindata=pickle.load(open( "{s}world_raindata_dm_lognorm_sigm_grp.p".format(s=pickdir), "rb" ) )
        world_dsd=pickle.load(open( "{s}world_dsd_dm_lognorm_sigm_grp.p".format(s=pickdir), "rb" ) )
        data_pos = pickle.load( open( '{s}data_pos_dm_grp.p'.format(s=pickdir), "rb" ) )
        data_size= pickle.load( open ( '{s}data_size_dm_grp.p'.format(s=pickdir), "rb") )
    
        self.distype = distype
        self.color = color
        self.exper = exper
        if self.exper != 'global':
            self.raindat,self.dsddat=self.read_raindata(sav=sav)
            self.eof1,self.eof2,self.eof3 = self.get_eof()
        self.extra = extra
        self.good_vars = good_vars

        self.labels = labels
        self.labelsunits = labelsunits
        self.thresh1 = thresh1
        self.thresh2 = thresh2
        self.thresh3 = thresh3
        self.thresh4 = thresh4
        self.locs = self.get_location()
        self.region = self.get_region()
        
#        self.rain, self.dsd,self.size = self.read_raindata(sav)
        if exper == 'global':
            self.rain = world_raindata
            self.dsd = world_dsd
            self.eof1 = world_dsd.get_eof(rank=1)
            self.eof2 =world_dsd.get_eof(rank=2)
            self.eof3 = world_dsd.get_eof(rank=3)
            self.pos = 0
            self.size = np.shape(self.rain['d00'])[0]
            print ('pos:size',self.pos,self.size)

        else:
            self.pos = data_pos[self.exper]
            self.size = data_size[self.exper]
            print ('pos:size',self.pos,self.size)
            self.dsd = world_dsd
            self.rain= self.sub_world(world_raindata,world_dsd,self.pos,self.size)
            self.datetime = self.get_apu_times()

#        self.dsd = dsddata
#        self.loc = loc
#        self.ctypes=['conv','strat','warm','ldlc']
        self.whppc1,self.whnpc1,self.whppc2,self.whnpc2,self.whnpc2c,self.whsnpc2,\
            self.whppc2c,self.whsppc2,self.whsppc1,self.whsnpc1,self.whall,self.whamb,\
            self.whambs = self.get_pcthresh(self.thresh1,self.thresh2,self.thresh3,self.thresh4,self.pos,self.size)
#        self.pc1,self.npc1,self.pc2,self.npc2,self.pc3,self.npc3 = self.get_pc()        

        self.cloudtypes()
        self.subrain = self.subset(self.thresh1,self.thresh2,self.thresh3,self.thresh4) 

    pass
    
    def get_apu_times(self):
        from datetime import datetime
        date_times = []
        print (self.raindat.keys())
        if 'datee' in self.raindat.keys():
            for dt_tuple in zip(self.raindat['datee'], ("{:06d}".format(int(i)) for i in self.raindat['timee'])):
                #print dt_tuple "{:06d}".format(int(
                date_times.append(datetime.strptime(','.join(dt_tuple), '%Y%m%d,%H%M%S'))
        elif 'doyy' in self.raindat.keys():
            if self.exper == 'twpice':
                for i,dt_tuple in enumerate(zip(self.raindat['doyy'],self.raindat['hhh'],self.raindat['minn'])):
                    base = datetime(self.raindat['year'][i],1,1)
                    date_times.append(base + timedelta(days = int(dt_tuple[0])-1, hours=int(dt_tuple[1]), minutes=int(dt_tuple[2])))
            else:
                yr=self.raindat['yyy'][0]
                base = datetime(yr,1,1)
                for dt_tuple in zip(self.raindat['doyy'],self.raindat['hhh'],self.raindat['mmm']):
                    date_times.append(base + timedelta(days = int(dt_tuple[0])-1, hours=int(dt_tuple[1]), minutes=int(dt_tuple[2])))
        else:
            print ('what is up with the date?')
        
        if self.exper == 'olympexapu':
            for i, d in enumerate(date_times):

        #print 'checking',d.year,d.month
                if (d.year == 2015 and d.month < 6):
        
                    date_times[i]= d.replace(year=2016)        
                    print ('fixing date',date_times[i])
                date_times[i]= date_times[i]+timedelta(days=-1)      

        
        return date_times

    def cloudtypes(self):
        self.ctypes=OrderedDict()
        self.ctypes['conv'] = {'color':'red',
                              'longname':'convective',
                              'pcnam':'Group 1 FULL'}
        self.ctypes['strat'] = {'color':'green',
                              'longname':'stratiform',
                              'pcnam':'Group 2 FULL'}
        self.ctypes['convs'] = {'color':'red',
                              'longname':'convective_sub',
                              'pcnam':'Group 1'}
        self.ctypes['strats'] = {'color':'green',
                              'longname':'stratiform_sub',
                              'pcnam':'Group 2'}

        self.ctypes['warm'] = {'color':'gold',
                              'longname':'warm',
                              'pcnam':'Group 3 FULL'}
        self.ctypes['ldlc'] = {'color':'dodgerblue',
                                'longname':'LDLC',
                                'pcnam':'Group 4 FULL'}
        self.ctypes['hstrat'] = {'color':'blue',
                              'longname':'heavy stratiform',
                              'pcnam':'Group 4'}

        self.ctypes['ice'] = {'color':'indigo',
                              'longname':'melted ice',
                              'pcnam':'Group 6'}
        self.ctypes['shallow'] = {'color':'gold',
                              'longname':'shallow, weak',
                              'pcnam':'Group 3'}

        self.ctypes['cc'] = {'color':'orange',
                              'longname':'Collision-coalescence',
                              'pcnam':'Group 5'}
        self.ctypes['all'] = {'color':'grey',
                             'longname':'allpoints',
                             'pcnam':'ALL'}
        self.ctypes['amb'] = {'color':'grey',
                             'longname':'ambiguous',
                             'pcnam':'AMB'}
        self.ctypes['ambs'] = {'color':'grey',
                             'longname':'ambiguous_sub',
                             'pcnam':'AMBS'}

    def get_pcthresh(self,thresh1,thresh2,thresh3,thresh4,pos,sz):

        condc = self.dsd.get_pc(rank =1) [pos:pos+sz] > thresh1
        conds = self.dsd.get_pc(rank =1 )[pos:pos+sz] < (-1.)*thresh2
        conda = self.dsd.get_pc(rank =1) [pos:pos+sz] > -999.

        condi = self.dsd.get_pc(rank =2 )[pos:pos+sz] < (-1.)*thresh4
        condw = self.dsd.get_pc(rank =2 )[pos:pos+sz] > thresh3

        condcc = np.logical_and(condc,condw)
        condci = np.logical_and(condc,condi)
    
        condc11 = np.logical_and(condc,~condcc)
        condc12 = np.logical_and(condc,~condci)
        condcon = np.logical_and(condc11,condc12)

        condshal = np.logical_and(condw,~condcc)
#        condshal = np.logical_and(conds1,~conds)
        condstrat = np.logical_and(conds,~condshal)
    
        condice = np.logical_and(condi,~condci)


        amb1 = np.logical_and(~condc,~conds)
        amb2 = np.logical_and(~condi,~condw)
        amb = np.logical_and(amb1,amb2)

        amb1s = np.logical_and(~condci,~condcc)
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
        return whppc1,whnpc1,whppc2,whnpc2,whnpc2c,whsnpc2,whppc2c,whsppc2,whsppc1,whsnpc1,whall,whamb,whambs


 #        cond1c = self.dsd.get_pc(rank =1) [pos:pos+sz] > -999.
#         cond1s = self.dsd.get_pc(rank =1 )[pos:pos+sz] < (-1.)*thresh2
#         cond2p = self.dsd.get_pc(rank =2 )[pos:pos+sz] > thresh1
#         cond2n = self.dsd.get_pc(rank =2 )[pos:pos+sz] < (-1.)*thresh1
#         conda = self.dsd.get_pc(rank =1) [pos:pos+sz] > -999.
#     
#     
#         whall = np.where(conda)
#         whppc1 = np.where(cond1c)
#         whnpc1 = np.where(cond1s)
#         whppc2 = np.where(cond2p)
#         whnpc2 = np.where(cond2n)
# 
#         wh_neg2c=np.logical_and(self.dsd.get_pc(rank =1)[pos:pos+sz] > thresh1,cond2n)
#         wh_neg2s=np.logical_and(self.dsd.get_pc(rank =1)[pos:pos+sz] <0,cond2n)
# 
#         wh_pos2c=np.logical_and(self.dsd.get_pc(rank =1)[pos:pos+sz] > thresh1,cond2p)
#         wh_pos2s=np.logical_and(self.dsd.get_pc(rank =1)[pos:pos+sz] <0,cond2p)
# 
#         whsppc1t = np.logical_and(cond1c,~wh_neg2c)
#         whsppc1s = np.logical_and(whsppc1t,~wh_pos2c)
#         whsnpc1t = np.logical_and(cond1s,~wh_neg2s)
#         whsnpc1s = np.logical_and(whsnpc1t,~wh_pos2s)
#         
# 
#         whnpc2c=np.where(wh_neg2c)
#         whnpc2s=np.where(wh_neg2s)
#         whppc2c=np.where(wh_pos2c)
#         whppc2s=np.where(wh_pos2s)
#         whsnpc1 =np.where(whsnpc1s)
#         whsppc1 =np.where(whsppc1s)

#        return whppc1,whnpc1,whppc2,whnpc2,whnpc2c,whnpc2s,whppc2c,whppc2s,whsppc1,whsnpc1,whall,whamb,whambs

    def sub_world(self,raindata,dsddata,pos,sz):
        subd = {}
        dsdd = {}
        for v in raindata.keys():
            subd[v] = raindata[v][pos:pos+sz]
        return subd


    def read_raindata(self,sav=1,pdir = 'scratch/PCA_NEW/SIGM_TEST/'):
        print (self.exper)
        if self.exper == 'sgp' or self.exper == 'fin' or self.exper == 'dar' or self.exper == 'man':
            raindata = pickle.load( open( '{s}{ex:}group_raindata.p'.format(ex=self.exper,s=pdir), "rb" ) )
            dsddata = pickle.load( open( '{s}{ex:}group_dsd_dm_lognorm_sigm.p'.format(ex=self.exper,s=pdir),"rb" ) )

            if self.exper == 'dar':
                self.exper = 'darwin'
        elif self.exper == 'lba':
            raindata = pickle.load( open( '{s}{ex:}2dvd_raindata.p'.format(ex=self.exper,s=pdir), "rb" ) )
            dsddata = pickle.load( open( '{s}{ex:}2dvd_dsd_dm_lognorm_sigm.p'.format(ex=self.exper,s=pdir),"rb" ) )
        elif self.exper =='sgpt':
            self.exper = 'sgp'
            raindata = pickle.load( open( '{s}{ex:}_raindata.p'.format(ex=self.exper,s=pdir), "rb" ) )
            dsddata = pickle.load( open( '{s}{ex:}_dsd_dm_lognorm_sigm.p'.format(ex=self.exper,s=pdir),"rb" ) )
        elif self.exper =='olympexaput':
            self.exper = 'olympexaput'
            raindata = pickle.load( open( '{s}olympexapu_raindata_dm_sub030508_sigm.p'.format(ex=self.exper,s=pdir), "rb" ) )
            dsddata = pickle.load( open( '{s}{ex:}_dsd_dm_lognorm_sigm.p'.format(ex=self.exper,s=pdir),"rb" ) )
        
        else:
            raindata = pickle.load( open( '{s}{ex:}_raindata.p'.format(ex=self.exper,s=pdir), "rb" ) )
            dsddata = pickle.load( open( '{s}{ex:}_dsd_dm_lognorm_sigm.p'.format(ex=self.exper,s=pdir),"rb" ) )

        print ('reading files', raindata.keys())

        return raindata,dsddata

#     def get_pc(self):
#         pc1=self.dsddat.get_pc(rank =1)
#         npc1=self.dsddat.get_pc(rank =1)
#         pc2=self.dsddat.get_pc(rank =2)
#         npc2=self.dsddat.get_pc(rank =2)
#         pc3=self.dsddat.get_pc(rank =3)
#         npc3=self.dsddat.get_pc(rank =3)
# 
#         return pc1,npc1,pc2,npc2,pc3,npc3


    def get_eof(self):
        eof1=self.dsddat.get_eof(rank =1)
        eof2=self.dsddat.get_eof(rank =2)
        eof3=self.dsddat.get_eof(rank =3)

        return eof1,eof2,eof3
        
    def linearize_dat(self,dat):
        return np.log10(dat)
    
    def subset(self, thresh1,thresh2,thresh3,thresh4):
        whppc1,whnpc1,whppc2,whnpc2,whnpc2c,whsnpc2,\
            whppc2c,whsppc2,whsppc1,whsnpc1,whall,whamb,whambs=self.get_pcthresh(thresh1,thresh2,thresh3,thresh4,self.pos,self.size)
        indlist = [whppc1,whnpc1,whppc2,whnpc2,whnpc2c,whsnpc2,whppc2c,whsppc2,whall,whamb,whsppc1,whsnpc1,whambs]
        namlist = ['conv','strat','warm','ldlc','ice','hstrat','cc','shallow','all','amb','convs','strats','ambs']
        subrain = {}
        for v in self.rain.keys():
            subrain[v] = {}
            for i in range(len(indlist)):
                subrain[v][namlist[i]] = self.rain[v][indlist[i]]
        return subrain
    
    def get_location(self):
        if self.exper == 'global':
            locs=0
        elif self.exper == 'ifloods':
            locs={"sn35":[-92.3653,42.1822],"sn36":[-92.2817,42.1255],"sn37":[-92.0717,41.9914],"sn38":[-92.8739,41.8603],"sn70":[-91.5417,41.6406]}
        elif self.exper == 'olympex' or self.exper == 'olympexapu':
            locs={"apu03":[-123.993,47.3599],"apu05":[-123.890,47.4596],"apu08":[-123.812,47.5135]}
#            locs={"sn35":[-123.993,47.360],"sn36":[-123.890,47.460],"sn38":[-123.812,479.514]}        elif self.exper == 'sgp':
        elif self.exper == 'sgp':
            locs={"sgp":[-97.485,36.605],"sn25":[-97.53222,36.6233],"sn35":[-97.4797,36.618],"sn36":[-97.4786,35.58138],"sn37":[-97.48083,36.633],"sn38":[-97.444,36.57833]}
#         elif self.exper == 'lpvex':
#             locs={"harmaja":[24.9749,60.1045],"emasalo":[25.6247,60.2037],"jarvenpaa":[25.0822,60.4846]}
#         elif self.exper == 'mc3e':
#             locs={"sn25":[-97.53222,36.6233],"sn35":[-97.4797,36.618],"sn36":[-97.4786,35.58138],"sn37":[-97.48083,36.633],"sn38":[-97.444,36.57833]}
        elif self.exper == 'iphex':
            locs={"sn25":[-82.0565,35.2266],"sn35":[-82.1707,35.2929],"sn36":[-82.3700,35.3728],"sn37":[-82.0948,35.5198],"sn38":[-82.0732,35.5858]}
        elif self.exper == 'lba':
            locs={"jwd":[-61.8489,-10.8748]}
        elif self.exper == 'dar' or self.exper == 'darwin':
            locs={"2dvd":[130.892,-12.425]}
        elif self.exper == 'twpice':
            locs={"jwd":[130.892,-12.425]}
        elif self.exper == 'man' or self.exper == 'MANNU':
            locs=man_locs={"man":[147.425,-2.06],"gan":[73.1501,-0.690368]}
#         elif self.exper == 'gan':
#             locs={}
        elif self.exper == 'fin' or self.exper == 'FIN':
            locs={"2dvd":[24.288,61.843],"harmaja":[24.9749,60.1045],"emasalo":[25.6247,60.2037],"jarvenpaa":[25.0822,60.4846]}
        elif self.exper == 'wallops' or self.exper == 'wallops15':
            locs={"sn25":[-75.4818,37.9448],"sn35":[-75.4668,37.9378],"sn36":[-75.4718,37.9348], "sn37":[-75.4738,37.9298],
              "sn38":[-75.4568,37.9388,],"sn70":[-75.4648,37.9448]}
        elif self.exper == 'oceanrain' or self.exper == 'oceanrainhigh':
            locs={"ODM":[self.raindat['lon'],self.raindat['lat']]}
        else:
            print ('No data for this experimeint')
        return locs
    
    def get_region(self):
        if self.exper == 'global':
            region =3

        elif self.exper == 'oceanrain' or self.exper == 'oceanrainhigh':
            region = 2
        else:
            lat = self.locs[self.locs.keys()[0]][1]
            lon = self.locs[self.locs.keys()[0]][0]
            #print lat
            if lat <=23.:
                region = 0
            elif lat > 23. and lat <= 45.:
                region = 1
            elif lat > 45.:
                region = 2
        return region
