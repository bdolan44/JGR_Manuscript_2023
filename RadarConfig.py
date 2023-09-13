# this houses the RadarConfig object which is basically just carrying
# along some options about names, colorbars, plotting limits, etc.

from __future__ import division


from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
import sys

plt.switch_backend('agg')
from copy import deepcopy
import matplotlib as mpl
from cmap import make_cmap
from ctables import *
from Config import Config # yaml stuff


class RadarConfig(object):
    
    def __init__(self, config_file=None, dz='DZ', zdr='DR', kdp='KD', ldr='LH', rho='RH', ph='FDP', hid='HID', vel='VR',rr='RR',wi='WI',
            d0='D0',temp='T', x='x', y='y', z='z', u='U', v='V', w='Wvar', lat='lat', lon='lon', cs='CS'):

                 # ******** first the polarimetric stuff *************
                 # establish some defaults
        self.dz_name = dz # reflecitiy name
        self.zdr_name = zdr # differential reflectivity name
        self.kdp_name = kdp # specific differential phase name
        self.ldr_name = ldr # linear depolarization ratio name
        self.rho_name = rho # rho hv/correlation coefficient name
        self.ph_name = ph # differential phase name
        self.rr_name = rr
        self.vel_name = vel # radial velocity name
        self.temp_name = temp
        self.hid_name = hid
        self.cs_name = cs
        self.wi_name = wi
        self.d0_name = d0

        # ********** next the cartesian coordinates **********
        self.x_name = x
        self.y_name = y
        self.z_name = z
        self.lat_name = lat
        self.lon_name = lon
    
        # ********** dual doppler names **********************
        self.u_name = u
        self.v_name = v
        self.w_name = w

        # now check for a file that can overwrite stuff

        if config_file is not None:
            self.read_config_file(config_file) # read in the config file




        self.species = np.array(['DZ','RN','CR','AG','WS','VI','LDG','HDG','HA','BD'])
        self.hid_colors = ['White', 'LightBlue', 'MediumBlue', 'Darkorange', 'LightPink', 'Cyan', 'DarkGray',\
            'Lime', 'Yellow', 'Red', 'DarkRed']

        self.cs_colors = ['#FFFFFF', 'DodgerBlue', 'Red', 'Khaki']
        self.cs_labels = ['', 'Strat', 'Conv', 'Mixed']


        self.pol_vars = np.array([self.dz_name, self.zdr_name, self.kdp_name, self.ldr_name, self.rho_name, self.hid_name])

    
        self.set_dbz_colorbar()
        self.set_hid_colorbar()
        self.set_cs_colorbar()
        self.set_zdr_colorbar()
        self.set_kdp_colorbar()
        self.set_phase_colorbar()


    #### Use your colormap


        #print 'self.plot params'

        self.plot_params = {
            self.dz_name: {'lims': [0,80], 'delta': 10, 'units': '(dBZ)', 'ticklabels': np.arange(0, 90, 10),
            'name': 'Z', 'longname': 'Reflectivity', 'cmap': self.temp_cmap, 'norm': self.normdbz,'ncname':'corrected_reflectivity'},

            self.zdr_name: {'lims': [-1, 4], 'delta': 1, 'units': '(dB)', 'ticklabels': np.arange(-1, 5, 1),
            'name': 'Z$_{DR}$', 'longname': 'Corrected Differential reflectivity', 'cmap': self.zdr_cmap, 
            'norm': self.zdr_norm,'ncname':'corrected_differential_reflectivity'},

            self.d0_name: {'lims': [-1, 4], 'delta': 1, 'units': '(mm)', 'ticklabels': np.arange(-1, 5, 1),
            'name': 'D$_{0}$', 'longname': 'D0', 'cmap': self.zdr_cmap, 
            'norm': self.zdr_norm,'ncname':'mean_drop_diameter'},

        #                'name': 'Z$_{DR}$', 'longname': 'Differntial reflectivity', 'cmap': plt.cm.jet, 
        #                'norm': colors.BoundaryNorm(np.arange(-1, 5, 0.5), plt.cm.jet.N)},

            self.kdp_name: {'lims': [-0.5,3], 'delta': 1, 'units': '(deg/km)', 'ticklabels': np.arange(-0.5, 4, 0.5),
            'name': 'K$_{dp}$', 'longname': 'Specific differential phase', 'cmap': self.kdp_cmap, 
            'norm': self.kdp_norm,'ncname':'specific_differential_phase'},

            self.rho_name: {'lims': [0.95,1.0], 'delta': 0.01, 'units': '', 'ticklabels': np.arange(0.95, 1.01, 0.01),
            'name': r'$\rho_{hv}$', 'longname': 'cross_correlation_coefficient_at_zero_lag', 'cmap': plt.cm.jet, 
            'norm': colors.Normalize(vmin=0.95, vmax=1.0),'ncname':'correlation_coefficient'},

            self.hid_name: {'lims': [0, len(self.species)+1], 'delta': 1, 'units': '', 
            'ticklabels': np.append('', self.species), 'name': '', 'longname': 'Hydrometeor classification', 
            'cmap': self.hid_cmap, 'norm': self.normhid,'ncname':'hydrometeor_identification'},

            self.cs_name: {'lims': [0, 4], 'delta': 1, 'units': '', 
            'ticklabels': self.cs_labels, 'name': '', 'longname': 'Convective/Stratiform', 
            'cmap': self.cs_cmap, 'norm': None,'ncname':'convective_stratiform_classification'},

            self.ph_name: {'lims': [-180, 180], 'delta': 30, 'units': '($^{\circ}$)', 
            'ticklabels': np.arange(-180, 210, 30), 'name': 'Phase', 'longname': 'Differential phase', 
            'cmap': self.ph_cmap, 'norm': self.ph_norm,'ncname':'differential_phase'},

            self.vel_name: {'lims': [-20, 20], 'delta': 5, 'units': '(m/s)', 
            'ticklabels': np.arange(-20, 25, 5), 'name': 'Velocity', 'longname': 'Radial Velocity', 
            'cmap': Carbone11, 'norm': colors.Normalize(vmin=-20, vmax=20),'ncname':'unfolded_radial_velocity'},

            self.rr_name: {'lims': [0, 120], 'delta': 5, 'units': '(mm/hr)', 
            'ticklabels': np.arange(0, 125, 5), 'name': 'RainRate', 'longname': 'Polarimetric Rain Rate', 
            'cmap': self.kdp_cmap, 'norm': colors.Normalize(vmin=0, vmax=120),'ncname':'rain_rate'},


            self.wi_name: {'lims': [0, 10], 'delta': 0.1, 'units': '(m/s)', 
            'ticklabels': np.arange(0, 10.1, 0.1), 'name': 'Spectrum Width', 'longname': 'Spectrum Width', 
            'cmap': Carbone17, 'norm': colors.Normalize(vmin=0, vmax=10),'ncname':'spectrum_width'},


            self.ldr_name: {'lims': [-35, -20], 'delta': 5, 'units': '(dB)', 
            'ticklabels': np.arange(-35, -20, 5), 'name': 'LDR', 'longname': 'Linear depolarization ratio', 
            'cmap': plt.cm.gist_rainbow_r, 'norm': colors.Normalize(vmin=-35, vmax=-20),'ncname':'linear_depolarization_ratio'},
                        }

        # # Now just set some defaults
        # self.lims = {self.dz_name: [0,80], self.zdr_name: [-1, 3.5], self.kdp_name: [-0.5, 3], self.ldr_name: [-35, -20], 
        #         self.rho_name: [0.95, 1.00], self.hid_name: [0, len(self.species)+1], 
        #                 self.cs_name: [0,4], self.ph_name: [-150, -80], self.vel_name: [-20, 20]}


        # self.delta = {self.dz_name: 10, self.zdr_name: 1, self.kdp_name: 1, self.ldr_name: 5, self.rho_name: 0.01, 
        #             self.hid_name: 1, self.cs_name: 1, self.ph_name: 30, self.vel_name: 5}


        # self.units = {self.dz_name: '(dBZ)', self.zdr_name: '(dB)', self.kdp_name: '($^{\circ}$/km)', self.ldr_name: '(dB)', 
        #                 self.rho_name: '', self.hid_name: '', 
        #                 self.cs_name: '', self.ph_name: '($^{\circ}$)', self.vel_name: '(m/s)'}


        # self.names = {self.dz_name: 'Z', self.zdr_name: 'Z$_{DR}$', self.kdp_name: 'K$_{dp}$', self.ldr_name: 'LDR', 
        #             self.rho_name: r'$\rho_{hv}$', self.hid_name: '', 
        #             self.cs_name: '', self.ph_name: 'Phase', self.vel_name: 'Velocity'}


        # self.longnames = {self.dz_name: 'Reflectivity', self.zdr_name: 'Differntial reflectivity', 
        #             self.kdp_name: 'Specific differential phase',\
        #         self.ldr_name: 'Linear depolarization ratio', self.rho_name: 'Correlation coefficient', 
        #         self.hid_name: 'Hydrometeor identification', 
        #         self.cs_name: 'Convective/Stratiform', self.ph_name: 'Differential phase', self.vel_name: 'Radial Velocity'}


        # self.cmaps = {self.dz_name: self.temp_cmap, self.zdr_name: zdr_cmap, self.kdp_name: kdp_cmap, 
        #         self.ldr_name: plt.cm.gist_rainbow_r, self.rho_name: plt.cm.jet, \
        #             self.hid_name: self.hid_cmap, self.cs_name: self.cs_cmap, self.ph_name: phase_cmap, 
        #             self.vel_name: ctables.Carbone11}


        # self.norms = {self.dz_name: self.normdbz, self.zdr_name: colors.Normalize(vmin=-1, vmax=4), 
        #             self.ldr_name: colors.Normalize(vmin=-35, vmax=-20), 
        #                 self.rho_name: colors.Normalize(vmin=0.95, vmax=1.0), self.hid_name: self.normhid, 
        #                 self.cs_name: None, 
        #                 self.ph_name: colors.Normalize(vmin=-180, vmax=180), 
        #                 self.vel_name: colors.Normalize(vmin=-20, vmax=20), 
        #                 self.kdp_name: colors.Normalize(vmin=-0.5, vmax=3.0)}


        # self.ticklabels = {self.dz_name: np.arange(0, 90, 10), self.zdr_name: np.arange(-1, 5, 1), self.kdp_name: np.arange(-0.5, 4.5, 1), 
        #         self.ldr_name: np.arange(-35, -15, 5), self.rho_name: np.arange(0.95, 1.01, 0.01), 
        #         self.hid_name: np.append('', self.species), 
        #         self.cs_name: self.cs_labels, self.ph_name: np.arange(-180, 210, 30), self.vel_name: np.arange(-20, 25, 5)}

    

#############################################################################################################

    def read_config_file(self, filename):
        # this will read in the yaml file and make the proper attributions
        self.cfg = Config(filename)

        for k in self.cfg.v.keys():
            setattr(self, k, self.cfg.v[k])

        pass


    def get_hid_indices(self, ind):
        # this can take in a list, or array of species names and convert back to the indices
        # which are actually used for the functions
        out = []
        ind = np.asarray(ind)
        for i in ind:
            i_ind = np.where(self.species==i)[0]
            if len(i_ind) == 1:
                out.append(i_ind[0]+1)

        return out


    def set_dbz_colorbar(self, color_list=None):
        if color_list is None:
            # just use the default here
            radarcbar = ['PeachPuff', 'Aqua', 'DodgerBlue', 'MediumBlue', 'Lime', \
                'LimeGreen', 'Green', 'Yellow', 'Orange', 'OrangeRed', \
                'Red', 'Crimson', 'Fuchsia', 'Indigo', 'DarkCyan', 'White']
        else:
            radarcbar = deepcopy(color_list)

        temp_cmap = colors.ListedColormap(radarcbar)
        self.temp_cmap = temp_cmap 
        self.boundsdbz = np.arange(0, 85, 5)
        self.normdbz = colors.BoundaryNorm(self.boundsdbz, self.temp_cmap.N)



#############################################################################################################


    def set_hid_colorbar(self, color_list=None):

        if color_list is None:
             hidcbar = deepcopy(self.hid_colors)

        else:
            hidcbar = deepcopy(color_list)
        self.hid_cmap = colors.ListedColormap(hidcbar)

        self.boundshid = np.arange(0,12)
        self.normhid = colors.BoundaryNorm(self.boundshid, self.hid_cmap.N)
#############################################################################################################

    def set_cs_colorbar(self, color_list=None):
        if color_list is None:
            cscbar = deepcopy(self.cs_colors)

        else:
            cscbar = deepcopy(color_list)

        self.cs_cmap = colors.ListedColormap(cscbar)
        self.cs_bounds = np.arange(0,5)
        self.cs_norm = colors.BoundaryNorm(self.cs_bounds, self.cs_cmap.N)

    def set_kdp_colorbar(self, color_list=None):

        kdp_ncolors = 8
        kdp_oldcmap = plt.cm.get_cmap("hot_r", kdp_ncolors) #generate a jet map with 10 values
        kdp_old_vals = kdp_oldcmap(np.arange(kdp_ncolors)) #extract those values as an array
        kdp_old_vals[0] = [0.65, 0.65, 0.65, 1]
        kdp_old_vals[1] = [0.92, 0.92, 0.6, 1]
        kdp_old_vals = np.insert(kdp_old_vals, 0, [0.2, 0.76, 0.86, 1], axis=0) #change the first value
        kdp_colors = np.delete(kdp_old_vals, -1, axis=0)
        kdp_colors = np.delete(kdp_colors, 4, axis=0)

    #    print kdp_colors


        # Maybe do the interpolation myself I guess??
        reds = np.array([_[0] for _ in kdp_colors])
        greens = np.array([_[1] for _ in kdp_colors])
        blues = np.array([_[2] for _ in kdp_colors])

        kdp_values = np.arange(-1, 6.25, 0.25)
        #print 'len of kdp vals: {}'.format(len(kdp_values))
        kdp_ncolors = len(kdp_values)

        kdp_colors_interp = []
        for i in range(kdp_ncolors):
            frac = (i+1)/kdp_ncolors
            #print frac
            interp_red = int(255*np.interp(frac, np.arange(len(kdp_colors))/len(kdp_colors), reds))
            interp_green = int(255*np.interp(frac, np.arange(len(kdp_colors))/len(kdp_colors), greens))
            interp_blue = int(255*np.interp(frac, np.arange(len(kdp_colors))/len(kdp_colors), blues))

            kdp_colors_interp.append((interp_red, interp_green, interp_blue))

        #print kdp_colors_interp

        # Could I just convert these to hex keys?
        kdp_colors_hex = ['#%02x%02x%02x'%_ for _ in kdp_colors_interp]
        
        #print "Kdp hex colors", kdp_colors_hex

    #################
        self.kdp_cmap = colors.ListedColormap(kdp_colors_hex)

        # define the bins and normalize
        self.kdp_norm = mpl.colors.BoundaryNorm(kdp_values, self.kdp_cmap.N)


    def set_zdr_colorbar(self, color_list=None):

        zdr_colors = [ (104, 243, 79), (134, 223, 120), 
                (154, 203, 150), (186, 186, 186), 
                (240, 240, 120), (255, 237, 60),
                (240, 180, 30), (255, 100, 0),
                (226, 42, 16), (214, 40, 116), 
                (224, 0, 216), (100, 0, 150), (100, 0, 250)]       

        # Maybe do the interpolation myself I guess??
        reds = np.array([_[0] for _ in zdr_colors])
        greens = np.array([_[1] for _ in zdr_colors])
        blues = np.array([_[2] for _ in zdr_colors])

        zdr_values = np.arange(-1, 4.2, 0.2)
        #print 'len of zdr vals: {}'.format(len(zdr_values))
        zdr_ncolors = len(zdr_values)

        zdr_colors_interp = []
        for i in range(zdr_ncolors):
            frac = (i+1)/zdr_ncolors
            #print frac
            interp_red = int(np.interp(frac, np.arange(len(zdr_colors))/len(zdr_colors), reds))
            interp_green = int(np.interp(frac, np.arange(len(zdr_colors))/len(zdr_colors), greens))
            interp_blue = int(np.interp(frac, np.arange(len(zdr_colors))/len(zdr_colors), blues))

            zdr_colors_interp.append((interp_red, interp_green, interp_blue))
        #print(zdr_colors_interp)
        #print('in zdr color hell')
        # Could I just convert these to hex keys?
#        zdr_colors_hex = ['#%02x%02x%02x'.format(i,i,i) for i in zdr_colors_interp]
####BD Changed this due to error reading 'ValueError: Invalid RGBA argument: '#%02x%02x%02x'
        zdr_colors_hex = ['#%02x%02x%02x'%_ for _ in zdr_colors_interp]
        #zdr_colors_hex = ['#%02x%02x%02x'%_ for _ in zdr_colors_interp]
        #print "zdr_colors_hex",zdr_colors_hex
        #print 'ZDR colors 4-6: {}'.format(zdr_colors_hex[4:7])
        zdr_colors_hex[4] = '#9a9a9a'
        zdr_colors_hex[5] = '#bababa'
        zdr_colors_hex[6] = '#dfdfca'
        #print 'ZDR colors 4-6: {}'.format(zdr_colors_hex[4:7])
# 

        #################
        self.zdr_cmap = colors.ListedColormap(zdr_colors_hex)

        # define the bins and normalize
        bounds = np.arange(-1, 4.2, 0.2)
        self.zdr_norm = mpl.colors.BoundaryNorm(bounds, self.zdr_cmap.N)


    def set_phase_colorbar(self, color_list=None):

        # phase colormap, which is a modified version of the existing RdYlBu colormap
        ph_ncolors = 5
        ph_oldcmap = plt.cm.get_cmap("RdYlBu_r", ph_ncolors) #generate a jet map with 10 values
        ph_old_vals = ph_oldcmap(np.arange(ph_ncolors)) #extract those values as an array
        #print 'phase green: {}'.format(ph_old_vals[1])
        ph_old_vals[1] = [0.2, 0.76, 0.86, 1]
        ph_old_vals[2] = [1.0, 1.0, 0.3, 1]

        #phase_cmap = mpl.colors.LinearSegmentedColormap.from_list("newphase", ph_old_vals)

        # Maybe do the interpolation myself I guess??
        reds = np.array([_[0] for _ in ph_old_vals])
        greens = np.array([_[1] for _ in ph_old_vals])
        blues = np.array([_[2] for _ in ph_old_vals])

        ph_values = np.arange(30, 120, 2)
        #print 'len of phase vals: {}'.format(len(ph_values))
        ph_ncolors = len(ph_values)

        ph_colors_interp = []
        for i in range(ph_ncolors):
            frac = (i+1)/ph_ncolors
            #print frac
            interp_red = int(255*np.interp(frac, np.arange(len(ph_old_vals))/len(ph_old_vals), reds))
            interp_green = int(255*np.interp(frac, np.arange(len(ph_old_vals))/len(ph_old_vals), greens))
            interp_blue = int(255*np.interp(frac, np.arange(len(ph_old_vals))/len(ph_old_vals), blues))

            ph_colors_interp.append((interp_red, interp_green, interp_blue))

        #print kdp_colors_interp

        # Could I just convert these to hex keys?
        ph_colors_hex = ['#%02x%02x%02x'%_ for _ in ph_colors_interp]
        #print ph_colors_hex

    #################
        self.ph_cmap = colors.ListedColormap(ph_colors_hex)

        # define the bins and normalize
        self.ph_norm = mpl.colors.BoundaryNorm(ph_values, self.ph_cmap.N)


    def HID_barplot_colorbar(self, figure, location = [0.9, 0.1, 0.03, 0.8]):

        scalarMap = plt.cm.ScalarMappable(norm=self.normhid,cmap=self.hid_cmap)
        axcb = figure.add_axes(location) # x pos, y pos, x width, y width
        cb = mpl.colorbar.ColorbarBase(axcb, cmap=self.hid_cmap, norm=self.normhid, boundaries=self.boundshid,\
                        orientation = 'vertical')
        cb.set_ticks(np.arange(0,11))
            # need to add a blank at the beginning of species to align labels correctly
        labs = np.concatenate((np.array(['']), np.array(self.species)))
        cb.set_ticklabels(labs)

        return cb


    def set_plot_param(self, var, key, value):
        "Reset limits of a variable for plotting, feed it a 2 element list"
        self.plot_params[var][key] = value
#############################################################################################################

    def get_plot_param(self, var, key):
        return self.plot_params[var][key]


    def add_plot_param(self, indict):
        self.plot_params[var] = indict
