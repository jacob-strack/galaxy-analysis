from starter import * 
from scipy import stats 
from matplotlib.colors import SymLogNorm 
from matplotlib.colors import LogNorm 

class profile: 
    def __init__(self, name, field, num_bins, x_bin = 'z', cen=[0.5,0.5,0.5], log_y = False): 
        self.name = name 
        self.num_bins = num_bins
        self.field = field 
        self.cen = cen 
        self.log_y = log_y
        self.x_bin = x_bin
        ds = yt.load(name + "/" + name) 
        ad = ds.all_data()
        if(self.x_bin == 'z'):  
            z = ad['z'].in_units('kpc')
            z -= ad.center[2].in_units('kpc') 
        field1 = ad[field]
        if(log_y == True): field1 = np.log10(field1)
        if(self.x_bin == 'z'): 
            self.hist = np.histogram2d(z, field1, bins = self.num_bins) 
        if(self.x_bin == 'r'): 
            r = np.sqrt((ad['x'] - ad.center[0])**2 + (ad['y'] - ad.center[1])**2 + (ad['z'] - ad.center[2])**2).in_units('kpc')
            self.hist = np.histogram2d(r, field1, bins = self.num_bins) 
        self.hist[0][self.hist[0] == 0] = np.nan
    def plot(self): 
        fig, ax = plt.subplots(1,1) 
        a = ax.imshow(self.hist[0].T, origin = 'lower', cmap = 'jet')
        if(self.x_bin == 'z'): 
            ax.set_xlabel("z [kpc]") 
        if(self.x_bin == 'r'): 
            ax.set_xlabel("r [kpc]")
        if(self.log_y == True): 
            ax.set_ylabel("log(" + self.field + ")")
        else:    
            ax.set_ylabel(self.field)
        ax.set_xticks(np.arange(0, len(self.hist[1]), 20),np.round(self.hist[1][0::20].d,0))
        if(self.log_y == True): 
            ax.set_yticks(np.arange(0, len(self.hist[2]), 10),np.round(self.hist[2][0::10],0))
        else:
            ax.set_yticks(np.arange(0, len(self.hist[2]), 10),np.round(self.hist[2][0::10].d,0))
        plt.colorbar(a, ax = ax)
        fig.savefig('frames/' + self.name + "_" + self.field + "_" + self.x_bin + ".png")

class volume_weighted_profile: 
    def __init__(self, name, field, num_bins, x_bin = 'z', cen=[0.5,0.5,0.5], log_y = False): 
        self.name = name 
        self.num_bins = num_bins
        self.field = field 
        self.cen = cen 
        self.log_y = log_y
        self.x_bin = x_bin
        ds = yt.load(name + "/" + name) 
        ad = ds.all_data()
        if(self.x_bin == 'z'):  
            z = ad['z'].in_units('kpc')
            z -= ad.center[2].in_units('kpc') 
        field1 = ad[field]
        if(log_y == True): field1 = np.log10(field1)
        if(self.x_bin == 'z'): 
            self.hist = stats.binned_statistic_2d(z, field1, ad['cell_volume'].in_units('code_length**3'), bins = self.num_bins, statistic = 'sum') 
        if(self.x_bin == 'r'): 
            r = np.sqrt((ad['x'] - ad.center[0])**2 + (ad['y'] - ad.center[1])**2 + (ad['z'] - ad.center[2])**2).in_units('kpc')
            self.hist = stats.binned_statistic_2d(r, field1, ad['cell_volume'].in_units('code_length**3'), bins = self.num_bins, statistic = 'sum') 
        self.hist[0][self.hist[0] == 0] = np.nan

    def plot(self): 
        fig, ax = plt.subplots(1,1) 
        a = ax.imshow(self.hist[0].T, origin = 'lower', cmap = 'jet')
        if(self.x_bin == 'z'): 
            ax.set_xlabel("z [kpc]") 
        if(self.x_bin == 'r'): 
            ax.set_xlabel("r [kpc]")
        if(self.log_y == True): 
            ax.set_ylabel("log(" + self.field + ")")
        else:    
            ax.set_ylabel(self.field)
        ax.set_xticks(np.arange(0, len(self.hist[1]), 20),np.round(self.hist[1][0::20],0))
        if(self.log_y == True): 
            ax.set_yticks(np.arange(0, len(self.hist[2]), 10),np.round(self.hist[2][0::10],0))
        else:
            ax.set_yticks(np.arange(0, len(self.hist[2]), 10),np.round(self.hist[2][0::10],0))
        plt.colorbar(a, ax = ax)
        fig.savefig('frames/' + self.name + "_" + self.field + "_" + self.x_bin + "_volume_weighted.png")

class mass_weighted_profile: 
    def __init__(self, name, field, num_bins, x_bin = 'z', cen=[0.5,0.5,0.5], log_y = False): 
        self.name = name 
        self.num_bins = num_bins
        self.field = field 
        self.cen = cen 
        self.log_y = log_y
        self.x_bin = x_bin
        ds = yt.load(name + "/" + name) 
        ad = ds.all_data()
        if(self.x_bin == 'z'):  
            z = ad['z'].in_units('kpc')
            z -= ad.center[2].in_units('kpc') 
        field1 = ad[field]
        if(log_y == True): field1 = np.log10(field1)
        if(self.x_bin == 'z'): 
            self.hist = stats.binned_statistic_2d(z, field1, ad['cell_mass']/ad['cell_mass'].sum(), bins = self.num_bins, statistic = 'sum') 
        if(self.x_bin == 'r'): 
            r = np.sqrt((ad['x'] - ad.center[0])**2 + (ad['y'] - ad.center[1])**2 + (ad['z'] - ad.center[2])**2).in_units('kpc')
            self.hist = stats.binned_statistic_2d(r, field1, ad['cell_mass']/ad['cell_mass'].sum(), bins = self.num_bins, statistic = 'sum') 
        self.hist[0][self.hist[0] == 0] = np.nan

    def plot(self): 
        fig, ax = plt.subplots(1,1) 
        a = ax.imshow(self.hist[0].T, origin = 'lower', cmap = 'jet')
        if(self.x_bin == 'z'): 
            ax.set_xlabel("z [kpc]") 
        if(self.x_bin == 'r'): 
            ax.set_xlabel("r [kpc]")
        if(self.log_y == True): 
            ax.set_ylabel("log(" + self.field + ")")
        else:    
            ax.set_ylabel(self.field)
        ax.set_xticks(np.arange(0, len(self.hist[1]), 20),np.round(self.hist[1][0::20],0))
        if(self.log_y == True): 
            ax.set_yticks(np.arange(0, len(self.hist[2]), 10),np.round(self.hist[2][0::10],0))
        else:
            ax.set_yticks(np.arange(0, len(self.hist[2]), 10),np.round(self.hist[2][0::10],0))
        plt.colorbar(a, ax = ax)
        fig.savefig('frames/' + self.name + "_" + self.field + "_" + self.x_bin + "_mass_weighted.png")

class profile_compare: 
    def __init__(self, name1, name2, name, field, num_bins, x_bin = 'z', log_y = False, max_x = None, sym_log_y = False, sym_lin_thresh = None, cen = [0.5,0.5,0.5]): 
        self.name = name
        self.name1 = name1
        self.name2 = name2 
        self.field = field
        self.num_bins = num_bins
        self.x_bin = x_bin
        self.log_y = log_y
        ds1 = yt.load(name1 + "/" + name + "/" + name)
        self.time = np.round(ds1.current_time.in_units('Myr'),1)
        ds2 = yt.load(name2 + "/" + name + "/" + name)
        ad1 = ds1.all_data()
        ad2 = ds2.all_data()
        if(self.x_bin == 'z'):  
            z1 = ad1['z'].in_units('kpc')
            z1 -= yt.YTQuantity(cen[2], ds1.units.code_length).in_units('kpc') 
            z2 = ad2['z'].in_units('kpc')
            z2 -= yt.YTQuantity(cen[2], ds2.units.code_length).in_units('kpc')
        field1 = ad1[field]
        self.units = str(field1.units)
        field2 = ad2[field]
        field_min = np.minimum(field1.min(), field2.min())
        field_max = np.maximum(field1.max(), field2.max())
        if(log_y == True):
            if(sym_log_y == True): 
                field1[field1 < -1*sym_lin_thresh] = -1*np.log10(np.abs(field1[field1 < -1*sym_lin_thresh]))
                field2[field2 < -1*sym_lin_thresh] = -1*np.log10(np.abs(field2[field2 < -1*sym_lin_thresh]))
                field1[field1 > sym_lin_thresh] = np.log10(field1[field1 > sym_lin_thresh])
                field2[field2 > sym_lin_thresh] = np.log10(field2[field2 > sym_lin_thresh])
            else:
                field1 = np.log10(field1)
                field2 = np.log10(field2)
        mask1 = ~np.isinf(field1)
        mask2 = ~np.isinf(field2)
        field1 = field1[mask1]
        field2 = field2[mask2]
        if(log_y == True): 
            field_min = np.minimum(field1.min(), field2.min())
            field_max = np.maximum(field1.max(), field2.max())
        if(self.x_bin == 'z'):
            z1 = z1[mask1]
            z2 = z2[mask2]
            if(max_x == None): 
                self.hist1 = stats.binned_statistic_2d(z1, field1, ad1["ones"][mask1], bins = self.num_bins,range = [[z1.min(),z1.max()],[field_min,field_max]],statistic = 'count') 
                self.hist2 = stats.binned_statistic_2d(z2, field2, ad2["ones"][mask2], bins = self.num_bins,range = [[z2.min(),z2.max()],[field_min,field_max]],statistic = 'count') 
            else:
                self.hist1 = stats.binned_statistic_2d(z1, field1, ad1["ones"][mask1], bins = self.num_bins,range = [[-max_x,max_x],[field_min, field_max]] ,statistic = 'count') 
                self.hist2 = stats.binned_statistic_2d(z2, field2, ad2["ones"][mask2], bins = self.num_bins,range = [[-max_x,max_x],[field_min,field_max]] ,statistic = 'count') 
        if(self.x_bin == 'r'): 
            r1 = np.sqrt((ad1['x'].in_units('kpc') - yt.YTQuantity(cen[0], ds1.units.code_length).in_units('kpc'))**2 + (ad1['y'].in_units('kpc') - yt.YTQuantity(cen[1], ds1.units.code_length).in_units('kpc'))**2 + (ad1['z'].in_units('kpc') - yt.YTQuantity(cen[2], ds1.units.code_length).in_units('kpc'))**2)
            r2 = np.sqrt((ad2['x'].in_units('kpc') - yt.YTQuantity(cen[0], ds2.units.code_length).in_units('kpc'))**2 + (ad2['y'].in_units('kpc') - yt.YTQuantity(cen[1], ds2.units.code_length).in_units('kpc'))**2 + (ad2['z'].in_units('kpc') - yt.YTQuantity(cen[2], ds2.units.code_length).in_units('kpc'))**2)
            r1 = r1[mask1]
            r2 = r2[mask2]
            if(max_x == None): 
                self.hist1 = stats.binned_statistic_2d(r1, field1, ad1['ones'][mask1], bins = self.num_bins, range=[[r1.min(),r1.max()],[field_min,field_max]],statistic = 'count') 
                self.hist2 = stats.binned_statistic_2d(r2, field2, ad2['ones'][mask2], bins = self.num_bins, range=[[r2.min(),r2.max()],[field_min,field_max]],statistic = 'count') 
            else:
                self.hist1 = stats.binned_statistic_2d(r1, field1, ad1['ones'][mask1], bins = self.num_bins,range = [[0, max_x],[field_min,field_max]], statistic = 'count') 
                self.hist2 = stats.binned_statistic_2d(r2, field2, ad2['ones'][mask2], bins = self.num_bins, range = [[0,max_x],[field_min, field_max]] ,statistic = 'count') 
        self.hist1[0][self.hist1[0] == 0] = np.nan
        self.hist2[0][self.hist2[0] == 0] = np.nan


    def plot(self): 
        fig, ax = plt.subplots(1,2)
        fig.set_size_inches(12,6)
        fig.suptitle('t = ' + str(self.time) + ' Myr')
        a = ax[0].imshow(self.hist1[0].T, origin = 'lower', cmap = 'jet')
        b = ax[1].imshow(self.hist2[0].T, origin = 'lower', cmap = 'jet')
        if(self.x_bin == 'z'): 
            ax[0].set_xlabel("z [kpc]") 
            ax[1].set_xlabel("z [kpc]") 
        if(self.x_bin == 'r'): 
            ax[0].set_xlabel("r [kpc]")
            ax[1].set_xlabel("r [kpc]")
        if(self.log_y == True): 
            ax[0].set_ylabel("log(" + self.field + ") [" + self.units +']' )
            ax[1].set_ylabel("log(" + self.field + ") [" + self.units +']' )
        else:    
            ax[0].set_ylabel(self.field + ' [' + self.units + ']')
            ax[1].set_ylabel(self.field + ' [' + self.units + ']')
        ax[0].set_xticks(np.arange(0, len(self.hist1[1]), 20),np.round(self.hist1[1][0::20],0))
        ax[1].set_xticks(np.arange(0, len(self.hist2[1]), 20),np.round(self.hist2[1][0::20],0))
        if(self.log_y == True): 
            ax[0].set_yticks(np.arange(0, len(self.hist1[2]), 10),np.round(self.hist1[2][0::10],0))
            ax[1].set_yticks(np.arange(0, len(self.hist2[2]), 10),np.round(self.hist2[2][0::10],0))
        else:
            ax[0].set_yticks(np.arange(0, len(self.hist1[2]), 10),np.round(self.hist1[2][0::10].d,0))
            ax[1].set_yticks(np.arange(0, len(self.hist2[2]), 10),np.round(self.hist2[2][0::10].d,0))
        ax[0].set_title(self.name1)
        ax[1].set_title(self.name2)
        plt.colorbar(a, ax = ax[0], shrink = 0.5)
        plt.colorbar(b, ax = ax[1], shrink = 0.5)
        fig.savefig(self.name1 + '/frames/' + self.name + "_" + self.field + "_" + self.x_bin + "_comp.png")
