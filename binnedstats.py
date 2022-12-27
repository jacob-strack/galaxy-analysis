from starter import * 
from scipy import stats
from matplotlib.colors import LogNorm
from decimal import Decimal
from mpl_toolkits.axes_grid1 import ImageGrid
#generalizing binnedstats to any field in dataset
class binnedstats: 
    def __init__(self, num, field):    
        self.num = num #number of dumps
        self.arr = []
        self.t_arr = []
        self.arr2 = []
        self.t_arr2 = []
        self.edges = 0
        for i in range (0, num):
            self.field = field
            ds_name = 'DD' + str(i).zfill(4)#filename used in loading dataset
            ds_2_name = '../7_level_tweaked_halo/DD' + str(i).zfill(4) + '/DD' + str(i).zfill(4)
            filename = ds_name + '/' + ds_name 
            ds = yt.load(filename) #loading dataset
            ds2 = yt.load(ds_2_name)
            sph = ds.sphere([0.5,0.5,0.5],(25, 'kpc')) #sphere chunk from frame
            sph2 = ds2.sphere([0.5,0.5,0.5], (25, 'kpc'))
            z = sph['z'].in_units('kpc') #array of z values on the sphere
            z2 = sph2['z'].in_units('kpc')
            z -= sph.center[2].in_units('kpc')
            z2 -= sph2.center[2].in_units('kpc')
            field1 = sph[field]
            field2 = sph2[field]
            self.bins = stats.binned_statistic(z, field1, bins = 16)
            self.bins2 = stats.binned_statistic(z2, field2, bins = 16)
            self.arr.append(self.bins[0])
            self.arr2.append(self.bins2[0])
            self.edges = self.bins[1]
            self.edges2 = self.bins2[1]
            self.t_arr.append(ds.current_time)
            self.t_arr2.append(ds2.current_time)
        self.arr = np.asarray(self.arr)
        self.arr2 = np.array(self.arr2)
        for i in range (0,len(self.edges)):
            self.edges[i] = np.format_float_scientific(float(self.edges[i]), 3)
        for i in range (0,len(self.edges2)):
            self.edges2[i] = np.format_float_scientific(float(self.edges2[i]), 3)
        self.t_arr = np.asarray(self.t_arr)
        self.t_arr2 = np.asarray(self.t_arr2)

    def plot(self, vmin, vmax):
        fig, ax = plt.subplots(2,1)
        fig.set_size_inches(12, 10)
        a = ax[0].imshow(self.arr.T, origin = 'lower', cmap = 'jet', norm=LogNorm(vmin = vmin, vmax = vmax))
        plt.colorbar(a, ax = ax[0])
        ax[0].set_xlabel('t [Myr]')
        ax[0].set_ylabel('z [kpc]')
        ax[0].title.set_text('<' + self.field + '>(z,t)')
        ax[0].set_xticks(np.arange(0, len(self.t_arr), 5),np.round(self.t_arr[0::5],2))
        ax[0].set_yticks(np.arange(0, len(self.edges), 5), self.edges[0::5])
        b = ax[1].imshow(self.arr2.T, origin = 'lower', cmap = 'jet',norm=LogNorm(vmin = vmin, vmax = vmax))
        plt.colorbar(b, ax = ax[1])
        ax[1].set_xlabel('t [Myr]')
        ax[1].set_ylabel('z [kpc]')
        ax[1].title.set_text('<' + self.field + '>(z,t) Halo')
        ax[1].set_xticks(np.arange(0, len(self.t_arr2), 5),np.round(self.t_arr2[0::5],2))
        ax[1].set_yticks(np.arange(0, len(self.edges2), 5), self.edges2[0::5])
        fig.savefig(self.field + 'comp.png')	
    
    def zprofile(self, timestep): 
        plt.close('all') 
        y_1 = self.arr[timestep, :]
        x_1 = np.linspace(self.edges[0], self.edges[-1], len(y_1))
        y_2 = self.arr2[timestep, :]
        plt.semilogy(x_1,y_1, label='<' + self.field + '>')
        x_2 = np.linspace(self.edges2[0], self.edges2[-1], (len(y_2)))
        plt.semilogy(x_2,y_2, label='<' + self.field + '> Halo')
        plt.xlabel('z [kpc]')
        plt.ylabel(self.field)
        plt.legend()
        plt.savefig("frames/zprofile" + self.field + ".png")


