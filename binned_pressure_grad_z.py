from starter import * 
from pressure_gradient_periodic import *
from scipy import stats 
from matplotlib.colors import SymLogNorm
class binnedpressurez: 
    def __init__(self, num): 
         self.num = num #number of dumps
         self.arr = []
         self.t_arr = []
         self.arr2 = []
         self.t_arr2 = []
         self.edges = 0
         for i in range (0, num):
             ds_name = 'DD' + str(i).zfill(4)#filename used in loading dataset
             ds_2_name = '../7_level_tweaked_halo/DD' + str(i).zfill(4) + '/DD' + str(i).zfill(4)
             filename = ds_name + '/' + ds_name 
             ds = yt.load(filename) #loading dataset
             ds2 = yt.load(ds_2_name)
             add_pressure_gradient(ds)
             add_pressure_gradient(ds2)
             sph = ds.sphere([0.5,0.5,0.5],(25, 'kpc')) #sphere chunk from frame
             sph2 = ds2.sphere([0.5,0.5,0.5], (25, 'kpc'))
             z = sph['z'].in_units('kpc') #array of z values on the sphere
             z2 = sph2['z'].in_units('kpc')
             z -= sph.center[2].in_units('kpc')
             z2 -= sph2.center[2].in_units('kpc')
             field1 = sph[yt_dpdz]
             field2 = sph2[yt_dpdz]
             self.field = 'dpdz'
             self.bins = stats.binned_statistic(z, field1, bins = 17)
             self.bins2 = stats.binned_statistic(z2, field2, bins = 17)
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
        a = ax[0].imshow(self.arr.T, origin = 'lower', cmap = 'jet', norm = SymLogNorm(linthresh = 10**-25, vmin = vmin,vmax = vmax,base = 10))
        plt.colorbar(a, ax = ax[0])
        ax[0].set_xlabel('t [Myr]')
        ax[0].set_ylabel('z [kpc]')
        ax[0].title.set_text('<' + self.field + '>(z,t)')
        ax[0].set_xticks(np.linspace(0, len(self.t_arr) - 1, 5))
        ax[0].set_yticks(np.linspace(0, len(self.edges) - 2, 5))
        ax[0].set_xticklabels([str((int(np.linspace(self.t_arr[0], self.t_arr[-1], 5)[i]) // 10) * 10) for i in range(5)])
        ax[0].set_yticklabels([str(int((np.linspace(self.edges[0], self.edges[-1], 5)[i]))) for i in range(5)])
        b = ax[1].imshow(self.arr2.T, origin = 'lower', cmap = 'jet', norm = SymLogNorm(linthresh = 10**-25, vmin = vmin, vmax = vmax, base = 10))
        plt.colorbar(b, ax = ax[1])
        ax[1].set_xlabel('t [Myr]')
        ax[1].set_ylabel('z [kpc]')
        ax[1].title.set_text('<' + self.field + '>(z,t) Halo')
        ax[1].set_xticks(np.linspace(0, len(self.t_arr2) - 1, 5))
        ax[1].set_yticks(np.linspace(0, len(self.edges2) - 2, 5))
        ax[1].set_xticklabels([str((int(np.linspace(self.t_arr2[0], self.t_arr2[-1], 5)[i]) // 10) * 10) for i in range(5)])
        ax[1].set_yticklabels([str(int((np.linspace(self.edges2[0], self.edges2[-1], 5)[i]))) for i in range(5)])
        fig.savefig(self.field + 'comp.png')	
         
    def single_plot(self, plot_num, vmin, vmax): 
        if plot_num == 0: 
            arr = self.arr.T
            t_arr = self.t_arr
            edges = self.edges
        elif plot_num == 1: 
            arr = self.arr2.T
            t_arr = self.t_arr2
            edges = self.edges2
        else:
            print("no \n")
            return
        ax = plt.subplot() 
        im = plt.imshow(arr, origin = 'lower', cmap = 'jet', norm = SymLogNorm(linthresh = 10**-25, vmin = vmin,vmax = vmax,base = 10))
        ax.set_xlabel('t [Myr]')
        ax.set_ylabel('z [kpc]')
        ax.title.set_text('<' + self.field + '>(z,t)')
        ax.set_xticks(np.linspace(0, len(t_arr) - 1, 5),[str((int(np.linspace(t_arr[0], t_arr[-1], 5)[i]) // 10) * 10) for i in range(5)])
        ax.set_yticks(np.linspace(0, len(edges) - 2, 5), [str(int((np.linspace(edges[0], edges[-1], 5)[i]))) for i in range(5)])
        plt.colorbar(im, ax = ax)
        plt.savefig(self.field + "single.png")
