from starter import *
from yt.data_objects.particle_filters import add_particle_filter 
import os

class sfr:
    def __init__(self, num): 
        self.num = num 
        add_particle_filter("formed_stars", function = self.formed_stars, filtered_type = "all", requires=['creation_time'])
        filename = 'DD' + str(self.num).zfill(4) + '/DD' + str(self.num).zfill(4)
        self.ds = yt.load(filename)
        self.ad = self.ds.all_data()
        self.ds.add_particle_filter("formed_stars")
        self.mass = self.ad["formed_stars", "particle_mass"].in_units("Msun")
        self.formation_time = self.ad["formed_stars", "creation_time"].in_units("yr")
        hist, bins = np.histogram(self.formation_time, bins = 1000, range = [0, self.ds.current_time.in_units("yr").to_value()])
        inds = np.digitize(self.formation_time, bins = bins)
        self.time = (bins[:-1] + bins[1:]) / 2
        self.sfr_arr = []
        for j in range(len(self.time)):
            self.sfr_arr.append(self.mass[inds == j + 1].sum() / (bins[j+1] - bins[j]))
        self.sfr_arr[0] = np.nan
        self.sfr_arr = np.asarray(self.sfr_arr)


    def formed_stars(self, pfilter, data):
         filter = data["all", "creation_time"] > 0
         return filter

    def plot(self):
        plt.plot(self.time/1e6, self.sfr_arr)
        plt.xlabel("t [Myr]")
        plt.ylabel(r"SFR [M$_\odot$ yr$yr^{-1}$]")
        dirname = os.getcwd().split('/')[-1]
        plt.savefig("frames/sfr_" + dirname + ".png")

