from starter import * 
from sfrvstime import sfr
import os

class halosfrcomp: 
    def __init__(self, num, directory): 
        self.no_halo = sfr(num)
        thisdir = os.getcwd()
        os.chdir("../" + directory)
        self.halo = sfr(num)
        os.chdir(thisdir)

    def plot(self): 
        plt.plot(self.no_halo.time / 1e6, self.no_halo.sfr_arr, label = 'no halo')
        plt.plot(self.halo.time / 1e6, self.halo.sfr_arr, label = 'halo')
        plt.title("SFR vs time")
        plt.xlabel("t [Myr]")
        plt.ylabel(r"SFR [M$_\odot$$yr^{-1}$]")
        plt.legend()
        dirname = os.getcwd().split('/')[-1]
        plt.savefig("/home/x-jstrack1/plots/halosfrcomp.png")
