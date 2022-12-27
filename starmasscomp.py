from starter import * 
import os
from starmassvstime import starmassvstime
class starmasscomp: 
    def __init__(self,num):
        self.num = num
        self.no_halo = starmassvstime(self.num)
        os.chdir("../7_level_tweaked_halo")
        self.halo = starmassvstime(self.num)
        os.chdir("../7_level_tweaked")
    
    def plot(self): 
        plt.plot(self.halo.t_arr, self.halo.sm_arr, label = "halo")
        plt.plot(self.no_halo.t_arr, self.no_halo.sm_arr, label = "no halo")
        plt.xlabel('t [Myr]')
        plt.ylabel(r'M$_[\odot]$')
        plt.title("Total Star Mass vs Time")
        plt.semilogy()
        plt.legend()
        plt.savefig('frames/starmasscomp.png')

