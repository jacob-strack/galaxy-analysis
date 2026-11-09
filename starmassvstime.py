from starter import * 
class starmassvstime: 
    def __init__(self, num, start_time = 0):
        self.num = num
        self.t_arr = []
        self.sm_arr = []
        yt.add_particle_filter("stars", function=stars, filtered_type='all', requires=["particle_type", "creation_time"])
        for i in range(start_time,num): 
            dsname = "DD" + str(i).zfill(4)
            filename = dsname + '/' + dsname
            ds = yt.load(filename)
            ds.add_particle_filter("stars")
            ad = ds.all_data()
            starmass = ad[("stars","particle_mass")].sum()
            time = ds.current_time
            self.sm_arr.append(starmass)
            self.t_arr.append(time)
        self.t_arr = np.array(self.t_arr)
        self.sm_arr = np.array(self.sm_arr)


    def plot(self): 
        plt.plot(self.t_arr, self.sm_arr)
        plt.xlabel('t [Myr]')
        plt.ylabel(r'M$_[\odot]$')
        plt.title("Total Star Mass vs Time")
        plt.savefig("frames/starmass.png")

def stars(pfilter, data, include_init = True):
    pt = data["all", "particle_type"] == 2
    ct = data["all", "creation_time"] > 0.0
    if include_init == False:
        fil = np.logical_and(pt, ct)a
    else: 
        fil = pt
    return fil


