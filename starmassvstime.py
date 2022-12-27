from starter import * 
class starmassvstime: 
	def __init__(self, num):
		self.num = num
		self.t_arr = []
		self.sm_arr = []
		for i in range(13,num): 
			dsname = "DD" + str(i).zfill(4)
			filename = dsname + '/' + dsname
			ds = yt.load(filename)
			ad = ds.all_data()
			starmass = ad["particle_mass"].sum()
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

