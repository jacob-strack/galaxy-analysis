from starter import *
from scipy import stats
from matplotlib.colors import LogNorm
class cumulativemass:
	def __init__(self, num): 
		self.num = num #number of data dumps
		self.arr = []
		self.arr2 = [] 
		self.t_arr = [] 
		self.t2_arr = []
		self.edges = 0
		self.edges2 = 0
		self.m = 0
		self.m2 = 0

	def run(self):
		for i in range(self.num): 
			ds_name = 'DD' + str(i).zfill(4)
			filename = ds_name + '/' + ds_name
			ds = yt.load(filename)
			ad = ds.all_data()
			sphere = ds.sphere([0.5,0.5,0.5], (25, 'kpc'))
			z = sphere['z']
			midplane = sphere.center[2].in_units('cm') 
			z -= midplane
			z = np.abs(z)
			z = z.in_units('kpc')
			sort1 =  np.argsort(z)
			self.m = np.cumsum(sphere['cell_mass'][sort1])
			print(self.m.in_units('M_sun'))
			bins = stats.binned_statistic(np.sort(z), self.m.in_units('M_sun'), bins = 16)
			self.edges = bins[1]
			self.t_arr.append(ds.current_time)
			self.arr.append(bins[0])

			#second dataset with halo
			ds2_name = '../7_level_tweaked_halo/' + filename
			ds2 = yt.load(ds2_name)
			sph2 = ds2.sphere([0.5, 0.5, 0.5], (25, 'kpc'))
			z2 = sph2['z']
			midplane2 = sph2.center[2].in_units('cm')
			z2 -= midplane2
			z2 = np.abs(z2)
			z2 = z2.in_units('kpc')
			sorted2 = np.argsort(z2)
			self.m2 = np.cumsum(sph2['cell_mass'][sorted2])
			bins2 = stats.binned_statistic(np.sort(z2), self.m2.in_units('M_sun'), bins = 16)
			self.edges2 = bins2[1]
			self.t2_arr.append(ds2.current_time)
			self.arr2.append(bins2[0])
		for i in range(0,len(self.edges)): 
			self.edges[i] = np.format_float_scientific(float(self.edges[i]),3)
		for i in range(0,len(self.edges2)): 
			self.edges2[i] = np.format_float_scientific(float(self.edges2[i]),3)
		self.t_arr = np.asarray(self.t_arr)
		self.t2_arr = np.asarray(self.t2_arr)
		self.arr = np.asarray(self.arr)
		self.arr2 = np.asarray(self.arr2)
		

	def plot(self, vmin, vmax):
		fig, ax = plt.subplots(2,1)
		fig.set_size_inches(12, 10)
		a = ax[0].imshow(self.arr.T, origin = 'lower', cmap = 'jet', norm=LogNorm(vmin = vmin, vmax = vmax))
		plt.colorbar(a, ax = ax[0])
		ax[0].set_xlabel('t [Myr]')
		ax[0].set_ylabel('z [kpc]')
		ax[0].title.set_text('<m>(z,t)')
		ax[0].set_xticks(np.arange(0, len(self.t_arr), 5),np.round(self.t_arr[0::5],2))
		ax[0].set_yticks(np.arange(0, len(self.edges), 5), self.edges[0::5])
		b = ax[1].imshow(self.arr2.T, origin = 'lower', cmap = 'jet',norm=LogNorm(vmin = vmin, vmax = vmax))
		plt.colorbar(b, ax = ax[1])
		ax[1].set_xlabel('t [Myr]')
		ax[1].set_ylabel('z [kpc]')
		ax[1].title.set_text('<m>(z,t) Halo')
		ax[1].set_xticks(np.arange(0, len(self.t2_arr), 5),np.round(self.t2_arr[0::5],2))
		ax[1].set_yticks(np.arange(0, len(self.edges2), 5), self.edges2[0::5])
		fig.savefig("frames/masscomp.png")	

	def massprofile_z(self, timestep): #the idea is right but i haven't tested this yet
		y_1 = self.arr[timestep, :]
		x_1 = np.linspace(0, self.edges[-1], len(y_1))
		y_2 = self.arr2[timestep, :]
		plt.semilogy(x_1,y_1)
		x_2 = np.linspace(0, self.edges2[-1], (len(y_2)))
		plt.semilogy(x_2,y_2)
		plt.xlabel('z [kpc]')
		plt.ylabel('Mass [M_Sun]')
		plt.savefig("frames/zmassprofile.png")
	


if "cmplot" not in dir(): 
	cmplot = cumulativemass(3)
	cmplot.run()

