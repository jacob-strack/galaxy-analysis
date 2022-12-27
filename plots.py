

def meanvarvstime(fieldname, dumpnum): 
	t_arr = [] 
	var_arr = [] 
	for i in range(0,num): 
	    filename ='DD' + str(i).zfill(4)
	    ds = yt.load(filename + '/' + filename)
	    ad = ds.all_data()
	    temp_var_arr = ad[var]
	    mean_var = np.mean(temp_var_arr)
	    t = ds.current_time
	    var_arr.append(mean_var)
	    t_arr.append(t)
	t_arr = np.array(t_arr)
	var_arr = np.array(var_arr)
	return plt.plot(t_arr, var_arr)

def starmassvstime(fieldname, dumpnum): 
	ts = yt.load("DD001?/DD????")
	storage = {}
	for store, ds in ts.piter(storage=storage):
	    ad = ds.all_data()
	    starmass = ad['particle_mass'].sum()
	    store.result = (ds.current_time.in_units("Myr"), starmass.in_units("Msun"))
	arr = np.array(list(storage.values()))
	return plt.semilogy(arr[:,0], arr[:,1])
	 
