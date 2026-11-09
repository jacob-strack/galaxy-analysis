from starter import *
import scipy.stats as sp
import matplotlib
import gc

def make_rotational_velocity(filename):
    ds = yt.load(filename)
    ds.add_field(name=("gas", "rotational_velocity"), function = rotational_velocity, sampling_type="local", units = "km/s")
    return ds

def rotational_velocity(field, data):
    #v_tan = sqrt(v_phi^2 + v_theta^2)
    #hard coded bs for now for center of galaxy, I could pass in a parameter
    #or read DiskPosition from output parameter file, but i'd need to retool where DiskPosition is defined
    x = data[("gas", "x")] - yt.YTQuantity(0.5 * data.ds.length_unit, "cm")
    y = data[("gas", "y")] - yt.YTQuantity(0.5 * data.ds.length_unit, "cm")
    z = data[("gas", "z")] - yt.YTQuantity(0.53 * data.ds.length_unit, "cm")
    print(np.shape(x))  
    theta = np.arctan(np.sqrt(x**2 + y**2) / z)
    phi = np.arctan(y/x)
    return np.sqrt((-data["velocity_x"]*np.sin(phi) + data["velocity_y"]*np.cos(phi))**2 + (data["velocity_x"]*np.cos(theta)*np.cos(phi) + data["velocity_y"]*np.cos(theta)*np.sin(phi) - data["velocity_z"]*np.sin(theta))**2).in_units("km/s")

def plot_rotation_curve(filename, center, color_val):
    ds = make_rotational_velocity(filename)
    ad = ds.all_data()
    r = np.sqrt((ad[("gas","x")] - yt.YTQuantity(center[0]*ds.length_unit, "cm")).in_units("kpc")**2 + (ad[("gas","y")] - yt.YTQuantity(center[1]*ds.length_unit, "cm")).in_units("kpc")**2) # + (ad[("gas","z")]  - yt.YTQuantity(center[2]*ds.length_unit, "cm")).in_units("kpc")**2)
    
    #only take points below abs(z) criterion (that's how i find the disk, could also do a density criterion)
    #bin by radius 
    hist = sp.binned_statistic(r[np.where(np.abs((ad[("gas", "z")] - yt.YTQuantity(center[2]*ds.length_unit,"cm")).in_units("kpc")) < 1)[0]], ad["rotational_velocity"][np.where(np.abs((ad[("gas","z")] - yt.YTQuantity(center[2]*ds.length_unit, "cm")).in_units("kpc")) < 1)[0]], bins = 100, range = (0,15))
    plt.plot(hist[1][:-1], hist[0], color=plt.cm.jet(color_val))
    plt.xlabel("r [kpc]")
    plt.ylabel("v [km/s]")
    plt.savefig("frames/rotation_curve.png")

def my_radius(field, data): 
    return np.sqrt((data[("gas","x")] - yt.YTQuantity(0.5*data.ds.length_unit, "cm")).in_units("kpc")**2 + (data[("gas","y")] - yt.YTQuantity(0.25*data.ds.length_unit, "cm")).in_units("kpc")**2 + (data[("gas","z")]  - yt.YTQuantity(0.5*data.ds.length_unit, "cm")).in_units("kpc")**2)

def make_my_radius(filename): 
    ds = yt.load(filename)
    ds.add_field(name=("gas","my_radius"), function = my_radius, sampling_type="local",units = "kpc")
    return ds

def turbulent_kinetic_energy(field, data): 
    rot_curve = data.get_field_parameter("rotation_curve", np.zeros(100)) 
    rot_curve_bins = data.get_field_parameter("rotation_curve_bins",np.zeros(100))
    x = data[("gas", "x")].in_units("cm") - yt.YTQuantity(0.5 * data.ds.length_unit, "cm")
    y = data[("gas", "y")].in_units("cm") - yt.YTQuantity(0.5 * data.ds.length_unit, "cm")
    z = data[("gas", "z")].in_units("cm") - yt.YTQuantity(0.53 * data.ds.length_unit, "cm")
    rot_vel = yt.ytarray(np.full(data["velocity_x"].shape,np.nan), data["velocity_x"].units)
    a = np.sqrt(data["x-velocity"]*data["x-velocity"] + data["y-velocity"]*data["y-velocity"] + data["z-velocity"]*data["z-velocity"])
    r_bins = np.digitize(r.in_units("kpc"), rot_curve_bins[:-1], right=False)
    r_bins[np.where(r_bins == 100)[0]] -= 1
    print(np.shape(r_bins))
    rot_vel = a - data.ds.arr(rot_curve[r_bins], "code_velocity") 
    rot_K = 0.5 * data["density"].in_units("g/cm**3") * (rot_vel.in_units("km/s"))* rot_vel.in_units("km/s")
    return rot_K
    

def make_turbulent_kinetic_energy(filename, axis, center = [0.5, 0.5, 0.5]):
    ds = yt.load(filename)
    ad = ds.all_data() 
    x = ad[("gas", "x")].in_units("cm") - yt.YTQuantity(0.5 * ds.length_unit, "cm")
    y = ad[("gas", "y")].in_units("cm") - yt.YTQuantity(0.5 * ds.length_unit, "cm")
    z = ad[("gas", "z")].in_units("cm") - yt.YTQuantity(0.53 * ds.length_unit, "cm")
    theta = np.arctan(np.sqrt(x**2 + y**2) / z)
    phi = np.arctan(y/x)
    r = np.sqrt(x**2 + y**2 + z**2)
    rot = np.sqrt((-1.0 * ad["velocity_x"]*np.sin(phi) + ad["velocity_y"]*np.cos(phi))**2 + (ad["velocity_x"]*np.cos(theta)*np.cos(phi) + ad["velocity_y"]*np.cos(theta)*np.sin(phi) - ad["velocity_z"]*np.sin(theta))**2)
    pts_mask = np.where(np.abs(z.in_units("pc")) < 20)
    hist = sp.binned_statistic(r[pts_mask[0]].in_units("kpc"), rot[pts_mask[0]], bins = 100, range = (0,50))
    ad.set_field_parameter("rotation_curve", hist[0])
    ad.set_field_parameter("rotation_curve_bins", hist[1])
    ds.add_field(name=("gas", "turbulent_kinetic_energy"), function = turbulent_kinetic_energy, sampling_type="cell", units = "g/cm**3 * km**2/s**2")
    ad["turbulent_kinetic_energy"]
    #prj = yt.ProjectionPlot(ds, axis, "turbulent_kinetic_energy",  center = center, data_source=ad)
    #prj.save()
    return ad["turbulent_kinetic_energy"]

def turbulent_kinetic_energy_r_profile(filename, center = [0.5, 0.5, 0.5]): 
    ds = yt.load(filename)
    ad = ds.all_data() 
    #compute rotation curve
    x = ad[("gas", "x")].in_units("cm") - yt.YTQuantity(0.5 * ds.length_unit, "cm")
    y = ad[("gas", "y")].in_units("cm") - yt.YTQuantity(0.5 * ds.length_unit, "cm")
    z = ad[("gas", "z")].in_units("cm") - yt.YTQuantity(0.53 * ds.length_unit, "cm")
    theta = np.arctan(np.sqrt(x**2 + y**2) / z)
    #rot_vel = yt.YTArray(np.full(ad["velocity_x"].shape,np.nan), ad["velocity_x"].units)
    phi = np.arctan(y/x)
    #was using just z, but that uses the center of cells, so big enough cells would fall outside the range 
    #using the edges might be a better thing to do 
    r = np.sqrt(x**2 + y**2 + z**2)
    a = np.sqrt(ad["x-velocity"]*ad["x-velocity"] + ad["y-velocity"]*ad["y-velocity"] + ad["z-velocity"]*ad["z-velocity"])
    rot = np.sqrt((-1.0 * ad["velocity_x"]*np.sin(phi) + ad["velocity_y"]*np.cos(phi))**2 + (ad["velocity_x"]*np.cos(theta)*np.cos(phi) + ad["velocity_y"]*np.cos(theta)*np.sin(phi) - ad["velocity_z"]*np.sin(theta))**2)
    lower_z = (z - ad["dz"]/2).in_units("kpc")
    upper_z = (z + ad["dz"]/2).in_units("kpc")
    pts_mask = np.where(np.logical_or(np.abs(lower_z) < 20, np.abs(upper_z) < 20))
    rot_curve, rot_curve_bins, vals = sp.binned_statistic(r[pts_mask[0]].in_units("kpc"), rot[pts_mask[0]], bins = 100, range = (0,50))
    r_bins = np.digitize(r.in_units("kpc"), rot_curve_bins[:-1], right=False)
    r_bins[np.where(r_bins == 100)[0]] -= 1
    print(len(np.where(np.isnan(rot_curve[r_bins]) == False)[0]))
    rot_vel = a - ds.arr(rot_curve[r_bins], "code_velocity") 
    rot_K = 0.5 * ad["density"].in_units("g/cm**3") * (rot_vel.in_units("km/s"))* rot_vel.in_units("km/s")
    print("rotK ",rot_K)
    #bin and plot values 
    arr, plt_bins, plt_vals = sp.binned_statistic(r.in_units("kpc"), rot_K, statistic=np.nanmean, bins=100, range=(0,50))
    return arr,plt_bins
    
def make_rot_K_plot(filename, center = [0.5, 0.5, 0.5]): 
    arr, plt_bins = turbulent_kinetic_energy_r_profile(filename, center)
    plt.semilogy(plt_bins[:-1], arr) 
    plt.xlabel("radius [kpc]") 
    plt.ylabel("Turbulent Kinetic Energy") 
    plt.savefig("KEvsr_profile.png")

def make_rot_K_plot_comp(filename, directory, center = [0.5, 0.5, 0.5]): 
    arr, plt_bins = turbulent_kinetic_energy_r_profile(filename, center)
    arr2, plt_bins2 = turbulent_kinetic_energy_r_profile(directory + "/" + filename, center) 
    plt.semilogy(plt_bins[:-1], arr, label = 'no halo') 
    plt.semilogy(plt_bins2[:-1], arr2, label = 'halo')
    plt.legend()
    plt.xlabel("radius [kpc]") 
    plt.ylabel("Turbulent Kinetic Energy") 
    plt.savefig("KEvsr_profilecomp.png")

def rotation_curve(filename, center = [0.5, 0.5,0.5]): 
    ds = yt.load(filename)
    ad = ds.all_data()
    x = ad[("gas", "x")].in_units("cm") - yt.YTQuantity(center[0] * ds.length_unit, "cm")
    y = ad[("gas", "y")].in_units("cm") - yt.YTQuantity(center[1] * ds.length_unit, "cm")
    z = ad[("gas", "z")].in_units("cm") - yt.YTQuantity(center[2] * ds.length_unit, "cm")
    theta = np.arctan(np.sqrt(x**2 + y**2) / z)
    phi = np.arctan(y/x)
    #was using just z, but that uses the center of cells, so big enough cells would fall outside the range 
    #using the edges might be a better thing to do 
    r = np.sqrt(x**2 + y**2 + z**2)
    rot = np.sqrt((-1.0 * ad["velocity_x"]*np.sin(phi) + ad["velocity_y"]*np.cos(phi))**2 + (ad["velocity_x"]*np.cos(theta)*np.cos(phi) + ad["velocity_y"]*np.cos(theta)*np.sin(phi) - ad["velocity_z"]*np.sin(theta))**2).in_units("km/s")
    lower_z = (z - ad["dz"]/2).in_units("kpc")
    upper_z = (z + ad["dz"]/2).in_units("kpc")
    pts_mask = np.where(np.logical_or(np.abs(lower_z) < 20, np.abs(upper_z) < 20))
    rot_curve, rot_curve_bins, vals = sp.binned_statistic(r[pts_mask[0]].in_units("kpc"), rot[pts_mask[0]], bins = 100, range = (0,25))
    del ds
    gc.collect()
    return rot_curve, rot_curve_bins 

def rotation_curve_comp(filename, directory, center = [0.5, 0.5, 0.5]): 
    arr, plt_bins = rotation_curve(filename, center)
    arr2, plt_bins2 = rotation_curve(directory + "/" + filename, center) 
    plt.semilogy(plt_bins[:-1], arr, label = 'no halo') 
    plt.semilogy(plt_bins2[:-1], arr2, label = 'halo')
    plt.legend()
    plt.xlabel("radius [kpc]") 
    plt.ylabel("Rotational Velocity [km/s]") 
    plt.savefig("rotation_curve_comp.png")

def rotation_curve_time_series(start_t, end_t, center = [0.5, 0.5, 0.5]): 
    jet = plt.get_cmap('jet')
    cNorm = matplotlib.colors.Normalize(vmin = 0, vmax = end_t-start_t) 
    scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=jet)
    for t in range(start_t, end_t): 
        filename = "DD" + str(t).zfill(4) + "/DD" + str(t).zfill(4)
        arr, bins = rotation_curve(filename, center)
        plt.semilogy(bins[:-1], arr, color = scalarMap.to_rgba(t)) 
    plt.xlabel("radius [kpc]") 
    plt.ylabel("rotation velocity")
    plt.savefig("rotation_curve_time_series.png")

def rotation_curve_comp_time(start_t, end_t, directory, center = [0.5, 0.5, 0.5]): 
    for t in range(start_t, end_t): 
        plt.close('all')
        dsname = "DD" + str(t).zfill(4)
        filename = dsname + "/" + dsname
        arr, plt_bins = rotation_curve(filename, center)
        arr2, plt_bins2 = rotation_curve(directory + "/" + filename, center) 
        plt.semilogy(plt_bins[:-1], arr, label = 'no halo') 
        plt.semilogy(plt_bins2[:-1], arr2, label = 'halo')
        del arr
        del arr2 
        del plt_bins
        del plt_bins2
        gc.collect()
        plt.legend()
        plt.xlabel("radius [kpc]") 
        plt.ylabel("Rotational Velocity [km/s]") 
        plt.savefig(dsname + "_rotation_curve_comp.png")
        plt.close('all')
