from starter import * 
import scipy.stats as sp 

def z_vel_vs_radius(filename, center = [0.5, 0.5, 0.5]): 
    ds = yt.load(filename) 
    ad = ds.all_data() 
    x = ad[("gas", "x")].in_units("kpc") - yt.YTQuantity(0.5 * ds.length_unit, "cm").in_units("kpc")    
    y = ad[("gas", "y")].in_units("kpc") - yt.YTQuantity(0.5 * ds.length_unit, "cm").in_units("kpc")
    radius = np.sqrt(x**2 + y**2) 
    z_vel, rad_bins, vals = sp.binned_statistic(radius, np.abs(ad["velocity_z"]), bins = 1000, range = (0, 50))
    return z_vel, rad_bins[:-1] 

def plot_z_vel_vs_radius(t, center = [0.5, 0.5, 0.5]): 
    filename = "DD" + str(t).zfill(4) + "/DD" + str(t).zfill(4)
    arr, bins = z_vel_vs_radius(filename, center) 
    plt.plot(bins, arr) 
    plt.xlabel("r_cyl [kpc]") 
    plt.ylabel("v_z [km/s]") 
    plt.savefig("DD" + str(t).zfill(4) + "_vz_vs_r.png")

def z_mass_flux_vs_radius(filename, center = [0.5, 0.5, 0.5]): 
    ds = yt.load(filename) 
    ad = ds.all_data() 
    x = ad[("gas", "x")].in_units("kpc") - yt.YTQuantity(0.5 * ds.length_unit, "cm").in_units("kpc")    
    y = ad[("gas", "y")].in_units("kpc") - yt.YTQuantity(0.5 * ds.length_unit, "cm").in_units("kpc")
    radius = np.sqrt(x**2 + y**2) 
    z, rad_bins, vals = sp.binned_statistic(radius, np.abs(ad["density"].in_units("Msun/cm**3") * ad["cell_volume"].in_units("cm**3") * ad["velocity_z"]), bins = 1000, range = (0, 50))
    return z, rad_bins[:-1] 

def plot_z_mass_flux_vs_radius(t, center = [0.5, 0.5, 0.5]): 
    filename = "DD" + str(t).zfill(4) + "/DD" + str(t).zfill(4)
    arr, bins = z_mass_flux_vs_radius(filename, center) 
    plt.plot(bins, arr) 
    plt.xlabel("r_cyl [kpc]") 
    plt.ylabel("mass flux z") 
    plt.savefig("DD" + str(t).zfill(4) + "_mass_flux_vs_r.png")
