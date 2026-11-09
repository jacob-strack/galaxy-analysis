from starter import * 
import scipy.stats as sp 

def mean_magnetic_field_profile(filename, component, center = [0.5, 0.5, 0.5]): 
    ds = yt.load(filename) 
    ad = ds.all_data() 
    x = ad[("gas", "x")].in_units("cm") - yt.YTQuantity(center[0] * ds.length_unit, "cm")
    y = ad[("gas", "y")].in_units("cm") - yt.YTQuantity(center[1] * ds.length_unit, "cm")
    z = ad[("gas", "z")].in_units("cm") - yt.YTQuantity(center[2] * ds.length_unit, "cm")
    radius = np.sqrt(x**2 + y**2 + z**2).in_units("kpc")
    print("radius", radius)
    field_name = "B" + component 
    field = ad[field_name].in_units("gauss")
    mean_field, mean_field_bins, vals = sp.binned_statistic(radius, field,bins = 100, range=(0, 50))
    print("mean field bins", mean_field_bins)
    return mean_field, mean_field_bins[:-1] 

def turbulent_magnetic_field_profile(filename, component, center = [0.5, 0.5, 0.5]): 
    mean_field, mean_field_bins = mean_magnetic_field_profile(filename, component, center) 
    ds = yt.load(filename) 
    ad = ds.all_data()
    x = ad[("gas", "x")].in_units("cm") - yt.YTQuantity(center[0] * ds.length_unit, "cm")
    y = ad[("gas", "y")].in_units("cm") - yt.YTQuantity(center[1] * ds.length_unit, "cm")
    z = ad[("gas", "z")].in_units("cm") - yt.YTQuantity(center[2] * ds.length_unit, "cm")
    radius = np.sqrt(x*x + y*y + z*z).in_units("kpc")
    field_name = "B" + component 
    field = ad[field_name].in_units("gauss")
    r_bins = np.digitize(radius, mean_field_bins)
    r_bins[np.where(r_bins == 100)[0]] -= 1 
    turb_field = field - yt.YTArray(mean_field[r_bins], "gauss") 
    mean_turb_field, mean_turb_field_bins, vals = sp.binned_statistic(radius, turb_field, bins = 100, range=(0,50))
    return mean_turb_field, mean_turb_field_bins[:-1]

def plot_turb_magnetic_field_profile(filename, component, center = [0.5, 0.5, 0.5]): 
    vals, bins = turbulent_magnetic_field_profile(filename, component, center)
    print(bins) 
    print(vals)
    plt.semilogy(bins, vals)
    plt.xlabel("r [kpc]") 
    plt.ylabel("b [G]") 
    plt.savefig("turb_magnetic_profile.png")

def plot_turb_magnetic_energy_profile(filename, center = [0.5, 0.5, 0.5]): 
    b_x, b_x_bins = turbulent_magnetic_field_profile(filename, "x", center)
    b_y, b_y_bins = turbulent_magnetic_field_profile(filename, "y", center)
    b_z, b_z_bins = turbulent_magnetic_field_profile(filename, "z", center)
    total_magnetic_energy = (b_x**2 + b_y**2 + b_z**2) / 2
    plt.semilogy(b_x_bins, total_magnetic_energy) 
    plt.xlabel("r [kpc]") 
    plt.ylabel("b^2 / 2")
    plt.savefig("turb_magnetic_energy_profile.png") 

def plot_turb_magnetic_energy_profile_dumber(timestep, center=[0.5,0.5,0.5]): 
    filename = "DD" + str(timestep).zfill(4) + "/DD" + str(timestep).zfill(4)
    ds = yt.load(filename) 
    ad = ds.all_data() 
    mean_Bx = ad["Bx"].mean().in_units("gauss")
    mean_By = ad["By"].mean().in_units("gauss")
    mean_Bz = ad["Bz"].mean().in_units("gauss")
    bx = ad["Bx"].in_units("gauss") - mean_Bx
    by = ad["By"].in_units("gauss") - mean_By
    bz = ad["Bz"].in_units("gauss") - mean_Bz
    turb_magnetic_energy = (bx**2 + by**2 + bz**2) / 2
    x = ad[("gas", "x")].in_units("cm") - yt.YTQuantity(center[0] * ds.length_unit, "cm")
    y = ad[("gas", "y")].in_units("cm") - yt.YTQuantity(center[1] * ds.length_unit, "cm")
    z = ad[("gas", "z")].in_units("cm") - yt.YTQuantity(center[2] * ds.length_unit, "cm")
    radius = np.sqrt(x*x + y*y + z*z).in_units("kpc")
    arr, bins, vals = sp.binned_statistic(radius, turb_magnetic_energy, bins = 100, range=(0,50))
    plt.semilogy(bins[:-1], arr) 
    plt.xlabel("r [kpc]") 
    plt.ylabel("b^2 / 2") 
    plt.savefig("DD" + str(timestep).zfill(4) + "_turbulent_magnetic_energy_r_profile.png")

    
