from starter import * 
import scipy.stats as sp 

def enclosed_mass_profile(filename, center = [0.5, 0.5, 0.5]): 
    ds = yt.load(filename) 
    ad = ds.all_data() 
    x = ad[("gas", "x")].in_units("kpc") - yt.YTQuantity(0.5 * ds.length_unit, "cm").in_units("kpc")    
    y = ad[("gas", "y")].in_units("kpc") - yt.YTQuantity(0.5 * ds.length_unit, "cm").in_units("kpc")
    radius = np.sqrt(x**2 + y**2) 
    m_enc, rad_bins, vals = sp.binned_statistic(radius, ad["density"]*ad["cell_volume"], bins = 100, statistic = "sum", range = (0,50)) 
    cumulative_mass = 0.0 
    for i in range(len(m_enc)): 
        cumulative_mass += m_enc[i] 
        m_enc[i] = cumulative_mass 
    return rad_bins[:-1], m_enc 

def plot_enclosed_mass_profile(t, center = [0.5, 0.5, 0.5]): 
    filename = "DD" + str(t).zfill(4) + "/DD" + str(t).zfill(4)
    bins, arr = enclosed_mass_profile(filename, center) 
    plt.plot(bins, arr) 
    plt.xlabel("r_cyl [kpc]") 
    plt.ylabel("Mass Enclosed") 
    plt.savefig("DD" + str(t).zfill(4) + "_enclosed_mass_profile.png") 
    
def mass_profile(filename, center = [0.5, 0.5, 0.5]): 
    ds = yt.load(filename) 
    ad = ds.all_data() 
    x = ad[("gas", "x")].in_units("kpc") - yt.YTQuantity(0.5 * ds.length_unit, "cm").in_units("kpc")    
    y = ad[("gas", "y")].in_units("kpc") - yt.YTQuantity(0.5 * ds.length_unit, "cm").in_units("kpc")
    radius = np.sqrt(x**2 + y**2) 
    m_enc, rad_bins, vals = sp.binned_statistic(radius, ad["density"]*ad["cell_volume"], bins = 100, statistic = "sum", range = (0,50))
    return rad_bins[:-1], m_enc 

def plot_mass_profile(t, center = [0.5, 0.5, 0.5]): 
    filename = "DD" + str(t).zfill(4) + "/DD" + str(t).zfill(4)
    bins, arr = mass_profile(filename, center) 
    plt.plot(bins, arr) 
    plt.xlabel("r_cyl [kpc]") 
    plt.ylabel("Mass") 
    plt.savefig("DD" + str(t).zfill(4) + "_mass_profile.png") 
