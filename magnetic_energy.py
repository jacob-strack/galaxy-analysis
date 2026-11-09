from starter import * 
import matplotlib
import scipy.stats as sp

def magnetic_energy_time_profile(start_time, end_time):
    field_arr = [] 
    time_arr = [] 
    for i in range(start_time, end_time): 
        filename = "DD" + str(i).zfill(4) + "/DD" + str(i).zfill(4)
        ds = yt.load(filename) 
        ad = ds.all_data() 
        field_arr.append((ad["magnetic_field_strength"].in_units("gauss")**2 / 2).mean()) 
        time_arr.append(ds.current_time.in_units("Myr")) 
    plt.semilogy(time_arr, field_arr) 
    plt.xlabel("time [Myr]") 
    plt.ylabel("Mean Magnetic Energy [G^2]")
    plt.savefig("magnetic_energy_time_profile.png")

def magnetic_energy_r_profile(t, center = [0.5, 0.5, 0.5]):
    ds = yt.load("DD" + str(t).zfill(4) + "/DD" + str(t).zfill(4))
    ad = ds.all_data()
    magnetic_energy = ad["magnetic_field_strength"].in_units("gauss")**2 / 2 
    x = ad[("gas", "x")].in_units("kpc") - yt.YTQuantity(center[0] * ds.length_unit, "cm").in_units("kpc")
    y = ad[("gas","y")].in_units("kpc") - yt.YTQuantity(center[1] * ds.length_unit, "cm").in_units("kpc")
    z = ad[("gas","z")].in_units("kpc") - yt.YTQuantity(center[2] * ds.length_unit, "cm").in_units("kpc")
    rad = np.sqrt(x**2 + y**2 + z**2) 
    arr, bins = sp.stats.binned_statistic(rad, magnetic_energy, bins = 1000) 
    plt.semilogy(bins[:-1], arr) 
    plt.xlabel("radius [kpc]") 
    plt.ylabel("Mean Magnetic Energy [G^2]")
    plt.savefig("magnetic_energy_time_profile" + str(t) + ".png")

def magnetic_energy_r_profile_time_series(start_t, end_t, center = [0.5, 0.5, 0.5]): 
    jet = plt.get_cmap('jet')
    cNorm = matplotlib.colors.Normalize(vmin = 0, vmax = end_t-start_t) 
    scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=jet)
    for t in range(start_t, end_t): 
        ds = yt.load("DD" + str(t).zfill(4) + "/DD" + str(t).zfill(4))
        ad = ds.all_data()
        magnetic_energy = ad["magnetic_field_strength"].in_units("gauss")**2 / 2 
        x = ad[("gas", "x")].in_units("kpc") - yt.YTQuantity(center[0] * ds.length_unit, "cm").in_units("kpc")
        y = ad[("gas","y")].in_units("kpc") - yt.YTQuantity(center[1] * ds.length_unit, "cm").in_units("kpc")
        z = ad[("gas","z")].in_units("kpc") - yt.YTQuantity(center[2] * ds.length_unit, "cm").in_units("kpc")
        rad = np.sqrt(x**2 + y**2 + z**2) 
        arr, bins, vals = sp.binned_statistic(rad, magnetic_energy, bins = 1000, range = (0,50)) 
        plt.semilogy(bins[:-1], arr, color = scalarMap.to_rgba(t)) 
    plt.xlabel("radius [kpc]") 
    plt.ylabel("Mean Magnetic Energy [G^2]")
    plt.savefig("magnetic_energy_time_series.png")

def magnetic_energy_r_profile_comp(t, directory, center = [0.5, 0.5, 0.5]):
    ds = yt.load("DD" + str(t).zfill(4) + "/DD" + str(t).zfill(4))
    ad = ds.all_data()
    ds2 = yt.load(directory + "/DD" + str(t).zfill(4) + "/DD" + str(t).zfill(4))
    ad2 = ds2.all_data()
    magnetic_energy = ad["magnetic_field_strength"].in_units("gauss")**2 / 2 
    magnetic_energy2 = ad2["magnetic_field_strength"].in_units("gauss")**2 / 2
    x = ad[("gas", "x")].in_units("kpc") - yt.YTQuantity(center[0] * ds.length_unit, "cm").in_units("kpc")
    y = ad[("gas","y")].in_units("kpc") - yt.YTQuantity(center[1] * ds.length_unit, "cm").in_units("kpc")
    z = ad[("gas","z")].in_units("kpc") - yt.YTQuantity(center[2] * ds.length_unit, "cm").in_units("kpc")
    rad = np.sqrt(x**2 + y**2 + z**2) 
    arr, bins = sp.stats.binned_statistic(rad, magnetic_energy, bins = 1000) 
    plt.semilogy(bins[:-1], arr, label = 'no_halo') 
    x2 = ad2[("gas", "x")].in_units("kpc") - yt.YTQuantity(center[0] * ds.length_unit, "cm").in_units("kpc")
    y2 = ad2[("gas","y")].in_units("kpc") - yt.YTQuantity(center[1] * ds.length_unit, "cm").in_units("kpc")
    z2 = ad2[("gas","z")].in_units("kpc") - yt.YTQuantity(center[2] * ds.length_unit, "cm").in_units("kpc")
    rad2 = np.sqrt(x2**2 + y2**2 + z2**2) 
    arr2, bins2 = sp.stats.binned_statistic(rad2, magnetic_energy2, bins = 1000) 
    plt.semilogy(bins2[:-1], arr2, label = 'halo') 
    plt.xlabel("radius [kpc]") 
    plt.ylabel("Mean Magnetic Energy [G^2]")
    plt.legend()
    plt.savefig("magnetic_energy_time_profile" + str(t) + ".png")
