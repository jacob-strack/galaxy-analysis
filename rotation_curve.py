from starter import *
import scipy.stats as sp

def make_rotational_velocity(filename):
    ds = yt.load(filename)
    ds.add_field(name=("gas", "rotational_velocity"), function = rotational_velocity, sampling_type="local", units = "km/s")
    return ds

def rotational_velocity(field, data):
    #v_tan = sqrt(v_phi^2 + v_theta^2)
    #hard coded bs for now for center of galaxy, I could pass in a parameter
    #or read DiskPosition from output parameter file, but i'd need to retool where DiskPosition is defined
    x = data[("gas", "x")] - yt.YTQuantity(0.25 * data.ds.length_unit, "cm")
    y = data[("gas", "y")] - yt.YTQuantity(0.25 * data.ds.length_unit, "cm")
    z = data[("gas", "z")] - yt.YTQuantity(0.25 * data.ds.length_unit, "cm")
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

