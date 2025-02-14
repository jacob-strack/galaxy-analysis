from starter import * 
from yt.data_objects.particle_filters import add_particle_filter 
import os 

def formed_stars(pfilter, data): 
    age = data.ds.current_time - data["all", "creation_time"]
    filter = age.in_units("Myr") <= 11
    return filter 

def make_sfr(filename): 
    yt.add_particle_filter("formed_stars", function = formed_stars, filtered_type = "all", requires=["creation_time"])
    ds = yt.load(filename) 
    ds.add_particle_filter("formed_stars") 
    deposited = ds.add_deposited_particle_field(("formed_stars", "particle_mass"), method = "sum")
    ds.add_field(name=("gas", "sfr_density"), function=sfr_density, sampling_type="local", units="Msun/cm**3/yr")
    print(deposited) 
    print(ds.r[deposited].in_units("Msun"))
    print(ds.r[deposited].max())
    return ds #probably wrong thing to do

def sfr_density(field, data): 
    return data["deposit", "formed_stars_sum_mass"].in_units("Msun") / (data["cell_volume"].in_units("cm**3") * (data.ds.quan(11.0, "Myr").in_units("yr"))) 

def ks_plot(filename, axis): 
    ds = make_sfr(filename)
    sfr_column = ds.proj("sfr_density", 0)
    column = ds.proj("density", 0)
    #flatten these projections and plot
    plt.scatter(column["density"], sfr_column["sfr_density"])
    plt.loglog()
    plt.savefig("test.png")

def ks_plot_all(axis, lower_t, upper_t): 
    for i in range(lower_t, upper_t): 
        filename = "DD" + str(i).zfill(4) + "/DD" + str(i).zfill(4)
        ds = make_sfr(filename)
        sfr_column = ds.proj("sfr_density", 0)
        column = ds.proj("density", 0)
        plt.scatter(column["density"], sfr_column["sfr_density"])
    plt.loglog()
    plt.savefig("test.png")
