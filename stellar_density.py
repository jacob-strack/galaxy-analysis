from starter import * 
from yt.data_objects.particle_filters import add_particle_filter
import os 

def stars(pfilter, data): 
    filter = data[pfilter.filtered_type, "particle_type"] == 2
    return filter

def make_stellar_density(filename): 
    yt.add_particle_filter("stars", function = stars, filtered_type = "all", requires = ["particle_type"]) 
    ds = yt.load(filename)
    ds.add_particle_filter("stars")
    deposited = ds.add_deposited_particle_field(("stars", "particle_mass"), method = "sum")
    ds.add_field(name=("gas", "stellar_density"), function = stellar_density, sampling_type="local", units = "Msun/pc**3")
    return ds 

def stellar_density(field,data): 
    return data["deposit", "stars_sum_mass"].in_units("Msun") / data["cell_volume"].in_units("pc**3")

def stellar_density_proj(name, axis = "z"):
    filename = name + "/" + name
    ds = make_stellar_density(filename)
    yt.ProjectionPlot(ds, axis, ("gas", "stellar_density"), width = (200.0,"kpc")).save("frames/" + name + "_stellar_density_" + axis + ".png") 

