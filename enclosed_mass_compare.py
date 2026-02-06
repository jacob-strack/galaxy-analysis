import yt 
import matplotlib.pyplot as plt 

def mass_enclosed_compare(ds, name, radius, center): 
    plt.close('all') 
    fig, ax = plt.subplots(figsize=(8,6))
    sph = ds.sphere(center, (radius, "code_length"))
    profiles = yt.create_profile(sph, ("index", "radius"), "MassEnclosed", weight_field="ones", extrema={("index", "radius"):(1e22,1e24)}, logs={"radius":True}) 
    profile_mass = yt.create_profile(sph, ("index", "radius"), "cell_mass", weight_field="ones", extrema={("index", "radius"):(1e22, 1e24)}, logs={"radius":True}, accululation = True)
    ax.loglog(profiles.x.in_units("kpc").d, profiles["MassEnclosed"].in_units("Msun").d)
    ax.loglog(profile_mass.x.in_units("kpc").d, profile_mass["total_mass"].d)
    ax.set_xlabel("r [kpc]")
    ax.set_ylabel("Mass Enclosed [" + str(sph["MassEnclosed"].units) + "]")
    plt.axvline(x = 206, color = "r", linestyle="--")
    plt.savefig("frames/" + name)

def make_total_mass(filename): 
    ds = yt.load(filename) 
    ds.add_field(name = ("gas","total_mass"), function = total_mass, sampling_type="local", units = "Msun")
    return ds 

def total_mass(field,data): 
    return (data["density"] * data["cell_volume"]).in_units("Msun") 

