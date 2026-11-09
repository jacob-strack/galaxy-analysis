import yt 
import matplotlib.pyplot as plt 

def make_profile(name, field, radius, center): 
    ds = yt.load(name + '/' + name) 
    sph = ds.sphere(center, (radius, "code_length"))
    profiles = yt.create_profile(sph, ("index", "radius"), field, weight_field="ones", extrema={("index", "radius"):(1e22, 1e24)}, logs={"radius":True})
    plot = yt.ProfilePlot.from_profiles(profiles)
    plot.set_unit(("index","radius"), "kpc")
    plot.save("frames/" + name) 

def make_profile_mp(ds, name, field, radius, center, power_law = False, power_law_index = 1): 
    plt.close('all')
    fig, ax = plt.subplots(figsize=(8,6))
    sph = ds.sphere(center, (radius, "code_length")) 
    profiles = yt.create_profile(sph, ("index", "radius"), field, weight_field="ones", extrema={("index", "radius"):(1e22,1e24)}, logs={"radius":True})
    ax.loglog(profiles.x.in_units("kpc").d, profiles[field].d, label=field)
    if(power_law):  
        ax.loglog(profiles.x.in_units("kpc").d,profiles[field].d[0] * (1 + (profiles.x.in_units("kpc").d)**(power_law_index)), linestyle='--', label = "$r^{2/3}$")
    ax.legend()
    ax.set_xlabel("r [kpc]")
    ax.set_ylabel(field + " [" + str(sph[field].units) + "]")
    plt.axvline(x = 206, color = "r", linestyle="--")
    plt.savefig("frames/" + name)
    

