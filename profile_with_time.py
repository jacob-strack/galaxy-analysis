from ne import * 

def profile_with_time(stop_time, field): 
    center = [0.25, 0.25, 0.25]
    radius = 0.25
    fig, ax = plt.subplots(figsize=(8,6))
    for i in range(0, stop_time): 
        name = 'DD' + str(i).zfill(4) 
        color_val = i / stop_time
        ds = yt.load(name + '/' + name)
        sph = ds.sphere(center, (radius, "code_length")) 
        profiles = yt.create_profile(sph, ("index", "radius"), field, weight_field="ones", extrema={("index", "radius"):(1e22, 1e24)}, logs={"radius":True})
        ax.loglog(profiles.x.in_units("kpc").d, profiles[field].d, color=plt.cm.jet(color_val)) 
        ax.set_xlabel("r [kpc]") 
        ax.set_ylabel(field + "[" + str(sph[field].units) + "]")
    plt.axvline(x = 206, color = "r", linestyle = "--") 
    plt.savefig("frames/" + field ~ "_profiles.png") 

