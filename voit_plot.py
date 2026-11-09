from ne import * 

center = [0.25, 0.25, 0.25]
radius = 1. #code 

ds_K = make_K("DD0000/DD0000") 
ds_kT = make_kT("DD0000/DD0000") 
ds_ne = make_ne("DD0000/DD0000")

sph_K = ds_K.sphere(center, (radius, "code_length")) 
sph_kT = ds_kT.sphere(center, (radius, "code_length"))
sph_ne = ds_ne.sphere(center, (radius, "code_length")) 

fig, ax = plt.subplots(2,2,sharex = 'col')
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)

#profile K 
field = "K"
profile = yt.create_profile(sph_K, ("index", "radius"), field, weight_field="ones", extrema={("index", "radius"):(1e22, 1e24)}, logs={"radius":True})
ax[0,0].loglog(profile.x.in_units("kpc").d, profile[field].d) 
ax[0,0].set_ylabel("K [KeV * cm**2]")

#profile kT
field = "kT"
profile = yt.create_profile(sph_kT, ("index", "radius"), field, weight_field="ones", extrema={("index", "radius"):(1e22, 1e24)}, logs={"radius":True})
ax[1,0].loglog(profile.x.in_units("kpc").d, profile[field].d) 
ax[1,0].set_ylabel("kT [keV]")

#profile ne 
field = "ne"
profile = yt.create_profile(sph_ne, ("index", "radius"), field, weight_field="ones", extrema={("index", "radius"):(1e22, 1e24)}, logs={"radius":True})
ax[0,1].loglog(profile.x.in_units("kpc").d, profile[field].d)
ax[0,1].yaxis.tick_right()
ax[0,1].set_ylabel("n_e [cm**-3]")
ax[0,1].yaxis.set_label_position("right")

#profile 2*ne*r
field = "ne"
profile = yt.create_profile(sph_ne, ("index", "radius"), field, weight_field="ones", extrema={("index", "radius"):(1e22, 1e24)}, logs={"radius":True})
ax[1,1].loglog(profile.x.in_units("kpc").d, 2 * profile.x.in_units("cm").d * profile[field].d)
ax[1,1].set_ylabel("2n_er [cm**2]")
ax[1,1].yaxis.tick_right()
ax[1,1].yaxis.set_label_position("right")

ax[0,0].axvline(x = 206, color = 'r', linestyle='--') 
ax[0,1].axvline(x = 206, color = 'r', linestyle='--') 
ax[1,0].axvline(x = 206, color = 'r', linestyle='--') 
ax[1,1].axvline(x = 206, color = 'r', linestyle='--') 

fig.savefig("frames/voit_comp.png") 





