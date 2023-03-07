from starter import * 
from mpl_toolkits.axes_grid1 import AxesGrid
import pcolormesh_helper as pch
profiles = []
labels = []
plot_specs = []
#halo_types = [0, 1, 2, 3, 4, 5, 6, 7,11]
halo_types = [2, 6]
fields = ["entropy", "pressure", "temperature", "density"]
fields_w_labels = ["entropy", "pressure", "temperature", "density"]
fig, axes = plt.subplots(2,3,figsize = (25, 12))
list_of_hist = []
list_of_edge = []
for i in halo_types:
    hist = {}
    edges = {}
    for j,field in enumerate(fields): 
        ds = yt.load("DD0000/HT" + str(i) + "0000")
        ax = axes.flatten()[j]
        sph = ds.sphere("c", (100, "kpc"))
        hist[field], edges[field], centers = pch.simple_profile(sph["radius"], sph[field], ax = ax, label = "HT" + str(i))
        ax.set(yscale = "log", xscale = "log", ylabel = field, xlabel = "radius")
        handles, labels = ax.get_legend_handles_labels()
        if j == 3:
            list_of_hist.append(hist)
            list_of_edge.append(edges)
for i in range(len(halo_types)):
    dTdr = []
    r = []
    dPdr = []
    for k in range(len(list_of_hist[i]["temperature"]) - 1): 
        print((list_of_hist[i]["temperature"][k+1] - list_of_hist[i]["temperature"][k]), list_of_edge[i]["temperature"][k + 1] - list_of_edge[i]["temperature"][k])
        dTdr.append((list_of_hist[i]["temperature"][k + 1] - list_of_hist[i]["temperature"][k]) / (list_of_edge[i]["temperature"][k+1] - list_of_edge[i]["temperature"][k]))
        r.append(list_of_edge[i]["temperature"][k+1] - list_of_edge[i]["temperature"][k])
        a =np.abs((list_of_hist[i]["pressure"][k + 1] - list_of_hist[i]["pressure"][k]) / (list_of_edge[i]["pressure"][k + 1] - list_of_edge[i]["pressure"][k]))
        dPdr.append((list_of_hist[i]["temperature"][k] / list_of_hist[i]["pressure"][k]) * a * 2 / 5)
    dPdr = np.asarray(dPdr)
    dTdr = np.asarray(dTdr)
    r = np.asarray(r)
    axes[1][2].plot(r, dTdr) 
    axes[1][1].plot(r, dPdr) 
    axes[1][2].set(yscale = "log", xscale = "log", ylabel = "dTdr", xlabel = "radius")
    axes[1][1].set(yscale = "log", xscale = "log", ylabel = "T / P * dPdr * 2 / 5", xlabel = "radius")
fig.legend(handles, labels, loc='upper center',ncol = len(halo_types))
fig.savefig("frames/2x2_EntopyPlot26.png") 
        
