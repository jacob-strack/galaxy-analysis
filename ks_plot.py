from starter import * 
from yt.data_objects.particle_filters import add_particle_filter 
import os 

def formed_stars(pfilter, data): 
    age = data.ds.current_time - data["all", "creation_time"]
    filter = age.in_units("code_time") <= data.ds["dtDataDump"]
    return filter 

def make_sfr(filename): 
    yt.add_particle_filter("formed_stars", function = formed_stars, filtered_type = "all", requires=["creation_time"])
    ds = yt.load(filename) 
    ds.add_particle_filter("formed_stars") 
    deposited = ds.add_deposited_particle_field(("formed_stars", "particle_mass"), method = "sum")
    ds.add_field(name=("gas", "sfr_density"), function=sfr_density, sampling_type="local", units="Msun/cm**3")
    print(deposited) 
    print(ds.r[deposited].in_units("Msun"))
    print(ds.r[deposited].max())
    return ds 

def sfr_density(field, data): 
    return data["deposit", "formed_stars_sum_mass"].in_units("Msun") / data["cell_volume"].in_units("pc**3") 

def ks_plot(filename, axis): 
    ds = make_sfr(filename)
    sfr_column = ds.proj("sfr_density", 0)
    column = ds.proj("density", 0)
    #flatten these projections and plot
    plt.scatter(column["density"].in_units("Msun/pc**2"), sfr_column["sfr_density"])
    plt.loglog()
    plt.savefig("test.png")

def ks_plot_all(axis, lower_t, num, lower_ind = 0, upper_ind = 0):
    plt.close('all')
    if lower_ind == 0 and upper_ind == 0: 
        lower_ind = lower_t 
        upper_ind = lower_t + num - 1
    ans_dens = yt.YTArray([], "Msun/pc**2")
    ans_sfr = yt.YTArray([], "Msun/pc**2/yr")
    ans_color = np.array([])
    yt.add_particle_filter("formed_stars", function = formed_stars, filtered_type = "all", requires=["creation_time"])
    for i in range(lower_t, lower_t + num):
        filename_lower = "DD" + str(i).zfill(4) + "/DD" + str(i).zfill(4)
        filename_upper = "DD" + str(i+1).zfill(4) + "/DD" + str(i+1).zfill(4)
        ds_lower = yt.load(filename_lower)
        ds_upper = yt.load(filename_upper)
        ds_lower.add_particle_filter("formed_stars")
        ds_upper.add_particle_filter("formed_stars") 
        deposited_lower = ds_lower.add_deposited_particle_field(("formed_stars", "particle_mass"), method = "sum")
        deposited_upper = ds_upper.add_deposited_particle_field(("formed_stars", "particle_mass"), method = "sum")
        ds_upper.add_field(name=("gas", "sfr_density"), function=sfr_density, sampling_type="local", units="Msun/cm**3")
        ds_lower.add_field(name=("gas", "sfr_density"), function=sfr_density, sampling_type="local", units="Msun/cm**3")
        upper_proj = ds_upper.proj("sfr_density", 2)
        lower_proj = ds_lower.proj("sfr_density", 2)
        upper_proj_fr = upper_proj.to_frb((1.0, "code_length"), 1024)
        lower_proj_fr = lower_proj.to_frb((1.0, "code_length"), 1024)
        sfr = (upper_proj_fr["sfr_density"] - lower_proj_fr["sfr_density"]) / (ds_upper.current_time - ds_lower.current_time).in_units("yr") 
        column = ds_upper.proj("density", 2)
        column_fr = column.to_frb((1.0, "code_length"), 1024)
        sfr_flat = sfr.in_units("Msun/pc**2/yr").flatten()
        dens_flat = column_fr["density"].in_units("Msun/pc**2").flatten()
        dens_pts = np.asarray(dens_flat)
        sfr_pts = np.asarray(sfr_flat) 
        sfr_pts = sfr_flat[np.where(sfr_flat > 0)]
        dens_pts = dens_flat[np.where(sfr_flat > 0)]
        ans_dens = np.append(ans_dens, dens_pts)
        ans_sfr = np.append(ans_sfr, sfr_pts)
        color_pts = np.ones_like(sfr_pts.d)*i
        ans_color = np.append(ans_color, color_pts)
    plt.scatter(ans_dens, ans_sfr, c = ans_color, cmap = 'jet', vmin = lower_t, vmax = lower_t + num - 1)
    #mask = (ans_color > lower_ind) * (ans_color <= upper_ind)
    sorted_dens_pts = np.sort((ans_dens))
    sorted_sfr_pts = (ans_sfr[np.argsort((ans_dens))])
    sorted_ans_color = ans_color[np.argsort(np.log10(ans_dens))]
    coefficients = np.polyfit(np.log10(sorted_dens_pts.d), np.log10(sorted_sfr_pts.d), 1)
    trendline = np.poly1d(coefficients)
    y_vals = trendline(np.log10(sorted_dens_pts.d))
    #plt.scatter(sorted_dens_pts, sorted_sfr_pts, c = sorted_ans_color, cmap = 'jet', vmin = lower_t, vmax = lower_t + num - 1)
    plt.plot(sorted_dens_pts.d, 10**y_vals, label = "n = " + str(np.round(coefficients[0],2)))
    plt.legend()
    plt.colorbar()
    print(coefficients)
    print(np.log10(sorted_dens_pts.d)) 
    print(np.log10(sorted_sfr_pts.d))
    print(np.log10(ans_dens))
    print(np.log10(ans_sfr))
    plt.ylabel(r"$\Sigma_{SFR}$")
    plt.xlabel(r"$\Sigma$")
    plt.loglog()
    plt.savefig("test.png")

def sfr_proj_plot(axis, lower_t, num): 
    for i in range(lower_t, lower_t + num): 
        yt.add_particle_filter("formed_stars", function = formed_stars, filtered_type = "all", requires=["creation_time"])
        filename = "DD" + str(i).zfill(4) + "/DD" + str(i).zfill(4)
        ds = yt.load(filename) 
        ds.add_particle_filter("formed_stars")
        deposited = ds.add_deposited_particle_field(("formed_stars", "particle_mass"), method="sum")
        ds.add_field(name=("gas", "sfr_density"), function=sfr_density, sampling_type="local", units="Msun/cm**3")
        prj = yt.ProjectionPlot(ds, axis, "sfr_density", width = (100,"kpc"))
        prj.save()
