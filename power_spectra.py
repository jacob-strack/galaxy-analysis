import yt 
import matplotlib.pyplot as plt 
import numpy as np
import matplotlib 

def power_spectra(timestep, level = 0): 
    #load dataset
    ds = yt.load("DD" + str(timestep).zfill(4) + "/DD" + str(timestep).zfill(4))
    ad = ds.all_data() 
    
    #use covering grid for uniform cells 
    low = [0.0, 0.0, 0.0]
    dims = ds.domain_dimensions * (2**level)
    cg = ds.covering_grid(level = level, left_edge = low, dims = dims) 
    width = ds.domain_dimensions * cg.dds.v
    
    #get magnetic field components of covering grid
    Bx = cg["Bx"].in_units("gauss")
    By = cg["By"].in_units("gauss")
    Bz = cg["Bz"].in_units("gauss")
    
    #ffts 
    bx_fft = np.fft.fftn(Bx)
    by_fft = np.fft.fftn(By)
    bz_fft = np.fft.fftn(Bz)
    kx = np.fft.fftfreq(dims[0], d=cg.dds[0].v)
    ky = np.fft.fftfreq(dims[1], d=cg.dds[1].v)
    kz = np.fft.fftfreq(dims[2], d=cg.dds[2].v)
    Kx, Ky, Kz = np.meshgrid(kx, ky, kz, indexing="ij")
    k = np.sqrt(Kx**2 + Ky**2 + Kz**2) 

    #power spectra components 
    res = (np.abs(bx_fft)**2 + np.abs(by_fft)**2 + np.abs(bz_fft)**2)
    kmin = np.min(1.0 / width) 
    kmax = np.min(0.5 * dims / width)
    k_vals = k.ravel() 
    power_vals = res.ravel()
    
    #make k bins 
    k_bins = np.arange(0.5, 100, 1.0) #was dims[0]. change to whatever max k you want, that makes sense
    k_center = (k_bins[:-1] + k_bins[1:])/2 
    power_spec, _ = np.histogram(k_vals, bins = k_bins, weights = power_vals) 
    
    #plot
    #plt.loglog(k_center, power_spec)
    #plt.savefig("spectrum.png") 
    return k_center, power_spec

def plot_time_power_spectra(start_time, end_time, level = 0): 
    jet = plt.get_cmap('jet')
    cNorm = matplotlib.colors.Normalize(vmin = 0, vmax = end_time-start_time) 
    scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=jet)
    for i in range(start_time, end_time): 
        kc, ps = power_spectra(i, level)
        plt.loglog(kc, ps, color = scalarMap.to_rgba(i))
    plt.colorbar()
    plt.savefig("timespec.png")

def plot_time_power_spectra(t, level = 0): 
    kc, ps = power_spectra(t, level)
    coefficients = np.polyfit(np.log10(kc), np.log10(ps), 1)
    trendline = np.poly1d(coefficients)
    y_vals = trendline(np.log10(kc))
    plt.loglog(kc, ps)
    plt.plot(kc, 10**y_vals, label = "n = " + str(np.round(coefficients[0],2)), linestyle = '--')
    plt.legend()
    plt.savefig("DD" + str(t).zfill(4) + "_Bspec.png")
