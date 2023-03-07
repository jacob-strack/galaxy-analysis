from starter import * 
def Steinwandel_2021(r): 
    rho_0 = 5e-26 #g/cm**3
    beta = 2/3
    r_s = 3.086e21  # 1 kpc in cm
    return rho_0 * (1 + r**2/r_s**2)**(-beta)

 
def Butsky_2017(r, z): 
    M_200 = 1.074e12 * 1.988e33 #g
    c = 10
    M_d = 4.297e10 * 1.988e33 #g
    r_d = 3.432 * 3.086e21 # cm
    z_d = 343.2 * 3.086e21 # cm
    rho_0 = M_d * .2 / (4 * np.pi * r_d**2 * z_d)
    R_200 = (3 * M_d / (4 * np.pi * 200 * rho_0))**1/3
    return rho_0 * np.exp(-r/r_d) * np.exp(-np.abs(z) / z_d)

def Wibking_2021(r, z): 
    R_0 = 3.43218 * 3.086e21 #cm 
    z_0 = 0.343218 * 3.086e21 #cm 
    M_d = 8.593e9 * 1.988e33 #g
    rho_0 = M_d / (4 * np.pi * R_0**2 * z_0) 
    return rho_0 * np.exp(-r/R_0) * np.exp(-np.abs(z) / z_0) 

def plot(): 
    plt.close('all')
    r = np.linspace(1, 10**23 ,num = 1000000) 
    plt.plot(r, Butsky_2017(0,r), label = 'Butsky 2017')
    plt.plot(r, Steinwandel_2021(r), label = 'Steinwandel 2021')
    plt.plot(r, Wibking_2021(0, r), label = 'Wibking 2021')
    plt.loglog()
    plt.legend()
    plt.ylabel('Density [g/cm^3]')
    plt.xlabel('z [cm]')
    plt.savefig("Disk_Density_Profiles.png")
