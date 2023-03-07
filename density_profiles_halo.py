from starter import * 
#convert stuff to pure cgs units
def NFW(r, rho_0, r_s):
    return rho_0 / (r/r_s * (1 + r/r_s)**2)
def Hernquist(r, rho_0, r_s): 
    return rho_0 / (r/r_s * (1 + r/r_s)**3)
def Rathjen_2021(r):
    rho_0 = 128.2 * 1.988435e33 / (3.086e21)**3 #solar mass / kpc**3 i think this is the right number 
    r_200 = 3.086e21 #cm did have a 200 * here but might have been busted
    return NFW(r, rho_0 , r_200 / 12 ) 
   
def Steinwandel_2021(r): 
    rho_0 = 5e-26 #g/cm**3
    beta = 2/3
    r_s = 3.086e21  # 1 kpc in cm
    return Hernquist(r, rho_0, r_s) 

def Pakmor_2013(r, z):#coming back to this, the problem uses dimensionless quantities, but i can't find the scaling factors used, maybe buried in another paper? 
    M_200 = 10**12 * 1.988e33 #g
    c = 7.2
    rho_0 = (200 / 3) * (c**3 / (np.log(1 + c) - c/(1+c))) * 8.67e-30 #g/cm**3
    r_200 = (M_200 / (4 * np.pi * rho_0  * (np.log(1 + c) - c/(1+c))))**(1/3) #cm
    r_s = r_200 / c #cm 
    ans = NFW(r, rho_0, r_s)
    return ans

def phi_ext(z): 
    z_star = 7.5607e20 #cm 
    sigma_star = 8.767e-3 #g/cm^2
    rho_dm = 4.329e-26 #g/cm^3
    R_0 = 2.4688e22 #cm
    G = 6.6743e-8 #cm^3/g/s^2
    return 2*np.pi*G*sigma_star*z_star*(np.sqrt(1 + z**2/z_star**2) - 1) + 2*np.pi*G*rho_dm*R_0**2*np.log(1 + z**2/R_0**2)
def phi_0(z): 
    sigma = 13 * 1.988e33 / (3.086e18)**2 #g/cm^2
    G = 6.6743e-8 #cm/g/s^2
    return phi_ext(z) + 2 * np.pi * sigma * G * np.abs(z)
def Kim_2017(z): 
    m_h = 1.673e-24 #g  
    rho_1 = 2.85 * m_h
    rho_2 = 10**-5 * rho_1
    sigma_1 = 700000 #cm/s
    sigma_2 = 10 * sigma_1
    ans = rho_1 * np.exp(-phi_0(z) / sigma_1**2) + rho_2 * np.exp(-phi_0(z) / sigma_2**2)
    return ans
def Butsky_2017(r, z): 
    M_200 = 1.074e12 * 1.988e33 #g
    c = 10
    M_d = 4.297e10 * 1.988e33 #g
    r_d = 3.432 * 3.086e21 # cm
    z_d = 343.2 * 3.086e21 # cm
    rho_0 = M_d * .2 / (4 * np.pi * r_d**2 * z_d)
    R_200 = (3 * M_d / (4 * np.pi * 200 * rho_0))**(1/3)
    return NFW(r,rho_0,R_200/c )
def Wibking_2021(r): 
    M_200 = 1.254e5 * 1.988e33 #g 
    R_200 = 205.5 * 3.086e21 
    c = 10
    R_s = R_200 / 10
    rho_0 = (3 * M_200 / (4 *np.pi * R_200**3))**(1/3)
    return NFW(r, rho_0, R_s)
def plot(): 
    plt.close('all')
    r = np.linspace(1, 10**23 ,num = 1000000) 
    plt.plot(r, Rathjen_2021(r), label = 'Rathjen 2021')
    plt.plot(r, Steinwandel_2021(r), label = 'Steinwandel 2021')
    #plt.plot(r, Pakmor_2013(10**-15, r), label = 'Pakmor 2013')
    plt.plot(r, Kim_2017(r), label = 'Kim 2017')
    plt.plot(r, Butsky_2017(r,r), label = 'Butsky 2017')
    plt.plot(r, Wibking_2021(r), label = 'Wibking 2021')
    plt.loglog()
    plt.legend()
    plt.ylabel('Density [M_s/kpc^3]')
    plt.xlabel('r [kpc]')
    plt.savefig("Halo_Density_Profiles.png")
