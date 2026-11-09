from starter import * 
from ne import * 

mu = 0.6 #mean molecular weight 
mh = 1.67e-24 
KEV_PER_ERG = 6.242e8 #conversion factor 
c = 10 #DM concentration 
M = 1e12 * 1.989e33 #g from solar masses 
rho_crit = 1.8788e-29*0.49
r_200 = (3 / (4*3.14159) * M / (200 * rho_crit))**(1./3.) 
G = 6.67e-8 
vc_max = 

def entropy_approx(r): 
    S_precip = (2 * mu * mh)**(1/3) * (r * 10 / 3)**(2/3) 
    return S_precip

def entropy_real(r): 
    S_precip = (2 * mu * mh)**(1/3) * r * 
