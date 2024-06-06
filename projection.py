import yt
from yt import YTArray
import healpy as hp
from yt.fields.api import ValidateParameter
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.colors as colors 
from mpi4py import MPI
import os 
import trident
import shutil
from yt.funcs import mylog
import math
import scipy as sp
#####################################################################


mylog.setLevel(40)  # only logs errors and critical problems

def _Q_integrand(field,data):
    theta = data['theta']
    phi = data['phi']
    B_x = data['magnetic_field_x']*np.cos(theta)*np.cos(phi) + data['magnetic_field_y']*np.cos(theta)*np.sin(phi) - data['magnetic_field_z']*np.sin(theta) 
    B_y = -data['magnetic_field_x']*np.sin(phi) + data['magnetic_field_y']*np.cos(phi)
    return data['density']*(B_x**2 - B_y**2)/(data['magnetic_field_strength'])**2 

def _U_integrand(field,data): 
    theta = data['theta']
    phi = data['phi']
    B_x = data['magnetic_field_x']*np.cos(theta)*np.cos(phi) + data['magnetic_field_y']*np.cos(theta)*np.sin(phi) - data['magnetic_field_z']*np.sin(theta) 
    B_y = -data['magnetic_field_x']*np.sin(phi) + data['magnetic_field_y']*np.cos(phi)
    return data['density']*2*B_x*B_y/(data['magnetic_field_strength'])**2

def _theta(field,data):
    domain_center = (data.ds.domain_right_edge - data.ds.domain_left_edge) / 2
    if data.has_field_parameter("obs_center"):
        obs_cent_arr = data.get_field_parameter("obs_center")
        center = domain_center + obs_cent_arr
    else: 
        center = domain_center
    z = data[("gas","z")] - center[2].in_units("cm")
    y = data[("gas","y")] - center[1].in_units("cm")
    x = data[("gas","x")] - center[0].in_units("cm")
    ans = YTArray(np.arccos(z/np.sqrt(x**2 + y**2+z**2)), 'dimensionless')
    return ans
def _phi(field,data):
    domain_center = (data.ds.domain_right_edge - data.ds.domain_left_edge) / 2
    if data.has_field_parameter("obs_center"):
        obs_cent_arr = data.get_field_parameter("obs_center")
        center = domain_center + obs_cent_arr
    else: 
        center = domain_center
    z = data[("gas","z")] - center[2].in_units("cm")
    y = data[("gas","y")] - center[1].in_units("cm")
    x = data[("gas","x")] - center[0].in_units("cm")
    original_shape = z.shape
    flat_x = x.flatten()
    flat_y = y.flatten()
    flat_z = z.flatten()
    res = np.arctan2(flat_y,flat_x)
    res = np.reshape(res, original_shape)
    ans = YTArray(res, 'dimensionless')
    return ans
def E_integrand(dsname,nside, comm, int_LOS=False):
    data = yt.load(dsname)
    Q_map = enzo_to_healpix(dsname, "Q_integrand",nside,comm,True)
    U_map = enzo_to_healpix(dsname, "U_integrand", nside,comm,True)
    I_map = enzo_to_healpix(dsname, "U_integrand", nside, comm, True) #I think this is fine?
    print(np.shape(Q_map))
    print(np.shape(U_map))
    print(np.shape(I_map))
    if comm.Get_rank() == 0: 
        alms = hp.sphtfunc.map2alm([I_map,Q_map,U_map])
        E_map = hp.sphtfunc.alm2map(alms[1],nside)
        return E_map
def B_integrand(dsname, nside, comm, int_LOS=False): 
    data = yt.load(dsname)
    Q_map = enzo_to_healpix(dsname, "Q_integrand",nside,comm,True)
    U_map = enzo_to_healpix(dsname, "U_integrand", nside,comm,True)
    I_map = enzo_to_healpix(dsname, "U_integrand", nside, comm, True)
    if comm.Get_rank() == 0:
        alms = hp.sphtfunc.map2alm([I_map,Q_map,U_map])
        B_map = hp.sphtfunc.alm2map(alms[2],nside)
        return B_map

yt.add_field("theta", function = _theta, sampling_type = "cell", units = 'dimensionless',validators=[ValidateParameter(["obs_center"])])
yt.add_field("phi", function = _phi, sampling_type = "cell", units = 'dimensionless',validators=[ValidateParameter(["obs_center"])])
yt.add_field("Q_integrand", function=_Q_integrand,sampling_type="cell", units = "g/cm**3")
yt.add_field("U_integrand", function=_U_integrand,sampling_type="cell", units = "g/cm**3")
def ray_initialize(ds, c, obs, ray_len): 
    center = ds.arr(c, 'code_length')
    center = center.in_units('kpc')
    observation_loc = ds.arr(obs,'kpc')
    observation_loc = observation_loc.in_units('kpc')
    length = ray_len
    ray_start = center - observation_loc
    return center, observation_loc, length, ray_start

def generate_dicts(fields,comm):
    dicts = [] 
    for field in fields: 
        name = str(field)
        new_dir = {'field':field, 'file_name':field, 'title':field, 'colormap':'jet'}
    dicts.append(new_dir)
    if comm.rank == 0:
        if os.path.exists('ray_files'):
            shutil.rmtree('ray_files')
            os.mkdir('ray_files')
        else:
            os.mkdir('ray_files')
    return dicts 

def generate_angles(pixels,comm): 
    dphi = 2.0*np.pi / (pixels)
    dtheta = np.pi / (pixels)
    phi,theta = np.mgrid[slice(-np.pi +dphi/2,np.pi+dphi/2,dphi),slice(-np.pi/2 + dtheta/2,np.pi/2+dtheta/2,dtheta)]
    unflat_shape = phi.shape
    phi = np.reshape(phi,-1)
    theta = np.reshape(theta,-1)
    dN = phi.size // comm.size
    print(phi,theta)
    start_index = comm.rank*dN 
    end_index = (comm.rank+1)*dN
    print("theta min/max generate_angles", theta.min(), theta.max(), phi.min(),phi.max())
    dN = theta.size // comm.size

    if (theta.size % comm.size == 0):
        if 1:
            print("Task", comm.rank, "no leftovers!")
    else:
        print("Task", comm.rank, "your array size and MPI tasks do not divide evenly!")
        comm.Abort(errorcode=123)
    return theta,phi,unflat_shape,start_index, end_index

def generate_ray(ds,field_list,theta, phi, R, ray_start,start_index,end_index,i,parameters):
    #dx = R*np.cos(theta)*np.cos(-phi)#did have negatives on phi, maybe not needed?
    #dy = R*np.cos(theta)*np.sin(-phi)
    #dz = R*np.sin(theta)
    dx = R*np.sin(theta)*np.cos(phi)
    dy = R*np.sin(theta)*np.sin(phi)
    dz = R*np.cos(theta)
    padded_num = '{:06d}'.format(i)
    rayfile = 'ray_files/ray'+str(ds)+'_'+padded_num+'.h5'
    delta = YTArray([dx,dy,dz],'kpc')
    ray_end = ray_start + delta
    lr = trident.LightRay(ds)
    ray = lr.make_light_ray(start_position=ray_start, end_position=ray_end,data_filename=rayfile, fields = field_list,field_parameters = parameters)
    ray_data = ray.all_data()
    return ray_data

def calculate_path(ray_data,remove_first_N_kpc=1): 
    path_length = ray_data['dl'].in_units('kpc').d
    path = np.zeros(len(path_length))
    for h in range(len(path_length) - 1): 
        dl = path_length[h]
        p = path[h]
        path[h+1] = dl + p
    for g in range(len(path)):
        if path[g] < remove_first_N_kpc:
            continue
        else:
            start = g
        break
    return start,path

def column_density(field_name,path,field,field_column,ray_data,start,i,int_path_weight=False):
    path_mod = path[start:]
    field_number_density_mod = ray_data[field_name][start:]
    field_density = ray_data[field_name]*ray_data['dl']
    field_density_mod = field_density[start:]
    field_column_density = sum(field_density_mod.d)
    if int_path_weight == False:
        field_column_density /= ray_data['dl'].sum()
    field_column = np.append(field_column, field_column_density)
    field[i] += field_column_density
    return field_column,field

def aitoff_projection(theta,phi,field,original_shape,comm,savename): 
    theta1 = np.reshape(theta,original_shape)
    print("reduce phi")
    phi = np.reshape(phi,original_shape)
    print("reduce_field")
    print(field,comm.rank)
    reduce_field = np.zeros_like(field)
    comm.Reduce(field,reduce_field,op=MPI.SUM)
    print("reshape field", reduce_field)
    #make projection plot on root grid
    if comm.Get_rank() == 0:
        print("make figure")
        #field = np.reshape(reduce_field, original_shape)
        print("project")
        print("pcolor")
        hp.visufunc.mollview(reduce_field)
        #plt.pcolor((phi-np.pi),(theta-(np.pi/2)),(theta-np.pi/2),cmap='Paired')
        #cbar=plt.colorbar(pad=0.02,shrink=0.55)
        plt.grid(True,alpha=0.5)
        plt.xticks(fontsize=7)
        plt.yticks(fontsize=7)
        print("saving figure")
        plt.savefig(savename + '.png')
        print("done save. closing")
        plt.close()
    return phi,theta,field

def create_aitoff_projection(dsn,field_name, pixels_per_dim,savename):
    comm = MPI.COMM_WORLD
    field_column = np.array([])
    ds = yt.load(dsn)
    dicts = generate_dicts(['magnetic_field_x','magnetic_field_y','magnetic_field_strength','theta','phi',field_name],comm)
    theta,phi,unflat_shape,start_index, end_index = generate_angles(pixels_per_dim,comm)
    field_arr = np.zeros_like(theta)
    c, loc, length, ray_start = ray_initialize(ds, [0.5,0.5,0.5], [8.,0.,0.],200.)
    for i in range(start_index, end_index):
        ray_data = generate_ray(ds,[field_name,'magnetic_field_x','magnetic_field_y','magnetic_field_z','magnetic_field_strength','density','theta','phi'],theta[i],phi[i],length,ray_start,start_index,end_index,i)
        start,path = calculate_path(ray_data)
        col,field = column_density(field_name,path,field_arr,field_column,ray_data,start,i)
    phi1,theta1,field = aitoff_projection(theta,phi,field,unflat_shape,comm,savename)
    return field

def QU_to_EB(Q,U): 
    Q_F = np.fft.fft2(Q)
    U_F = np.fft.fft2(U)
    k_x = np.fft.fftfreq()
    k_y = np.fft.fftfreq()
    theta = np.arccos(k_x/np.sqrt(k_x**2 + k_y**2))
    res = (Q_F + j*U_F)*np.exp(2*j*theta)
    E_F = np.real(res) 
    B_F = np.imag(res) 
    E = np.fft.ifft2(E_F)
    B = np.fft.ifft2(B_F)
    return E,B

def EB_Projection(dsname,field_name,nside,savename,comm): 
    #fuction to take in a ds name, make frb, use healpix to make EB projection
    ds = yt.load(dsname)
    dicts = generate_dicts(['magnetic_field_x','magnetic_field_y','magnetic_field_strength','theta','phi',field_name],comm)
    theta,phi,unflat_shape,start_index,end_index = generate_angles_healpix(2**20, comm)
    print("theta",theta,theta.shape,type(theta))
    print("phi",phi,phi.shape,type(phi))
    field_column = np.array([])
    dicts = generate_dicts(['magnetic_field_x','magnetic_field_y','magnetic_field_strength',field_name],comm)
    field_arr = np.zeros_like(theta)
    c, loc, length, ray_start = ray_initialize(ds, [0.5,0.5,0.5], [8.,0.,0.],200.)
    for i in range(start_index, end_index):
        ray_data = generate_ray(ds,[field_name,'magnetic_field_x','magnetic_field_y','magnetic_field_z','magnetic_field_strength','density'],phi[i],theta[i],length,ray_start,start_index,end_index,i,"obs_center")
        start,path = calculate_path(ray_data)
        col,field_arr = column_density(field_name,path,field_arr,field_column,ray_data,start,i)
    print("Done with loop")
    phi1,theta1,field = aitoff_projection(theta,phi,field_arr,unflat_shape,comm,savename)

def plot_phi(dsname,field_name,nside,savename,comm,int_LOS=False): 
    ds = yt.load(dsname)
    theta,phi,unflat_shape,start_index,end_index = generate_angles_healpix(nside, comm)
    print("theta", theta.max(), theta.min())
    field_column = np.array([])
    field_arr = np.zeros_like(theta)
    c, loc, length, ray_start = ray_initialize(ds, [0.5,0.5,0.5], [8.,0.,0.],200.)
    ad = ds.all_data()
    ad.set_field_parameter("obs_center", YTArray([8.,0.,0.],"kpc"))
    for i in range(start_index, end_index):
        ray_data = generate_ray(ds,[field_name,'magnetic_field_x','magnetic_field_y','magnetic_field_z','magnetic_field_strength','density'],theta[i],phi[i],length,ray_start,start_index,end_index,i,{"obs_center":YTArray([8.,0.,0.],"kpc")})
        start,path = calculate_path(ray_data)
        col,field_arr = column_density(field_name,path,field_arr,field_column,ray_data,start,i,int_LOS)
    reduce_field = np.zeros_like(phi)
    comm.Reduce(field_arr,reduce_field,op=MPI.SUM)
    if comm.Get_rank() == 0:
        hp.mollview(reduce_field)
        hp.visufunc.graticule()
        plt.savefig(savename + ".png")
        print("theta", theta.max(), theta.min(),"phi", phi.max(),phi.min())

def enzo_to_healpix(dsname, field_name, nside, comm, int_LOS=False):
    ds = yt.load(dsname)
    theta, phi, unflat_shape, start_index, end_index = generate_angles_healpix(nside,comm)
    field_column = np.array([])
    field_arr = np.zeros_like(theta)
    c, loc, length, ray_start = ray_initialize(ds, [0.5,0.5,0.5], [8.,0.,0.],200.)
    ad = ds.all_data()
    ad.set_field_parameter("obs_center", YTArray([8.,0.,0.], "kpc"))
    for i in range(start_index,end_index): 
        ray_data = generate_ray(ds,[field_name,'magnetic_field_x','magnetic_field_y','magnetic_field_z','magnetic_field_strength','density'],theta[i],phi[i],length,ray_start,start_index,end_index,i,{"obs_center":YTArray([8.,0.,0.],"kpc")})
        start,path = calculate_path(ray_data)
        col,field_arr = column_density(field_name,path,field_arr,field_column,ray_data,start,i,int_LOS)
    red = np.zeros_like(phi)
    comm.Reduce(field_arr,red,op=MPI.SUM)
    print(red)
    if comm.Get_rank() == 0:
        return red

def healpix_generate_angles(nside): 
    npix = hp.nside2npix(nside)
    arr = np.arange(npix)
    angles = hp.pix2ang(nside, arr,lonlat = False)
    return angles[1], angles[0]

def generate_angles_healpix(nside,comm): 
    phi,theta = healpix_generate_angles(nside)
    unflat_shape = phi.shape
    dN = phi.size // comm.size
    start_index = comm.rank*dN 
    end_index = (comm.rank+1)*dN
    dN = theta.size // comm.size
    if (theta.size % comm.size == 0):
        if 1:
            print("Task", comm.rank, "no leftovers!")
    else:
        print("Task", comm.rank, "your array size and MPI tasks do not divide evenly!")
        comm.Abort(errorcode=123)
    return theta,phi,unflat_shape,start_index, end_index

def fast_ray_caster(ds, field_name,nside, comm): 
    ad = ds.all_data()
    ad.set_field_parameter("obs_center", YTArray([0.,0.,0.], "kpc"))
    phi, theta = healpix_generate_angles(nside)
    phi_edges = np.unique(phi)
    theta_edges = np.unique(theta)
    phi_edges = phi
    theta_edges = theta
    theta_arr = ad["theta"] 
    phi_arr = ad["phi"] 
    field_arr = ad[field_name] 
    indices = hp.ang2pix(nside, theta_arr.d, phi_arr.d)
    hist = sp.stats.binned_statistic(indices, field_arr,bins = hp.nside2npix(nside),statistic = "mean")
    hp.mollview(hist[0])
    hp.graticule()
    plt.savefig("fasttest" + field_name + ".png")
    return hist[0]

comm = MPI.COMM_WORLD
#plot_phi("DD0100/DD0100","E_integrand",4,'Etest',comm,int_LOS=True)
res = fast_ray_caster(yt.load("DD0100/DD0100"),"x",512, comm)
res = fast_ray_caster(yt.load("DD0100/DD0100"),"y",64, comm)
res = fast_ray_caster(yt.load("DD0100/DD0100"),"z",64, comm)
res = fast_ray_caster(yt.load("DD0100/DD0100"),"theta",8, comm)
res = fast_ray_caster(yt.load("DD0100/DD0100"),"phi",64, comm)
