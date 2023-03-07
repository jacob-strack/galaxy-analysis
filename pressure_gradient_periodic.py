from starter import *

def grad(data, fieldname, direction):
    iM1 = slice(None, -2)
    iP1 = slice(2,None)
    all = slice(1,-1)
    all_all=tuple([all]*3)
    dxi=1./(2*data.dds )
    out = np.zeros_like(data[fieldname]*dxi[0])
    Left = [all]*3
    Right = [all]*3
    Right[direction] = iP1
    Left[direction] = iM1
    Left=tuple(Left); Right=tuple(Right)
    out[all_all] = (data[fieldname][ Right ]- data[fieldname][ Left]) *dxi[direction]
    return out    

yt_dpdx = ("gas", "dpdx") 
yt_dpdy = ("gas", "dpdy")
yt_dpdz = ("gas", "dpdz") 

def add_pressure_gradient(obj):
    def _pressure_gradient_x(field, data): 
        out = data.ds.arr(grad(data, yt_pressure, 0), 'dyne/cm**3')
        return out
    obj.add_field(name = yt_dpdx,function = _pressure_gradient_x,validators = yt.ValidateSpatial(5, yt_pressure), units='dyne/cm**3', sampling_type = 'cell')
    def _pressure_gradient_y(field, data): 
        out = data.ds.arr(grad(data, yt_pressure, 1), 'dyne/cm**3')
        return out
    obj.add_field(yt_dpdy,function = _pressure_gradient_y,validators = yt.ValidateSpatial(5, yt_pressure), units='dyne/cm**3', sampling_type = 'cell')
    def _pressure_gradient_z(field, data): 
        out = data.ds.arr(grad(data, yt_pressure, 2), 'dyne/cm**3')
        return out
    obj.add_field(name = yt_dpdz,function = _pressure_gradient_z,validators = yt.ValidateSpatial(5, yt_pressure), units='dyne/cm**3', sampling_type = 'cell')

test_ds = yt.load("DD0001/DD0001")
add_pressure_gradient(test_ds)
ad = test_ds.all_data()
a = ad[yt_dpdx]
b = ad[yt_dpdy]
c = ad[yt_dpdz]
