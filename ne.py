from starter import * 

def make_ne(filename): 
    ds = yt.load(filename)
    ds.add_field(name = ("gas", "ne"), function = ne, sampling_type="local", units="1/cm**3")
    return ds 

def ne(field,data): 
    return data["Electron_Density"].in_units("kg/cm**3") / yt.YTQuantity(1836.15267 * 9.11e-31, "kg") #the electron density field in enzo isn't the actual electron density

def make_kT(filename): 
    ds = yt.load(filename) 
    ds.add_field(name = ("gas", "kT"), function = kT, sampling_type="local", units="keV") 
    return ds 

def kT(field,data): 
    return (data["temperature"].in_units("K") * yt.YTQuantity(8.617333e-5, "eV/K")).in_units("keV")

def make_K(filename): 
    ds = yt.load(filename)
    ds.add_field(name = ("gas", "K"), function = K, sampling_type="local", units="keV * cm**2")
    return ds 

def K(field,data): 
    return (yt.YTQuantity(8.617333e-5, "eV/K") * data["temperature"].in_units("K") * ((data["Electron_Density"]/yt.YTQuantity(1836.15267, "dimensionless")).in_units("kg/cm**3") / yt.YTQuantity(9.11e-31, "kg"))**(-2/3)).in_units("keV * cm**2")

def make_n(filename): 
    ds = yt.load(filename)
    ds.add_field(name=("gas", "n"), function = n, sampling_type="local", units="1/cm**3") 
    return ds 

def n(field,data): 
    return (data["Electron_Density"].in_units("kg / cm**3") / yt.YTQuantity(1836.15267 * 9.11e-31, "kg") + data["HI_Density"].in_units("kg / cm**3") / yt.YTQuantity(1.67e-27, "kg") + data["HII_Density"].in_units("kg/cm**3") / yt.YTQuantity(1.67e-27, "kg") + data["HeI_Density"].in_units("kg / cm**3") / yt.YTQuantity(6.64e-27, "kg") + data["HeII_Density"].in_units("kg/cm**3") / yt.YTQuantity(6.64e-27, "kg") + data["HeIII_Density"].in_units("kg/cm**3") / yt.YTQuantity(6.64e-27, "kg")) 

def make_ni(filename): 
    ds = yt.load(filename)
    ds.add_field(name=("gas", "ni"), function = ni, sampling_type="local", units="1/cm**3")
    return ds

def ni(field, data): 
    return (data["HII_Density"].in_units("kg / cm**3") / yt.YTQuantity(1.67e-27, "kg") + data["HeII_Density"].in_units("kg/cm**3") / yt.YTQuantity(6.64e-27, "kg") + data["HeIII_Density"].in_units("kg/cm**3") / yt.YTQuantity(6.64e-27, "kg"))

def make_nHe(filename): 
    ds = yt.load(filename) 
    ds.add_field(name=("gas", "nHe"), function = nHe, sampling_type="local", units="1/cm**3") 
    return ds 

def nHe(field, data): 
    return (data["HeI_Density"] + data["HeII_Density"] + data["HeIII_Density"]).in_units("kg/cm**3") / yt.YTQuantity(6.64e-27, "kg")

def make_nH(filename): 
    ds = yt.load(filename)
    ds.add_field(name=("gas", "nH"), function = nH, sampling_type="local", units="1/cm**3")
    return ds 

def nH(field, data): 
    return (data["HI_Density"] + data["HII_Density"]).in_units("kg/cm**3") / yt.YTQuantity(1.67e-27, "kg")
