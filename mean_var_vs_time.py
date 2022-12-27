import yt
import numpy as np 
import matplotlib.pyplot as plt
var = input("Plotting Variable:") #input a string which corresponds to a dataset variable 
t_arr = [] 
var_arr = [] 
for i in range(0,20): 
    filename ='DD' + str(i).zfill(4)
    ds = yt.load(filename + '/' + filename)
    ad = ds.all_data()
    temp_var_arr = ad[var]
    mean_var = np.mean(temp_var_arr)
    t = ds.current_time
    var_arr.append(mean_var)
    t_arr.append(t)
t_arr = np.array(t_arr)
var_arr = np.array(var_arr)
plt.plot(t_arr, var_arr)
plt.savefig('frames/meanvstime.png') 

