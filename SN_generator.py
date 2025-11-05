import yt 
import numpy as np 

num = 150 

box_size = 2 * 3.147e20 #1 kpc in cm 

t_max = 2e14 #max time for star birth 

v_max = 1e6 #cm/s, max component of v

with open("random_SN.txt", "w") as file: 
    file.write(str(num) + "\n") #AMReX expects first line to contain number of particles contained in file 
    for i in range(num): 
        x = (np.random.random() - 0.5) * box_size
        y = (np.random.random() - 0.5) * box_size
        z = (np.random.random() - 0.5) * box_size
        m = 1.988409871e38
        v_x = (np.random.random() - 0.5) * 2 * v_max
        v_y = (np.random.random() - 0.5) * 2 * v_max
        v_z = (np.random.random() - 0.5) * 2 * v_max
        birth_time = np.random.random() * t_max * 0 + 1
        death_time = birth_time + np.random.random()*3e13
        lum = 0. 
        E_alpha = 0. 
        E_beta = 0. 
        E_gamma = 0. 
        to_write = [x, y, z, m, v_x, v_y, v_z, birth_time, death_time, lum, E_alpha, E_beta, E_gamma] 
        for j in to_write: 
            file.write(str(j) + " ")
        file.write("\n")
