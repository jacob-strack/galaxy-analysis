import yt 
import numpy as np 
import os 

#script that plots a bunch of density, temperature, magnetic field, etc. images and stores them in the right place automatically

fields = ['density', 'magnetic_field_strength', 'Temperature'] 

t = 0
axis = "z"
'''
while True:
    print(t)
    try:
        dirname = 'DD' + str(t).zfill(4)
        print("trying to load ds")
        ds = yt.load(dirname + '/' + 'DD' + str(t).zfill(4))
        print("loaded ds")
    except: 
        break
    for field in fields: 
        proj = yt.ProjectionPlot(ds, axis, field, weight_field = "ones")
        filename = dirname + field + axis + ".png"
        proj.save('frames/' + filename)
        proj_zoom = yt.ProjectionPlot(ds, axis, field, width = (100,'kpc'), weight_field="ones")
        filename_zoom = dirname + field + axis + "zoom.png"
        proj_zoom.save('frames/' + filename_zoom)
        print("saving plot ...")
        curr_dir = os.getcwd() 
        os.chdir(curr_dir + '/frames')
        if not os.path.exists(field):
            os.mkdir(field) 
        if not os.path.exists(field + "_zoom"):
            os.mkdir(field + "_zoom")
        os.rename(filename, field + '/' + filename)
        os.rename(filename_zoom, field + '_zoom/' + filename_zoom)
        os.chdir(curr_dir)
    t += 1
    '''
t_max = 101
part_type = 1
print("trying particle plots . . .")
for t in range(0,t_max):
    print(t)
    try:
        field = "density"
        dirname = 'DD' + str(t).zfill(4)
        ds = yt.load(dirname + '/DD' + str(t).zfill(4))
        filename = dirname + field + axis + "_particles.png"
        proj = yt.ProjectionPlot(ds, axis, "density", width = (200,'kpc'), weight_field = "ones")
        proj.annotate_particles((200,'kpc'),ptype = part_type)
        proj.save('frames/' + filename)
    except:
        print("No particles found!")

