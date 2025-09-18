import yt 
import numpy as np 
import sys

t_max = 100 
axis = 'z' 

for t in range(0,100): 
    try: 
        print(len(sys.argv))
        if(len(sys.argv) == 4): 
            field = sys.argv[1]
            axis = sys.argv[2]
            c = sys.argv[3]
            print(axis)
            ds = yt.load('DD' + str(t).zfill(4) + '/DD' + str(t).zfill(4))
            proj = yt.ProjectionPlot(ds, axis, field, width = (200,"kpc"), weight_field = "ones", center=[c, c, c])
            proj.annotate_timestamp(corner = "lower_right", redshift = False, draw_inset_box = True)
            proj.save('frames/DD' + str(t).zfill(4) + '_Projection_' + axis + '_' + field + '.png')
        if(len(sys.argv) == 3): 
            field = sys.argv[1]
            axis = sys.argv[2]
            print(axis)
            ds = yt.load('DD' + str(t).zfill(4) + '/DD' + str(t).zfill(4))
            proj = yt.ProjectionPlot(ds, axis, field, width = (200,"kpc"), weight_field = "ones")
            proj.annotate_timestamp(corner = "lower_right", redshift = False, draw_inset_box = True)
            proj.save('frames/DD' + str(t).zfill(4) + '_Projection_' + axis + '_' + field + '.png')
        if(len(sys.argv) == 2):
            field = sys.argv[1]
            ds = yt.load('DD' + str(t).zfill(4) + '/DD' + str(t).zfill(4))
            proj = yt.ProjectionPlot(ds, axis, field, width = (200,"kpc"), weight_field = "ones")
            proj.annotate_timestamp(corner = "lower_right", redshift = False, draw_inset_box = True)
            proj.save('frames/DD' + str(t).zfill(4) + '_Projection_' + axis + '_' + field + '.png')
        if(len(sys.argv) == 1): 
            ds = yt.load('DD' + str(t).zfill(4) + '/DD' + str(t).zfill(4))
            proj = yt.ProjectionPlot(ds, axis, "density", width = (200,"kpc"), weight_field = "ones")
            proj.annotate_timestamp(corner = "lower_right", redshift = False, draw_inset_box = True)
            proj.save('frames/DD' + str(t).zfill(4) + '_Projection_' + axis + '_density.png')
    except: 
        print("No data dump for ",t)
