import yt 

directory = "/scratch3/jstrack/star2"
frame = 17
filename = "%s/DD%04d/data%04d"%(directory,frame,frame)
cpufile = "%s.cpu0000"%(filename)


def read_stuff(filename): 
	fptr = h5py.File(filename,"r") 
        output = dict()
	try: 
		for field in fptr: 
			output[field] = fptr[field][()]
	except:
		raise
	finally:
		fptr.close()	
	return output

