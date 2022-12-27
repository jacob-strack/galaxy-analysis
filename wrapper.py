from starter import * 
import sys 
from plotscripts import *


#Assume sysarg layout of plotscript, fieldname, dumpnumber
def plot_wrapper():
	plot = sys.argv[0](sys.argv[1], sys.argv[2])
	return plot
