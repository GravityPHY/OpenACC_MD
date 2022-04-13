import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import math
plt.rcParams.update({'font.size':14})

def read_dat(address):
    file=open(address,'r')
    N=[]
    time=[]
    first_line = file.readline()
    for line in file:
        nums=line.split()
        N.append(int(nums[0]))
        time.append(float(nums[1]))
    return N, time

def loglog_plot(title,x,y,):
    fig, ax = plt.subplots(1, 1,figsize=(9, 8))
    ax.set_title(title)
    ax.plot(x, y, 'bo', )
    fig.savefig(title+'.png')