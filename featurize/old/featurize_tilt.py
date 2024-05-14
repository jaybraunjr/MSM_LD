
import matplotlib.pyplot as plt
# import matplotlib as mpl
import pyemma
from pyemma.util.contexts import settings
import mdtraj
import itertools
from itertools import combinations
import pandas as pdnp
import numpy as np
from pyemma.coordinates import source
import mdtraj as md
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.core.groups import AtomGroup
import math



def magnitude(a,b):
    v = a-b
    v=v*v
    sum = np.sum(v)
    return math.sqrt(sum)

def solve(a,b):

    N = len(a)
    sum = 0

    for j in range(N):
        mag = magnitude(a[j].position,b[j].position)
        c_theta = (a[j].position[2]-b[j].position[2])/mag
        arccos = np.arccos(c_theta)
        degrees = np.degrees(arccos)
        x_norm = degrees - 90
        sum += x_norm
    sum /= N
    return sum



def helix_tilt(a,b,u):
    
    mylist =[]

    for ts in u.trajectory:
        # calling function to solve order param in each timestep
        t = u.trajectory.time/1000
        y = solve(a,b)
        mylist.append([y])
    return np.array(mylist)