
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





def dist_resids(traj:traj):
    centers1 = a.transform(traj)  # yields ndarray
    centers2 = b.transform(traj)  # yields ndarray
    cen1 = centers1[:,-1]
    cen2 = centers2[:,-1]

    xyz = np.hstack((centers1, centers2))
    traj = md.Trajectory(xyz.reshape(-1, 2, 3), topology=None)
    
    cont = md.compute_distances(traj,atom_pairs=[[0, 1]], periodic=False)
    print(np.shape(cont))
    return cont
        
### compute distance between 11 and 27
from pyemma.coordinates.data.featurization.misc import GroupCOMFeature
c_1a = GroupCOMFeature(feat.topology, [list(range(157,180))])
c_2a = GroupCOMFeature(feat.topology, [list(range(437,443))])
# c3 = GroupCOMFeature(feat.topology, [list(range(90,468))])

def feature_function_trp_gly(traj: traj):
    centers1 = c_1a.transform(traj)  # yields ndarray
    centers2 = c_2a.transform(traj)  # yields ndarray
#     print(centers1[:,-1])
    cen1 = centers1[:,-1]
    cen2 = centers2[:,-1]
#     centers3 = c3.transform(traj)  # yields ndarray

    xyz = np.hstack((centers1, centers2))
    traj = md.Trajectory(xyz.reshape(-1, 2, 3), topology=None)
    # this has shape (n_frames, 1)
    print(traj)
    
    
    cont = md.compute_distances(traj,atom_pairs=[[0, 1]], periodic=False)
    print(np.shape(cont))
    return cont

########### COM, accounting only for z-dim

from pyemma.coordinates.data.featurization.misc import GroupCOMFeature
c1 = GroupCOMFeature(feat.topology, [list(range(325,367))])
c2 = GroupCOMFeature(feat.topology, [list(range(628,33617))])
# c3 = GroupCOMFeature(feat.topology, [list(range(90,468))])

def feature_function_com(traj: traj):
    centers1 = c1.transform(traj)  # yields ndarray
    centers2 = c2.transform(traj)  # yields ndarray

#     print(centers1[:,-1])
    cen1 = centers1[:,-1]
    cen2 = centers2[:,-1]
    sub=(cen2-cen1)
    reshaped=(sub.reshape(-1,1))
    print(np.shape(reshaped))
    j=np.absolute(reshaped)
    k=j.astype('float32')
#     dt=j.dtype(np.float32)
    return k

### z-com of tyrosine on end (37)

c1a = GroupCOMFeature(feat.topology, [list(range(604,627))])
c2a = GroupCOMFeature(feat.topology, [list(range(628,33617))])
# c3 = GroupCOMFeature(feat.topology, [list(range(90,468))])

def feature_function_com2(traj: traj):
    centers1 = c1a.transform(traj)  # yields ndarray
    centers2 = c2a.transform(traj)  # yields ndarray

#     print(centers1[:,-1])
    cen1 = centers1[:,-1]
    cen2 = centers2[:,-1]
    sub=(cen2-cen1)
    reshaped=(sub.reshape(-1,1))
    print(np.shape(reshaped))
    j=np.absolute(reshaped)
    k=j.astype('float32')
#     dt=j.dtype(np.float32)
    return k




atoms=[]
for number in range(355,367):
    atoms.append(number)
    
def contacts_arg(traj: traj):

    cutoff=0.3
    j=md.compute_neighbors(traj,cutoff=cutoff,query_indices=atoms)
    ls=[]
    for arr in j:
        lst=arr.tolist()
        sort = [x for x in lst if x>628 and x<36617]
        ls.append(len(sort))
    
    arr=np.array(ls)
    reshaped=(arr.reshape(-1,1))
    print(np.shape(reshaped))
    j=np.absolute(reshaped)
    k=j.astype('float32')
    return k


#### Phe

atoms_phe=[]
for number in range(15,24):
    atoms_phe.append(number)

def contacts_phe(traj: traj):

    cutoff=0.35
    j=md.compute_neighbors(traj,cutoff=cutoff,query_indices=atoms_phe)
    ls=[]
    for arr in j:
        lst=arr.tolist()
        sort = [x for x in lst if x>628 and x<36617]
        ls.append(len(sort))
    
    arr=np.array(ls)
    reshaped=(arr.reshape(-1,1))
    j=np.absolute(reshaped)
    k=j.astype('float32')
    return k

######## contacts leu21

atoms_leu21=[]
for number in range(332,340):
    atoms_leu21.append(number)

def contacts_leu21(traj: traj):

    cutoff=0.35
    j=md.compute_neighbors(traj,cutoff=cutoff,query_indices=atoms_leu21)
    ls=[]
    for arr in j:
        lst=arr.tolist()
        sort = [x for x in lst if x>628 and x<36617]
        ls.append(len(sort))
    
    arr=np.array(ls)
    reshaped=(arr.reshape(-1,1))
    j=np.absolute(reshaped)
    k=j.astype('float32')
    return k

###### leu 36

atoms_leu36=[]
for number in range(594,600):
    atoms_leu36.append(number)


def contacts_leu36(traj: traj):

    cutoff=0.35
    j=md.compute_neighbors(traj,cutoff=cutoff,query_indices=atoms_leu36)
    ls=[]
    for arr in j:
        lst=arr.tolist()
        sort = [x for x in lst if x>628 and x<36617]
        ls.append(len(sort))
    
    arr=np.array(ls)
    reshaped=(arr.reshape(-1,1))
    j=np.absolute(reshaped)
    k=j.astype('float32')
    return k


###### lys33

atoms_lys33=[]
for number in range(548,551):
    atoms_lys33.append(number)


def contacts_lys33(traj: traj):

    cutoff=0.35
    j=md.compute_neighbors(traj,cutoff=cutoff,query_indices=atoms_lys33)
    ls=[]
    for arr in j:
        lst=arr.tolist()
        sort = [x for x in lst if x>628 and x<36617]
        ls.append(len(sort))
    
    arr=np.array(ls)
    reshaped=(arr.reshape(-1,1))
    j=np.absolute(reshaped)
    k=j.astype('float32')
    return k

### compute distance between 11 and 31
from pyemma.coordinates.data.featurization.misc import GroupCOMFeature
c_1 = GroupCOMFeature(feat.topology, [list(range(157,180))])
c_2 = GroupCOMFeature(feat.topology, [list(range(496,512))])


def feature_function_trp_gln(traj: traj):
    centers1 = c_1.transform(traj)  # yields ndarray
    centers2 = c_2.transform(traj)  # yields ndarray

    cen1 = centers1[:,-1]
    cen2 = centers2[:,-1]

    xyz = np.hstack((centers1, centers2))
    traj = md.Trajectory(xyz.reshape(-1, 2, 3), topology=None)


    cont = md.compute_distances(traj,atom_pairs=[[0, 1]], periodic=False)
    return cont
    
### compute distance between 8 and 31
from pyemma.coordinates.data.featurization.misc import GroupCOMFeature
c_1a = GroupCOMFeature(feat.topology, [list(range(112,118))])
c_2a = GroupCOMFeature(feat.topology, [list(range(496,512))])