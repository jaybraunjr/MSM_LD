
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

#### this is inefficient, should find a way to make this simple

# traj = md.load(../)



def f1(traj,c1,c2):
    centers1 = c1.transform(traj)
    bulk = c2.transform(traj)
    cen1 = centers1[:,-1]
    cen_b = bulk[:,-1]
    sub1=(cen1-cen_b)
    reshaped1=(sub1.reshape(-1,1))

    j1=np.absolute(reshaped1)
    k1=j1.astype('float32')
    return 1/k1

# def f2(traj,c1,c2):
#     centers1 = c1.transform(traj)
#     bulk = c2.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# # def f3(traj:traj):
#     centers1 = c3.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f4(traj:traj):
#     centers1 = c4.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f5(traj:traj):
#     centers1 = c5.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f6(traj:traj):
#     centers1 = c6.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f7(traj:traj):
#     centers1 = c7.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f8(traj:traj):
#     centers1 = c8.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f9(traj:traj):
#     centers1 = c9.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f10(traj:traj):
#     centers1 = c10.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f11(traj:traj):
#     centers1 = c11.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f12(traj:traj):
#     centers1 = c12.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f13(traj:traj):
#     centers1 = c13.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f14(traj:traj):
#     centers1 = c14.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f15(traj:traj):
#     centers1 = c15.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f16(traj:traj):
#     centers1 = c16.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f17(traj:traj):
#     centers1 = c17.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f18(traj:traj):
#     centers1 = c18.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f19(traj:traj):
#     centers1 = c19.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f20(traj:traj):
#     centers1 = c20.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f21(traj:traj):
#     centers1 = c21.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f22(traj:traj):
#     centers1 = c22.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f23(traj:traj):
#     centers1 = c23.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f24(traj:traj):
#     centers1 = c24.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f25(traj:traj):
#     centers1 = c25.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f26(traj:traj):
#     centers1 = c26.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f27(traj:traj):
#     centers1 = c27.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f28(traj:traj):
#     centers1 = c28.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f29(traj:traj):
#     centers1 = c29.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f30(traj:traj):
#     centers1 = c30.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f31(traj:traj):
#     centers1 = c31.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f32(traj:traj):
#     centers1 = c32.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f33(traj:traj):
#     centers1 = c33.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f34(traj:traj):
#     centers1 = c34.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f35(traj:traj):
#     centers1 = c35.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1

# def f36(traj:traj):
#     centers1 = c36.transform(traj)
#     bulk = b.transform(traj)
#     cen1 = centers1[:,-1]
#     cen_b = bulk[:,-1]
#     sub1=(cen1-cen_b)
#     reshaped1=(sub1.reshape(-1,1))

#     j1=np.absolute(reshaped1)
#     k1=j1.astype('float32')
#     return k1