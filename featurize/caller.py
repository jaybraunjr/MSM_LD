### this will call on features

import MDAnalysis as mda
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pyemma
from pyemma.util.contexts import settings
import mdtraj as md
import itertools
from itertools import combinations
import pandas as pd
from pyemma.coordinates import source


import numpy as np
import itertools


def feature(input):
	feat = pyemma.coordinates.featurizer(input)
	return feat


#### here we will use MDAnalysis to create features as arrays. This goes against the standard of pyemma (mdtraj), but same results. 
### we just have to modify a few things

# def concat(a,b,c,d,e,f,g,h,i,j):
# 	cont = [np.concatenate((a,b,c,d,e,f,g,h,i,j),axis=1)]
# 	return cont

def make_arr(tic):
	ls = tic.tolist()
	ls1 = np.array(ls)
	return(ls1)


# def save_reader(inputs,name):
# 	for i in range(len(inputs)):
# 	    with open(f'reader_{i:03d}.npy','wb') as handle:           
# 	        np.save(handle, inputs[i])
def save_reader(inputs, name_prefix):
    
    for i in range(len(inputs)):
        with open(f'{name_prefix}_{i:03d}.npy','wb') as handle:
            np.save(handle, inputs[i])
"""
Now have to make functions that include reshaping the arrays and what not to be ready for tICA.

- Takes in data from featurizing
- Reshapes and returns
- Saves for future use/tICA

"""


## the below function will always change depending on what the action is

ls1=[]
def run_func(ls,feat_list,feature_function,in_sel2):
    for traj in ls:
        print(traj)
        for feat in feat_list:
            print(feat)
            ls1.append(feature_function(traj,traj.select_atoms(feat),traj.select_atoms(in_sel2)))
    return ls1



""" 
The below functions are used to reshape
"""


class ChunkProcessor:
    def __init__(self, lst, chunk_size):
        self.lst = lst
        self.chunk_size = chunk_size

    def split_list_into_chunks(self, lst, chunk_size):
        return [lst[i:i+chunk_size] for i in range(0, len(lst), chunk_size)]
    
    def fubar(self, a, chunk_size):
        return [self.split_list_into_chunks(sub, chunk_size) for sub in a]
    


class ReturnInputs:
    def iteration(self, array):
        arrr = np.array(array)
        unested1 = [list(itertools.chain(*sub)) for sub in arrr]
        arr_ = [np.array(unested1)]

        return(arr_)

    def get_array(self, array, chunk):
        ls=[]
        for stuff in array[chunk]:
            ls.append(self.iteration(stuff))
        return ls
    
    def return_inputs(self, array, chunk, feat_len):
        r1 = self.get_array(array,chunk)
        shape = np.shape(r1)
        print(shape)
        arr = np.array(r1)
        reshape = arr.reshape(feat_len,shape[3])
        stuff2 = [np.transpose(reshape)]
        return(stuff2)



def get_inputs(m, splits, feat_len):
    processors = [ReturnInputs() for i in range(len(splits))]
    inputs = [processor.return_inputs(m, i, feat_len)[0] for i, processor in enumerate(processors)]
    return inputs

def get_fubar_output(splits, feat_len):
    processor = ChunkProcessor(splits, feat_len)
    return processor.fubar(splits, feat_len)
