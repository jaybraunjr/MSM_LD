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
import os

def feature(input):
    feat = pyemma.coordinates.featurizer(input)
    return feat

def make_arr(tic):
    ls = tic.tolist()
    ls1 = np.array(ls)
    return(ls1)

def save_reader(inputs, name_prefix):
    for i in range(len(inputs)):
        with open(f'{name_prefix}_{i:03d}.npy','wb') as handle:
            np.save(handle, inputs[i])

def run_func(ls,feat_list,feature_function,in_sel2):
    ls1=[]
    for traj in ls:
        print(traj)
        for feat in feat_list:
            print(feat)
            ls1.append(feature_function(traj,traj.select_atoms(feat),traj.select_atoms(in_sel2)))
    return ls1

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


def make(data_inp):
    splits = np.array_split(data_inp, len(files))
    processor_ = caller.ChunkProcessor(splits, int(len(data_inp)/length) )
    m1 = processor_.fubar(splits, int(len(data_inp)/(len(files))))
    return(m1)
    
def combine(data_inp):
    # Create a list of ReturnInputs instances, one for each trajectory
    processors = [caller.ReturnInputs() for _ in range(len(files))]
    # Reshape the processed results using the ReturnInputs instances
    inputs = [processor.return_inputs(make(data_inp), i, int(len(data_inp)/length)) for i, processor in enumerate(processors)]
    inputs = [x[0] for x in inputs]
    return(inputs)




