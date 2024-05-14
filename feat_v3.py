import numpy as np
import MDAnalysis as mda
import mdtraj as md
import pyemma
from pyemma.coordinates.data.featurization.misc import GroupCOMFeature

import featurize.featurize as feat
import featurize.caller as caller


def calculate_dists(top_file, bil, a, b, c, d, e, f, g, h, i, j, k, l, m, n, nlist):
    bil = GroupCOMFeature(top_file.topology,bil)
    a = GroupCOMFeature(top_file.topology, a)
    b = GroupCOMFeature(top_file.topology, b)
    c = GroupCOMFeature(top_file.topology, c)
    d = GroupCOMFeature(top_file.topology, d)
    e = GroupCOMFeature(top_file.topology, e)
    f = GroupCOMFeature(top_file.topology, f)
    g = GroupCOMFeature(top_file.topology, g)
    h = GroupCOMFeature(top_file.topology, h)
    i = GroupCOMFeature(top_file.topology, i)
    j = GroupCOMFeature(top_file.topology, j)
    k = GroupCOMFeature(top_file.topology, k)
    l = GroupCOMFeature(top_file.topology, l)
    m = GroupCOMFeature(top_file.topology, m)
    n = GroupCOMFeature(top_file.topology, n)
    group_features_dict = {'a': a, 'b': b, 'c': c, 'd': d, 'e': e, 'f': f, 'g': g, 'h': h, 'i': i, 'j': j, 'k': k, 'l': l, 'm': m, 'n': n}
    d_bil = feat.dist()
    ls2 = []
    for traj in nlist:
        for name, feature in group_features_dict.items():
            ls2.append(d_bil.dist_bil(traj, bil, feature))
    return ls2

def get_inputs(m, splits, feat_len):
    processors = [caller.ReturnInputs() for i in range(len(splits))]
    inputs = [processor.return_inputs(m, i, feat_len)[0] for i, processor in enumerate(processors)]
    return inputs

def get_fubar_output(splits, feat_len):
    processor = caller.ChunkProcessor(splits, feat_len)
    return processor.fubar(splits, feat_len)

def make_arrlist(input_cont, dims):
    arrs = [np.array(input_cont[i]).reshape(*dims[i]) for i in range(len(input_cont))]
    return arrs


