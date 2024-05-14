import itertools
import numpy as np
import pandas as pd
import MDAnalysis as mda
import matplotlib as mpl
print('testing')
import matplotlib.pyplot as plt
import mdtraj as md
import pyemma
from pyemma.coordinates import source
from pyemma.util.contexts import settings
from pyemma.coordinates.data.featurization.misc import GroupCOMFeature

import sys
sys.path.insert(0, '../')
import featurize.featurize_v2 as feat_
import featurize.caller_v2 as caller
import os

path = '../../adaptive_sampling_trajs/send/'
files = [os.path.join(path, f) for f in os.listdir(path) if f.endswith('.xtc')]

for file in files:
    print(file)
mda_list = []
for item in files:
    u_ = mda.Universe('10.gro', item)
    mda_list.append(u_)
    print(item)

# Define lipid and protein selections
lipid_selection = '(resname POPC DOPE SAPI)'
prot_selections = [f"(resid {i}) and (not backbone)" for i in range(1, 36)]

# Initialize contact analysis class
contact_analysis = feat_.cont()

# Run analysis for each trajectory and each protein selection
mda_conts = []
for universe in mda_list:
    for prot_sel in prot_selections:
        print(prot_sel)
        result = contact_analysis.contacts_MDA(universe, universe.select_atoms(prot_sel), universe.select_atoms(lipid_selection))
        mda_conts.append(result)


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

# Combine the data from each dataset for each trajectory
result_list = []
for i in range(len(files)):
    arr = combine(results_coord)[i]
    result_list2.append(arr)

caller.save_reader(result_list, 'contacts_v1')
