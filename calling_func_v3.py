import feat_v3 as f_v2
import featurize.caller as caller
import MDAnalysis as mda
import numpy as np
import pyemma
import mdtraj as md
from featurize.featurize_v2 import cont, dist

class TrajectoryAnalysis:
    def __init__(self, top_file, traj_files):
        self.top_file = top_file
        self.traj_files = traj_files

    def load_md_and_mda_list(self):
        return [(md.load(item, top=self.top_file), mda.Universe(self.top_file, item)) for item in self.traj_files]

    def calculate_contacts(self, lipid_selection, prot_selections, universes):
        contact_analysis = cont()
        return [contact_analysis.contacts_MDA(universe, universe.select_atoms(prot_sel), universe.select_atoms(lipid_selection))
                for universe in universes for prot_sel in prot_selections]

    def calculate_penetration_depths(self, md_list, protein_ca_indices, lipid_p_indices):
        dist_instance = dist()
        penetration_depths = []
        for traj in md_list:
            penetration_depth = dist_instance.compute_penetration_depth(traj, protein_ca_indices, lipid_p_indices)
            penetration_depths.append(penetration_depth)
        print(np.shape(penetration_depths))
        print(np.shape(penetration_depths[0]))
        return penetration_depths

    def combine_results(self, data_output, inputs1, inputs2):
        return [np.concatenate((arr1, arr2, arr3), axis=1) for arr1, arr2, arr3 in zip(data_output, inputs1, inputs2)]

    # def calculate_dists(self, feat_md, md_list):
    #     a = [list(range(1, 367))]
    #     bil = [list(range(628, 33617))]
    #     b = [list(range(1, 2))]
    #     c = [list(range(626, 627))]
    #     d = [list(range(128, 129))]
    #     e = [list(range(228, 229))]
    #     f = [list(range(328, 339))]
    #     g = [list(range(428, 429))]
    #     h = [list(range(528, 529))]
    #     i = [list(range(600, 601))]
    #     j = [list(range(28, 29))]
    #     k = [list(range(439, 440))]
    #     l = [list(range(462, 463))]
    #     m = [list(range(498, 499))]
    #     n = [list(range(534, 535))]

    #     return f_v2.calculate_dists(feat_md, bil, a, b, c, d, e, f, g, h, i, j, k, l, m, n, md_list)

    def run_analysis(self):
        # Load both MD and MDA lists and unzip them into separate lists
        md_list, mda_list = zip(*self.load_md_and_mda_list())
        
        # Generate protein selections based on residue IDs
        prot_selections = [f"(resid {i}) and (not backbone)" for i in range(1, 36)]
        
        # Define the lipid selection string
        lipid_selection = 'resname POPC DOPE SAPI'
        
        # Calculate contacts between lipids and proteins
        contacts = self.calculate_contacts(lipid_selection, prot_selections, mda_list)

        # Create a featurizer object for the topology file
        feat_md = pyemma.coordinates.featurizer(self.top_file)

            # Calculate penetration depths for each trajectory in md_list
        traj = md.load('traj_dat/w1.xtc', top='traj_dat/10.gro')
        protein_ca_indices = traj.topology.select('type C and protein')
        lipid_p_indices = traj.topology.select('name P and resname POPC')
        inputs2 = self.calculate_penetration_depths(md_list, protein_ca_indices, lipid_p_indices)

        # dist_instance = dist()
        # penetration_depths = []
        # for traj in md_list:
        #     print(traj)
        #     penetration_depth = dist_instance.compute_penetration_depth(traj, protein_ca_indices, lipid_p_indices)
        #     print(penetration_depth)
        #     penetration_depths.append(penetration_depth)
        
        # # Calculate distances between groups of atoms
        # dists = self.calculate_dists(feat_md, md_list)
        
        # Split the contacts and distances lists into equal parts
        # This changes a lot depending on how many features we want to use
        # splits1, splits2 = np.array_split(contacts, len(md_list)), np.array_split(dists, len(md_list))
        splits1 = np.array_split(contacts, len(md_list))

        # Get FUBAR outputs for contacts and distances
        # This changes a lot as well
        # m1, m2 = f_v2.get_fubar_output(splits1, 35), f_v2.get_fubar_output(splits2, 14)
        m1 = f_v2.get_fubar_output(splits1, 35)
        
        # Get inputs for contacts and distances
        # This changes alot depending on the number of features you want to use, and which features we use
        # inputs1, inputs2 = f_v2.get_inputs(m1, splits1, 35), f_v2.get_inputs(m2, splits2, 14)
        inputs1 = f_v2.get_inputs(m1, splits1, 35)

        # Read the PyEMMA data output
        data_output = self.pyemma_reader(self.top_file, self.traj_files)
        
        # Combine the results from PyEMMA data output, contacts, and distances
        # We change this up a lot depending on what we want to use
        result_list = self.combine_results(data_output, inputs1, inputs2)
        
        # Get the number of frames in each MD trajectory
        n_frames = [len(item) for item in md_list]
        
        # Create a list of tuples with the number of frames and the second dimension of the results
        t_list = [(frames, result_list[0].shape[1]) for frames in n_frames]
        
        # Create an array list with the appropriate dimensions
        arrlist = f_v2.make_arrlist(result_list, dims=t_list)
        
        # Print the shapes of the array list and the first element of the array list
        print(np.shape(arrlist))
        print(np.shape(arrlist[0]))
        
        # Save the array list to a file
        caller.save_reader(arrlist, 'arrlist')



    def pyemma_reader(self, top, files):
        feat = pyemma.coordinates.featurizer(top)
        feat.add_backbone_torsions()
        reader = pyemma.coordinates.source(files, features=feat)
        return reader.get_output()

def main():
    top_file = 'traj_dat/10.gro'
#     files = ['traj_dat/whole.1.xtc', 'traj_dat/whole.3.xtc', 'traj_dat/whole.2.xtc']
    files =['../trajdat/dat2/rep1_tot.xtc','../trajdat/dat2/rep2_tot.xtc','../trajdat/dat2/rep3_tot.xtc',\
       '../trajdat/dat2/rep4.xtc','../trajdat/dat2/rep6_tot2.xtc','../trajdat/dat/1_tot.xtc','../trajdat/dat/3_tot.xtc',\
        '../trajdat/dat/5_tot.xtc','../trajdat/dat/9_tot.xtc']
    traj_analysis = TrajectoryAnalysis(top_file, files)
    traj_analysis.run_analysis()

if __name__ == '__main__':
    main()
