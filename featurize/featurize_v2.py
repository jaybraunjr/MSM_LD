import matplotlib.pyplot as plt
import numpy as np
import pyemma
from pyemma.util.contexts import settings
import mdtraj as md
from itertools import combinations
from pyemma.coordinates import source
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.core.groups import AtomGroup
import math
from MDAnalysis.analysis import contacts
import MDAnalysis as mda


class dist:
    """
    This class provides methods for calculating distances and penetration depths
    in molecular dynamics simulations.
    """

    def min_dist(self, mda_list):
        # Dictionary to store results for each trajectory
        all_min_distances = {}
        for u in mda_list:
            # Select protein residues
            protein_residues = u.select_atoms('protein').residues
            lipids = u.select_atoms('resname POPC DOPE SAPI')
            min_distances = np.zeros((len(protein_residues), len(u.trajectory)))
            for ts in u.trajectory:
                # For each protein residue
                for i, residue in enumerate(protein_residues):
                    distances = np.min(mda.lib.distances.distance_array(residue.atoms.positions, lipids.positions, box=ts.dimensions))
                    min_distances[i, ts.frame] = distances
            all_min_distances[u.trajectory.filename] = min_distances/10

        trajectories_list = [value.T for value in all_min_distances.values()]

        return trajectories_list

    def dist_bil(self, traj, c1, c2):
        """
        Calculate the absolute distance between the last atoms of two components
        in a trajectory.
        
        Parameters:
        - traj: MD trajectory object.
        - c1, c2: Components for distance calculation.
        
        Returns:
        - j1: The absolute distances between the last atoms of components c1 and c2.
        """
        centers1 = c1.transform(traj)
        bulk = c2.transform(traj)
        cen1 = centers1[:, -1]
        cen_b = bulk[:, -1]
        sub1 = (cen1 - cen_b)
        reshaped1 = (sub1.reshape(-1, 1))

        j1 = np.absolute(reshaped1)
        return j1

    def dist_resids(self, traj, c1, c2):
        """
        Calculate distances between residues of two components in a trajectory.
        
        Parameters:
        - traj: MD trajectory object.
        - c1, c2: Components for distance calculation.
        
        Returns:
        - cont: The distances between residues of c1 and c2.
        """
        centers1 = c1.transform(traj)
        centers2 = c2.transform(traj)
        cen1 = centers1[:, -1]
        cen2 = centers2[:, -1]

        xyz = np.hstack((centers1, centers2))
        traj = md.Trajectory(xyz.reshape(-1, 2, 3), topology=None)

        cont = md.compute_distances(traj, atom_pairs=[[0, 1]], periodic=False)
        print(np.shape(cont))
        return cont

    def compute_penetration_depth(self, traj, protein_ca_indices, lipid_p_indices):

        """
        Compute the penetration depth of protein CA atoms into a lipid bilayer.
        
        Parameters:
        - traj: MD trajectory object.
        - protein_ca_indices: Indices of CA atoms in the protein.
        - lipid_p_indices: Indices of P atoms in the lipid.
        
        Returns:
        - penetration_depth: The penetration depth of protein CA atoms into the lipid bilayer.
        """
        lipid_p_positions = traj.xyz[:, lipid_p_indices]

        lipid_p_indices_top = lipid_p_indices[lipid_p_positions[..., 2].mean(axis=0) >= traj.unitcell_lengths[:, 2].mean() / 2]
        lipid_p_indices_bottom = lipid_p_indices[lipid_p_positions[..., 2].mean(axis=0) < traj.unitcell_lengths[:, 2].mean() / 2]

        lipid_p_com_top = md.compute_center_of_mass(traj.atom_slice(lipid_p_indices_top))
        lipid_p_com_bottom = md.compute_center_of_mass(traj.atom_slice(lipid_p_indices_bottom))

        protein_ca_positions = traj.xyz[:, protein_ca_indices]

        penetration_depth = np.empty((traj.n_frames, len(protein_ca_indices)))

        for frame in range(traj.n_frames):
            dist_top = np.linalg.norm(protein_ca_positions[frame] - lipid_p_com_top[frame], axis=-1)
            dist_bottom = np.linalg.norm(protein_ca_positions[frame] - lipid_p_com_bottom[frame], axis=-1)
            binding_to_top = dist_top < dist_bottom

            penetration_depth[frame] = np.where(binding_to_top, protein_ca_positions[frame, :, 2] - lipid_p_com_top[frame, 2],
                                                protein_ca_positions[frame, :, 2] - lipid_p_com_bottom[frame, 2])

            box_half_z = traj.unitcell_lengths[frame, 2] / 2

            for i, depth in enumerate(penetration_depth[frame]):
                if depth > 0 and protein_ca_positions[frame, i, 2] < lipid_p_com_bottom[frame, 2] and not binding_to_top[i]:
                    penetration_depth[frame, i] = -depth
                elif depth < 0 and protein_ca_positions[frame, i, 2] < box_half_z:
                    penetration_depth[frame, i] = np.abs(depth)

                if depth > 0 and protein_ca_positions[frame, i, 2] < box_half_z and not binding_to_top[i]:
                    penetration_depth[frame, i] = -depth

        return penetration_depth


class cont:

    """
    This class provides methods for analyzing contacts within molecular dynamics simulations.
    """
    
    def contacts_mdtraj(self, traj, cutoff, query_indices):

        """
        Calculate the number of contacts within a cutoff distance using mdtraj.

        """
        j = md.compute_neighbors(traj, cutoff=cutoff, query_indices=query_indices)
        ls = []
        for arr in j:
            lst = arr.tolist()
            sort = [x for x in lst if x > 628 and x < 36617]
            ls.append(len(sort))
        
        arr = np.array(ls)
        reshaped = (arr.reshape(-1, 1))
        j = np.absolute(reshaped)
        k = j.astype('float32')
        return k

    def contacts_MDA(self, u, a, b, radius=3.5):

        """
        Use MDAnalysis to calculate the number of contacts between two atom groups.

        """
        timeseries = []

        for ts in u.trajectory:
            dist = contacts.distance_array(a.positions, b.positions)
            n_contacts = contacts.contact_matrix(dist, radius).sum()
            timeseries.append(n_contacts)
            arr = np.array(timeseries)
        return np.transpose([arr])

    def calculate_contacts(self, mda_list,protein_selection, lipid_selection):

        """Calculate number of contacts for each residue in a trajectory."""
        # Select protein residues
        ls=[]
        for u in mda_list:
            protein_residues = u.select_atoms(protein_selection).residues
            lipids = u.select_atoms(lipid_selection)
            # Cutoff distance for contact (in Angstroms)
            cutoff = 3.5

            contacts = np.zeros((len(u.trajectory), len(protein_residues)))

            # Loop over each frame
            for ts in u.trajectory:
                # For each protein residue
                for i, residue in enumerate(protein_residues):
                    # Calculate distances between the residue and all lipids
                    distances = mda.lib.distances.distance_array(residue.atoms.positions, lipids.positions, box=ts.dimensions)
                    # Count the number of lipids that have any atom within the cutoff distance to the residue
                    num_contacts = np.sum(np.any(distances < cutoff, axis=1))
                    contacts[ts.frame, i] = num_contacts
            ls.append(contacts)
        return ls

    def run_coord(self, u_, ag1, ag2, nn=6, mm=12, d0=0, r0=2.5, density=False, b=0, e=None, skip=1):

        """
        Calculate coordination numbers for atom groups in an MDAnalysis Universe.
        
        Parameters:
        - u_: MDAnalysis Universe object.
        - ag1, ag2: Atom groups for which coordination numbers are calculated.
        - nn, mm, d0, r0: Parameters for the coordination number calculation.
        - density: If True, normalize the coordination numbers by the product of the numbers of atoms in ag1 and ag2.
        - b, e, skip: Start, end, and skip values for trajectory frames.
        
        Returns:
        - Coordination numbers for ag1 and ag2.
        """
        times = []
        coords = []
        for i, ts in enumerate(u_.trajectory[b:e:skip]):
            times.append(i)
            d = distance_array(ag1.positions, ag2.positions, box=u_.dimensions)
            d[d < d0] = 1
            D = (d - d0) / r0
            sij = (1 - np.power(D, nn)) / (1 - np.power(D, mm))
            coords.append(np.sum(sij))
        if density:
            coords = np.array(coords)
            coords /= (ag1.n_atoms * ag2.n_atoms)
        return np.transpose([coords])



class HelixAnalysis:
    """
    A class dedicated to analyzing helical structures within molecular dynamics simulations.
    """

    def __init__(self, top_file):
        """
        Initializes the HelixAnalysis with a given MDAnalysis Universe object.

        Parameters:
        - universe: The MDAnalysis Universe object containing the trajectory.
        """
        self.universe = None
        self.top_file = top_file

    def calculate_tilt_angle(self, start_residue, end_residue):
        """
        Compute the tilt angle of a helix with respect to the x-y plane.

        Parameters:
        - start_residue: The starting residue number of the helix.
        - end_residue: The ending residue number of the helix.

        Returns:
        - angle_with_xy: The tilt angle of the helix with respect to the x-y plane, adjusted for orientation.
        """
        start_atom = self.universe.select_atoms(f"resid {start_residue} and name CA").center_of_mass()
        end_atom = self.universe.select_atoms(f"resid {end_residue} and name CA").center_of_mass()

        helix_vector = end_atom - start_atom
        helix_vector /= np.linalg.norm(helix_vector)
        dot_product = np.dot(helix_vector, [0, 0, 1])
        angle_with_z = np.arccos(dot_product) * 180.0 / np.pi
        angle_with_xy = 90 - angle_with_z

        center_atom_z = self.universe.select_atoms(f"resid {start_residue} and name CA").center_of_mass()[2]
        box_center_z = self.universe.dimensions[2] / 2.0
        if center_atom_z < box_center_z:
            angle_with_xy = -angle_with_xy

        return angle_with_xy


    def analyze_multiple_files(self, traj_files, helix_ranges):
        """
        Analyzes multiple trajectories for specified helix ranges.

        Parameters:
        - traj_files: List of trajectory file paths.
        - helix_ranges: List of tuples specifying start and end residues for each helix to analyze.
        """
        trajectories_list = []

        for traj_file in traj_files:
            u = mda.Universe(self.top_file, traj_file)
            self.universe = u
            angles = []
            for ts in u.trajectory:
                angles_ = [self.calculate_tilt_angle(start_residue, end_residue) for start_residue, end_residue in helix_ranges]
                angles.append(tuple(angles_))
            trajectories_list.append(angles)
        return trajectories_list

