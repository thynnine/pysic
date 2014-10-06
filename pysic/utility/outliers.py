#! /usr/bin/env python
# coding=utf-8
"""
Detection of irregular atomic neighborhoods by bond distance and angle analysis.
"""

# Ville Parkkinen, Teemu Hynninen

import numpy as np
from math     import *
from scipy    import spatial
from warnings import catch_warnings, simplefilter
import ase
import pysic


class Angle(object):
    __slots__ = ['center_index', 'type1', 'type2', 'type3', 'value']

    def __init__(self, center_index, type1, type2, type3, angle):
        self.center_index                  = center_index
        self.type1, self.type2, self.type3 = np.sort([type1, type2, type3])
        self.value                         = angle



class Distance(object):
    __slots__ = ['primary_index', 'type1', 'type2', 'value']

    def __init__(self, primary_index, type1, type2, distance):
        self.primary_index     = primary_index
        self.type1, self.type2 = np.sort([type1, type2])
        self.value             = distance



class Distribution:
    def __init__(self):
        self.items = []


class BondType:
    def __init__(self, elements, cutoff):
        self.elements = elements
        self.cutoff = cutoff



def angle(A, O, B):
        """Return the angle between vectors OA and OB

        Parameters O, A, B: coordinates in 3d-space
        """

        OA  = A - O
        OB  = B - O

        div = sqrt(np.power(OA, 2).sum()) * sqrt(np.power(OB, 2).sum())
        if div == 0: return None
        else       : s = (OA * OB).sum() / div
        if s >  1  : return 0   # due to floating-point arithmetic
        if s < -1  : return pi  # due to floating-point arithmetic

        return acos(s) 
    
def vec_angle(A, B):
        """Return the angle between vectors A and B

        Parameters A, B: coordinates in 3d-space
        """

        OA  = A
        OB  = B

        div = sqrt(np.power(OA, 2).sum()) * sqrt(np.power(OB, 2).sum())
        if div == 0: return None
        else       : s = (OA * OB).sum() / div
        if s >  1  : return 0   # due to floating-point arithmetic
        if s < -1  : return pi  # due to floating-point arithmetic

        return acos(s) 

class Structure:
    def __init__(self, atoms):
        self.system = atoms
        self.nbl = None
        self.bond_list = []


    
    def add_bond(self, elements, cutoff):
        """Adds a bond between atoms of the given elements, up to the cutoff separation.
        """
        self.bond_list.append(BondType(elements, cutoff))
        
    def get_bond_length(self, elements):
        """Returns the stored bond length for the given pair of elements.
        """
        rev = elements[::-1]
        for b in self.bond_list:
            if b.elements == elements:
                return b.cutoff
            elif b.elements == rev:
                return b.cutoff
    
    def create_neighbor_lists(self):
        """Creates neighbor lists for the structure.
        """
        skin = pysic.calculator.FastNeighborList.neighbor_marginal;

        # create a dummy calculator
        dummy = pysic.Pysic()
        old_calc = self.system.get_calculator()
        
        # Note that the list finds all neighbors within the maximum interaction
        # radius of the particular atom.
        for bond in self.bond_list:
            pot = pysic.Potential('LJ', cutoff=bond.cutoff-skin, symbols=bond.elements)
            dummy.add_potential(pot)
        self.system.set_calculator(dummy)
        
        # It is important to manually initialize the core since no actual calculations are carried out.
        dummy.set_core()

        # get the list and access its contents
        self.nbl = dummy.get_neighbor_list()
        
        # restore the original calculator, if there was one
        self.system.set_calculator(old_calc)
        
    
    def get_neighbors(self, index):
        """Returns the indices of neighboring atoms for the given atom. (Sorted according to distance.)
        """
        neighbors, offsets = self.nbl.get_neighbors(index, self.system, True)
        return neighbors
    
    def get_separations(self, index):
        """Returns the separation vectors from a given atom to its neighbors. (Sorted according to distance.)
        """
        return self.nbl.get_neighbor_separations(index, self.system, True)
        
    def get_distances(self, index):
        """Returns the distances between a given atom and its neighbors. (Sorted.)
        """
        return self.nbl.get_neighbor_distances(index, self.system, True)
        
    def get_all_angles(self):
        """Returns a list of all 3-atom angles as Angle objects.
        """
        
        angles = []
        symbs = self.system.get_chemical_symbols()
        for atom in self.system:
            nbors = self.get_neighbors(atom.index)
            vecs = self.get_separations(atom.index)
            ds = self.get_distances(atom.index)
            centre = atom.symbol
                
            for n1 in range(len(nbors)):
                for n2 in range(len(nbors)):
                    if(n1 > n2):
                        
                        type1 = symbs[nbors[n1]]
                        type2 = symbs[nbors[n2]]
                        
                        if(ds[n1] < self.get_bond_length([centre,type1]) and
                            ds[n2] < self.get_bond_length([centre,type2])):
                        
                            new_angle = vec_angle(vecs[n1],vecs[n2])
                            angles.append(Angle(atom.index,
                                                atom.symbol,
                                                type1,
                                                type2,
                                                new_angle))
        return angles
                
                
    def get_all_distances(self):
        """Returns a list of all 2-atom distances as Distance objects.
        """
        
        dists = []
        symbs = self.system.get_chemical_symbols()
        for atom in self.system:
            
            nbors = self.get_neighbors(atom.index)
            ds = self.get_distances(atom.index)
            type1 = atom.symbol
            
            for n1 in range(len(nbors)):
                type2 = symbs[nbors[n1]]
                        
                if(ds[n1] < self.get_bond_length([type1,type2])):
                    dists.append(Distance(atom.index,
                                          type1,
                                          type2,
                                          ds[n1]))
                                      
        return dists
    
    
def get_distributions(angles, distances, radii):
    """Return observed log-distributions of angles and distances

    Distributions are obtained for all combinations of observed atom types.
    They are histogram-based, using the breakpoints specified by parameters
    'angle.grid' and 'dist.grid'. The distributions are represented by
    vectors giving the log-probability of each "bin". These are not normalized
    by the "bin" widths. This is not a problem as these would cancel out later.
    """

    angle_distribs = {}
    dist_distribs  = {}

    min_angle, max_angle = np.inf, 0
    for angle in angles:
        if angle.value < min_angle: min_angle = angle.value
        if angle.value > max_angle: max_angle = angle.value
    
    for angle in angles:
        label = angle.type1 + '-' + angle.type2 + '-' + angle.type3
        if not label in angle_distribs.keys():
            angle_distribs[label] = Distribution()
        angle_distribs[label].items.append(angle.value)

    for label in angle_distribs:
        a, b = np.histogram(angle_distribs[label].items,
                            bins=100, range=[0.99*min_angle, 1.01*max_angle])
        a = a.astype(np.float)
        angle_distribs[label].grid = b
        with catch_warnings(): # ignore divbyzero-whining
            simplefilter('ignore')
            angle_distribs[label].distribution = np.log(a / sum(a))
        del(angle_distribs[label].items)



    min_dist , max_dist  = np.inf, 0
    for dist in distances:
        if dist.value < min_dist: min_dist = dist.value
        if dist.value > max_dist: max_dist = dist.value

    for dist in distances:
        label = dist.type1 + '-' + dist.type2
        if not label in dist_distribs.keys():
            dist_distribs[label] = Distribution()
        dist_distribs[label].items.append(dist.value)

    for label in dist_distribs:
        a, b = np.histogram(dist_distribs[label].items,
                            bins=100, range=[0.99*min_dist, 1.01*max_dist])
        a = a.astype(np.float)
        dist_distribs[label].grid = b
        with catch_warnings(): # ignore divbyzero-whining
            simplefilter('ignore')
            dist_distribs[label].distribution = np.log(a / sum(a))
        del(dist_distribs[label].items)

    return angle_distribs, dist_distribs



def get_log_likelihoods(angles, distances,angle_distribs,dist_distribs,n_atoms):
    """Calculate log-likelihoods of each atom

    All observed angles and distances are considered independent.
    Parameters 'angles' and .distances' supplie the observations,
    'angle_distribs' and 'dist_distribs'' the histogram-based distributions
    and 'angle.grid' & 'dist.grid' the breakpoints of these histograms.
    Parameter 'n_atoms' gives the number of atoms in the original data.
    """

    angles    = np.array(angles)
    distances = np.array(distances)

    a_logls   = np.tile(-np.inf, n_atoms)
    d_logls   = np.tile(-np.inf, n_atoms)

    for i in range(n_atoms): 
        # All angles/distances originating from atom i
        angleset = angles[np.array(map(lambda x: x.center_index, angles)) == i]
        distset  = distances[np.array(
                               map(lambda x: x.primary_index, distances)) == i]
        if np.size(angleset) > 0:
            a_logls[i] = 0
            for angle in angleset:
                label = angle.type1 + '-' + angle.type2 + '-' + angle.type3
                slots = np.histogram(angle.value,
                                     bins=angle_distribs[label].grid)[0]
                with catch_warnings(): # ignore nan-whining
                    simplefilter('ignore')
                    a_logls[i] += np.nansum(angle_distribs[label].distribution *
                                            slots)
            a_logls[i] = a_logls[i] / np.size(angleset)

            d_logls[i] = 0
            for dist in distset:
                label = dist.type1 + '-' + dist.type2
                slots = np.histogram(dist.value,
                                     bins=dist_distribs[label].grid)[0]
                with catch_warnings(): # ignore nan-whining
                    simplefilter('ignore')
                    d_logls[i] += np.nansum(dist_distribs[label].distribution *
                                            slots)
            d_logls[i] = d_logls[i] / np.size(distset)

    return a_logls, d_logls



def write_to_file(filename, boxsize, atoms, coordinates, a_logls, d_logls):
    """Write original data + results into a file

    Creates/overwrites an xyz-like file with two additional columns. These
    additional columns will contain the log-likelihood contributions from
    (1) angles and (2) distances.
    """

    f = open(filename, 'w')

    f.write('     {}\n'.format(np.shape(atoms)[0]))
    f.write('{} {} {}\n'.format(boxsize[0], boxsize[1], boxsize[2]))

    n = 10
    for i in range(np.shape(atoms)[0]):
        f.write(' {0:6}{1:14}{2:14}{3:14}{4:14}{5:10}\n'.format(
                  atoms[i],
                  str(coordinates[i, 0])[0:n],
                  str(coordinates[i, 1])[0:n],
                  str(coordinates[i, 2])[0:n],
                  str(a_logls[i])[0:n],
                  str(d_logls[i])[0:n]))
    f.close()
            

