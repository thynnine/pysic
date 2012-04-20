"""Contains necessary and supporting auxiliary functions for Pysic."""

try:
    import matplotlib.pyplot as plt
except:
    print "error importing matplotlib - the plotting tools are not available"
    
from math import sqrt, pi, sin, cos, exp, fabs, floor

from ase.calculators.neighborlist import NeighborList
import numpy as np


neighbor_marginal = 0.5

codec_s2i = {'1' : -1,
             '2' : -2,
             '3' : -3,
             '4' : -4,
             '5' : -5,
             '6' : -6,
             '7' : -7,
             '8' : -8,
             '9' : -9,
             '0' : -10,
             'a' : 1,
         'b' : 2, 
         'c' : 3, 
         'd' : 4, 
         'e' : 5, 
         'f' : 6, 
         'g' : 7, 
         'h' : 8, 
         'i' : 9, 
         'j' : 10,
         'k' : 11,
         'l' : 12,
         'm' : 13,
         'n' : 14,
         'o' : 15,
         'p' : 16,
         'q' : 17,
         'r' : 18,
         's' : 19,
         't' : 20,
         'u' : 21,
         'v' : 22,
         'w' : 23,
         'x' : 24,
         'y' : 25,
         'z' : 26,
         'A' : 101,
         'B' : 102,
         'C' : 103,
         'D' : 104,
         'E' : 105,
         'F' : 106,
         'G' : 107,
         'H' : 108,
         'I' : 109,
         'J' : 110,
         'K' : 111,
         'L' : 112,
         'M' : 113,
         'N' : 114,
         'O' : 115,
         'P' : 116,
         'Q' : 117,
         'R' : 118,
         'S' : 119,
         'T' : 120,
         'U' : 121,
         'V' : 122,
         'W' : 123,
         'X' : 124,
         'Y' : 125,
         'Z' : 126,
         '_' : 201 }

codec_i2s = {-1 : '1',
             -2 : '2',
             -3 : '3',
             -4 : '4',
             -5 : '5',
             -6 : '6',
             -7 : '7',
             -8 : '8',
             -9 : '9',
             -10 : '0',
             1 : 'a',
         2 : 'b',
         3 : 'c',
         4 : 'd',
         5 : 'e',
         6 : 'f',
         7 : 'g',
         8 : 'h',
         9 : 'i',
         10 : 'j',
         11 : 'k',
         12 : 'l',
         13 : 'm',
         14 : 'n',
         15 : 'o',
         16 : 'p',
         17 : 'q',
         18 : 'r',
         19 : 's',
         20 : 't',
         21 : 'u',
         22 : 'v',
         23 : 'w',
         24 : 'x',
         25 : 'y',
         26 : 'z',
         101 : 'A',
         102 : 'B',
         103 : 'C',
         104 : 'D',
         105 : 'E',
         106 : 'F',
         107 : 'G',
         108 : 'H',
         109 : 'I',
         110 : 'J',
         111 : 'K',
         112 : 'L',
         113 : 'M',
         114 : 'N',
         115 : 'O',
         116 : 'P',
         117 : 'Q',
         118 : 'R',
         119 : 'S',
         120 : 'T',
         121 : 'U',
         122 : 'V',
         123 : 'W',
         124 : 'X',
         125 : 'Y',
         126 : 'Z',
         201 : '_' }



class Cell:
    """Cell describing the simulation volume of a subvolume.
        
        The Cell object is extensively used by the :class:`~pysic_utility.FastNeighborList`
        for dividing the simulation volume and locating atom neighbors. 
        It can also be used by the user for coordinate manipulation. Note however,
        that ASE does not use on this class, as it is part of Pysic. The class is
        merely a tool for examining the geometry.
        
        Parameters:
        
        vector1: list of doubles
            3-vector specifying the first vector spanning the cell :math:`\mathbf{v}_1`
        vector2: list of doubles
            3-vector specifying the second vector spanning the cell :math:`\mathbf{v}_2`
        vector3: list of doubles
            3-vector specifying the third vector spanning the cell :math:`\mathbf{v}_3`
        pbc: list of logicals
            three logic switches for specifying periodic boundaries - True denotes periodicity
        """
    
    def __init__(self,vector1,vector2,vector3,pbc):
        self.vector1 = np.array( vector1 )
        self.vector2 = np.array( vector2 )
        self.vector3 = np.array( vector3 )
        self.matrix = np.array( [vector1,vector2,vector3] ).transpose()
        self.inverse = np.linalg.inv(self.matrix)
        self.pbc = pbc
    
    
    def get_relative_coordinates(self,coordinates):
        """Returns the coordinates of the given atom in fractional coordinates.
            
            The absolute position of the atom is given by multiplying the cell vectors
            by the fractional coordinates.
            
            Parameters:
            
            coordinates: numpy double 3-vector
                the absolute coordinates 
            """
        return self.inverse.dot(coordinates)
    
    
    def get_absolute_coordinates(self,fractional):
        """For the given fractional coordinates, returns the absolute coordinates.
            
            The absolute coordinates are the cell vectors multiplied by the fractional
            coordinates.
            
            Parameters:
            
            coordinates: numpy double 3-vector
                the fractional coordinates 
            """
        return self.matrix.dot(fractional)
    
    def get_wrapped_coordinates(self,coordinates):
        """Wraps the coordinates of the given atom inside the simulation cell.
            
            This method return the equivalent position (with respect to the periodic boundaries)
            of the atom inside the cell.
            
            For instance, if the cell is spanned by vectors of length 10.0 in 
            directions of :math:`x`, :math:`y`, and :math:`z`,
            an the coordinates [-1.0, 12.0, 3.0] wrap to [9.0, 2.0, 3.0].
            
            Parameters:
            
            coordinates: numpy double 3-vector
                the absolute coordinates 
            """
        rel_coord = self.get_relative_coordinates(coordinates)
        new_coord = np.array([0.0,0.0,0.0])
        for ind in range(3):
            if self.pbc[ind]:
                new_coord[ind] = rel_coord[ind] % 1.0
            else:
                new_coord[ind] = rel_coord[ind]

        return self.get_absolute_coordinates(new_coord)
    
    
    def estimate_splits(self,max_cutoff):
        """Estimates how many times can the supercell be divided in the directions of the spanning vectors to ensure neighboring subcells contain all points at the given maximum cutoff.
            
            The :class:`~pysic_utility.FastNeighborList` requires a division of the
            system in subcells.
            """
        splits = np.array([0,0,0])
        lengths = np.array([0.0,0.0,0.0]) 
        normal = np.cross(self.vector2,self.vector3)
        lengths[0] = fabs(np.dot(self.vector1,normal)/sqrt(np.dot(normal,normal)))
        normal = np.cross(self.vector1,self.vector3)
        lengths[1] = fabs(np.dot(self.vector2,normal)/sqrt(np.dot(normal,normal)))
        normal = np.cross(self.vector1,self.vector2)
        lengths[2] = fabs(np.dot(self.vector3,normal)/sqrt(np.dot(normal,normal)))
    
        for i in range(3):
            splits[i] = int( lengths[i]/max_cutoff )
    
        return splits
    
    
    def split(self,splits):
        """Split the cell in subcells according to the given number of divisions.
            
            The argument 'splits' should be a list of three integers determining how many
            times the cell is split. For instance, if splits = [3,3,5], the cell is divided in
            3*3*5 = 45 subcells: 3 cells along the first two cell vectors and 5 along the third.
            
            The Cell itself is not changed, but an array 'subcells' is created, containing
            the subcells which are Cell instances themselves. These cells will contain additional
            data arrays 'neighbors' and 'offsets'. These are 3-dimensional arrays with each dimension
            running from -1 to 1. The neighbors array contains references to the neighboring subcell
            Cell instances.
            The offsets contain coordinate offsets with respect to the periodic boundaries. In other words,
            if a subcell is at the border of the original Cell, it will have neighbors at the other side
            of the cell due to periodic boundary conditions. But from the point of view of the subcell,
            the neighboring cell is not on the other side of the master cell, but a periodic image of that
            cell. Therefore, any coordinates in the the subcell to which the neighbors array refers to must
            in fact be shifted by a vector of the master cell. The offsets list contains the multipliers
            for the cell vectors to make these shifts.
            
            Example in 2D for simplicity:
            
            split [3,4]
            
            creates subcells
            
            (0,3) (1,3) (2,3)
            (0,2) (1,2) (2,2)
            (0,1) (1,1) (2,1)
            (0,0) (1,0) (2,0)
            
            subcell (0,3) will have the neighbors 
            (2,0) (0,0) (1,0)
            (2,3) (0,3) (1,3)
            (2,2) (0,2) (1,2)
            
            and offsets
            [-1,1] [0,1] [0,1]
            [-1,0] [0,0] [0,0]
            [-1,0] [0,0] [0,0]
            
            Note that the central 'neighbor' is the cell itself.
            
            If a boundary is not periodic, extra subcells with indices -1 and split
            are created to pad the simulation cell. These will contain the atoms that
            are outside the simulation cell.
            """
        new_vector1 = self.vector1/splits[0]
        new_vector2 = self.vector2/splits[1]
        new_vector3 = self.vector3/splits[2]        
        
        self.subcells = [[[None for i in range(-1,splits[2]+1)] for j in range(-1,splits[1]+1)] for k in range(-1,splits[0]+1) ]
        
        for i in range(-1,splits[0]+1):
            for j in range(-1,splits[1]+1):
                for k in range(-1,splits[2]+1):
                    self.subcells[i][j][k] = Cell(new_vector1,new_vector2,new_vector3,
                                                  [False,False,False])
                    self.subcells[i][j][k].subcell_indices = (i,j,k)
                    self.subcells[i][j][k].neighbors = [[[None for ic in range(-1,2)] for jc in range(-1,2)] for kc in range(-1,2)]
                    self.subcells[i][j][k].offsets = [[[None for ic in range(-1,2)] for jc in range(-1,2)] for kc in range(-1,2)]
                    self.subcells[i][j][k].include_nbor = [[[True for ic in range(-1,2)] for jc in range(-1,2)] for kc in range(-1,2)]
        
        
        for i in range(-1,splits[0]+1):
            for j in range(-1,splits[1]+1):
                for k in range(-1,splits[2]+1):
                    
                    for i_nbor in range(-1,2):
                        for j_nbor in range(-1,2):
                            for k_nbor in range(-1,2):
                                
                                i_nbor_index = (i+i_nbor)
                                j_nbor_index = (j+j_nbor)
                                k_nbor_index = (k+k_nbor)
                                                                                                    
                                offsets = [0]*3
                                
                                if self.pbc[0]:
                                    
                                    # with periodicity, record if we 'loop' around the 
                                    # periodic boundary
                                    if i_nbor_index == splits[0]:
                                        i_nbor_index -= splits[0]
                                        offsets[0] = 1
                                    elif i_nbor_index == -1:
                                        i_nbor_index += splits[0]
                                        offsets[0] = -1
                                        
                                else:
                                    
                                    # If the boundary is not periodic, make a note
                                    # if the neighbor index points to a non-existing
                                    # subcell.
                                    if i_nbor_index < -1 or i_nbor_index > splits[0]:
                                        self.subcells[i][j][k].include_nbor[i_nbor][j_nbor][k_nbor] = False
                                            
                    
                                if self.pbc[1]:
                                    
                                    if j_nbor_index == splits[1]:
                                        j_nbor_index -= splits[1]
                                        offsets[1] = 1
                                    elif j_nbor_index == -1:
                                        j_nbor_index += splits[1]
                                        offsets[1] = -1
                                        
                                else:
                                    if j_nbor_index < -1 or j_nbor_index > splits[1]:
                                        self.subcells[i][j][k].include_nbor[i_nbor][j_nbor][k_nbor] = False
                                        
                                if self.pbc[2]:
                                    
                                    if k_nbor_index == splits[2]:
                                        k_nbor_index -= splits[2]
                                        offsets[2] = 1
                                    elif k_nbor_index == -1:
                                        k_nbor_index += splits[2]
                                        offsets[2] = -1
                                    
                                else:
                                    if k_nbor_index < -1 or k_nbor_index > splits[2]:
                                        self.subcells[i][j][k].include_nbor[i_nbor][j_nbor][k_nbor] = False
                                
                                self.subcells[i][j][k].neighbors[i_nbor][j_nbor][k_nbor] = \
                                    self.subcells[i_nbor_index][j_nbor_index][k_nbor_index]
                                self.subcells[i][j][k].offsets[i_nbor][j_nbor][k_nbor] = \
                                    np.array( offsets )
    
    
    
    
    
    def find_subcell_atoms(self,atoms):
        """Allocates atoms to subcells.
            
            Before running this, make sure to first run the method :meth:`~pysic_utility.Cell.split`.
            
            This method searches for each atom the subcell in which the atom
            is. The subcell :class:`~pysic_utility.Cell` instances are given lists of the Atom instances
            in them, named ``atoms``. The Atom instances are given a data field ``subcell_indices``
            to denote the subcell in which they are. (This is just a triplet of integers, not
            a reference to the actual Cell.)
            """
        sub_vector1 = self.subcells[0][0][0].vector1
        sub_vector2 = self.subcells[0][0][0].vector2
        sub_vector3 = self.subcells[0][0][0].vector3
        split1 = len(self.subcells)-2
        split2 = len(self.subcells[0])-2
        split3 = len(self.subcells[0][0])-2
        rel_length1 = 1.0/split1
        rel_length2 = 1.0/split2
        rel_length3 = 1.0/split3
        self.atom_subcell_indices = [ None for i in range(len(atoms)) ]

        
        for i in range(-1,split1+1):
            for j in range(-1,split2+1):
                for k in range(-1,split3+1):
                    self.subcells[i][j][k].atoms = []
        
        for atom in atoms:
            coord = self.get_relative_coordinates(self.get_wrapped_coordinates(atom.position))
            
            if self.pbc[0]:
                sub_i = int(floor(coord[0]*split1)) % split1
            else:
                sub_i = min(split1+1,max(-1,int(floor(coord[0]*split1))))

            if self.pbc[1]:
                sub_j = int(floor(coord[1]*split2)) % split2
            else:
                sub_j = min(split2+1,max(-1,int(floor(coord[1]*split2))))
            
            if self.pbc[2]:
                sub_k = int(floor(coord[2]*split3)) % split3
            else:
                sub_k = min(split3+1,max(-1,int(floor(coord[2]*split3))))
                    
            self.subcells[sub_i][sub_j][sub_k].atoms.append(atom)
            self.atom_subcell_indices[atom.index] = (sub_i,sub_j,sub_k)
                             

    def get_atom_subcell(self,index):
        return self.atom_subcell_indices[index]
    
            
    def get_distance(self,atom1,atom2,offsets=None):
        """Calculates the distance between two atoms.
            
            Offsets are multipliers for the cell vectors to be added to the
            plain separation vector r1-r2 between the atoms.
            """
        vec = self.get_separation(atom1,atom2,offsets)
        return sqrt(vec.dot(vec))

                
    def get_separation(self, atom1, atom2, offsets=None):
        """Returns the separation vector between two atoms, r1-r2.
            
            Offsets are multipliers for the cell vectors to be added to the
            plain separation vector r1-r2 between the atoms.
            """
        p1 = self.get_wrapped_coordinates(atom1.get_position())
        p2 = self.get_wrapped_coordinates(atom2.get_position())
        if offsets == None:
            return p1-p2
        else:
            return p1-p2 - self.matrix.dot(offsets.transpose())


class FastNeighborList(NeighborList):
    """ASE has a neighbor list class built in, but its implementation is
        currently inefficient, and building of the list is an :math:`O(n^2)`
        operation. This neighbor list class overrides the 
        :meth:`~pysic_utility.FastNeighborList.build` method with
        an :math:`O(n)` time routine. The fast routine is based on 
            
            - spacial partition of the simulation cell
            - locating the containing subcell for each particle
            - locating the neighbors for each atom by only examining the nearby subcells
        
        The way cutoffs are handled is also somewhat different to the original
        ASE list. In ASE, the distances for two atoms are compared against
        the sum of the individual cutoffs + neighbor list skin. This list, however,
        searches for the neighbors of each atom at a distance of the cutoff of the
        given atom only, plus skin.
        """
        
    def __init__(self, cutoffs, skin=neighbor_marginal):
        NeighborList.__init__(self, 
                                cutoffs=cutoffs, 
                                skin=skin, 
                                sorted=False, 
                                self_interaction=False,
                                bothways=True)
        self.cell = None
        self.cellmatrix = None
    
    
    def build(self,atoms):
        self.positions = atoms.get_positions()
        self.pbc = atoms.get_pbc()
        
        if len(self.cutoffs) > 0:
            rcmax = self.cutoffs.max()
        else:
            rcmax = 1.0

        if self.cellmatrix == None or self.cellmatrix != atoms.get_cell():
            self.cellmatrix = atoms.get_cell()
            self.cell = Cell(self.cellmatrix[0],
                             self.cellmatrix[1],
                             self.cellmatrix[2],
                             self.pbc)
            
            # split the cell in subcells
            splits = self.cell.estimate_splits(rcmax)
            self.cell.split(splits)
            self.cell.find_subcell_atoms(atoms)
    
                            
                            
        neighbors = [[] for a in range(len(atoms))]
        displacements = [[] for a in range(len(atoms))]
    
        for atom in atoms:
            print "atom ", atom.index
            sub_ind = self.cell.atom_subcell_indices[atom.index]
            subcell = self.cell.subcells[sub_ind[0]][sub_ind[1]][sub_ind[2]]
                            
                            
            for i_nbor in range(-1,2):
                for j_nbor in range(-1,2):
                    for k_nbor in range(-1,2):
                            
                        nbor_cell = subcell.neighbors[i_nbor][j_nbor][k_nbor]
                        nbor_offset = subcell.offsets[i_nbor][j_nbor][k_nbor]
                        nbor_include = subcell.include_nbor[i_nbor][j_nbor][k_nbor]
                        
                        if(nbor_include):

                            for atom2 in nbor_cell.atoms:
                                """
                                if atom.index == 0:
                                    print "***", i_nbor, j_nbor, k_nbor
                                    print nbor_cell.subcell_indices
                                    print nbor_offset
                                    print nbor_include
                                    print atom.index, atom2.index
                                    print atom.position, atom2.position
                                    print self.cell.get_relative_coordinates(atom.position), self.cell.get_relative_coordinates(atom2.position)
                                    print self.cell.get_wrapped_coordinates(atom.position), self.cell.get_wrapped_coordinates(atom2.position)
                                    print self.cell.get_relative_coordinates(self.cell.get_wrapped_coordinates(atom.position)), self.cell.get_relative_coordinates(self.cell.get_wrapped_coordinates(atom2.position))
                                    print self.cell.get_absolute_coordinates(nbor_offset)
                                    print self.cell.get_separation(atom,atom2,nbor_offset)
                                    print self.cell.get_distance(atom,atom2,nbor_offset)
                                """
                                if(atom.index < atom2.index):
                                    d = self.cell.get_distance(atom,atom2,nbor_offset)
                                    if d < self.cutoffs[atom.index]+self.skin:
                                        neighbors[atom.index].append(atom2.index)
                                        displacements[atom.index].append(np.array(nbor_offset))
                                    if d < self.cutoffs[atom2.index]+self.skin:
                                        neighbors[atom2.index].append(atom.index)
                                        displacements[atom2.index].append(np.array(nbor_offset))

                            
        self.nupdates += 1
        self.neighbors = np.array(neighbors)
        self.displacements = np.array(displacements)
    
    
                                
    # the old build method
    def brute_force_build(self, atoms):
        """Build the list."""
        self.positions = atoms.get_positions()
        self.pbc = atoms.get_pbc()
        self.cell = atoms.get_cell()
        if len(self.cutoffs) > 0:
            rcmax = self.cutoffs.max()
        else:
            rcmax = 0.0
        
        icell = np.linalg.inv(self.cell)
        scaled = np.dot(self.positions, icell)
        scaled0 = scaled.copy()
        
        N = []
        for i in range(3):
            if self.pbc[i]:
                scaled0[:, i] %= 1.0
                v = icell[:, i]
                h = 1 / sqrt(np.dot(v, v))
                n =  int(2 * rcmax / h) + 1
            else:
                n = 0
            N.append(n)
        
        offsets = np.empty((len(atoms), 3), int)
        (scaled0 - scaled).round(out=offsets)
        positions0 = np.dot(scaled0, self.cell)
        natoms = len(atoms)
        indices = np.arange(natoms)
        
        self.nneighbors = 0
        self.npbcneighbors = 0
        self.neighbors = [np.empty(0, int) for a in range(natoms)]
        self.displacements = [np.empty((0, 3), int) for a in range(natoms)]
        for n1 in range(0, N[0] + 1):
            for n2 in range(-N[1], N[1] + 1):
                for n3 in range(-N[2], N[2] + 1):
                    if n1 == 0 and (n2 < 0 or n2 == 0 and n3 < 0):
                        continue
                    displacement = np.dot((n1, n2, n3), self.cell)
                    for a in range(natoms):
                        d = positions0 + displacement - positions0[a]
                        i = indices[(d**2).sum(1) <
                                    (self.cutoffs + self.cutoffs[a])**2]
                        if n1 == 0 and n2 == 0 and n3 == 0:
                            if self.self_interaction:
                                i = i[i >= a]
                            else:
                                i = i[i > a]
                        self.nneighbors += len(i)
                        self.neighbors[a] = np.concatenate(
                                                           (self.neighbors[a], i))
                        disp = np.empty((len(i), 3), int)
                        disp[:] = (n1, n2, n3)
                        disp += offsets[i] - offsets[a]
                        self.npbcneighbors += disp.any(1).sum()
                        self.displacements[a] = np.concatenate(
                                                               (self.displacements[a], disp))
        
        if self.bothways:
            neighbors2 = [[] for a in range(natoms)]
            displacements2 = [[] for a in range(natoms)]
            for a in range(natoms):
                for b, disp in zip(self.neighbors[a], self.displacements[a]):
                    neighbors2[b].append(a)
                    displacements2[b].append(-disp)
            for a in range(natoms):
                self.neighbors[a] = np.concatenate((self.neighbors[a],
                                                    neighbors2[a]))
                self.displacements[a] = np.array(list(self.displacements[a]) +
                                                 displacements2[a])
        
        if self.sorted:
            for a, i in enumerate(self.neighbors):
                mask = (i < a)
                if mask.any():
                    j = i[mask]
                    offsets = self.displacements[a][mask]
                    for b, offset in zip(j, offsets):
                        self.neighbors[b] = np.concatenate(
                                                           (self.neighbors[b], [a]))
                        self.displacements[b] = np.concatenate(
                                                               (self.displacements[b], [-offset]))
                    mask = np.logical_not(mask)
                    self.neighbors[a] = self.neighbors[a][mask]
                    self.displacements[a] = self.displacements[a][mask]
        
        self.nupdates += 1



def char2int(char_in):
    """Codes a single character to an integer.
    """
    output = codec_s2i.get(char_in)
    if output == None:
        return 0
    else:
        return output

def int2char(int_in):
    """Decodes an integer to a single character.
    """
    output = codec_i2s.get(int_in)
    if output == None:
        return ' '
    else:
        return output

def str2ints(string_in,target_length=0):
    """Codes a string to a list of integers.

    Turns a string to a list of integers for f2py interfacing.
    If required, the length of the list can be specified and trailing spaces
    will be added to the end.
    """
    ints_out = []
    for char in string_in:
        ints_out.append(char2int(char))
    while len(ints_out) < target_length:
        ints_out.append(char2int(' '))
    return ints_out

def ints2str(ints_in):
    """Decodes a list of integers to a string.
    """
    string_out = ""
    for number in ints_in:
        string_out += int2char(number)
    return string_out


def expand_symbols_table(symbol_list,type=None):
    """Creates a table of symbols for a BondOrderParameters object.
            
            The syntax for defining the targets of bond order factors is precise
            but somewhat cumbersome due to the large number of permutations one gets
            when the number of bodies increases. Oftentimes one does not need such
            fine control over all the parameters since many of them have the same
            numerical values. Therefore it is convenient to be able to define
            the targets in a more compact way.
            
            This method generates the detailed target tables from compact syntax.
            By default, the method takes a list of list and multiplies each list
            with the others (note the call for a static method)::
            
             >>> pysic.BondOrderParameters.expand_symbols_table([  'Si',
             ...                                                  ['O', 'C'],
             ...                                                  ['H', 'O'] ])
             [['Si', 'O', 'H'],
              ['Si', 'C', 'H'],
              ['Si', 'O', 'O'],
              ['Si', 'C', 'O']]
            
            Other custom types of formatting can be defined with the type parameter.
            
            For type 'triplet', the target list is created for triplets A-B-C from an input list of the form::
            
             ['A', 'B', 'C']
            
            Remember that in the symbol table accepted by the BondOrderParameters, one needs to define
            the B-A and B-C bonds separately and so B appears as the first symbol in the output and the other
            two appear as second and third (both cases)::
            
             [['B', 'A', 'C'],
              ['B', 'C', 'A']]
            
            However, for an A-B-A triplet, the A-B bond should only be defined once to prevent double counting.
            Like the default function, also here several triplets can be defined at once::
            
             >>> pysic.BondOrderParameters.expand_symbols_table([ ['H', 'O'],
             ...                                                   'Si',
             ...                                                  ['O', 'C'] ],
             ...                                                type='triplet')
             [['Si', 'H', 'O'],
              ['Si', 'O', 'H'],
              ['Si', 'H', 'C'],
              ['Si', 'C', 'H'],
              ['Si', 'O', 'O'],
              ['Si', 'O', 'C'],
              ['Si', 'C', 'O']]
            
            
            Parameters:
            
            symbol_list: list of strings
            list to be expanded to a table
            type: string
            specifies a custom way of generating the table
            """
    if not isinstance(symbol_list,list):
        return [symbol_list]
    n_slots = len(symbol_list)
    n_subslots = []
    for sub in symbol_list:
        if isinstance(sub,list):
            n_subslots.append(len(sub))
        else:
            n_subslots.append(1)
        
    table = []
    if(type == None):
        slot_indices = n_slots*[0]
        while (slot_indices[n_slots-1] < n_subslots[n_slots-1]):
            row = []
            for i in range(n_slots):
                if isinstance(symbol_list[i],list):
                    symb = symbol_list[i][slot_indices[i]]
                else:
                    symb = symbol_list[i]
                row.append(symb)
            table.append(row)
            slot_indices[0] += 1
            for i in range(n_slots):
                if slot_indices[i] >= n_subslots[i]:
                    if(i < n_slots-1):
                        slot_indices[i] = 0
                        slot_indices[i+1] += 1
                    else:
                        pass
    elif(type == 'triplet'):
        if n_slots != 3:
            return None
        for i in range(n_subslots[1]):
            for j in range(n_subslots[0]):
                for k in range(n_subslots[2]):
                    
                    if isinstance(symbol_list[1],list):
                        symb1 = symbol_list[1][i]
                    else:
                        symb1 = symbol_list[1]
                    if isinstance(symbol_list[0],list):
                        symb0 = symbol_list[0][j]
                    else:
                        symb0 = symbol_list[0]
                    if isinstance(symbol_list[2],list):
                        symb2 = symbol_list[2][k]
                    else:
                        symb2 = symbol_list[2]
                        
                    table.append([symb1,symb0,symb2])
                    if symb0 != symb2:
                        table.append([symb1,symb2,symb0])
        
        
    return table




def plot_energy_on_line(index,system,direction=None,length=None,steps=100,start=None,end=None,lims=[-1e10,1e10]):
    """Plots the energy of the system as a function of the position of a single particle.

    The method probes the system by moving a single particle on a line
    and recording the energy. A plot is drawn. Also a tuple containing arrays of
    the distance traveled and the recorded energies is returned.

    After the operation is complete, the initial structure is restored.

    .. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html

    Parameters:
    
    index: integer
        index of the particle to be moved
    system: `ASE Atoms`_ object
        the structure to be explored
    direction: double 3-vector
        the direction where the atom is moved
    length: double
        the distance moved
    steps: integer
        number of points (taken uniformly on the movement path) for measuring the energy
    start: double 3-vector (array or list)
        starting point for the trajectory - if not specified, the position of the particle in 'system' is used
    end: double 3-vector (array or list)
        end point for the trajectory - alternative for direction and length (will override them)
    lims: double 2-vector (array or list)
        lower and upper truncation limits - if a recorded value is smaller than the lower limit or larger than the upper, it is replaced by the corresponding truncation value
    """
    x = np.array([0.0]*(steps+1))
    y = np.array([0.0]*(steps+1))

    orig_pos = system.get_positions()[index]
    if start == None:
        start = [ system[index].x,
                  system[index].y,
                  system[index].z ]
    if end != None:
        direction = np.array(end)-np.array(start)
        length = sqrt( direction.dot(direction) )
    if direction == None:
        print "Either specify direction and length or an end point."
        return
    if length == None:
        length = sqrt( direction.dot(direction) )
        
    dire = np.array(direction)
    unit = dire / sqrt( dire.dot(dire) )
    delta = unit * (length+0.0)/steps
    dx = (length+0.0)/steps
    
    for i in range(steps+1):
        xval = i*dx
        new_pos = start + i * delta
        system[index].x = new_pos[0]
        system[index].y = new_pos[1]
        system[index].z = new_pos[2]
        x[i] = xval
        value = system.get_potential_energy()
        y[i] = max(min(value,lims[1]),lims[0])

    system[index].x = orig_pos[0]
    system[index].y = orig_pos[1]
    system[index].z = orig_pos[2]
    system.get_calculator().update_core_coordinates()

    plt.plot(x,y)
    return (x,y)

def plot_abs_force_on_line(index,system,direction=None,length=None,steps=100,start=None,end=None,lims=[-1e10,1e10]):
    """Plots the absolute value of the force on a particle as a function of the position.

    The method probes the system by moving a single particle on a line
    and recording the force. A plot is drawn. Also a tuple containing arrays of
    the distance traveled and the recorded forces is returned.

    After the operation is complete, the initial structure is restored.

    .. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html

    Parameters:
    
    index: integer
        index of the particle to be moved
    system: `ASE Atoms`_ object
        the structure to be explored
    direction: double 3-vector
        the direction where the atom is moved
    length: double
        the distance moved
    steps: integer
        number of points (taken uniformly on the movement path) for measuring the force
    start: double 3-vector (array or list)
        starting point for the trajectory - if not specified, the position of the particle in 'system' is used
    end: double 3-vector (array or list)
        end point for the trajectory - alternative for direction and length (will override them)
    lims: double 2-vector (array or list)
        lower and upper truncation limits - if a recorded value is smaller than the lower limit or larger than the upper, it is replaced by the corresponding truncation value
    """
    x = np.array([0.0]*(steps+1))
    y = np.array([0.0]*(steps+1))

    orig_pos = system.get_positions()[index]
    if start == None:
        start = [ system[index].x,
                  system[index].y,
                  system[index].z ]
    if end != None:
        direction = np.array(end)-np.array(start)
        length = sqrt( direction.dot(direction) )
    if direction == None:
        print "Either specify direction and length or an end point."
        return
    if length == None:
        length = sqrt( direction.dot(direction) )
        
    dire = np.array(direction)
    unit = dire / sqrt( dire.dot(dire) )
    delta = unit * (length+0.0)/steps
    dx = (length+0.0)/steps
    
    for i in range(steps+1):
        xval = i*dx
        new_pos = start + i * delta
        system[index].x = new_pos[0]
        system[index].y = new_pos[1]
        system[index].z = new_pos[2]
        x[i] = xval
        value = system.get_forces()[index]
        y[i] = max(min(sqrt(value.dot(value)),lims[1]),lims[0])

    system[index].x = orig_pos[0]
    system[index].y = orig_pos[1]
    system[index].z = orig_pos[2]
    system.get_calculator().update_core_coordinates()

    plt.plot(x,y)
    return (x,y)

def plot_tangent_force_on_line(index,system,direction=None,length=None,steps=100,start=None,end=None,lims=[-1e10,1e10]):
    """Plots the tangential force on a particle as a function of the position.

    The method probes the system by moving a single particle on a line
    and recording the force tangent. A plot is drawn. Also a tuple containing arrays of
    the distance traveled and the recorded forces is returned.

    After the operation is complete, the initial structure is restored.

    .. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html

    Parameters:
    
    index: integer
        index of the particle to be moved
    system: `ASE Atoms`_ object
        the structure to be explored
    direction: double 3-vector
        the direction where the atom is moved
    length: double
        the distance moved
    steps: integer
        number of points (taken uniformly on the movement path) for measuring the energy
    start: double 3-vector (array or list)
        starting point for the trajectory - if not specified, the position of the particle in 'system' is used
    end: double 3-vector (array or list)
        end point for the trajectory - alternative for direction and length (will override them)
    lims: double 2-vector (array or list)
        lower and upper truncation limits - if a recorded value is smaller than the lower limit or larger than the upper, it is replaced by the corresponding truncation value
    """
    x = np.array([0.0]*(steps+1))
    y = np.array([0.0]*(steps+1))

    orig_pos = system.get_positions()[index]
    if start == None:
        start = [ system[index].x,
                  system[index].y,
                  system[index].z ]
    if end != None:
        direction = np.array(end)-np.array(start)
        length = sqrt( direction.dot(direction) )
    if direction == None:
        print "Either specify direction and length or an end point."
        return
    if length == None:
        length = sqrt( direction.dot(direction) )
        
    dire = np.array(direction)
    unit = dire / sqrt( dire.dot(dire) )
    delta = unit * (length+0.0)/steps
    dx = (length+0.0)/steps
    
    for i in range(steps+1):
        xval = i*dx
        new_pos = start + i * delta
        system[index].x = new_pos[0]
        system[index].y = new_pos[1]
        system[index].z = new_pos[2]
        x[i] = xval
        value = system.get_forces()[index]
        y[i] = max(min(value.dot(unit),lims[1]),lims[0])

    system[index].x = orig_pos[0]
    system[index].y = orig_pos[1]
    system[index].z = orig_pos[2]
    system.get_calculator().update_core_coordinates()

    plt.plot(x,y)
    return (x,y)

def plot_energy_on_plane(index,system,directions,lengths,steps=100,start=None,lims=[-1e10,1e10]):
    """Plots the energy of the system as a function of the position of a particle.

    The method probes the system by moving a single particle on a plane
    and recording the energy. A contour plot is drawn. Also a tuple containing arrays of
    the distances traveled on the plane and the recorded energies is returned.

    After the operation is complete, the initial structure is restored.

    .. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html

    Parameters:
    
    index: integer
        index of the particle to be moved
    system: `ASE Atoms`_ object
        the structure to be explored
    directions: double 2x3-matrix
        The directions where the atom is moved - i.e., the vectors defining the plane. If the second vector is not perpendicular to the first, the normal component is automatically used instead.
    lengths: double 2-vector
        the distances moved in the given directions
    steps: integer
        number of points (taken uniformly on the movement plane) for measuring the energy
    start: double 3-vector (array or list)
        starting point for the trajectory, i.e., a corner for the plane to be probed - if not specified, the position of the particle in 'system' is used
    lims: double 2-vector (array or list)
        lower and upper truncation limits - if a recorded value is smaller than the lower limit or larger than the upper, it is replaced by the corresponding truncation value
    """
    x = np.array([0.0]*(steps+1))
    y = np.array([0.0]*(steps+1))
    z = np.array([[0.0]*(steps+1)]*(steps+1))

    orig_pos = system.get_positions()[index]
    if start == None:
        start = [ system[index].x,
                  system[index].y,
                  system[index].z ]
    dir1 = np.array(directions[0])
    dir2 = np.array(directions[1])
    dir2 = dir2 - dir2 * dir2.dot( dir1 ) # take the perpendicular component
    unit1 = dir1 / sqrt( dir1.dot(dir1) )
    unit2 = dir2 / sqrt( dir2.dot(dir2) )
    unit_perp = np.array( [[unit1[1]*unit2[2] - unit1[2]*unit2[1]],
                           [unit1[2]*unit2[0] - unit1[0]*unit2[2]],
                           [unit1[0]*unit2[1] - unit1[1]*unit2[0]]] )
    delta1 = unit1 * (lengths[0]+0.0)/steps
    delta2 = unit2 * (lengths[1]+0.0)/steps
    dx = (lengths[0]+0.0)/steps
    dy = (lengths[1]+0.0)/steps
    
    for i in range(steps+1):
        xval = i*dx
        x[i] = xval       
        for j in range(steps+1):
            yval = j*dx
            y[j] = yval
            new_pos = start + i * delta1 + j * delta2
            system[index].x = new_pos[0]
            system[index].y = new_pos[1]
            system[index].z = new_pos[2]
            value = system.get_potential_energy()
            z[i,j] = max(min(value,lims[1]),lims[0])

    system[index].x = orig_pos[0]
    system[index].y = orig_pos[1]
    system[index].z = orig_pos[2]
    system.get_calculator().update_core_coordinates()

    plt.contourf(x,y,z)
    return (x,y,z)

def plot_abs_force_on_plane(index,system,directions,lengths,steps=100,start=None,lims=[-1e10,1e10]):
    """Plots the absolute value of force on a particle as a function of the position.

    The method probes the system by moving a single particle on a plane
    and recording the force. A contour plot is drawn. Also a tuple containing arrays of
    the distances traveled on the plane and the recorded forces is returned.

    After the operation is complete, the initial structure is restored.

    .. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html

    Parameters:
    
    index: integer
        index of the particle to be moved
    system: `ASE Atoms`_ object
        the structure to be explored
    directions: double 2x3-matrix
        The directions where the atom is moved - i.e., the vectors defining the plane. If the second vector is not perpendicular to the first, the normal component is automatically used instead.
    lengths: double 2-vector
        the distances moved in the given directions
    steps: integer
        number of points (taken uniformly on the movement plane) for measuring the force
    start: double 3-vector (array or list)
        starting point for the trajectory, i.e., a corner for the plane to be probed - if not specified, the position of the particle in 'system' is used
    lims: double 2-vector (array or list)
        lower and upper truncation limits - if a recorded value is smaller than the lower limit or larger than the upper, it is replaced by the corresponding truncation value
    """
    x = np.array([0.0]*(steps+1))
    y = np.array([0.0]*(steps+1))
    z = np.array([[0.0]*(steps+1)]*(steps+1))

    orig_pos = system.get_positions()[index]
    if start == None:
        start = [ system[index].x,
                  system[index].y,
                  system[index].z ]
    dir1 = np.array(directions[0])
    dir2 = np.array(directions[1])
    dir2 = dir2 - dir2 * dir2.dot( dir1 ) # take the perpendicular component
    unit1 = dir1 / sqrt( dir1.dot(dir1) )
    unit2 = dir2 / sqrt( dir2.dot(dir2) )
    unit_perp = np.array( [[unit1[1]*unit2[2] - unit1[2]*unit2[1]],
                           [unit1[2]*unit2[0] - unit1[0]*unit2[2]],
                           [unit1[0]*unit2[1] - unit1[1]*unit2[0]]] )
    delta1 = unit1 * (lengths[0]+0.0)/steps
    delta2 = unit2 * (lengths[1]+0.0)/steps
    dx = (lengths[0]+0.0)/steps
    dy = (lengths[1]+0.0)/steps
    
    for i in range(steps+1):
        xval = i*dx
        x[i] = xval       
        for j in range(steps+1):
            yval = j*dx
            y[j] = yval
            new_pos = start + i * delta1 + j * delta2
            system[index].x = new_pos[0]
            system[index].y = new_pos[1]
            system[index].z = new_pos[2]
            value = system.get_forces()[index]
            z[i,j] = max(min(sqrt( value.dot(value) ),lims[1]),lims[0])

    system[index].x = orig_pos[0]
    system[index].y = orig_pos[1]
    system[index].z = orig_pos[2]
    system.get_calculator().update_core_coordinates()

    plt.contourf(x,y,z)
    return (x,y,z)

def plot_tangent_force_on_plane(index,system,directions,lengths,steps=100,start=None,lims=[-1e10,1e10]):
    """Plots the absolute value of the tangent component of force on a particle as a function of the position.

    The method probes the system by moving a single particle on a plane
    and recording the force. The force is projected on the same plane, and the absolute
    value of the projection is calculated.
    A contour plot is drawn. Also a tuple containing arrays of
    the distances traveled on the plane and the recorded forces is returned.

    After the operation is complete, the initial structure is restored.

    .. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html

    Parameters:
    
    index: integer
        index of the particle to be moved
    system: `ASE Atoms`_ object
        the structure to be explored
    directions: double 2x3-matrix
        The directions where the atom is moved - i.e., the vectors defining the plane. If the second vector is not perpendicular to the first, the normal component is automatically used instead.
    lengths: double 2-vector
        the distances moved in the given directions
    steps: integer
        number of points (taken uniformly on the movement plane) for measuring the force
    start: double 3-vector (array or list)
        starting point for the trajectory, i.e., a corner for the plane to be probed - if not specified, the position of the particle in 'system' is used
    lims: double 2-vector (array or list)
        lower and upper truncation limits - if a recorded value is smaller than the lower limit or larger than the upper, it is replaced by the corresponding truncation value
    """
    x = np.array([0.0]*(steps+1))
    y = np.array([0.0]*(steps+1))
    z = np.array([[0.0]*(steps+1)]*(steps+1))

    orig_pos = system.get_positions()[index]
    if start == None:
        start = [ system[index].x,
                  system[index].y,
                  system[index].z ]
    dir1 = np.array(directions[0])
    dir2 = np.array(directions[1])
    dir2 = dir2 - dir2 * dir2.dot( dir1 ) # take the perpendicular component
    unit1 = dir1 / sqrt( dir1.dot(dir1) )
    unit2 = dir2 / sqrt( dir2.dot(dir2) )
    unit_perp = np.array( [[unit1[1]*unit2[2] - unit1[2]*unit2[1]],
                           [unit1[2]*unit2[0] - unit1[0]*unit2[2]],
                           [unit1[0]*unit2[1] - unit1[1]*unit2[0]]] )
    delta1 = unit1 * (lengths[0]+0.0)/steps
    delta2 = unit2 * (lengths[1]+0.0)/steps
    dx = (lengths[0]+0.0)/steps
    dy = (lengths[1]+0.0)/steps
    
    for i in range(steps+1):
        xval = i*dx
        x[i] = xval       
        for j in range(steps+1):
            yval = j*dx
            y[j] = yval
            new_pos = start + i * delta1 + j * delta2
            system[index].x = new_pos[0]
            system[index].y = new_pos[1]
            system[index].z = new_pos[2]
            value = system.get_forces()[index]
            tangent = value - value.dot(unit_perp)
            z[i,j] = max(min(sqrt( tangent.dot(tangent) ),lims[1]),lims[0])

    system[index].x = orig_pos[0]
    system[index].y = orig_pos[1]
    system[index].z = orig_pos[2]
    system.get_calculator().update_core_coordinates()

    plt.contourf(x,y,z)
    return (x,y,z)



def plot_force_component_on_plane(index,system,directions,lengths,component,steps=100,start=None,lims=[-1e10,1e10]):
    """Plots the projected component of force on a particle as a function of the position.

    The method probes the system by moving a single particle on a plane
    and recording the force. The component of the force projected on a
    given vector is recorded.
    A contour plot is drawn. Also a tuple containing arrays of
    the distances traveled on the plane and the recorded forces is returned.

    After the operation is complete, the initial structure is restored.

    .. _ASE Atoms: https://wiki.fysik.dtu.dk/ase/ase/atoms.html

    Parameters:
    
    index: integer
        index of the particle to be moved
    system: `ASE Atoms`_ object
        the structure to be explored
    directions: double 2x3-matrix
        The directions where the atom is moved - i.e., the vectors defining the plane. If the second vector is not perpendicular to the first, the normal component is automatically used instead.
    lengths: double 2-vector
        the distances moved in the given directions
    component: double 3-vector
        the direction on which the force is projected - e.g., if component is [1,0,0], the x-component is recorded
    steps: integer
        number of points (taken uniformly on the movement plane) for measuring the force
    start: double 3-vector (array or list)
        starting point for the trajectory, i.e., a corner for the plane to be probed - if not specified, the position of the particle in 'system' is used
    lims: double 2-vector (array or list)
        lower and upper truncation limits - if a recorded value is smaller than the lower limit or larger than the upper, it is replaced by the corresponding truncation value
    """
    x = np.array([0.0]*(steps+1))
    y = np.array([0.0]*(steps+1))
    z = np.array([[0.0]*(steps+1)]*(steps+1))

    orig_pos = system.get_positions()[index]
    if start == None:
        start = [ system[index].x,
                  system[index].y,
                  system[index].z ]
    dir1 = np.array(directions[0])
    dir2 = np.array(directions[1])
    dir2 = dir2 - dir2 * dir2.dot( dir1 ) # take the perpendicular component
    unit1 = dir1 / sqrt( dir1.dot(dir1) )
    unit2 = dir2 / sqrt( dir2.dot(dir2) )
    unit_perp = np.array( [[unit1[1]*unit2[2] - unit1[2]*unit2[1]],
                           [unit1[2]*unit2[0] - unit1[0]*unit2[2]],
                           [unit1[0]*unit2[1] - unit1[1]*unit2[0]]] )
    delta1 = unit1 * (lengths[0]+0.0)/steps
    delta2 = unit2 * (lengths[1]+0.0)/steps
    dx = (lengths[0]+0.0)/steps
    dy = (lengths[1]+0.0)/steps
    component = np.array(component)
    unit_comp = component / sqrt( component.dot(component) )
    
    for i in range(steps+1):
        xval = i*dx
        x[i] = xval       
        for j in range(steps+1):
            yval = j*dx
            y[j] = yval
            new_pos = start + i * delta1 + j * delta2
            system[index].x = new_pos[0]
            system[index].y = new_pos[1]
            system[index].z = new_pos[2]
            value = system.get_forces()[index]
            z[i,j] = max(min(value.dot(component),lims[1]),lims[0])

    system[index].x = orig_pos[0]
    system[index].y = orig_pos[1]
    system[index].z = orig_pos[2]
    system.get_calculator().update_core_coordinates()

    plt.contourf(x,y,z)
    return (x,y,z)
