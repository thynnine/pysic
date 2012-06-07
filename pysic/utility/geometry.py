#! /usr/bin/env python

from math import sqrt, pi, sin, cos, exp, fabs, floor
import numpy as np


class Cell:
    """Cell describing the simulation volume of a subvolume.
         
        This class can be used by the user for coordinate manipulation. Note however,
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
            
            fractional: numpy double 3-vector
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
    
    
                
    def get_distance(self,atom1,atom2,offsets=None):
        """Calculates the distance between two atoms.
            
            Offsets are multipliers for the cell vectors to be added to the
            plain separation vector r1-r2 between the atoms.
            
            Parameters:
             
            atom1: ASE Atoms object
                first atom
            atom2: ASE Atoms object
                second atom
            offsets: Numpy integer 3-vector
                the periodic boundary offsets
            """
        vec = self.get_separation(atom1,atom2,offsets)
        return sqrt(vec.dot(vec))

                
    def get_separation(self, atom1, atom2, offsets=None):
        """Returns the separation vector between two atoms, r1-r2.
            
            Offsets are multipliers for the cell vectors to be added to the
            plain separation vector r1-r2 between the atoms.
            
            Parameters:
            
            atom1: ASE Atoms object
                first atom
            atom2: ASE Atoms object
                second atom
            offsets: Numpy integer 3-vector
                the periodic boundary offsets
            """
        p1 = self.get_wrapped_coordinates(atom1.get_position())
        p2 = self.get_wrapped_coordinates(atom2.get_position())
        if offsets == None:
            return p1-p2
        else:
            return p1-p2 - self.matrix.dot(offsets.transpose())


