#! /usr/bin/env python

try:
    import h5py
except:
    print "h5py not found, archiving tools not available"
import numpy as np
from ase import Atoms

class Archive:
    """A class representing a data archive in the hdf5 format.
    
    The archive is built on h5py. See its 
    `documentation <http://docs.h5py.org/en/latest/quick.html>`_
    to access further functionality.

    Parameters:
    
    filename: string
        path to the file storing the data
    """
    def __init__(self, filename):
        self.file = h5py.File(filename)
        self.active_group = self.file
        self.active_dataset = None


    def get_contents(self):
        return self.active_group.keys()

    def list(self):
        for name in self.active_group:
            print name

    def state(self):
        print "current group: "+self.active_group.name
        try:
            print "current dataset: "+self.active_dataset.name
        except:
            print "no dataset accessed"

    def group(self):
        print "current group: "+self.active_group.name
        meta = self.get_group_metadata()
        print "attributes:"
        for attribute in meta:
            print "  "+attribute+" = "+meta[attribute]

        return self.active_group


    def data(self):
        print "current dataset: "+self.active_dataset.name
        meta = self.get_dataset_metadata()
        print "attributes:"
        for attribute in meta:
            print "  "+attribute+" = "+meta[attribute]

        return self.active_dataset

    def get_group_name(self):
        return self.active_group.name

    def get_group(self):
        return self.active_group

    def create_subgroup(self, name):
        new_group = self.active_group.create_group(name)
        #self.active_group = new_group

    def delete_subgroup(self, name):
        del self.active_group[name]

    def move_to_group(self, name):
        self.active_group = self.active_group[name]

    def move_to_parent_group(self):
        self.active_group = self.active_group.parent
        self.active_dataset = None

    def move_to_root_group(self):
        self.active_group = self.file
        self.active_dataset = None

    def delete_group(self):
        delname = self.active_group.name.split('/')[-1]
        self.move_to_parent_group()
        self.delete_subgroup(delname)
        self.active_dataset = None

    def move(source, destination):
        self.active_group.move(source, destination)


    def get_dataset_name(self):
        return self.active_dataset.name

    def get_dataset(self):
        return self.active_dataset

    def access_dataset(self, name):
        self.active_dataset = self.active_group[name]
        return self.active_dataset

    def create_dataset(self, name, data, **kwds):
        dataset = self.active_group.create_dataset(name, data=data, **kwds)
        self.active_dataset = dataset

    def delete_dataset(self, name):
        del self.active_group[name]

    def delete_current_dataset(self):
        del self.file[self.active_dataset.name]



    def get_group_metadata(self):
        return self.active_group.attrs

    def get_dataset_metadata(self):
        return self.active_dataset.attrs

    def add_group_metadata(self, attribute, metadata):
        self.active_group.attrs[attribute] = metadata

    def add_dataset_metadata(self, attribute, metadata):
        self.active_dataset.attrs[attribute] = metadata

    def remove_group_metadata(self, attribute):
        del self.active_group.attrs[attribute]

    def remove_dataset_metadata(self, attribute):
        del self.active_dataset.attrs[attribute]



    def store_system(self, name, atoms, calculate=True):
        """Stores the information of an atomic system.
        
        The routine creates a group with the given name and
        stores the atomic information there as named datasets::
        
         'atomic numbers'
         'positions'
         'charges'
         'magnetic moments'
         'tags'
         'cell'
         'pbc'
         'momenta'
         'potential energy'
         'forces'
         'stress'
         'electronegativites'
         
        The data requiring calculation is stored only if a calculator
        is attached to the atoms, electronegativities only if the
        calculator is :class:`~pysic.calculator.Pysic`.
        
        Parameters:
        
        name: string
            the name of the group for storing the data
        atoms: Atoms object
            the structure to be stored
        calculate: boolean
            if True, the data requiring calculations will not be stored
        """
        
        self.create_subgroup(name)
        self.move_to_group(name)

        self.create_dataset('atomic numbers', atoms.get_atomic_numbers())
        self.create_dataset('positions', atoms.get_positions())
        self.create_dataset('charges', atoms.get_initial_charges())
        self.create_dataset('magnetic moments', atoms.get_initial_magnetic_moments())
        self.create_dataset('tags', atoms.get_tags())
        self.create_dataset('cell', atoms.get_cell())
        self.create_dataset('pbc', atoms.get_pbc())
        self.create_dataset('momenta', atoms.get_momenta())

        if calculate:
            try:
                self.create_dataset('potential energy', atoms.get_potential_energy())
            except:
                pass
            try:
                self.create_dataset('forces', atoms.get_forces())
            except:
                pass
            try:
                self.create_dataset('stress', atoms.get_stress())
            except:
                pass
            try:
                self.create_dataset('electronegativities', atoms.get_calculator().get_electronegativities())
            except:
                pass
        self.move_to_parent_group()
        

    def get_system(self, name):
        """Restore a stored system.
        
        The function returns a tuple with the data::
        
          system, energy, forces, stress, electronegativities
        
        Parameters:
        
        name: string
            name of the group storing the system
        """

        self.move_to_group(name)

        numbers = self.access_dataset('atomic numbers')[:]
        system = Atoms('H'*len(numbers))
        system.set_atomic_numbers(numbers)
        system.set_positions(np.transpose(self.access_dataset('positions')[:,:]))
        system.set_initial_charges(self.access_dataset('charges')[:])
        try:
            system.set_initial_magnetic_moments(self.access_dataset(np.transpose('magnetic moments')[:,:]))
        except:
            pass
        
        system.set_tags(self.access_dataset('tags')[:])
        system.set_cell(self.access_dataset('cell')[:])
        system.set_pbc(self.access_dataset('pbc')[:])
        system.set_momenta(self.access_dataset('momenta')[:,:])

        try:
            ene = self.access_dataset('potential energy')[()]
        except:
            ene = None
        try:
            forces = self.access_dataset('forces')[()]
        except:
            forces = None
        try:
            stress = self.access_dataset('stress')[()]
        except:
            stress = None

        try:
            enega = self.access_dataset('electronegativities')[:]
        except:
            enega = None

        self.move_to_parent_group()
        return (system, ene, forces, stress, enega)






