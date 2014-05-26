#! /usr/bin/env python
#
# This test uses GPAW for the primary system and pysic for the secondary.
# The embedding scheme is set at mechanical embedding with hydrogen links.
#
#===============================================================================
from ase import Atoms, Atom
from gpaw import GPAW
from ase.io import write

for charge in [-1, 0 ,1]:

    # Prepare the system
    atom = Atoms('H', cell=[6.,6.,6.], pbc=False, charges=[charge])
    atom.center()

    # Initialize calculators
    gpaw_calc = GPAW(hund=True)
    atom.set_calculator(gpaw_calc)
    energy = atom.get_potential_energy()

    # Write wave functions to gpw file
    #gpaw_calc.write('H.gpw', mode='all')
    
    # Generate cube-files of the orbitals.
    wf = gpaw_calc.get_pseudo_wave_function(band=0)
    write('H%d.cube' % charge, atom, data=wf)


