#! /usr/bin/env python
"""Provides a function for calculating the Bader charges for a given atomic
system with the given calculator."""

from distutils import spawn
from pysic.utility.error import error
import numpy as np
import os
from ase.parallel import rank, barrier
from ase.units import Bohr
from ase.io import write, bader
import subprocess
import shutil


def get_bader_charges(atoms, calc, charge_source="all-electron", gridrefinement=4):
        """This function uses an external Bader charge calculator from
        http://theory.cm.utexas.edu/henkelman/code/bader/. This tool is
        provided also in pysic/tools. Before using this function the bader
        executable directory has to be added to PATH.

        Parameters:
            atoms: ASE Atoms
                The structure from which we want to calculate the charges from.
            calc: ASE calculator
            charge_source: string
                Indicates the electron density that is used in charge calculation.
                Can be "pseudo" or "all-electron".
            gridrefinement: int
                The factor by which the calculation grid is densified in charge
                calculation.

        Returns: numpy array of the atomic charges
        """
        # First check that the bader executable is in PATH
        if spawn.find_executable("bader") is None:
            error((
                "Cannot find the \"bader\" executable in PATH. The bader "
                "executable is provided in the pysic/tools folder, or it can be "
                "downloaded from http://theory.cm.utexas.edu/henkelman/code/bader/. "
                "Ensure that the executable is named \"bader\", place it in any "
                "directory you want and then add that directory to your system"
                "PATH."))

        atoms_copy = atoms.copy()
        calc.set_atoms(atoms_copy)

        if charge_source == "pseudo":
            try:
                density = np.array(calc.get_pseudo_density())
            except AttributeError:
                error("The calculator doesn't provide pseudo density.")

        if charge_source == "all-electron":
            try:
                density = np.array(calc.get_all_electron_density(gridrefinement=gridrefinement))
            except AttributeError:
                error("The calculator doesn't provide all electron density.")

        wrk_dir = os.getcwd()+"/.BADERTEMP"
        dir_created = False

        # Write the density in bader supported units and format
        if rank == 0:

            # Create temporary folder for calculations
            if not os.path.exists(wrk_dir):
                os.makedirs(wrk_dir)
                dir_created = True
            else:
                error("Tried to create a temporary folder in " + wrk_dir + ", but the folder already existed. Please remove it manually first.")

            rho = density * Bohr**3
            write(wrk_dir + '/electron_density.cube', atoms, data=rho)

            # Run the bader executable in terminal. The bader executable included
            # int pysic/tools has to be in the PATH/PYTHONPATH
            command = "cd " + wrk_dir + "; bader electron_density.cube"
            subprocess.check_output(command, shell=True)
            #os.system("gnome-terminal --disable-factory -e '"+command+"'")

        # Wait for the main process to write the file
        barrier()

        # ASE provides an existing function for attaching the charges to the
        # atoms (safe because using a copy). Although we don't want to actually
        # attach the charges to anything, we use this function and extract the
        # charges later.
        bader.attach_charges(atoms_copy, wrk_dir + "/ACF.dat")

        # The call for charges was changed between
        # ASE 3.6 and 3.7
        try:
            bader_charges = np.array(atoms_copy.get_initial_charges())
        except:
            bader_charges = np.array(atoms_copy.get_charges())

        # Remove the temporary files
        if rank == 0:
            if dir_created:
                shutil.rmtree(wrk_dir)

        return bader_charges
