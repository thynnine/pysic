#! /usr/bin/env python
"""Used to make an ASE Atoms solution system from the given solvent and
solute.
"""

import numpy as np


def make_solution(solute, solvent, grid, lattice_constants, safe_distance, padding=5):
    """Surrounds the given solute with solvent.

    Surrounds the given solute with a grid of solvent molecules. A periodic
    grid of solvent molecules are created and the solute is placed in the
    center of this grid. The solute molecules that are too close to the solute
    are then removed.

    Parameters:
        solute: ASE Atoms object
            The solute system.
        solvent: ASE Atoms object
            The solvent molecules.
        grid: tuple
            Contains the amount of solvent atoms in each of the three
            dimensions.
        lattice_constants: tuple
            Contains the lattice constants for each of the three dimensions.
        safe_distance: float
            The minimum distance that should be between the solute and solvent
            atoms. Solvent molecule is not placed in a certain grid site if the
            distance is smaller than this.
        padding: float
            The distance left between the solution and the cell boundary.

    Returns: tuple
        The solution as an ASE Atoms object and the number of solvent molecules created
    """
    a = lattice_constants[0]
    b = lattice_constants[1]
    c = lattice_constants[2]

    cell = ((grid[0]-1)*a+2*padding, (grid[1]-1)*b+2*padding, (grid[2]-1)*c+2*padding)
    solute_copy = solute.copy()
    solute_copy.set_cell(cell)
    solute_copy.center()

    solution = solute_copy.copy()

    n_solvents = 0
    for x in range(grid[0]):
        for y in range(grid[1]):
            for z in range(grid[2]):
                position = np.array((x*a+padding, y*b+padding, z*c+padding), dtype=float)
                i_solvent = solvent.copy()
                cm = i_solvent.get_center_of_mass()
                translation = position-cm
                i_solvent.translate(translation)
                safe = True
                for solvent_atom in i_solvent:
                    if not safe:
                        break
                    for solute_atom in solute_copy:
                        solvent_r = solvent_atom.position
                        solute_r = solute_atom.position
                        distance = np.linalg.norm(solvent_r-solute_r)
                        if distance < safe_distance:
                            safe = False
                            break
                if safe:
                    solution += i_solvent
                    n_solvents += 1

    return (solution, n_solvents)
