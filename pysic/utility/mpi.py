#! /usr/bin/env python
"""Defines MPI safe routines for printing, file access etc.
    """
import pysic.pysic_fortran as pf
import os


def finish_mpi():
    """Terminates the MPI framework.
        
        If the Fortran core is compiled in MPI mode, :mod:`~pysic` will automatically
        initialize MPI upon being imported. This method terminates the MPI.
        """
    pf.pysic_interface.finish_mpi()


def sync_mpi():
    """Calls MPI barrier from the Fortran core MPI framework."""
    
    pf.pysic_interface.sync_mpi()


def get_number_of_cpus():
    """Gets the number of cpus from the Fortran MPI.
        """
    return pf.pysic_interface.get_number_of_cpus()


def get_cpu_id():
    """Gets the cpu ID from the Fortran MPI.
        """
    return pf.pysic_interface.get_cpu_id()


def cpu_id():
    """Gets the cpu ID from the Fortran MPI.

        Equivalent to :func:`~pysic.utility.mpi.get_cpu_id`.
        """
    return get_cpu_id()


def is_master():
    """Returns ``True`` for the cpu with ID 0.
        """
    return cpu_id() == 0


def mkdir(dir):
    """Creates a new directory only with the master cpu.
        
    Parameters:
        
    dir: string
        the name of the directory
        """
    if is_master():
        try:
            os.mkdir(dir)
        except:
            mprint("Failed to create directory "+dir)


def cd(dir):
    """Changes to the given directory on all cpus.
        
        Parameters:
        
        dir: string
            the name of the directory
        """
    try:
        os.chdir(dir)
    except:
        mprint("Failed to change to directory "+dir)


def mprint(string):
    """Prints the string to stdout only from the master cpu.
        
    Parameters:
     
    string: string
        the string to be written        
        """
    
    if is_master():
        print string


def write_file(lines, filename):
    """Writes the given lines to a text file so that only the master cpu writes.
        
    Parameters:
        
    lines: list of strings
        lines to be written
        
    filename: string
        the name of the file
    """
    if is_master():
        f = file(filename,'w')
        string = ""
        for line in lines:
            string += line+"\n"
        f.write(string[:-2])
        f.close()


def append_file(lines, filename):
    """Appends the given lines to a text file so that only the master cpu writes.

    Parameters:
        
    lines: list of strings
        lines to be written
    
    filename: string
        the name of the file
    """
    if is_master():
        f = file(filename,'a')
        string = ""
        for line in lines:
            string += line+"\n"
        f.write(string[:-2])
        f.close()


def mpi_barrier():
    """Calls MPI barrier from the Fortran core MPI framework.
    
    This is equivalent to :func:`~pysic.utility.mpi.sync_mpi`
    """
    pf.pysic_interface.sync_mpi()




