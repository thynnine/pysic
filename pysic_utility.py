"""Contains necessary and supporting auxiliary functions for Pysic."""

try:
    import matplotlib.pyplot as plt
except:
    print "error importing matplotlib - the plotting tools are not available"
    
from math import sqrt, pi, sin, cos, exp
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
