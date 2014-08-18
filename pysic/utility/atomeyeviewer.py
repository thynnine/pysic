#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Used for visualizing ASE Atoms with AtomEye3."""

from ase.io import write
from pysic.utility.error import style_message
import os
from multiprocessing import Process
import shutil
from ase.parallel import rank
from distutils import spawn


#==============================================================================
class AtomEyeViewer(object):

    """Class for viewing ASE Atoms with AtomEye3:
    http://li.mit.edu/Archive/Graphics/A3/A3.html.

    The viewer is MPI compatible, and will only save files or display
    visualizations on the rank 0 process.

    This is an extra utility class, which relies on an existing installation of
    the AtomEye3 viewer. You have to add the AtomEye3 executable to your and
    PATH in order for this class to use it. By default this class will look for
    an executable named 'atomeye3' in that folder. If you wish to use a
    different name, use the 'atomeyecmd' parameter.

    You don't need to provide a wrk_dir if you wish to only view some
    configuration. If you wish to save any pictures or configuration files,
    provide the name of an existing folder as 'wrk_dir' parameter, or later
    on use the set_wrk_dir() method.

    Setup: download the AtomEye3 viewer, place the executable to a folder add
    that appropriate folder to and PATH (on Linux add the lines::

        export PATH=/path/to/AtomEye3Folder:$PATH

    to .bashrc or .profile), and possible provide the name of the executable as
    argument to this class.
    """
    def warn(self, message):
        """Prints a warning message to console.

        Parameters:
            message: string
                The message to print.
        """
        print style_message("ATOMEYE3VIEWER", message)

    def __init__(
            self,
            atoms=None,
            wrk_dir=None,
            subdirectory="atomeyeviewer",
            atomeyecmd='atomeye3'):
        """Initializes the viewer.

        Parameters:
            atoms: `ASE Atoms`_ object
                The structure to visualize.
            wrk_dir: string
                The working wrk_dir, current working wrk_dir by default.
            initial_script: string
                Path to an AtomEye3 script file. Loaded when opening AtomEye3 window.
            atomeyecmd: string
                Name of the AtomEye3 executable, 'atomeye3' by default.
        """
        self.atomeyecmd = atomeyecmd
        self.atoms = atoms
        self.new_folder_created = None
        self.frame_counter = 0
        self.colors = None
        self.radii = None

        # First check that the atomeye3 executable is in PATH
        if spawn.find_executable("atomeye3") is None:
            self.warn((
                "Cannot find the \"atomeye3\" executable in PATH. This "
                "executable is provided in the pysic/tools folder, or it can "
                "be downloaded from "
                "http://li.mit.edu/Archive/Graphics/A3/A3.html#download. "
                "Ensure that the executable is named appropriately (by default "
                "the name is \"atomeye3\", but you can specify any name in the "
                "constructor), place it in any directory you want and then add "
                "that directory to your system" "PATH."))

        # If no wrk_dir is provided, use the current working directory
        if wrk_dir is None:
            self.set_subdirectory(os.getcwd(), ".ATOMEYEVIEWERTEMP")
            self.explicit_dir_given = False
            self.warn(("Working directory not set. Working in temporary directory" + self.subdir + "\n" +
                       "In order to save any files, please provide a working directory."))
        else:
            self.set_subdirectory(wrk_dir, subdirectory)
            self.explicit_dir_given = True

        # Setup the initial conf script
        self.conf_lines = {
            "set n->bond_mode": [1],
            "set n->atom_r_ratio": [0.75],
            "redraw": [],
            "resize": [512, 512]}

    def __del__(self):
        """If a directory was not provided and a new folder was created, then
        delete it.
        """
        if rank == 0:
            if not self.explicit_dir_given and self.folder_created:
                if os.path.isdir(self.subdir):
                    shutil.rmtree(self.subdir)

    def set_atoms(self, atoms):
        """ Set the ASE Atoms that are visualized.
        """
        self.atoms = atoms

    def set_subdirectory(self, wrk_dir, subdirectory):
        """Set the working wrk_dir.

        Parameters:
            wrk_dir: string
                Path of the working wrk_dir.
        """
        if rank == 0:
            if not os.path.exists(wrk_dir):
                    self.warn("Please provide an existing wrk_dir")
                    return
            self.wrk_dir = wrk_dir
            self.subdir = wrk_dir+'/'+subdirectory
            self.cfg_dir = self.subdir + "/cfgs"
            self.jpg_dir = self.subdir + "/jpgs"
            self.animation_script = self.cfg_dir+"/scr_anim"
            self.jpg_script = self.cfg_dir+"/jpg_scr"
            self.startup_script = self.subdir+"/startup_script"
            self.conf_script = self.subdir+"/conf_script"
            self.viewer_cfg = self.subdir+"/view.cfg"

            # Create subfolder if necessary
            if not os.path.exists(self.subdir):
                os.makedirs(self.subdir)
                self.folder_created = True
            else:
                self.folder_created = False

            # Create cfg folder if necessary
            if not os.path.exists(self.cfg_dir):
                os.makedirs(self.cfg_dir)

            # Create jpg folder if necessary
            if not os.path.exists(self.jpg_dir):
                os.makedirs(self.jpg_dir)

            self.temporary_session = False

    def save_cfg(self, name):
        """Save the current structure as a .cfg file to the /cfgs folder under
        the working wrk_dir.

        Parameter:
            name: string
                Name of the file. The .cfg postfix is automatically appended.
        """
        if rank == 0:
            write(self.cfg_dir+'/'+name+".cfg", self.atoms)

    def save_cfg_frame(self):
        """A convenience function for saving structure as .cfg file during
        molecular dynamics or structure optimization. An internal counter is
        used for keeping count of the frame number. Example usage::

        atoms = Atoms('H2O')
        dyn = VelocityVerlet(atoms, 3*units.fs)
        visuals = AtomEyeViewer(atoms, "/folder")
        dyn.attach(visuals.save_cfg_frame, interval=1)
        dyn.run(10)

        """
        if rank == 0:
            self.frame_counter += 1
            write(self.cfg_dir+'/'+str(self.frame_counter)+".cfg", self.atoms)

    def call_terminal(self, command):
        """Opens the terminal and issues the given command.

        Parameters:
            command: string
                The command to run.
        """
        os.system("gnome-terminal -e '"+command+"'")

    def view(self, width=512, height=512):
        """Opens the AtomEye3 window in a separate process through terminal and
        displays the structure.
        """

        self.set_size(width, height)

        if rank == 0:
            # Add the conf script to startup
            with open(self.startup_script, 'w') as myfile:
                myfile.write("load_script "+self.conf_script+"\n")

            self.write_usr_file(self.subdir, "view.usr")
            self.write_conf_script()
            write(self.viewer_cfg, self.atoms)

            command = self.atomeyecmd+" "+self.viewer_cfg+" -f="+self.startup_script

            # Open the viewer in another process
            p = Process(target=self.call_terminal, args=(command,))
            p.start()

    def set_size(self, width, height):
        self.conf_lines["resize"] = [width, height]

    def view_series(self):
        """Used to view a series of .cfg files created with save_cfg_frame()
        method. Opens the first (1.cfg) frame in the /cfgs folder in a AtomEye3
        window. The next frame can be loaded with Delete key, and the previous
        with Insert-key.
        """
        if rank == 0:
            # Check that the 1.cfg file exists
            if not os.path.exists(self.cfg_dir+"/1.cfg"):
                self.warn("The first frame called 1.cfg does not exist")
                return

            # Write the conf script
            self.write_usr_file(self.subdir, "view.usr")
            self.write_conf_script()

            # add the conf script to startup
            with open(self.startup_script, 'w') as myfile:
                myfile.write("load_script "+self.conf_script+"\n")

            command = self.atomeyecmd+" "+self.cfg_dir+"/1.cfg -f="+self.startup_script

            # Open the viewer in another process
            p = Process(target=self.call_terminal, args=(command,))
            p.start()

    def save_jpg_series(self, quality=90, resolution=(256, 256)):
        """A convenience function for creating a series of .jpg files in the
        /jpgs folder. Creates jpg images of all the of .cfg files created with
        the function save_cfg_frame().

        Parameters:
            quality: int
                The image quality, between 0-100, 90 by default.
            resolution: tuple
                The image resolution. 256x256 by default.
        """
        if rank == 0:
            # Check that the 1.cfg file exists
            if not os.path.exists(self.cfg_dir+"/1.cfg"):
                self.warn("The first frame called 1.cfg does not exist")
                return

            # Write the scr_anim file
            with open(self.animation_script, 'w') as myfile:
                myfile.write(str(quality)+'\n')
                for i in range(self.frame_counter):
                    myfile.write(self.cfg_dir+'/'+str(i+1)+".cfg "+self.jpg_dir+'/'+str(i+1)+".jpg\n")

            # Write the jpg script
            with open(self.jpg_script, 'w') as myfile:
                myfile.write("resize "+str(resolution[0])+' '+str(resolution[1])+'\n')
                myfile.write("script_animate "+self.animation_script+'\n')
                myfile.write("quit\n")

            # Add the jpg script to startup
            with open(self.startup_script, 'w') as myfile:
                myfile.write("load_script "+self.conf_script+"\n")
                myfile.write("load_script "+self.jpg_script+"\n")

            # Run atomeye to create the jpegs
            command = self.atomeyecmd+" "+self.cfg_dir+'/'+str(self.frame_counter)+".cfg -f="+self.startup_script
            self.call_terminal(command)

    #---------------------------------------------------------------------------
    # Visualization settings

    def set_colors(self, colors):
        """Sets the colors of the atoms.

        Parameters:
            colors: list of tuples
                A list containing a tuple of three floats for each atom. The
                color tuple specifies the RGB color of the atom. The RGB values
                are specified in the interval 0.0-1.0. Example (0, 1, 0) =
                Green.
        """
        if not (len(self.atoms) == len(colors)):
            self.warn("The color list length does not match the number of atoms", 3)
            return
        self.colors = colors

    def set_radii(self, radii):
        """Sets the radii of the atoms.

        Parameters:
            radii: list of floats
                A list of atomic radius in units of Å
        """
        if not (len(self.atoms) == len(radii)):
            self.warn("The radii list length does not match the number of atoms", 3)
            return
        self.radii = radii

    def write_usr_file(self, directory, name):
        if rank == 0:
            """Writes the .usr file containing extra visualization parameters.
            """
            # Write the .usr-file if some visualization properties have been set
            if self.colors is not None or self.radii is not None:
                with open(directory+'/'+name, 'w') as myfile:
                    for i_atom, atom in enumerate(self.atoms):
                        if self.colors is not None:
                            for i_channel, channel in enumerate(self.colors[i_atom]):
                                myfile.write(str(channel))
                                if i_channel is not 2:
                                    myfile.write(' ')
                            if self.radii is not None:
                                myfile.write(' ')
                                myfile.write(str(self.radii[i_atom]))
                        myfile.write('\n')

                # Add load line to conf_script
                self.conf_lines["load_atom_color"] = [directory+'/'+name]

    def write_conf_script(self):
        if rank == 0:
            """Writes the conf script according to the inforamtion in the
            dictionary self.conf_lines.
            """
            with open(self.conf_script, 'w') as myfile:
                for key, value in self.conf_lines.iteritems():
                    myfile.write(key)
                    for item in value:
                        myfile.write(' ')
                        myfile.write(str(item))
                    myfile.write('\n')
