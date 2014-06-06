#! /usr/bin/env python
from ase.io import write
from pysic.utility.error import *
import os
from multiprocessing import Process
from shutil import copyfile

#==============================================================================
class AtomEyeViewer(object):

    """Class for viewing atomistic configurations with AtomEye3:
    http://li.mit.edu/Archive/Graphics/A3/A3.html.
    
    This is an extra utility class, which relies on an existing installation of
    the AtomEye3 viewer. You have to add the AtomEye3 executable to your
    PYTHONPATH in order for this class to use it. By default this class will
    look for an executable named 'A' in that folder. If you wish to use a
    different name, use the 'atomeyecmd' parameter.

    You don't need to provide a directory if you wish to only view some
    configuration. If you wish to save any pictures or configuration files,
    provide the name of an existing folder as 'directory' parameter, or later
    on use the set_directory() method.
    
    Setup: download the AtomEye3 viewer, place the executable to a folder
    called /home/user/atomeye, add that appropriate folder to PYTHONPATH and
    PATH (on Linux add the lines::
    
        export PYTHONPATH=/path/to/AtomEye3Folder:$PYTHONPATH
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
        print "\n=============== AtomEyeViewer ===============\n" + message +'\n'

    def __init__(
            self,
            atoms=None,
            directory=None,
            initial_script=None,
            atomeyecmd='A'
            ):
        """Initializes the viewer.

        Parameters:
            atoms: `ASE Atoms`_ object
                The structure to visualize.
            directory: string
                The working directory, current working directory by default.
            initial_script: string
                Path to an AtomEye3 script file. Loaded when opening AtomEye3 window.
            atomeyecmd: string
                Name of the AtomEye3 executable, 'A' by default.
        """
        self.atomeyecmd = atomeyecmd
        self.atoms = atoms

        # If no directory is provided, use the current working directory
        if directory == None:
            self.set_directory(os.getcwd())
            self.warn("using directory: " + self.directory)
        else:
            self.set_directory(directory)

        self.startup_script = self.directory+"/startup_script"

        # Setup loading of initial script
        self.initial_script = initial_script
        if self.initial_script is None:
            setup = ("set n->bond_mode 1\n"
                     "set n->atom_r_ratio 0.75\n"
                     "resize 512 512\n"
                     "redraw\n"
                     )
            with open(self.directory+"/conf_script", "w") as myfile:
                myfile.write(setup)
        else:
            copyfile(initial_script, self.directory+"/conf_script")

        self.conf_script = self.directory+"/conf_script"
        self.animation_script = None
        self.jpg_script = None
        self.viewer_cfg = self.directory+"/view.cfg"
        self.frame_counter = 0
        self.folder_created = False

    def __del__(self):
        """Cleanup all temporary and unnecessary script and .cfg files."""
        if os.path.isfile(self.startup_script):
            os.remove(self.startup_script)
        if os.path.isfile(self.conf_script):
            os.remove(self.conf_script)
        if os.path.isfile(self.viewer_cfg):
            os.remove(self.viewer_cfg)

    def set_directory(self, directory):
        """Set the working directory.
        
        Parameters:
            directory: string
                Path of the working directory.
        """
        if not os.path.exists(directory):
                self.warn("Please provide an existing directory")
                return
        self.directory = directory

    def save_cfg(self, name=None):
        """Save the current structure as a .cfg file to the /cfgs folder under
        the working directory. 

        Parameter:
            name: string
                Name of the file. The .cfg postfix is automatically appended.
        """
        # Create cfg folder is necessary
        if self.cfg_dir is None:
            self.cfg_dir = self.directory + "/cfgs"
            if not os.path.exists(self.cfg_dir):
                os.makedirs(self.cfg_dir)

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
        # Create cfg folder is necessary
        if self.cfg_dir is None:
            self.cfg_dir = self.directory + "/cfgs"
            if not os.path.exists(self.cfg_dir):
                os.makedirs(self.cfg_dir)

        self.frame_counter += 1
        write(self.cfg_dir+'/'+str(self.frame_counter)+".cfg", self.atoms)

    def call_terminal(self, command):
        """Opens the terminal and issues the given command.
        
        Parameters:
            command: string
                The command to run.
        """
        os.system("gnome-terminal -e '"+command+"'")

    def view(self):
        """Opens the AtomEye3 window in a separate process through terminal and
        displays the structure.
        """
        # Add the conf script to startup
        with open(self.startup_script, 'w') as myfile:
            myfile.write("load_script "+self.conf_script+"\n")

        write(self.viewer_cfg, self.atoms)
        command = self.atomeyecmd+" "+self.directory+"/view.cfg -f="+self.startup_script

        # Open the viewer in another process
        p = Process(target=self.call_terminal, args=(command,))
        p.start()

    def view_series(self):
        """Used to view a series of .cfg files created with save_cfg_frame()
        method. Opens the first (1.cfg) frame in the /cfgs folder in a AtomEye3
        window. The next frame can be loaded with Delete key, and the previous
        with Insert-key.
        """
        # Check that the 1.cfg file exists
        if not os.path.exists(self.cfg_dir+"/1.cfg"):
            self.warn("The first frame called 1.cfg does not exist")
            return

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
        # Create jpg folder if necessary
        if self.cfg_dir is None:
            self.jpg_dir = self.directory + "/jpgs"
            if not os.path.exists(self.jpg_dir):
                os.makedirs(self.jpg_dir)

        # Check that the 1.cfg file exists
        if not os.path.exists(self.cfg_dir+"/1.cfg"):
            self.warn("The first frame called 1.cfg does not exist")
            return

        # Write the scr_anim file
        self.animation_script = self.cfg_dir+"/scr_anim"
        with open(self.animation_script, 'w') as myfile:
            myfile.write(str(quality)+'\n')
            for i in range(self.frame_counter):
                myfile.write(self.cfg_dir+'/'+str(i+1)+".cfg "+self.jpg_dir+'/'+str(i+1)+".jpg\n")

        # Write the jpg script
        self.jpg_script = self.cfg_dir+"/jpg_scr"
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


