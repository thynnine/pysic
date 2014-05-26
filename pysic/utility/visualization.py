#! /usr/bin/env python
from ase.io import write
from pysic.utility.error import *
import os
from shutil import copyfile

#==============================================================================
class AtomEyeViewer(object):

    """Class for viewing atomistic configurations with AtomEye3:
    http://li.mit.edu/Archive/Graphics/A3/A3.html.
    
    This is an extra utility class, which relies on an existing installation of
    the AtomEye3 viewer. You have to add the AtomEye3 executable to your
    PYTHONPATH in order for this class to use it. You will have to provide the
    name of this executable as an argument.
    
    Setup: download the AtomEye3 viewer, place the executable to a folder
    called /home/user/atomeye, add that folder to PYTHONPATH (on Linux add the
    line: export PYTHONPATH=$HOME/Dropbox/SIN/AtomEye3:$PYTHONPATH to .bashrc
    or .profile), provide the name of the executable as argument to this class.
    """

    def __init__(self, atomeyecmd, atoms=None, directory=None, initial_script=None):
        """@todo: to be defined1.

        :atoms: @todo

        """
        self.atomeyecmd = atomeyecmd
        self.atoms = atoms
        self.counter = 0
        self.directory = directory
        self.startup_script = self.directory+"/startup_script"
        if directory is not None:
            self.set_directory(directory)

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

    def set_directory(self, directory):
        """ """
        if not os.path.exists(directory):
                warn("Please provide an existing directory", 2)
                return
        self.cfg_dir = self.directory + "/cfgs"
        self.jpg_dir = self.directory + "/jpgs"
        if not os.path.exists(self.cfg_dir):
            os.makedirs(self.cfg_dir)
        if not os.path.exists(self.jpg_dir):
            os.makedirs(self.jpg_dir)

    def save_cfg(self, name=None):
        """@todo: Docstring for function.

        :arg1: @todo
        :returns: @todo

        """
        self.counter += 1
        write(self.cfg_dir+'/'+str(self.counter)+".cfg", self.atoms)

    def view(self):
        """@todo: docstring for create_jpgs.
        :returns: @todo

        """
        # add the conf script to startup
        with open(self.startup_script, 'w') as myfile:
            myfile.write("load_script "+self.conf_script+"\n")

        write(self.directory+"/temp.cfg", self.atoms)
        command = self.atomeyecmd+" "+self.directory+"/temp.cfg -f="+self.startup_script
        os.system("gnome-terminal -e '"+command+"'")
        os.remove(self.directory+"/temp.cfg")
        
    def view_series(self):
        """@todo: docstring for create_jpgs.
        :returns: @todo

        """
        # add the conf script to startup
        with open(self.startup_script, 'w') as myfile:
            myfile.write("load_script "+self.conf_script+"\n")

        command = self.atomeyecmd+" "+self.cfg_dir+"/1.cfg -f="+self.startup_script
        os.system("gnome-terminal -e '"+command+"'")

    def create_jpgs(self, quality=90, resolution=[256, 256]):
        """@todo: Docstring for function.

        :arg1: @todo
        :returns: @todo

        """
        # Write the scr_anim file
        self.animation_script = self.cfg_dir+"/scr_anim"
        with open(self.animation_script, 'w') as myfile:
            myfile.write(str(quality)+'\n')
            for i in range(self.counter):
                myfile.write(self.cfg_dir+'/'+str(i+1)+".cfg "+self.jpg_dir+'/'+str(i+1)+".jpg\n")

        # Write the jpg script
        self.jpg_script = self.cfg_dir+"/jpg_scr"
        with open(self.jpg_script, 'w') as myfile:
            myfile.write("script_animate "+self.animation_script+"\n")
            myfile.write("quit\n")

        # Add the jpg script to startup
        with open(self.startup_script, 'w') as myfile:
            myfile.write("load_script "+self.conf_script+"\n")
            myfile.write("load_script "+self.jpg_script+"\n")

        # Run atomeye to create the jpegs
        command = self.atomeyecmd+" "+self.cfg_dir+'/'+str(self.counter)+".cfg -f="+self.startup_script
        os.system("gnome-terminal -e '"+command+"'")

