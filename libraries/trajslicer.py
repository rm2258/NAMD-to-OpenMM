"""
This module is for slicing trajectory files
and outputting it to a new trajectory file.

Made by Ramon Mendoza, ramendoza@uchicago.edu
Last update Februrary 2, 2022
"""

import MDAnalysis as mda
import numpy as np


class TrajSlicer:

    def __init__(self, trajectory_file, psf_file, output_file):
        self.trajectory_file = str(trajectory_file)
        self.psf_file = str(psf_file)
        self.output_file = str(output_file)
        self.trajectory = self.load_trajectory()
        txt = str(self.trajectory.trajectory).split()
        self.last_frame = int(txt[3])
        self.total_atoms = int(txt[6])
        print('Use the help() method for class use')

    @staticmethod
    def help():
        help_message = 'class trajslicer.TrajSlicer(trajectory_file, psf_file, output_file)\n' + \
            '    This module is for slicing trajectory files into smaller trajectories\n' + \
            '    and outputting it to a new trajectory file.\n\n' + \
            '    __init__(self, trajectory_file, psf_file, output_file)\n\n' + \
            '        Parameters:\n' + \
            '            - trajectory_file (str) - DCD trajectory file name\n' + \
            '            - psf_file (str) - Protein Structure File (PSF) file name\n' + \
            '            - output_file (str) - DCD trajectory output file name\n\n' + \
            '        Attributes\n' + \
            '            - trajectory_file - DCD trajectory file name\n' + \
            '            - psf_file - Protein Structure File (PSF) file name\n' + \
            '            - output_file - DCD trajectory output file name\n' + \
            '            - trajectory - An MDAnalysis Universe object\n' + \
            '            - last_frame - Last frame of the trajectory file\n' + \
            '            - total_atoms - Total atoms in the trajectory file\n\n' + \
            '    Methods\n' + \
            '        help()\n\n' + \
            '            Print a help message\n\n' + \
            '            Returns:\n\n' + \
            '                - A help message\n\n' + \
            '        load_trajectory(self)\n\n' + \
            '            Makes an **MDAnalysis** Universe object; invoked upon initializaion\n\n' + \
            '            Returns:\n\n' + \
            '                - An MDAnalysis.Universe object\n\n' + \
            '        traj_slicer(self, first=0, last=last_frame, stride=1)\n\n' + \
            '            Slices the trajectory file from first to last and writes a new trajectory file\n\n' + \
            '            Parameters:\n\n' + \
            '                - first (int) - First frame to begin slicing; inclusive\n' + \
            '                - last (int) - Last frame to end slicing; inclusive\n' + \
            '                - stride (int) - Take every n-th frame from first to last\n\n' + \
            '            Returns:\n\n' + \
            '                - A new DCD trajectory file\n\n'
        print(help_message)

    def load_trajectory(self):
        traj = mda.Universe(self.psf_file, self.trajectory_file)
        return traj

    def traj_slicer(self, first=0, last=0, stride=1):
        if last == 0:
            last = self.last_frame
        else:
            last = last
        print('Outputting new trajectory to: {}\nStart at frame: {}\nFinish at frame: {}'.format(self.output_file, first, last))
        with mda.Writer(self.output_file, self.total_atoms) as w:
            i = 0
            for ts in self.trajectory.trajectory:
                txt = str(ts).split()
                current_frame = int(txt[2])
                i %= stride
                if not i and first <= current_frame <= last:
                    w.write(self.trajectory.atoms)
                i += 1