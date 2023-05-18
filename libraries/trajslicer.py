"""
This module is for slicing trajectory files
and outputting it to a new trajectory file.

Made by Ramon Mendoza, ramendoza@uchicago.edu
Last update Februrary 2, 2022
"""

import MDAnalysis as mda
import numpy as np


class TrajSlicer:
    """
    class trajslicer.TrajSlicer(trajectory_file, psf_file, output_file)
        This module is for slicing trajectory files into smaller trajectories
        and outputting it to a new trajectory file.
        __init__(self, trajectory_file, psf_file, output_file)
        Parameters: trajectory_file (str) - DCD trajectory file name
                    psf_file (str) - Protein Structure File (PSF) file name
                    output_file (str) - DCD trajectory output file name
        Attributes: trajectory_file - DCD trajectory file name
                    psf_file - Protein Structure File (PSF) file name
                    output_file - DCD trajectory output file name
                    trajectory - An MDAnalysis Universe object
                    last_frame - Last frame of the trajectory file
                    total_atoms - Total atoms in the trajectory file
        Methods:
            load_trajectory(self)
                Makes an MDAnalysis Universe object; invoked upon initializaion
                Returns: An MDAnalysis.Universe object
            traj_slicer(self, first=0, last=last_frame, stride=1)
                Slices the trajectory file from first to last and writes a new trajectory file
                Parameters: first (int) - First frame to begin slicing; inclusive
                            last (int) - Last frame to end slicing; inclusive
                            stride (int) - Take every n-th frame from first to last
                Returns: A new DCD trajectory file
    """
    
    def __init__(self, trajectory_file, psf_file, output_file):
        self.trajectory_file = str(trajectory_file)
        self.psf_file = str(psf_file)
        self.output_file = str(output_file)
        self.trajectory = self.load_trajectory()
        txt = str(self.trajectory.trajectory).split()
        self.last_frame = int(txt[3])
        self.total_atoms = int(txt[6])

    def load_trajectory(self):
        """Returns an MDAnalysis.Universe object
        """
        traj = mda.Universe(self.psf_file, self.trajectory_file)
        return traj

    def traj_slicer(self, first=0, last=0, stride=1):
        """Returns a new DCD trajectory file from first to last taken stride frames
        """
        if last == 0:
            last = self.last_frame
        else:
            last = last
        print('Outputting new trajectory to: ' 
              + '{}\nStart at frame: {}\nFinish at frame: {}'.format(self.output_file, first, 
                                                                     last))
        with mda.Writer(self.output_file, self.total_atoms) as w:
            i = 0
            for ts in self.trajectory.trajectory:
                txt = str(ts).split()
                current_frame = int(txt[2])
                i %= stride
                if not i and first <= current_frame <= last:
                    w.write(self.trajectory.atoms)
                i += 1
                
                