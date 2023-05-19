"""
This module is for adding common and not so common restraints to OpenMM simulations using 
contants in 'NAMD'-like units. Some examples include dihedral, COM, distance, and Positional 
restraints.

Made by Ramon Mendoza, ramendoza@uchicago.edu
Last update May 19, 2023
"""

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

import MDAnalysis as mda
import numpy as np

class OpenMMRestratints:
    """A container of common and not so common restraints and biases that can be used in 
    an OpenMM MD simulation.
    """
    
    def __int__(self, system):
        """Initialize the OpenMMRestraints module for the current simulation system
        Paramters: system - the OpenMM system object
        """
        self.system = system
    
    def dihedralRestraints(self, force_constant, dihedral_file):
        """Returns an OpenMM CustomTorsionForce with a flat-bottom restraint to dihedrals.
        
        Parameters: force_constant - force constant in units of kcal/mol/degree^2
                    dihedral_file - text file containing 6 columns; columns 1-4 (zero-based)
                                    index of atoms involved in the dihedral (order matters!), 
                                    column 5-6 the position of the lower and upper walls
        """
        k = force_constant * kilocalorie_per_mole/degree**2 # force constant in units of 
                                                            # kJ/mol/rad^2
        
        flatbottomRestraint = CustomTorsionForce("select(a, a, b) * (k/2)" \
                                                 "*(theta-select(a, theta0, theta1))^2;" \
                                                 "a=step(theta-theta0); b=step(theta1-theta)")
        
        flatbottomRestraint.setName('DihedralRestraint')
        flatbottomRestraint.addGlobalParameter('k', k)
        flatbottomRestraint.addPerTorsionParameter('theta0')
        flatbottomRestraint.addPerTorsionParameter('theta1')
        
        atomIndices = np.genfromtxt(dihedral_file, usecols=(0, 1, 2, 3))
        wallPositions = np.genfromtxt(dihedral_file, usecols=(4, 5))
        
        # need to check is sorting is working as expected
        wallPositions = np.sort(wallPositions) # don't trust that user will be consistent as to
                                               # which number from columns 5 and 6 will be the
                                               # lower or upper wall position!
        
        # Need to check if zip is doing what I think it will be doing!
        for atomSet, wallPosition in zip(atomIndices, wallPositions):
            flatbottomRestraint.addTorsion(atomSet, [wallPosition])
            
        self.system.addForce(flatbottomRestraint) 
    
    