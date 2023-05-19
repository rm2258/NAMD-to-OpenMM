"""
This module is for reporting energy output in a NAMD fashion during the course of a MD 
simulation

Made by Ramon Mendoza, ramendoza@uchicago.edu
Last update May 18, 2023
"""

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import mdtraj as md
import numpy as np

class EnergyReproter(object):
    """A custom 'NAMD'-like energy reporter for OpenMM
    """
    def __init__(self, file, reportInterval, output=False):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval
        self._output = output
        
    def __del__(self):
        self._out.close()
    
    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, False, True, False)
    
    def report(self, simulation, state):
        system = simulation.context.getSystem()
        step = simulation.currentStep
        self._out.write('Step: {}\n'.format(step))
        print("Step: ", step) if self._output else 0
        for i, f in enumerate(system.getForces()):
            state = simulation.context.getState(getEnergy=True, groups={i})
            print(f.getName() + ":\t", 
                  state.getPotentialEnergy().in_units_of(kilocalories_per_mole)) if self._output else 0
            self._out.write(f.getName() + ":\t" 
                            + str(state.getPotentialEnergy() \
                                  .in_units_of(kilocalories_per_mole))
                            + '\n')
        print("Potential energy: " 
              + str(simulation.context.getState(getEnergy=True) \
                    .getPotentialEnergy().in_units_of(kilocalories_per_mole)), 
              '\n') if self._output else 0
        self._out.write("Potential energy: " 
                        + str(simulation.context.getState(getEnergy=True) \
                              .getPotentialEnergy().in_units_of(kilocalories_per_mole))
                        + '\n\n')
        
        