"""
This module is for computing 'NAMD'-like energies in OpenMM

Made by Ramon Mendoza, ramendoza@uchicago.edu
Last update Februrary 2, 2022
"""

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import mdtraj as md
import numpy as np


class NamdToOpenmmTools:
    """
    class namdtoopenmm.NamdToOpenmmTools(trajectory_file, psf_file, parameter_files, 
                                         temperature, pbc)
        This module is for computing 'NAMD'-like energies in OpenMM
        __init__(self, trajectory_file, psf_file, parameter_files, temperature, pbc)
        Parameters: trajectory_file (str) - File name of DCD trajectory file
                    psf_file (str) - File name of PSF file
                    parameter_files (list) - list of parameter files to load
                    temperature (float) - simulation temperature in kelvin\
                    pbc (tuple) - a, b, c for periodic boundary conditions in angstroms
        Attributes: trajectory_file - Trajectory file to load
                    psf_file - Protein Structure File (psf) to load
                    parameter_files - A list of parameter files to load
                    temperature - The simulation temperature in kelvin
                    pbc - a, b, c for periodic boundary conditions in angstroms
                    time_conversions - Time conversion constants
                    trajectory - Loaded trajectory file as an mdtraj.Trajectory object
                    system - An OpenMM system object
                    simulation - Instantiate an OpenMM simulation object
        Methods
            load_trajectory(self)
                Load a trajectory file
                Returns: An mdtraj.Trajectory object
                
            load_system(self, parameter_files, temperature, pbc)
                Makes an OpenMM simulation object
                Parameters: parameter_files (list) - list of parameter files to load
                            temperature (float) - simulation temperature in kelvin
                            pbc (tuple) - a, b, c for periodic boundary conditions in angstroms
                Returns: system, simulation
                         system - An OpenMM system object
                         simulation - An OpenMM simulation object
            
            compute_energies(self, dcdfreq, time='femtoseconds', print_output=False)
                Computes 'NAMD'-like potential energies in OpenMM
                    Parameters: dcdfreq (int) - DCD frequency output
                                time (str) - Desired unit of time
                                print_output (bool) - Print energy breakdown output
                    Returns: potential_energies, output
                             potential_energies (numpy array) - Potential energies
                             output (list) - List of strings of energy breakdown output
    """

    def __init__(self, trajectory_file, psf_file, parameter_files, temperature, pbc):
        self.trajectory_file = str(trajectory_file)
        self.psf_file = str(psf_file)
        self.parameter_files = list(parameter_files)
        self.temperature = float(temperature)
        self.pbc = tuple(pbc)
        self.time_conversions = {'femtoseconds': 2, 'picoseconds': 0.002, 'nanoseconds': 2 
                                 * 10 ** (-6), 'microseconds': 2 * 10 ** (-9)}
        self.trajectory = self.load_trajectory()
        self.system, self.simulation = self.load_system(self.parameter_files, self.temperature, 
                                                        self.pbc)


    def load_trajectory(self):
        """Returns an mdtraj.Trajectory object from a trajectory file
        """
        traj = md.load_dcd(self.trajectory_file, self.psf_file)
        return traj

    def load_system(self, parameter_files, temperature, pbc):
        """Returns the OpenMM system and simulation object
        """
        psf = CharmmPsfFile(self.psf_file)
        params = CharmmParameterSet(*parameter_files)
        psf.setBox(pbc[0] * angstroms, pbc[1] * angstroms, pbc[2] * angstroms)
        system = psf.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=1.2 
                                  * nanometers, switchDistance=1.0 * nanometers,
                                  constraints=HBonds, ewaldErrorTolerance=0.0005)
        pressure_controller = MonteCarloMembraneBarostat
        system.addForce(pressure_controller(1.01325 * bar, 0.0 * bar * nanometers, temperature 
                                            * kelvin, pressure_controller.XYIsotropic,
                                            pressure_controller.ZFree, 100))
        for i, f in enumerate(system.getForces()):
            f.setForceGroup(i)
        integrator = LangevinIntegrator(temperature * kelvin, 1 / picosecond, 0.002 
                                        * picoseconds)
        _cudaPlatform = False
        _numPlatforms = Platform.getNumPlatforms()
        for i in range(_numPlatforms):
            _namePlatform = Platform.getPlatform(i).getName()
            if _namePlatform is 'CUDA':
                _cudaPlatform = True
                break
        
        if _cudaPlatform:
            platform = Platform.getPlatformByName('CUDA')
            prop = {'Precision': 'single'}
            simulation = Simulation(psf.topology, system, integrator, platform, prop)
        else:
            platform = Platform.getPlatformByName('CPU') # The cpu platform also has no. of
                                                         # threads property could add this
            simulation = Simulation(psf.topology, system, integrator)
        return system, simulation

    def compute_energies(self, dcdfreq, time='femtoseconds', print_output=False):
        """Returns the energies computed in OpenMM in a 'NAMD'-like output
        """
        potential_energies = np.array([])
        output = []
        if time not in self.time_conversions:
            time_error = 'The desired time is not included in time\n' + \
                         'conversion either use allowed time (i.e.,\n' + \
                         'femtoseconds, picoseconds, nanoseconds,\n' + \
                         'microseconds) or add a new time conversion\n' + \
                         'with add_time_conversion method.'  # To do add method add_time_conversion
            return print(time_error)
        for j in range(self.trajectory.n_frames):
            self.simulation.context.setPositions(self.trajectory.openmm_positions(j))
            self.simulation.context.setPeriodicBoxVectors(self.trajectory.openmm_boxes(j)[0],
                                                          self.trajectory.openmm_boxes(j)[1],
                                                          self.trajectory.openmm_boxes(j)[2])
            output.append('Time:{0:> 45.4f} {1}'.format(dcdfreq * self.time_conversions[time] * (j + 1), time))
            for i, f in enumerate(self.system.getForces()):
                state = self.simulation.context.getState(getEnergy=True, groups={i})
                output.append('{0}:{2:> {1}.4f}'.format(f.getName(),
                                                        49 - len(f.getName()),
                                                        state.getPotentialEnergy().value_in_unit(
                                                            kilocalories_per_mole)))
            potential_energy = self.simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
                kilocalories_per_mole)
            potential_energies = np.append(potential_energies, round(potential_energy, 4))
            output.append("Total Potential energy:{:> 27.4f}\n".format(potential_energy))
            if print_output:
                print(*output[12 * j:12 * (j + 1)], sep='\n')
        return potential_energies, output
