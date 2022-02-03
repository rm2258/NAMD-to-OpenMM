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

    def __init__(self, trajectory_file, psf_file, parameter_files, temperature, pbc):
        self.trajectory_file = str(trajectory_file)
        self.psf_file = str(psf_file)
        self.parameter_files = list(parameter_files)
        self.temperature = float(temperature)
        self.pbc = tuple(pbc)
        self.time_conversions = {'femtoseconds': 2, 'picoseconds': 0.002, 'nanoseconds': 2 * 10 ** (-6),
                                 'microseconds': 2 * 10 ** (-9)}
        self.trajectory = self.load_trajectory()
        self.system, self.simulation = self.load_system(self.parameter_files, self.temperature, self.pbc)
        print('Use the help() method for class use')

    @staticmethod
    def help():
        help_message = 'class namdtoopenmm.NamdToOpenmmTools(trajectory_file, psf_file, parameter_files, temperature,' \
                       ' pbc)\n' + \
                       '  This module is for computing \'NAMD\'-like energies in OpenMM\n\n' + \
                       '  __init__(self, trajectory_file, psf_file, parameter_files, temperature, pbc)\n' + \
                       '    Parameters: trajectory_file (str) - File name of DCD trajectory file\n' + \
                       '                psf_file (str) - File name of PSF file\n' + \
                       '                parameter_files (list) - list of parameter files to load\n' + \
                       '                temperature (float) - simulation temperature in kelvin\n' + \
                       '                pbc (tuple) - a, b, c for periodic boundary conditions in angstroms\n\n' + \
                       '  Attributes\n' + \
                       '    trajectory_file\n' + \
                       '      Trajectory file to load\n\n' + \
                       '    psf_file\n' + \
                       '      Protein Structure File (psf) to load\n\n' + \
                       '    parameter_files\n' + \
                       '      A list of parameter files to load\n\n' + \
                       '    temperature\n' + \
                       '      The simulation temperature in \n\n' + \
                       '    pbc\n' + \
                       '      a, b, c for periodic boundary conditions in angstroms\n\n' + \
                       '    time_conversions\n' + \
                       '      Time conversion constants\n\n' + \
                       '    trajectory\n' + \
                       '      Loaded trajectory file as an mdtraj.Trajectory object\n\n' + \
                       '   system\n' + \
                       '     An OpenMM system object\n\n' + \
                       '    simulation\n' + \
                       '      Instantiate an OpenMM simulation object\n\n' + \
                       '  Methods\n' + \
                       '    help()\n' + \
                       '      Print a help message\n\n' + \
                       '      Returns: A help message\n\n' + \
                       '    load_trajectory(self)\n' + \
                       '      Load a trajectory file\n\n' + \
                       '      Returns: An mdtraj.Trajectory object\n\n' + \
                       '    load_system(self, parameter_files, temperature, pbc)\n' + \
                       '      Makes an OpenMM simulation object\n\n' + \
                       '      Parameters: parameter_files (list) - list of parameter files to load\n' + \
                       '                  temperature (float) - simulation temperature in kelvin\n' + \
                       '                  pbc (tuple) - a, b, c for periodic boundary conditions in angstroms\n' + \
                       '      Returns: system, simulation\n' + \
                       '               system - An OpenMM system object\n' + \
                       '               simulation - An OpenMM simulation object\n\n' + \
                       '    compute_energies(self, dcdfreq, time=\'femtoseconds\', print_output=False)\n' + \
                       '      Computes \'NAMD\'-like potential energies in OpenMM\n\n' + \
                       '      Parameters: dcdfreq (int) - DCD frequency output\n' + \
                       '                  time (str) - Desired unit of time\n' + \
                       '                  print_output (bool) - Print energy breakdown output\n' + \
                       '      Returns: potential_energies, output\n' + \
                       '               potential_energies (numpy array) - Potential energies\n' + \
                       '               output (list) - List of strings of energy breakdown output\n\n'
        print(help_message)

    def load_trajectory(self):
        traj = md.load_dcd(self.trajectory_file, self.psf_file)
        return traj

    def load_system(self, parameter_files, temperature, pbc):
        psf = CharmmPsfFile(self.psf_file)
        params = CharmmParameterSet(*parameter_files)
        psf.setBox(pbc[0] * angstroms, pbc[1] * angstroms, pbc[2] * angstroms)
        system = psf.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=1.2 * nanometers,
                                  switchDistance=1.0 * nanometers, constraints=HBonds, ewaldErrorTolerance=0.0005)
        pressure_controller = MonteCarloMembraneBarostat
        system.addForce(pressure_controller(1.01325 * bar, 0.0 * bar * nanometers, temperature * kelvin,
                                            pressure_controller.XYIsotropic, pressure_controller.ZFree, 100))
        for i, f in enumerate(system.getForces()):
            f.setForceGroup(i)
        integrator = LangevinIntegrator(temperature * kelvin, 1 / picosecond, 0.002 * picoseconds)
        try:
            platform = Platform.getPlatformByName('CUDA')
            prop = {'Precision': 'single'}
            simulation = Simulation(psf.topology, system, integrator, platform, prop)
        except OpenMMException:
            platform = Platform.getPlatformByName('CPU')
            simulation = Simulation(psf.topology, system, integrator)
        return system, simulation

    def compute_energies(self, dcdfreq, time='femtoseconds', print_output=False):
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
