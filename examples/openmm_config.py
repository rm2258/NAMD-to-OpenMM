"""
JOB DESCRIPTION

Production run with PME, Constant Pressure and Temperature.
Demonstrates how to perform MD in OpenMM

"""

from sys import stdout, exit, stderr

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

# Load CHARMM files
print("Loading parameters")
psf = CharmmPsfFile('example.psf')          # psf file
pdb = PDBFile('example.pdb')                # coordinate file
outputName = 'example'                      # output name
saveXML = False                             # output a state or checkpoint file?

temp = 320                                  # Temperature (K)

restartFile = ''                     # Restart file to continue a simulation 
                                     # in either xml or chk format

#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

params = CharmmParameterSet('par_all36_prot.prm', 'top_all36_prot.rtf',
                            'par_all36_lipid.prm', 'top_all36_lipid.rtf',
                            'par_all36_carb.prm', 'top_all36_carb.rtf',
			    'par_all36_cgenff.prm', 'top_all36_cgenff.rtf',
			    'par_all36_na.prm', 'top_all36_na.rtf',
			    'toppar_water_ions_namd.str')                                          # Parameter files

psf.setBox(90.3397026816*angstroms, 91.3840923081*angstroms, 83.3302382086*angstroms)              # Box dimensions for PBC

# Build openmm system
system = psf.createSystem(params, nonbondedMethod=PME, 
                          nonbondedCutoff=1.2*nanometers, switchDistance=1.0*nanometers,
			  constraints=HBonds, ewaldErrorTolerance=0.0005)
system.addForce(MonteCarloMembraneBarostat(1.01325*bar, 0.0*bar*nanometers,
                                           temp*kelvin, MonteCarloMembraneBarostat.XYIsotropic,
		     			   MonteCarloMembraneBarostat.ZFree, 100))                 # Constant Pressure Control (variable volume)
                                                                                                   # 1.01325 bar -> 1 atm
integrator = LangevinIntegrator(temp*kelvin, 1/picosecond,
                                0.002*picoseconds)                                                 # Constant Temperature Control
				                                                                   # Time-step (ps)

# Use GPU
platform = Platform.getPlatformByName('CUDA')
prop = {'Precision': 'single'}                                                                     # Precision mode can be single, mixed, or double
                                                                                                   # double give most accurate calculations, but
												   # at a performance cost
# Build simulation context
simulation = Simulation(psf.topology, system, integrator, platform, prop)
simulation.context.setPositions(pdb.positions)

if 'xml' in restartFile:                     # Continue a job from a saved state
    simulation.loadState(restartFile)
if 'chk' in restartFile:                     # Continue a job from a checkpoint
    simulation.loadCheckpoint(restartFile)

# Calculate initial system energy
print("\nInitial system energy: " + str(simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalories_per_mole)))

# Energy minimization
if False:
    mini_nsteps = 100
    print("\nEnergy minimization: %s steps" % mini_nstep)
    simulation.minimizeEnergy(tolerance=100.0*kilojoule/mole, maxIterations=mini_nstep)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalories_per_mole))

# Generate initial velocities
if False:
    print("\nGenerate initial velocities")
    simulation.context.setVelocitiesToTemperature(temp)

# Production
if False:
    nstep = 300000000                          # 500000 steps = 1 ns
    dcdfreq = 5000                             # 5000 steps = every 0.01 ns
    outputEnergies = 5000
    print("\nMD run: %s steps" % nstep)
    simulation.reporters.append(DCDReporter(outputName+'.dcd', dcdfreq))
    simulation.reporters.append(StateDataReporter(sys.stdout, outputEnergies, step=True, time=True,
                                potentialEnergy=True, temperature=True, progress=True,
				remainingTime=True, speed=True, totalSteps=nstep, separator='\t'))
    simulation.step(nstep)

# Write restart file
if saveXML:
    simulation.saveState(outputName+'.xml')
else:
    simulation.saveCheckpoint(outputName+'.chk')


