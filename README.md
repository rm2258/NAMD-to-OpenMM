# Contents:

## 1.      Summary
## 2.      Documentation
### 2.1.     NamdToOpenmmTools
####  2.1.1.    Attributes
####  2.1.2.    Methods
####  2.1.3.    Constraints
####  2.1.4.    Hard-coded parameters
### 2.2.     NamdEnergyExtractor
####  2.2.1.    Attributes
####  2.2.2.    Methods
### 2.3 namdenergycomputer.tcl
#### 2.3.1 Parameters
## 3.      Common Errors
## 4.      Examples
## 5.      Useful links


## 1.      Summary
An all-in-one repository to get you up in running with **NAMD**, **CHARMM**, or **OpenMM** MD simulations.
Including useful functionalities for swapping between MD packages consisting of simulation
configuration files, and functionalities for calculating '**NAMD**'-like potential energies in
**NAMD** or **OpenMM** from DCD trajectory files and extracting energies from **NAMD** output files. 


## 2.      Documentation
### 2.1.    `NamdToOpenmmTools`
`class namdtoopenmm.NamdToOpenmmTools(trajectory_file, psf_file, parameter_files, temperature, pbc)`

This module is for computing '**NAMD**'-like energies in **OpenMM**

`__init__(self, trajectory_file, psf_file, parameter_files, temperature, pbc)`

**Parameters**: 
- **trajectory_file** (str) - DCD trajectory file name
- **psf_file** (str) - Protein Structure File (PSF) file name
- **parameter_files** (list) - list of parameter files to load
- **temperature** (float) - simulation temperature in kelvin
- **pbc** (tuple) - a, b, c for periodic boundary conditions in angstroms


####  2.1.1.    Attributes
- **trajectory_file**: 
Trajectory file to load

- **psf_file**:
Protein Structure File (psf) to load

- **parameter_files**:
A list of parameter files to load

- **temperature**:
The simulation temperature in kelvin

- **pbc**:
a, b, c for periodic boundary conditions in angstroms

- **time_conversions**:
Time conversion constants

- **trajectory**:
Loaded trajectory file as an mdtraj.Trajectory object

- **system**:
An **OpenMM** system object

- **simulation**:
An **OpenMM** simulation object


####  2.1.2.    Methods
`help()`

Print a help message

**Returns**:
- A help message

`load_trajectory(self)`

Load a trajectory file

**Returns**: 
- An mdtraj.Trajectory object

`load_system(self, parameter_files, temperature, pbc)`

Makes an **OpenMM** simulation object

**Parameters**: 
- **parameter_files** (list) - list of parameter files to load
- **temperature** (float) - simulation temperature in kelvin
- **pbc** (tuple) - a, b, c for periodic boundary conditions in angstroms
**Returns**: 
- **system**, **simulation**
- **system** - An **OpenMM** system object
- **simulation** - An **OpenMM** simulation object

`compute_energies(self, dcdfreq, time='femtoseconds', print_output=False)`

Computes '**NAMD**'-like potential energies in **OpenMM**

**Parameters**: 
- **dcdfreq** (int) - DCD frequency output
- **time** (str) - Desired unit of time
- **print_output** (bool) - Print energy breakdown output
**Returns**: 
- **potential_energies**, **output**
- **potential_energies** (numpy array) - Potential energies
- **output** (list) - List of strings of energy breakdown output


####  2.1.3.    Constraints
Keep in mind that the module was designed to use Protein Structure File (psf)
and the **CHARMM** force field. Therefore, the library at its current state will
only work for '**CHARMM**'-like files.

If a GPU is available in the system then it is used for the calculations.

Additionally, at its current state the library is designed to compute
energies from NPT simulations.


####  2.1.4.    Hard-coded parameters
The following parameters are hard-coded to emulate common **NAMD** and **CHARMM**
force field parameters. 
- Bonds involving hydrogens are constrained
- Timestep is hard coded to 2 fs
- Pressure is hard coded to 1 atm
- Long-range electrostatics are computed with PME
- Nonbonded cutoff is set to 12 angstroms
- vdw are smoothly switched off with a potential-based switching function
  set to 10 angstroms.

### 2.2.    `NamdEnergyExtractor`
`class namdextractor.NamdEnergyExtractor(namd_output_file)`

This module is for extracting energy data from **NAMD** output files

`__init__(self, namd_output_file)`

**Parameters**: 
- **namd_output_file** (str) - Name of **NAMD** log file


####  2.2.1.    Attributes
- **namd_output_file**:
**NAMD** log file

- **data_from_file**:
Dictionary where keys are the ETITLE headers and values are the list energies

- **keys**:
List of ETITLE headers


####  2.2.2.    Methods
`help()`

Print a help message

**Returns**: 
- A help message

`etitle(self)`

Get ETITLE from **NAMD** log file

**Returns**: 
- A dictionary where keys are the ETITLE entries and values are
  empty lists

`energy_extractor(self)`

Extract energies from **NAMD** log file

**Returns**: 
- A numpy array containing potential energies from the **NAMD** log file

### 2.3 namdenergycomputer.tcl
This module is for computing potential energies from 
trajectory files in **NAMD**.

#### 2.3.1 Parameters
- **ts** (int) - The number steps between each trajectory frame, and is eqiuvalent to **NAMD**'s dcdfreq.

- **energy_precision** (int) - The number of digits to appear to the right of the decimal point in the dedicated output file if it is written.

- **trajectory_file** (str) - The file name of the trajectory file

- **print_output** (bool) - Should a dedicated output file be written.

- **output_file** (str) - The file name for the energy output



## 3.      Common Errors
**ModuleNotFoundError** occurs when python library is not in the immediate working
directory. Fixes include copying the modules into your working directory,
defining the full path to the modules, or adding the full path to the
modules to your **PATH** environment variables.


## 4.      Examples

Examples are located in **examples/** directory which includes **NamdToOpenmmTools**,
**NamdEnergyExtractor**, **namdenergycomputer.tcl** examples, and example simulation configuration files.


## 5.      Useful links

**VMD/NAMD**: 
- [**VMD** Download](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)
- [**NAMD** Download](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD)
- [**VMD/NAMD** Tutorials](https://www.ks.uiuc.edu/Training/Tutorials/)
- [**NAMD** User's Guide](https://www.ks.uiuc.edu/Research/namd/2.14/ug/)
- [**NAMD** forums](http://www.ks.uiuc.edu/Research/namd/mailing_list/)

**CHARMM**: 
- [**CHARMM** Download](https://academiccharmm.org/program)
- [**CHARMM** Documentation](https://academiccharmm.org/documentation)

**OpenMM**: 
- [**OpenMM** User's Guide](http://docs.openmm.org/latest/userguide/)
- [**OpenMM** Documentation](http://docs.openmm.org/latest/api-python/)
