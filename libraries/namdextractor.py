"""
This module is for computing 'NAMD'-like energies in OpenMM

Made by Ramon Mendoza, ramendoza@uchicago.edu
Last update January 10, 2022
"""

import numpy as np


class NamdEnergyExtractor:
    """
    class namdextractor.NamdEnergyExtractor(namd_output_file)
                         This module is for extracting energy data from NAMD output files
                         __init__(self, namd_output_file)
                         Parameters: 
                             namd_output_file (str) - Name of NAMD log file
                         Attributes:
                             namd_output_file: NAMD log file
                             data: Dictionary where keys are the ETITLE entries and values are the energies
                             keys: List of ETITLE headers
                         Methods
                             etitle(self)
                                 Get ETITLE from NAMD log file
                                 Returns: A dictionary where keys are the ETITLE entries and 
                                          values are empty lists
                             energy_extractor(self)
                                 Extract energies from NAMD log file
                                 Returns: A numpy array containing potential energies from NAMD 
                                          log file
    """

    def __init__(self, namd_output_file):
        self.namd_output_file = str(namd_output_file)
        # data = etitle() # To use in case log file is different in the future
        self.data_from_file = {'TS': [], 'BOND': [], 'ANGLE': [], 'DIHED': [], 'IMPRP': [], 'ELECT': [], 'VDW': [],
                               'BOUNDARY': [], 'MISC': [], 'KINETIC': [], 'TOTAL': [], 'TEMP': [], 'POTENTIAL': [],
                               'TOTAL3': [], 'TEMPAVG': [], 'PRESSURE': [], 'GPRESSURE': [], 'VOLUME': [],
                               'PRESSAVG': [], 'GPRESSAVG': []}
        self.keys = ['ETITLE:', 'TS', 'BOND', 'ANGLE', 'DIHED', 'IMPRP', 'ELECT', 'VDW', 'BOUNDARY', 'MISC', 'KINETIC',
                     'TOTAL', 'TEMP', 'POTENTIAL', 'TOTAL3', 'TEMPAVG', 'PRESSURE', 'GPRESSURE', 'VOLUME', 'PRESSAVG',
                     'GPRESSAVG']

    def etitle(self):
        """Returns the headers of a log file
        """
        with open(self.namd_output_file, 'r') as f:
            lines = f.readlines()
            first_etitle = True
            for line in lines:
                if first_etitle and 'ETITLE' in line:
                    keys = line.split()
                    data = {keys[k]: [] for k in range(1, len(keys))}
                    first_etitle = False
            f.close()
        return data

    def energy_extractor(self):
        """Returns the energies for each of the headers
        """
        with open(self.namd_output_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'ENERGY:' in line and 'Info:' not in line:
                    entries = line.split()
                    for i in range(1, len(entries)):
                        self.data_from_file[self.keys[i]].append(float(entries[i]))
            f.close()
        datanp = np.array(self.data_from_file['POTENTIAL'])
        print('Note: Energy at index 0 in NAMD is the starting energy and is redundant if data is from a '
              'previous simulation')
        return datanp
