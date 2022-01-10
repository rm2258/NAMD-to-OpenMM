"""
This module is for computing 'NAMD'-like energies in OpenMM

Made by Ramon Mendoza, ramendoza@uchicago.edu
Last update January 10, 2022
"""

import numpy as np


class NamdEnergyExtractor:
    """
    Add a step to time converter in this class as well
    will need to import matplotlib.pyplot as plt
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
        print('Use the help() method for class use')

    @staticmethod
    def help():
        help_message = 'class namdextractor.NamdEnergyExtractor(namd_output_file)\n' + \
                       '  This module is for extracting energy data from NAMD output files\n\n' + \
                       '  __init__(self, namd_output_file)\n' + \
                       '    Parameters: namd_output_file (str) - Name of NAMD log file\n\n' + \
                       '  Attributes\n' + \
                       '    namd_output_file\n' + \
                       '      NAMD log file\n\n' + \
                       '    data\n' + \
                       '      Dictionary where keys are the ETITLE entries and values are the energies\n\n' + \
                       '   keys\n' + \
                       '     List of ETITLE headers\n\n' + \
                       '  Methods\n' + \
                       '    help()\n' + \
                       '      Print a help message\n\n' + \
                       '      Returns: A help message\n\n' + \
                       '    etitle(self)\n' + \
                       '      Get ETITLE from NAMD log file\n\n' + \
                       '      Returns: A dictionary where keys are the ETITLE entries and values are\n' + \
                       '               empty lists\n\n' + \
                       '    energy_extractor(self)\n' + \
                       '      Extract energies from NAMD log file\n\n' + \
                       '      Returns: A numpy array containing potential energies from NAMD log file'
        print(help_message)  # Decide if you want to keep this as return

    def etitle(self):
        """
        Use in case log file formatting Changes
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
