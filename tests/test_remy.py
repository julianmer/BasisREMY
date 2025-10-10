####################################################################################################
#                                          test_remy.py                                            #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 09/10/25                                                                                #
#                                                                                                  #
# Purpose: Test script for the metadata extraction provided by REMY.                               #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# own
from core.basisremy import BasisREMY

if __name__ == "__main__":

    # define test cases
    tests = {

        # Big GABA test sets
        'REMY_tests/Big_GABA_P1P_S01': {
            'import_fpath': './example_data/BigGABA_P1P_S01/S01_PRESS_35_act.SPAR',
            'export_fpath': './output/',
            'plot': True,
        },

        'REMY_tests/Big_GABA/G1_P': {   # NOT WORKING YET
            'import_fpath': './example_data/BigGABA_G1P_S01/S01_PRESS_35.7',
            'export_fpath': './output/',
            'plot': True,
        },

        'REMY_tests/Big_GABA/S1_P': {   # NOT WORKING YET
            'import_fpath': './example_data/BigGABA_S1P_S01/S01_PRESS_35.dat',
            'export_fpath': './output/',
            'plot': True,
        },

        # REMY test sets
        'REMY_tests/Dataset_00': {  # TODO: fix issue with bandwidth not being an int
            'import_fpath': './example_data/REMY_tests/Dataset_00_Bruker_14T_STEAM_08/Dataset_00_Bruker_14T_STEAM_08_method',
            'export_fpath': './output/',
            'userParams': {},
            'plot': True,
        },

        'REMY_tests/Dataset_01': {
            'import_fpath': './example_data/REMY_tests/Dataset_01_Bruker_14T_STEAM_09/Dataset_01_Bruker_14T_STEAM_09_method',
            'export_fpath': './output/',
            'userParams': {},
            'plot': True,
        },
        
        'REMY_tests/Dataset_02': {  # NOT WORKING YET
            'import_fpath': './example_data/REMY_tests/Dataset_02_GE_3T_pFile_MRSI_P13312/Dataset_02_GE_3T_pFile_MRSI_P13312.7',
            'export_fpath': './output/',
            'userParams': {},
            'plot': True,
        },

        'REMY_tests/Dataset_04': {
            'import_fpath': './example_data/REMY_tests/Dataset_04_GE_3T_pFile_P23040/Dataset_04_GE_3T_pFile_P23040.7',
            'export_fpath': './output/',
            'userParams': {},
            'plot': True,
        },
    }

    # initialize BasisREMY and run tests
    br = BasisREMY(backend='LCModel')

    for key, value in tests.items():
        print(f"Running test: {key}")
        br.run(**value)
