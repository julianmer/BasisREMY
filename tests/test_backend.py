####################################################################################################
#                                         test_backend.py                                          #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 13/10/25                                                                                #
#                                                                                                  #
# Purpose: Test script for running the basisREMY backend on various example datasets to simulate   #
#          the basis sets.                                                                         #
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

    # configuration
    config = {
        'backend': 'LCModel',   # the backend to use for the simulation

        'run_config': {
            'import_fpath': None,   # the file to import
            'export_fpath': './output/',   # the output file path
            'method': None,       # the type of raw MRS data (if None, auto-detect)
            'userParams': {},    # user parameters to override or add to the REMY parameters
            'optionalParams': {},   # optional parameters for the backend
            'plot': False,       # whether to plot the simulated basis set
        },

        'metab_config': {   # configuration for the metabolites to include in the basis set
            'Ala': False,
            'Asc': False,
            'Asp': False,
            'Ch': False,
            'Cit': False,
            'Cr': True,
            'EtOH': False,
            'GABA': False,
            'GPC': False,
            'GSH': False,
            'Glc': False,
            'Gln': False,
            'Glu': False,
            'Gly': False,
            'H2O': False,
            'Ins': False,
            'Lac': False,
            'Lip': False,
            'NAA': True,
            'NAAG': False,
            'PCh': False,
            'PCr': False,
            'PE': False,
            'Phenyl': False,
            'Ref0ppm': False,
            'Scyllo': False,
            'Ser': False,
            'Tau': False,
            'Tyros': False,
        },
    }

    # define test cases
    tests = {

        # Big GABA test sets
        'REMY_tests/Big_GABA_P1P_S01': {
            'import_fpath': './example_data/BigGABA_P1P_S01/S01_PRESS_35_act.SPAR',
        },

        'REMY_tests/Big_GABA/G1_P': {
            'import_fpath': './example_data/BigGABA_G1P_S01/S01_PRESS_35.7',
        },

        'REMY_tests/Big_GABA/S1_P': {
            'import_fpath': './example_data/BigGABA_S1P_S01/S01_PRESS_35.dat',
        },
    }

    # add all REMY tests
    remy_test_dir = './example_data/REMY_tests/'
    for test in sorted(os.listdir(remy_test_dir)):
        if test.startswith('Dataset_'):
            # go through all files in the directory and find the raw data file
            test_dir = os.path.join(remy_test_dir, test)
            for file in os.listdir(test_dir):
                if file.lower().endswith(('.spar', '.7', '.dat', '.ima', '.rda', '.nii', '.nii.gz', '_method')):
                    tests[test] = {'import_fpath': os.path.join(test_dir, file)}
                    break

    # initialize BasisREMY
    br = BasisREMY(backend=config['backend'])

    # metabs
    br.backend.mandatory_params['Metabolites'] = [key for key, value in config['metab_config'].items() if value]

    # run tests
    successes = []
    for key, value in tests.items():

        params = br.runREMY(**value)
        params, opt = br.backend.parseREMY(params)   # parser will return None for unsupported sequences
                                                     # of the given backend

        print(f"Running test: {key} with backend {config['backend']} "
              f"and sequence {params['Sequence']}")

        try:
            if config['backend'] == 'LCModel':
                if params['Sequence'] in ['Spin Echo', 'PRESS', 'STEAM', 'LASER']:
                    config['run_config'].update(value)
                    basis, params = br.run(**config['run_config'])
                    print("Success!")
                    successes.append(key)

            elif config['backend'] == 'sLaserSim':
                if params['Sequence'] == 'sLASER':
                    config['run_config'].update(value)
                    basis, params = br.run(**config['run_config'])
                    print("Success!")
                    successes.append(key)

            else:
                print("Currently no other backends are supported. "
                      "Please add them to the test script.")

        except Exception as e:
            print(f"Test {key} with backend {config['backend']} failed with error: {e}")

    print(f"Successful tests: {successes}")

    print("All tests completed.")