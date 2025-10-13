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

    # initialize BasisREMY and run tests
    br = BasisREMY()

    check = ['Protocol', 'NumberOfDatapoints', 'SpectralWidth', 'B0', 'TE']
    filled = 0
    failed = []
    for key, value in tests.items():
        print(f"Running test: {key}")

        params = br.runREMY(**value)
        print(params)

        # compute percentage of fields filled for lcmodel backend
        # needed are: sequence, samples, bandwidth, bfield, te
        fill = sum([1 for k in check if params.get(k) not in [None, '']])
        print(f"Filled {fill}/{len(check)} mandatory fields ({(fill/len(check))*100:.1f}%)")

        filled += fill
        if fill < len(check):
            failed.append(key)

    print(f"Overall filled {filled}/{len(tests)*len(check)} mandatory fields "
          f"({(filled/(len(tests)*len(check)))*100:.1f}%)")

    print(f"Failed tests: {failed}")

    print("Tests complete.")


