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
    tests = {}

    # quick test cases
    tests.update({

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
    })

    # add all REMY tests
    remy_test_dir = './example_data/REMY_tests/'
    for test in sorted(os.listdir(remy_test_dir)):
        if test.startswith('Dataset_'):
            # go through all files in the directory and find the raw data file
            test_dir = os.path.join(remy_test_dir, test)
            for file in os.listdir(test_dir):
                if file.lower().endswith(('.spar', '.7', '.dat', '.ima', '.rda', '.nii', '.nii.gz', 'method')):
                    tests[test] = {'import_fpath': os.path.join(test_dir, file)}
                    break

    # add all Big GABA test sets
    biggaba_test_dir = './example_data/BigGABA/'
    for site in sorted(os.listdir(biggaba_test_dir)):
        if site in ['G4_P', 'G6_P', 'G8_P', 'S1_P', 'S5_P', 'S6_P', 'P1_P', 'P3_P', 'P6_P']:
            site_dir = os.path.join(biggaba_test_dir, site)
            for subject in sorted(os.listdir(site_dir)):
                subject_dir = os.path.join(site_dir, subject)
                for file in os.listdir(subject_dir):
                    if file.lower().endswith(('.spar', '.7', '.dat', '.ima', '.rda', '.nii', '.nii.gz', 'method')):
                        tests[f'BigGABA_{site}_{subject}'] = {'import_fpath': os.path.join(subject_dir, file)}
                        break

    # add all spec2nii test sets
    spec2nii_test_dir = './example_data/spec2nii_tests/'
    for vendor in sorted(os.listdir(spec2nii_test_dir)):
        if vendor.lower() not in ['ge']: continue#['ge', 'philips', 'bruker', 'siemens']: continue
        vendor_dir = os.path.join(spec2nii_test_dir, vendor)
        for study in sorted(os.listdir(vendor_dir)):
            study_dir = os.path.join(vendor_dir, study)
            if os.path.isdir(study_dir):
                for file in os.listdir(study_dir):
                    if os.path.isfile(os.path.join(study_dir, file)):
                        if file.lower().endswith(('.spar', '.7', '.dat', '.ima', '.rda', '.nii', '.nii.gz', 'method')):
                            if file == 'T1.nii.gz': continue
                            tests[f'spec2nii_{vendor}_{study}_{file}'] = {'import_fpath': os.path.join(study_dir, file)}
                            # break
                    elif os.path.isdir(os.path.join(study_dir, file)):
                        subdir = os.path.join(study_dir, file)
                        for subfile in os.listdir(subdir):
                            if os.path.isfile(os.path.join(subdir, subfile)):
                                if subfile == 'T1.nii.gz': continue
                                if subfile.lower().endswith(('.spar', '.7', '.dat', '.ima', '.rda', '.nii', '.nii.gz', 'method')):
                                    tests[f'spec2nii_{vendor}_{study}_{file}_{subfile}'] = {'import_fpath': os.path.join(subdir, subfile)}
                                    # break
                            elif os.path.isdir(os.path.join(subdir, subfile)):
                                subsubdir = os.path.join(subdir, subfile)
                                for subsubfile in os.listdir(subsubdir):
                                    if subsubfile == 'T1.nii.gz': continue
                                    if subsubfile.lower().endswith(('.spar', '.7', '.dat', '.ima', '.rda', '.nii', '.nii.gz', 'method')):
                                        tests[f'spec2nii_{vendor}_{study}_{file}_{subfile}_{subsubfile}'] = {'import_fpath': os.path.join(subsubdir, subsubfile)}
                                        # break

    print(tests.keys())

    # initialize BasisREMY and run tests
    br = BasisREMY()

    # check = ['Protocol', 'NumberOfDatapoints', 'SpectralWidth', 'B0', 'TE']   # check REMY parsing
    check = ['Sequence', 'Samples', 'Bandwidth', 'Bfield', 'TE']   # check backend parsing
    filled = 0
    failed = []
    for key, value in tests.items():
        print(f"Running test: {key}")

        try:
            params = br.runREMY(**value)
            params, opt = br.backend.parseREMY(params)  # parser will return None for unsupported sequences
                                                        # of the given backend
            print(params)
        except Exception as e:
            print(f"  Test {key} failed with error: {e}")
            params = {}

        # compute percentage of fields filled for lcmodel backend
        # needed are: sequence, samples, bandwidth, bfield, te
        fill = sum([1 for k in check if params.get(k) not in [None, '']])
        print(f"Filled {fill}/{len(check)} mandatory fields ({(fill/len(check))*100:.1f}%)")

        filled += fill
        if fill < len(check):
            failed.append(key)

    print(f"Overall filled {filled}/{len(tests)*len(check)} mandatory fields "
          f"({(filled/(len(tests)*len(check)))*100:.1f}%)")

    print(f"Total tests: {len(tests)}, success rate: {(len(tests)-len(failed))/len(tests)*100:.1f}%")
    print(f"Failed tests: {failed}")

    print("Tests complete.")


