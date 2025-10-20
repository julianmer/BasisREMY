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

        'test_vendor_file': '.dat',  # filter based on file type (None for all)

        'run_config': {
            'import_fpath': None,   # the file to import
            'export_fpath': './output/',   # the output file path
            'method': None,       # the type of raw MRS data (if None, auto-detect)
            'userParams': {'Sequence': 'PRESS', 'Samples': 8096},    # user parameters to override or add to the REMY parameters
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

    # filetypes
    filetype = config['test_vendor_file'] if config['test_vendor_file'] is not None \
        else ('.spar', '.7', '.dat', '.ima', '.rda', '.nii', '.nii.gz', 'method')

    # define test cases
    tests = {}

    # quick test cases
    # tests.update({
    #
    #     # Big GABA test sets
    #     'REMY_tests/Big_GABA_P1P_S01': {
    #         'import_fpath': './example_data/BigGABA_P1P_S01/S01_PRESS_35_act.SPAR',
    #     },
    #
    #     'REMY_tests/Big_GABA/G1_P': {
    #         'import_fpath': './example_data/BigGABA_G1P_S01/S01_PRESS_35.7',
    #     },
    #
    #     'REMY_tests/Big_GABA/S1_P': {
    #         'import_fpath': './example_data/BigGABA_S1P_S01/S01_PRESS_35.dat',
    #     },
    # })

    # # add all REMY tests
    # remy_test_dir = './example_data/REMY_tests/'
    # for test in sorted(os.listdir(remy_test_dir)):
    #     if test.startswith(('Dataset_00', 'Dataset_01', 'Dataset_11', 'Dataset_12', 'Dataset_13', 'Dataset_28')):  # Filter on PRESS, STEAM, LASER, Spin Echo
    #         # go through all files in the directory and find the raw data file
    #         test_dir = os.path.join(remy_test_dir, test)
    #         for file in os.listdir(test_dir):
    #             if file.lower().endswith(filetype):
    #                 tests[test] = {'import_fpath': os.path.join(test_dir, file)}
    #                 break

    # # add all NIfTI REMY tests
    # remy_test_dir = './example_data/REMY_tests/Datasets_nifti/'
    # # go through all files in the directory and find the raw data file
    # for file in os.listdir(remy_test_dir):
    #     if file.lower().endswith(filetype):
    #         tests[file] = {'import_fpath': os.path.join(remy_test_dir, file)}

    # add (all) Big GABA test sets
    biggaba_test_dir = './example_data/BigGABA/'
    for site in sorted(os.listdir(biggaba_test_dir)):
        if site in ['G4_P', 'G6_P', 'G8_P', 'S1_P', 'S5_P', 'S6_P', 'P1_P', 'P3_P', 'P6_P']:
            site_dir = os.path.join(biggaba_test_dir, site)
            for subject in sorted(os.listdir(site_dir)):
                subject_dir = os.path.join(site_dir, subject)
                for file in os.listdir(subject_dir):
                    if file.lower().endswith(filetype):
                        tests[f'BigGABA_{site}_{subject}'] = {'import_fpath': os.path.join(subject_dir, file)}
                        break

    # # add all spec2nii test sets
    # spec2nii_test_dir = './example_data/spec2nii_tests/'
    # for vendor in sorted(os.listdir(spec2nii_test_dir)):
    #     vendor_dir = os.path.join(spec2nii_test_dir, vendor)
    #     for study in sorted(os.listdir(vendor_dir)):
    #         study_dir = os.path.join(vendor_dir, study)
    #         if os.path.isdir(study_dir):
    #             for file in os.listdir(study_dir):
    #                 if os.path.isfile(os.path.join(study_dir, file)):
    #                     if file.lower().endswith(filetype):
    #                         if file == 'T1.nii.gz': continue
    #                         tests[f'spec2nii_{vendor}_{study}_{file}'] = {'import_fpath': os.path.join(study_dir, file)}
    #                         # break
    #                 elif os.path.isdir(os.path.join(study_dir, file)):
    #                     subdir = os.path.join(study_dir, file)
    #                     for subfile in os.listdir(subdir):
    #                         if os.path.isfile(os.path.join(subdir, subfile)):
    #                             if subfile == 'T1.nii.gz': continue
    #                             if subfile.lower().endswith(filetype):
    #                                 tests[f'spec2nii_{vendor}_{study}_{file}_{subfile}'] = {'import_fpath': os.path.join(subdir, subfile)}
    #                                 # break
    #                         elif os.path.isdir(os.path.join(subdir, subfile)):
    #                             subsubdir = os.path.join(subdir, subfile)
    #                             for subsubfile in os.listdir(subsubdir):
    #                                 if subsubfile == 'T1.nii.gz': continue
    #                                 if subsubfile.lower().endswith(filetype):
    #                                     tests[f'spec2nii_{vendor}_{study}_{file}_{subfile}_{subsubfile}'] = {'import_fpath': os.path.join(subsubdir, subsubfile)}
    #                                     # break


    # initialize BasisREMY
    br = BasisREMY(backend=config['backend'])

    # metabs
    br.backend.mandatory_params['Metabolites'] = [key for key, value in config['metab_config'].items() if value]

    # run tests
    successes, failures = [], []
    check = ['Sequence', 'Samples', 'Bandwidth', 'Bfield', 'TE']   # check backend parsing
    results = {c: {'Success': 0, 'Fail': 0} for c in check}
    for key, value in tests.items():

        try:
            params = br.runREMY(**value)
        except Exception as e:
            print(f"Test {key} failed during REMY run or parsing with error: {e}")
            params = {}
        params, opt = br.backend.parseREMY(params)   # parser will return None for unsupported sequences
                                                     # of the given backend

        # count filled parameters
        for c in check:
            if params[c] not in [None, '']:
                results[c]['Success'] += 1
            else:
                results[c]['Fail'] += 1

        # update the params with the input from users
        params.update(config['run_config']['userParams'])
        opt.update(config['run_config']['optionalParams'])

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
            failures.append(key)

    print("Test summary:")
    print(f"All tests ({len(tests)}): {list(tests.keys())}")
    print(f"Successful tests ({len(successes)}/{len(tests)}): {successes}")
    print(f"Failed tests ({len(failures)}/{len(tests)}): {failures}")

    print("Parameter filling results:")
    for c in check:
        print(f"  {c}: {results[c]['Success']} filled, {results[c]['Fail']} missing, "
              f"({(results[c]['Success']/(results[c]['Success']+results[c]['Fail']))*100:.1f}%)")

    print("All tests completed.")