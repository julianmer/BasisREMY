####################################################################################################
#                                      basic_usage.py                                              #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                   #
#                                                                                                  #
# Created: 18/02/26                                                                                #
#                                                                                                  #
# Purpose: Example showing how to use BasisREMY without the GUI.                                  #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from core.basisremy import BasisREMY

if __name__ == "__main__":

    # =============================================================================================
    # Initialize BasisREMY
    # =============================================================================================

    br = BasisREMY()


    # =============================================================================================
    # Load MRS data file and extract parameters
    # =============================================================================================
    # Supported formats: .spar (Philips), .7 (GE), .dat/.rda/.ima (Siemens),
    #                    method (Bruker), .nii/.nii.gz (NIfTI)

    import_fpath = './example_data/BigGABA_P1P_S01/S01_PRESS_35_act.SPAR'

    print(f"Processing file: {import_fpath}")
    params = br.runREMY(import_fpath=import_fpath)


    # =============================================================================================
    # Set backend and parse parameters
    # =============================================================================================
    # Available backends:
    #   - LCModel  : Fast ideal pulse simulations (PRESS, STEAM, LASER)
    #   - sLaserSim: Realistic Bloch simulations (sLASER, semi-LASER)

    br.set_backend('LCModel')
    parsed_params, opt = br.backend.parseREMY(params)

    # Handle empty parameters
    def get_param(key, default):
        value = parsed_params.get(key, default)
        if value in [None, '', 'nan']:
            return default
        return value


    # =============================================================================================
    # Configure simulation parameters
    # =============================================================================================

    simulation_params = {
        'Sequence': get_param('Sequence', 'PRESS'),
        'Samples': get_param('Samples', 2048),
        'Bandwidth': get_param('Bandwidth', 2000),
        'Bfield': get_param('Bfield', 3.0),
        'Linewidth': 1,
        'TE': get_param('TE', 35),
        'TE2': 0,
        'Add Ref.': 'No',
        'Make .raw': 'Yes',
        'Output Path': './output/my_basis_set',
        'Center Freq': get_param('Center Freq', 127736713),
        'Metabolites': ['NAA', 'Cr', 'PCr', 'Glu', 'Gln', 'Ins', 'GABA', 'GSH', 'Lac', 'Tau'],
    }

    print(f"Simulating basis set with the following parameters:")
    for key, value in simulation_params.items():
        if key != 'Metabolites':
            print(f"{key}: {value}")
    print(f"Metabolites: {simulation_params['Metabolites']}")


    # =============================================================================================
    # Run simulation
    # =============================================================================================

    basis_set = br.backend.run_simulation(simulation_params)

    print(f"\nSimulation complete.")
    print(f"Generated {len(basis_set)} metabolite spectra")
    print(f"Output location: {simulation_params['Output Path']}")

