####################################################################################################
#                                      basic_usage.py                                              #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 18/02/26                                                                                #
#                                                                                                  #
# Purpose: Example showing how to use BasisREMY without the GUI.                                   #
#                                                                                                  #
#          Run from the repository root:  python examples/basic_usage.py                           #
#          (needs an Octave runtime — Docker or a local Octave — for the simulation step).         #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from basisremy.core.basisremy import BasisREMY
from basisremy.core.exporters import export

if __name__ == "__main__":

    # =============================================================================================
    # Initialize BasisREMY
    # =============================================================================================

    br = BasisREMY()
    print(f"Available backends: {br.available_backends}")


    # =============================================================================================
    # Load MRS data file and extract parameters
    # =============================================================================================
    # Supported formats: .spar (Philips), .7 (GE), .dat/.rda/.ima (Siemens),
    #                    method (Bruker), .nii/.nii.gz (NIfTI)

    import_fpath = './example_data/BigGABA_P1P_S01/S01_PRESS_35_act.SPAR'

    print(f"Processing file: {import_fpath}")
    br.runREMY(import_fpath=import_fpath)


    # =============================================================================================
    # Choose a backend
    # =============================================================================================
    # set_backend() re-parses the extracted parameters for the chosen backend and
    # fills its mandatory_params (Samples, Bandwidth, Bfield, TE, Center Freq, ...).

    br.set_backend('FidaIdeal')   # fast ideal-pulse FID-A simulation (SE / PRESS / STEAM / LASER)


    # =============================================================================================
    # Configure simulation parameters
    # =============================================================================================
    # Start from what REMY extracted and override / complete by hand.

    params = dict(br.backend.mandatory_params)
    params.update({
        'Sequence': 'PRESS',
        'Linewidth': 1,
        'Metabolites': ['NAA', 'Cr', 'PCr', 'Glu', 'Gln', 'Ins', 'GABA', 'GSH', 'Lac', 'Tau'],
    })

    missing = [k for k, v in params.items() if v in (None, '')]
    if missing:
        raise SystemExit(f"Parameters not found in the data, please set them by hand: {missing}")

    print("Simulating basis set with the following parameters:")
    for key, value in params.items():
        print(f"  {key}: {value}")


    # =============================================================================================
    # Run simulation
    # =============================================================================================
    # Returns { metabolite -> complex FID }. Octave (Docker or local) is started on first use.

    basis = br.backend.run_simulation(params)
    print(f"\nSimulation complete: {len(basis)} metabolite spectra")


    # =============================================================================================
    # Export
    # =============================================================================================
    # Formats: lcmodel_basis, lcmodel_raw, jmrui_txt, fsl_json, osprey_mat, fida_mat,
    #          inspector_mat, profit_mat, marss_mat, mrscloud_mat, spinwizard
    # A *_sidecar.json with the parameters and provenance is written next to every export.

    out = export(basis, './output/PRESS_TE35_3T', 'lcmodel_basis', params)
    print(f"Exported to: {out}")
