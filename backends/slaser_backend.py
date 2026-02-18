####################################################################################################
#                                         slaser_backend.py                                        #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 08/10/25                                                                                #
#                                                                                                  #
# Purpose: Defines the sLaserBackend class for simulating MRS basis sets using the                 #
#          sLASER sequence with the sLASER_makebasisset_function.m function.                       #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import os
import numpy as np

# own
from backends.base import Backend


#**************************************************************************************************#
#                                          sLaserBackend                                           #
#**************************************************************************************************#
#                                                                                                  #
# Implements the basis set simulation backend for the sLASER sequence.                             #
#                                                                                                  #
#**************************************************************************************************#
class sLaserBackend(Backend):
    def __init__(self):
        super().__init__()
        self.name = 'sLaserSim'

        # Mark that this backend requires Octave
        self.requires_octave = True
        # Octave will be initialized lazily when needed


        # define possible metabolites
        self.metabs = {
            'Ala': False,
            'Asc': True,
            'Asp': False,
            'Bet': False,
            'Ch': False,
            'Cit': False,
            'Cr': True,
            'EtOH': False,
            'GABA': True,
            'GABA_gov': False,
            'GABA_govind': False,
            'GPC': True,
            'GSH': True,
            'GSH_v2': False,
            'Glc': True,
            'Gln': True,
            'Glu': True,
            'Gly': True,
            'H2O': False,
            'Ins': True,
            'Lac': True,
            'NAA': True,
            'NAAG': True,
            'PCh': True,
            'PCr': True,
            'PE': True,
            'Phenyl': False,
            'Ref0ppm': False,
            'Scyllo': True,
            'Ser': False,
            'Tau': True,
            'Tau_govind': False,
            'Tyros': False,
            'bHB': False,
            'bHG': False,
        }

        # dropdown options
        self.dropdown = {
            'System': ['Philips', 'Siemens'],
            'Sequence': ['sLASER'],
            'Make .raw': ['Yes', 'No'],
        }

        # add file selection fields
        self.file_selection = ['Path to Pulse']

        # define dictionary of mandatory parameters
        self.mandatory_params = {
            "System": None,
            "Sequence": "sLASER",
            "Basis Name": "test.basis",
            "B1max": 22.,
            "Flip Angle": 180.,
            "RefTp": 4.5008,   # duration of the refocusing pulse
            "Samples": None,
            "Bandwidth": None,
            "Linewidth": 1.,
            "Bfield": None,

            # TODO: see that REMY gets the right values
            "thkX": 2.,  # in cm
            "thkY": 2.,
            "fovX": 3.,  # in cm   (if not found, default to +1 slice thickness)
            "fovY": 3.,

            "nX": 64.,
            "nY": 64.,

            "TE": None,
            "Center Freq": None,
            "Metabolites": [key for key, value in self.metabs.items() if value],

            "Tau 1": 15.,   # fake timing
            "Tau 2": 13.,

            "Path to Pulse": None,
            "Output Path": None,
            'Make .raw': 'No',   # TODO: fix io_writelcmraw in sLASER_makebasisset_function
        }


        # define dictionary of optional parameters
        self.optional_params = {
            'Nucleus': None,
            'TR': None,
        }

    def parseREMY(self, MRSinMRS):
        # extract as much information as possible from the MRSinMRS dict

        mandatory = {
            'System': self.parseSystem(MRSinMRS.get('Manufacturer', None)),
            'Sequence': self.parseProtocol(MRSinMRS.get('Protocol', None)),
            'B1max': None,  # TODO: find way to get from REMY? or set literature guided default
            'Flip Angle': MRSinMRS.get('ExcitationFlipAngle', None),  # TODO: find way to get from REMY? or set literature guided default
            'RefTp': None,   # duration of the refocusing pulse
            'Samples': MRSinMRS.get('NumberOfDatapoints', None),
            'Bandwidth': MRSinMRS.get('SpectralWidth', None),
            'Linewidth': None,   # TODO: find how to handle best...
            'Bfield': MRSinMRS.get('B0', None),

            'thkX': MRSinMRS.get('LeftRightSize', None),  # TODO: check for correctness
            'thkY': MRSinMRS.get('AnteriorPosteriorSize', None),

            'fovX': MRSinMRS.get('LeftRightSize', None),  # TODO: maybe get from VOI?
            'fovY': MRSinMRS.get('AnteriorPosteriorSize', None),

            'nX': None,
            'nY': None,

            'TE': MRSinMRS.get('TE', None),
            'Center Freq': MRSinMRS.get('Center Freq', None),

            'Tau 1': None,   # TODO: find way to get from REMY? or set literature guided default
            'Tau 2': None,
        }

        optional = {
            'Nucleus': MRSinMRS.get('Nucleus', None),
            'TR': MRSinMRS.get('TR', None),
            'Model': MRSinMRS.get('Model', None),
            'SoftwareVersion': MRSinMRS.get('SoftwareVersion', None),
            'BodyPart': MRSinMRS.get('BodyPart', None),
            'VOI': MRSinMRS.get('VOI', None),
            'CranioCaudalSize': MRSinMRS.get('CranioCaudalSize', None),
            'NumberOfAverages': MRSinMRS.get('NumberOfAverages', None),
            'WaterSuppression': MRSinMRS.get('WaterSuppression', None),
        }
        return mandatory, optional

    def parseProtocol(self, protocol):
        # backend only supports sLASER sequences for now

        if protocol is None:
            return None

        if 'mega' in protocol.lower():  # TODO: better check for editing
            print("Warning: LCModelBackend does not support MEGA sequences. "
                  "Ignoring MEGA part of the protocol.")

        if 'slaser' in protocol.lower():
            return 'sLASER'
        else:
            print("Warning: sLaserBackend only supports sLASER sequences. ")
            return None

    def parseSystem(self, system):
        # backend only supports Philips and Siemens systems for now
        if system is None:
            return None

        if 'philips' in system.lower():
            return 'Philips'
        elif 'siemens' in system.lower():
            return 'Siemens'
        else:
            print("Warning: sLaserBackend only supports Philips and Siemens systems. ")
            return None

    def parse2fidA(self, params):
        # convert parameters to fidA format if needed
        params['Make .basis'] = params['Make .basis'][0].lower()
        params['Make .raw'] = params['Make .raw'][0].lower()
        params['Display'] = params['Display'][0].lower()
        return params

    def setup_octave_paths(self):
        """Setup Octave paths for FID-A and sLASER toolboxes."""
        if self.octave is None:
            raise RuntimeError("Octave not initialized. Call initialize_octave() first.")

        self.octave.eval("warning('off', 'all');")
        self.octave.addpath('./externals/fidA/inputOutput/')
        self.octave.addpath('./externals/fidA/processingTools/')
        self.octave.addpath('./externals/fidA/simulationTools/')
        self.octave.addpath('./externals/jbss/')
        self.octave.addpath(self.octave.genpath('./adapters/'))

    def run_simulation(self, params, progress_callback=None):
        # Initialize Octave if not already done
        if self.octave is None:
            print("Initializing Octave runtime...")
            self.initialize_octave(prefer_docker=True)

        # Always setup paths (in case octave was initialized but paths weren't set)
        self.setup_octave_paths()

        # create the output directory if it does not exist
        if not os.path.exists(params['Output Path']):
            os.makedirs(params['Output Path'])

        # Convert output path to relative path for Docker compatibility
        output_path = params['Output Path']
        if os.path.isabs(output_path):
            # Convert absolute path to relative from current directory
            try:
                output_path = os.path.relpath(output_path)
            except ValueError:
                # If on different drive (Windows), keep absolute but remove drive letter
                output_path = output_path.replace('\\', '/')

        # Convert pulse path to relative path for Docker compatibility
        pulse_path = params['Path to Pulse']
        if pulse_path and os.path.isabs(pulse_path):
            # Convert absolute path to relative from current directory
            try:
                pulse_path = os.path.relpath(pulse_path)
            except ValueError:
                # If on different drive (Windows), keep absolute but remove drive letter
                pulse_path = pulse_path.replace('\\', '/')

        # Update params with converted pulse path
        params['Path to Pulse'] = pulse_path

        # fixed parameters
        params.update({
            "Curfolder": "./externals/jbss/",
            "Path to FIA-A": "./externals/fidA/",
            "Path to Spin System": "./externals/jbss/my_mets/my_spinSystem.mat",
            "Display": 'No',  # 0 (no display) <-> 1 (display)
            "Make .raw": 'No',  # TODO: fix io_writelcmraw in sLASER_makebasisset_function
            "Make .basis": 'Yes',
        })
        params = self.parse2fidA(params)   # convert parameters to fidA format if needed

        # define wrapper for octave function
        def sLASER_makebasisset_function(curfolder, pathtofida, system,
                                         seq_name, basis_name, B1max, flip_angle, refTp,
                                         Npts, sw, lw, Bfield, thkX, thkY, fovX, fovY, nX, nY, te,
                                         centreFreq, metab_list, tau1, tau2, path_to_pulse,
                                         path_to_save, path_to_spin_system, display,
                                         make_basis, make_raw):
            # Enable verbose output for Docker Octave to help debug issues
            verbose = hasattr(self.octave, 'verbose') and self.octave.verbose

            results = self.octave.feval('sLASER_makebasisset_function', curfolder, pathtofida,
                                        system, seq_name, basis_name, B1max, flip_angle, refTp,
                                        Npts, sw, lw, Bfield, thkX, thkY, fovX, fovY, nX, nY, te,
                                        centreFreq, metab_list, tau1, tau2, path_to_pulse,
                                        path_to_save, path_to_spin_system, display,
                                        make_basis, make_raw, verbose=verbose)
            return metab_list, results

        # prepare tasks for each metabolite - use converted output_path
        tasks = [(params['Curfolder'], params['Path to FIA-A'], params['System'],
                  params['Sequence'], params['Basis Name'], params['B1max'],
                  params['Flip Angle'], params['RefTp'], params['Samples'],
                  params['Bandwidth'], params['Linewidth'], params['Bfield'],
                  params['thkX'], params['thkY'], params['fovX'], params['fovY'],
                  params['nX'], params['nY'], params['TE'],
                  params['Center Freq'], params['Metabolites'], params['Tau 1'], params['Tau 2'],
                  params['Path to Pulse'], output_path,
                  params['Path to Spin System'], params['Display'], params['Make .basis'],
                  params['Make .raw'])]

        basis_set = {}
        total_steps = len(tasks)
        for task_idx, task in enumerate(tasks):
            metab_list, outputs = sLASER_makebasisset_function(*task)

            # Debug: Show what we got from Octave
            if hasattr(self.octave, 'verbose') and self.octave.verbose:
                print(f"\nDebug: Received from Octave:")
                print(f"  metab_list type: {type(metab_list)}, value: {metab_list}")
                print(f"  outputs type: {type(outputs)}")
                print(f"  outputs shape: {outputs.shape if hasattr(outputs, 'shape') else 'N/A'}")
                print(f"  outputs dtype: {outputs.dtype if hasattr(outputs, 'dtype') else 'N/A'}")
                if hasattr(outputs, '__len__'):
                    print(f"  outputs length: {len(outputs)}")

            # Handle metabolite list - might be a single item or array
            if not isinstance(metab_list, (list, tuple)):
                if hasattr(metab_list, '__iter__') and not isinstance(metab_list, str):
                    metab_list = list(metab_list)
                else:
                    metab_list = [metab_list]

            # outputs is a MATLAB cell array loaded by scipy
            # It might be 2D array, object array, or mat_struct depending on squeeze_me setting
            for met_idx, metab_name in enumerate(metab_list):
                try:
                    if hasattr(self.octave, 'verbose') and self.octave.verbose:
                        print(f"\n  Processing metabolite {met_idx}: {metab_name}")

                    # Access the output struct for this metabolite
                    # Check if it has ndim attribute (numpy array) or if it's a mat_struct
                    if hasattr(outputs, 'ndim'):
                        # It's a numpy array
                        if outputs.ndim == 2:
                            output_struct = outputs[0, met_idx]
                            if hasattr(self.octave, 'verbose') and self.octave.verbose:
                                print(f"    Accessed as outputs[0, {met_idx}] (2D array)")
                        elif outputs.ndim == 1:
                            output_struct = outputs[met_idx]
                            if hasattr(self.octave, 'verbose') and self.octave.verbose:
                                print(f"    Accessed as outputs[{met_idx}] (1D array)")
                        else:
                            output_struct = outputs
                            if hasattr(self.octave, 'verbose') and self.octave.verbose:
                                print(f"    Accessed as outputs (scalar)")
                    else:
                        # It's a mat_struct or similar object - treat as scalar
                        output_struct = outputs
                        if hasattr(self.octave, 'verbose') and self.octave.verbose:
                            print(f"    Accessed as outputs (mat_struct/scalar, single metabolite)")

                    if hasattr(self.octave, 'verbose') and self.octave.verbose:
                        print(f"    output_struct type: {type(output_struct)}")
                        print(f"    output_struct dtype: {output_struct.dtype if hasattr(output_struct, 'dtype') else 'N/A'}")
                        if hasattr(output_struct, '_fieldnames'):
                            print(f"    Available fields: {output_struct._fieldnames}")

                    # Extract fids - it's a struct field
                    if hasattr(output_struct, 'fids'):
                        fids_data = output_struct.fids
                        if hasattr(self.octave, 'verbose') and self.octave.verbose:
                            print(f"    ✓ Accessed via attribute: output_struct.fids")
                    elif isinstance(output_struct, dict):
                        fids_data = output_struct['fids']
                        if hasattr(self.octave, 'verbose') and self.octave.verbose:
                            print(f"    ✓ Accessed via dict: output_struct['fids']")
                    else:
                        # Try as item access
                        fids_data = output_struct['fids']
                        if hasattr(self.octave, 'verbose') and self.octave.verbose:
                            print(f"    ✓ Accessed via item: output_struct['fids']")

                    # Debug output
                    if hasattr(self.octave, 'verbose') and self.octave.verbose:
                        print(f"    fids_data type: {type(fids_data)}")
                        if hasattr(fids_data, 'shape'):
                            print(f"    fids_data shape: {fids_data.shape}")
                        if hasattr(fids_data, 'dtype'):
                            print(f"    fids_data dtype: {fids_data.dtype}")

                    # Convert to proper numpy array
                    if not isinstance(fids_data, np.ndarray):
                        if hasattr(fids_data, 'tolist'):
                            fids_data = np.array(fids_data, dtype=complex)
                        elif isinstance(fids_data, (list, tuple)):
                            fids_data = np.array(fids_data, dtype=complex)
                        else:
                            fids_data = np.asarray(fids_data, dtype=complex)

                    # Ensure dtype is complex
                    if fids_data.dtype not in (complex, np.complex64, np.complex128):
                        # If it's a structured array or object array, extract the actual data
                        if fids_data.dtype == object:
                            if hasattr(self.octave, 'verbose') and self.octave.verbose:
                                print(f"    ⚠️  fids_data is object dtype, unwrapping...")
                            # Try to extract first element if it's an object array wrapping the real data
                            if fids_data.size == 1:
                                fids_data = np.asarray(fids_data.item(), dtype=complex)
                                if hasattr(self.octave, 'verbose') and self.octave.verbose:
                                    print(f"    ✓ Unwrapped single object to shape: {fids_data.shape}")
                            else:
                                fids_data = np.array([complex(x) for x in fids_data.flat])
                                if hasattr(self.octave, 'verbose') and self.octave.verbose:
                                    print(f"    ✓ Converted object array to complex")
                        else:
                            fids_data = fids_data.astype(complex)

                    # Flatten if multidimensional
                    if fids_data.ndim > 1:
                        if hasattr(self.octave, 'verbose') and self.octave.verbose:
                            print(f"    Flattening from {fids_data.shape} to 1D")
                        fids_data = fids_data.flatten()

                    # Final validation
                    if fids_data.size == 0:
                        raise ValueError(f"Empty fids data for {metab_name}")

                    if hasattr(self.octave, 'verbose') and self.octave.verbose:
                        print(f"    ✓ Final fids_data: shape={fids_data.shape}, dtype={fids_data.dtype}")

                    basis_set[metab_name] = fids_data

                except Exception as e:
                    print(f"\n✗ Error extracting {metab_name} data: {e}")
                    print(f"  Output struct type: {type(output_struct)}")
                    if hasattr(output_struct, '__dict__'):
                        print(f"  Available attributes: {dir(output_struct)}")
                    import traceback
                    traceback.print_exc()
                    raise

            if progress_callback:
                progress_callback(task_idx + 1, total_steps)
        return basis_set

