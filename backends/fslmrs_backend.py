####################################################################################################
#                                      fslmrs_backend.py                                           #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 18/02/26                                                                                #
#                                                                                                  #
# Purpose: FSL-MRS backend for basis set simulation using density matrix methods.                  #
#          Uses FSL-MRS Python package (no Octave required).                                       #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import os
import sys
import json
import numpy as np
from backends.base import Backend

# denmatsim lives at externals/fsl_mrs/fsl_mrs/denmatsim/
# add its parent to sys.path so 'from denmatsim import ...' works
_denmatsim_parent = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', 'externals', 'fsl_mrs', 'fsl_mrs')
)
if _denmatsim_parent not in sys.path:
    sys.path.insert(0, _denmatsim_parent)


#**************************************************************************************************#
#                                          FSL-MRS Backend                                         #
#**************************************************************************************************#
#                                                                                                  #
# Implements basis set simulation using FSL-MRS density matrix simulations.                        #
# Advantages:                                                                                      #
#   - Pure Python (no Octave/MATLAB required)                                                      #
#   - Quantum-mechanically accurate simulations                                                    #
#   - Supports custom pulse sequences via JSON                                                     #
#   - Parallel processing built-in                                                                 #
#                                                                                                  #
#**************************************************************************************************#
class FSLMRSBackend(Backend):
    def __init__(self):
        super().__init__()

        self.name = 'FSL-MRS'
        self.requires_octave = False  # Pure Python!

        # Mode support (overrides base class)
        self.modes = ['Simple', 'Template', 'Custom']
        self.current_mode = 'Simple'

        # Supported sequences
        # - PRESS, STEAM: Have both ideal and template options
        # - Others: Ideal pulses only (no templates yet)
        # - Custom: User provides complete file
        self.supported_sequences = [
            'PRESS', 'STEAM', 'LASER', 'sLASER',
            'MEGA-PRESS', 'HERMES', 'HERCULES', 'MEGA-sLASER',
            'Custom'
        ]

        # Predefined sequence files available in externals/fsl_mrs
        # These contain REAL pulse shapes for specific parameters!
        # Using them with different parameters may be inaccurate
        self.predefined_sequences = {
            'PRESS_7T': {
                'path': 'examplePRESS.json',
                'description': 'PRESS at 7T with real pulse shapes',
                'B0': 7.0,
                'TE': None,
                'Rx_Points': 4096,
                'Rx_SW': 6000,
                'notes': 'Real Siemens 7T PRESS pulse shapes'
            },
            'STEAM_7T_11ms': {
                'path': 'example.json',
                'description': '11ms STEAM at 7T with real pulse shapes',
                'B0': 6.98,
                'TE': 11,
                'Rx_Points': 8192,
                'Rx_SW': 6000,
                'notes': 'FMRIB 7T STEAM sequence with real pulse shapes'
            },
        }

        # Simplified mapping for sequence selection
        # Maps generic sequence names to specific predefined files
        self.sequence_to_predefined = {
            'PRESS': 'PRESS_7T',
            'STEAM': 'STEAM_7T_11ms',
        }

        # Default metabolites
        # TODO: Add macromolecules sim parts
        self.default_metabolites = [
            'Ala', 'Asc', 'Asp', 'Cr', 'GABA', 'GPC', 'GSH', 'Glc', 'Gln', 'Glu', 'Gly',
            'Ins', 'Lac', 'NAA', 'NAAG', 'PCh', 'PCr', 'PE', 'Scyllo', 'Tau',
            'Ace', 'Acn', 'Ala', 'Asc', 'Asp', 'Bet', 'bHB', 'bHG', 'Cit', 'Cr',
            'Cystat', 'EtOH', 'GABA', 'Glc', 'Gln', 'Glu', 'Gly', 'GPC', 'GSH',
            'Gua', 'HCar', 'Hist', 'Hypotau', 'Ins', 'Lac', 'Lip09', 'Lip13a',
            'Lip13b', 'Lip20', 'NAA', 'NAAG', 'Oac', 'PCh', 'PCr', 'PE', 'Phenyl',
            'Pyruv', 'Scyllo', 'Ser', 'Suc', 'Tau', 'Thr', 'Tyros'
        ]

        # Metabolite selection for UI
        self.metabs = {name: False for name in self.default_metabolites}

        # Dropdown options (Mode NOT here - handled by base class)
        self.dropdown = {
            'Sequence': [
                'PRESS', 'STEAM', 'LASER', 'sLASER',
                'MEGA-PRESS', 'HERMES', 'HERCULES', 'MEGA-sLASER'
            ],
            'Template File': [info['description'] for info in self.predefined_sequences.values()],
            'Output Format': ['LCModel RAW', 'JSON'],
        }

        # File selection fields
        self.file_selection = ['Custom Sequence']

        # Mandatory parameters (Mode NOT here - handled by base class)
        self.mandatory_params = {
            'Sequence': 'PRESS',
            'Samples': 2048,
            'Bandwidth': 2000,
            'Bfield': 3.0,
            'TE': 35,
            'Nucleus': '1H',
            'Center Freq': 123.2,
            'Output Path': './output',
            'Metabolites': [],
            'Output Format': 'LCModel RAW',
        }

        # Optional parameters
        self.optional_params = {
            'TM': 10,
            'Template File': None,
            'Edit Frequency': 1.9,
            'Linewidth': 2.0,
            'Add Reference': False,
            'Auto Phase': None,
            'Parallel': True,
            'Custom Sequence': None,
        }

    def get_params_for_mode(self, mode=None):
        """
        Return parameters to display in GUI for the given mode.
        Overrides base class.
        """
        if mode is None:
            mode = self.current_mode

        # Common parameters shown in every mode
        common = {
            'Output Path': self.mandatory_params['Output Path'],
            'Metabolites': self.mandatory_params['Metabolites'],
            'Output Format': self.mandatory_params['Output Format'],
        }

        if mode == 'Simple':
            return {
                'Sequence': self.mandatory_params['Sequence'],
                'Bfield': self.mandatory_params['Bfield'],
                'TE': self.mandatory_params['TE'],
                'TM': self.optional_params['TM'],
                'Samples': self.mandatory_params['Samples'],
                'Bandwidth': self.mandatory_params['Bandwidth'],
                'Edit Frequency': self.optional_params['Edit Frequency'],
                **common,
            }

        elif mode == 'Template':
            return {
                'Template File': self.optional_params['Template File'],
                'Samples': self.mandatory_params['Samples'],
                'Bandwidth': self.mandatory_params['Bandwidth'],
                'Linewidth': self.optional_params['Linewidth'],
                **common,
            }

        elif mode == 'Custom':
            return {
                'Custom Sequence': self.optional_params['Custom Sequence'],
                'Samples': self.mandatory_params['Samples'],
                'Bandwidth': self.mandatory_params['Bandwidth'],
                **common,
            }

        return dict(self.mandatory_params)

    def parseREMY(self, MRSinMRS):
        """
        Parse REMY output to FSL-MRS backend parameters

        Args:
            MRSinMRS: Dictionary of parameters extracted by REMY

        Returns:
            tuple: (mandatory_params_dict, optional_params_dict)
        """
        params = {}
        opt = {}

        # Required parameters
        params['Samples'] = MRSinMRS.get('NumberOfDatapoints', MRSinMRS.get('Samples', 2048))
        params['Bandwidth'] = MRSinMRS.get('SpectralWidth', MRSinMRS.get('Bandwidth', 2000))
        params['Bfield'] = MRSinMRS.get('B0', MRSinMRS.get('Bfield', 3.0))
        params['TE'] = MRSinMRS.get('TE', 35)
        params['Nucleus'] = MRSinMRS.get('Nucleus', '1H')

        # Calculate center frequency from field strength if not provided
        center_freq = MRSinMRS.get('Center Freq', None)
        if center_freq is None:
            # Calculate for 1H at given field strength
            # gamma_1H = 42.577 MHz/T
            params['Center Freq'] = 42.577 * float(params['Bfield'])
        else:
            params['Center Freq'] = center_freq

        # Sequence detection
        sequence = self.parseProtocol(MRSinMRS.get('Protocol', ''))
        if sequence and sequence in self.supported_sequences:
            params['Sequence'] = sequence
        else:
            params['Sequence'] = None  # Will need manual selection

        params['Output Path'] = './output'
        params['Metabolites'] = []
        params['Output Format'] = 'LCModel RAW'

        # Optional parameters
        opt['Linewidth'] = 2.0
        opt['Add Reference'] = False
        opt['Auto Phase'] = None
        opt['Parallel'] = True
        opt['Custom Sequence'] = None

        return params, opt

    def show_predefined_sequences(self):
        """
        Display information about available predefined sequence files
        """
        print("\n" + "="*80)
        print("Available Predefined Sequences (with REAL pulse shapes)")
        print("="*80)

        for key, info in self.predefined_sequences.items():
            print(f"\n{key}:")
            print(f"  Description: {info['description']}")
            print(f"  Field Strength: {info['B0']} T")
            print(f"  TE: {info['TE']} ms" if info['TE'] else "  TE: Depends on delays")
            print(f"  Acquisition: {info['Rx_Points']} points @ {info['Rx_SW']} Hz")
            print(f"  Notes: {info['notes']}")
            print(f"  File: {info['path']}")

        print("\n" + "="*80)
        print("⚠️  These files contain real RF pulse waveforms!")
        print("   Using them with significantly different parameters may be inaccurate.")
        print("   For other field strengths/TEs, provide a custom sequence JSON file.")
        print("="*80 + "\n")

    def parseProtocol(self, protocol):
        """
        Parse sequence name from protocol string

        Args:
            protocol: Protocol string from scanner

        Returns:
            str: Standardized sequence name or None
        """
        protocol_upper = protocol.upper()

        if 'MEGA' in protocol_upper:
            return 'MEGA-PRESS'
        elif 'HERMES' in protocol_upper:
            return 'HERMES'
        elif 'PRESS' in protocol_upper:
            return 'PRESS'
        elif 'STEAM' in protocol_upper:
            return 'STEAM'
        elif 'SLASER' in protocol_upper:
            return 'sLASER'
        elif 'LASER' in protocol_upper:
            return 'LASER'

        return None

    def _generate_sequence_json(self, params):
        """
        Generate FSL-MRS sequence JSON with IDEAL PULSES (FID-A style)

        Uses instantaneous rotation operators (~10 μs) for all pulses.
        These are mathematically rigorous and standard in NMR simulation.

        For accurate simulations with real pulse shapes:
          - Use Template Mode (predefined files)
          - Use Custom Mode (provide your own JSON)

        Args:
            params: Dictionary of simulation parameters

        Returns:
            dict: Sequence definition in FSL-MRS JSON format with ideal pulses
        """
        sequence = params['Sequence']
        te = params['TE']
        bandwidth = params['Bandwidth']
        samples = params['Samples']
        bfield = params['Bfield']

        print(f"Generating IDEAL pulse sequence: {sequence}")
        print("  Using instantaneous rotation operators (FID-A style)")
        print("  Perfect flip angles, no realistic pulse effects")
        print("  For accurate simulations, use Template or Custom mode")

        # Base sequence structure for denmatsim
        seq_def = {
            'sequenceName': f'{sequence}_ideal',
            'description': f'Ideal {sequence} with instantaneous pulses',
            'B0': bfield,
            'centralShift': 4.65,  # ppm - typical for 1H MRS
            'Rx_Points': samples,
            'Rx_SW': bandwidth,
            'Rx_LW': 2.0,
            'Rx_Phase': 0.0,
            'x': [-15, 15],
            'y': [-15, 15],
            'z': [-15, 15],
            'resolution': [8, 8, 8],
            'RFUnits': 'Hz',
            'GradUnits': 'mT',
            'spaceUnits': 'mm',
        }

        # Ideal pulse duration (essentially instantaneous)
        ideal_pulse_duration = 0.00001  # 10 microseconds

        # For ideal pulses: flip_angle = 2π * B1_Hz * duration
        # So B1_Hz = flip_angle_rad / (2π * duration)
        import math
        amp_90 = (math.pi / 2) / (2 * math.pi * ideal_pulse_duration)   # ~25000 Hz
        amp_180 = math.pi / (2 * math.pi * ideal_pulse_duration)         # ~50000 Hz

        # Define sequences
        if sequence == 'PRESS':
            # PRESS: 90° - delay - 180° - delay - 180° - delay - ACQ
            # 3 RF pulses = 3 delays, 3 rephaseAreas, 3 CoherenceFilter
            te_half = te / 2000  # convert ms to seconds
            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                ],
                'delays': [te_half, te_half, te_half],
                'rephaseAreas': [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                'CoherenceFilter': [-1, 1, -1],
            })

        elif sequence == 'STEAM':
            # STEAM: 90° - TE/2 - 90° - TM - 90° - TE/2 - ACQ
            # 3 RF pulses = 3 delays, 3 rephaseAreas, 3 CoherenceFilter
            tm = params.get('TM', 10)
            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [3.14159], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                ],
                'delays': [te/2000, tm/1000, te/2000],
                'rephaseAreas': [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                'CoherenceFilter': [1, 0, -1],
                'Rx_Phase': 1.5708,
            })

        elif sequence == 'LASER':
            # LASER: 90° + 3 pairs of 180° AFP = 7 RF pulses
            # 7 RF = 7 delays, 7 rephaseAreas, 7 CoherenceFilter
            d = te / 6000  # TE/6 in seconds
            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                ],
                'delays': [d, d, d, d, d, d, d],
                'rephaseAreas': [[0, 0, 0]] * 7,
                'CoherenceFilter': [-1, 1, -1, 1, -1, 1, -1],
            })

        elif sequence == 'sLASER':
            # sLASER: 90° + 2 pairs of 180° AFP = 5 RF pulses
            # 5 RF = 5 delays, 5 rephaseAreas, 5 CoherenceFilter
            d = te / 4000  # TE/4 in seconds
            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                ],
                'delays': [d, d, d, d, d],
                'rephaseAreas': [[0, 0, 0]] * 5,
                'CoherenceFilter': [-1, 1, -1, 1, -1],
            })

        elif sequence == 'MEGA-PRESS':
            # MEGA-PRESS with ideal editing pulses
            # Default Siemens timing
            edit_freq = params.get('Edit_Frequency', 1.9)  # ppm (for GABA)
            edit_freq_hz = edit_freq * bfield * 42.577  # Convert to Hz

            # Siemens timing (ms)
            t1 = 4.545
            t2 = 12.7025
            t3 = 21.7975
            t4 = 12.7025
            t5 = 17.2526

            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},  # 90° excite
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},  # 180° refocus
                    {'time': ideal_pulse_duration, 'frequencyOffset': edit_freq_hz, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [0], 'grad': [0, 0, 0]},  # EDIT @ 1.9ppm
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},  # 180° refocus
                    {'time': ideal_pulse_duration, 'frequencyOffset': edit_freq_hz, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [0], 'grad': [0, 0, 0]},  # EDIT @ 1.9ppm
                ],
                'delays': [t1/1000, t2/1000, t3/1000, t4/1000, t5/1000],
                'rephaseAreas': [[0, 0, 0]] * 5,
                'CoherenceFilter': [-1, 1, -1, 1, -1],
            })
            print(f"  MEGA-PRESS editing frequency: {edit_freq} ppm ({edit_freq_hz:.1f} Hz)")

        elif sequence == 'HERMES':
            # HERMES - will need multiple sub-spectra
            # This generates one sub-spectrum (edit both GABA and GSH)
            gaba_freq_hz = 1.9 * bfield * 42.577
            gsh_freq_hz = 4.56 * bfield * 42.577

            # Use MEGA-PRESS timing
            t1, t2, t3, t4, t5 = 4.545, 12.7025, 21.7975, 12.7025, 17.2526

            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': gaba_freq_hz, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [0], 'grad': [0, 0, 0]},  # Edit GABA
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': gsh_freq_hz, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [0], 'grad': [0, 0, 0]},  # Edit GSH
                ],
                'delays': [t1/1000, t2/1000, t3/1000, t4/1000, t5/1000],
                'rephaseAreas': [[0, 0, 0]] * 5,
                'CoherenceFilter': [-1, 1, -1, 1, -1],
            })
            print(f"  HERMES editing: GABA @ 1.9 ppm, GSH @ 4.56 ppm")

        elif sequence == 'HERCULES':
            # HERCULES - extension of HERMES
            print(f"  HERCULES: Multi-metabolite editing")
            print(f"  Using HERMES framework with optimized frequencies")
            # Similar to HERMES but with additional editing targets
            gaba_freq_hz = 1.9 * bfield * 42.577
            glu_freq_hz = 2.3 * bfield * 42.577

            t1, t2, t3, t4, t5 = 4.545, 12.7025, 21.7975, 12.7025, 17.2526

            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': gaba_freq_hz, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [0], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': glu_freq_hz, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [0], 'grad': [0, 0, 0]},
                ],
                'delays': [t1/1000, t2/1000, t3/1000, t4/1000, t5/1000],
                'rephaseAreas': [[0, 0, 0]] * 5,
                'CoherenceFilter': [-1, 1, -1, 1, -1],
            })

        elif sequence == 'MEGA-sLASER':
            # MEGA-sLASER: Combination of MEGA editing with sLASER localization
            edit_freq = params.get('Edit_Frequency', 1.9)
            edit_freq_hz = edit_freq * bfield * 42.577

            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},  # 90° excite
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},  # AFP pair 1a
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},  # AFP pair 1b
                    {'time': ideal_pulse_duration, 'frequencyOffset': edit_freq_hz, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [0], 'grad': [0, 0, 0]},  # EDIT
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},  # AFP pair 2a
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},  # AFP pair 2b
                    {'time': ideal_pulse_duration, 'frequencyOffset': edit_freq_hz, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [0], 'grad': [0, 0, 0]},  # EDIT
                ],
                'delays': [te/7000] * 7,
                'rephaseAreas': [[0, 0, 0]] * 7,
                'CoherenceFilter': [-1, 1, -1, 1, -1, 1, -1],
            })
            print(f"  MEGA-sLASER: sLASER localization + MEGA editing @ {edit_freq} ppm")

        else:
            # Generic single pulse excitation for unknown sequences
            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                ],
                'delays': [te / 1000.0],
                'rephaseAreas': [[0, 0, 0]],
                'CoherenceFilter': [-1],
            })

        return seq_def

    def run_simulation(self, params, progress_callback=None):
        """
        Run FSL-MRS basis set simulation using denmatsim

        Args:
            params: Dictionary of simulation parameters
            progress_callback: Optional callback function(current, total)

        Returns:
            dict: {metabolite_name: FID_array}
        """
        # Import denmatsim (path already set up at module level)
        try:
            from denmatsim import simseq, utils as simutils
        except ImportError as e:
            raise RuntimeError(
                f"denmatsim not found at {_denmatsim_parent}/denmatsim/\n"
                f"Error: {e}\n\n"
                f"Run: git submodule update --init --recursive"
            )
        print("✓ denmatsim imported successfully")

        print(f"\n{'='*80}")
        print(f"Running FSL-MRS simulation (denmatsim)")
        print(f"{'='*80}")
        print(f"Sequence: {params['Sequence']}")
        print(f"TE: {params['TE']} ms")
        print(f"Field strength: {params['Bfield']} T")
        print(f"Metabolites: {len(params['Metabolites'])}")

        # Create output directory
        output_path = params['Output Path']
        os.makedirs(output_path, exist_ok=True)

        # Get or generate sequence JSON
        if params.get('Custom Sequence') and os.path.exists(params['Custom Sequence']):
            # User provided custom sequence file
            print(f"Using custom sequence: {params['Custom Sequence']}")
            with open(params['Custom Sequence'], 'r') as f:
                seq_params = json.load(f)

        elif params['Sequence'] in self.sequence_to_predefined:
            # Try to use predefined sequence file
            predefined_key = self.sequence_to_predefined[params['Sequence']]
            predefined_info = self.predefined_sequences[predefined_key]

            # Check if user's parameters match the predefined file
            param_match = True
            warnings = []

            # Check B0
            if abs(params['Bfield'] - predefined_info['B0']) > 0.5:
                param_match = False
                warnings.append(f"Field strength mismatch: requested {params['Bfield']}T, predefined is {predefined_info['B0']}T")

            # Check TE (if specified in predefined)
            if predefined_info['TE'] is not None:
                if abs(params['TE'] - predefined_info['TE']) > 5:  # Allow 5ms tolerance
                    param_match = False
                    warnings.append(f"TE mismatch: requested {params['TE']}ms, predefined is {predefined_info['TE']}ms")

            # If parameters match well enough, use predefined file
            if param_match:
                seq_rel_path = predefined_info['path']
                # JSON example files live inside denmatsim itself
                denmatsim_path = os.path.abspath(os.path.join(
                    os.path.dirname(__file__), '..', 'externals', 'fsl_mrs', 'fsl_mrs', 'denmatsim'
                ))
                seq_file_path = os.path.join(denmatsim_path, seq_rel_path)

                if os.path.exists(seq_file_path):
                    print(f"✓ Using predefined {params['Sequence']} sequence: {predefined_info['description']}")
                    print(f"  File: {seq_rel_path}")
                    print(f"  Parameters: B0={predefined_info['B0']}T, TE={predefined_info['TE']}ms, Points={predefined_info['Rx_Points']}, BW={predefined_info['Rx_SW']}Hz")
                    print(f"  Note: This file contains REAL pulse shapes - highly accurate!")

                    with open(seq_file_path, 'r') as f:
                        seq_params = json.load(f)

                    # Update ONLY acquisition parameters (not pulse shapes!)
                    # User can override Rx_Points and Rx_SW for their specific needs
                    if params['Samples'] != predefined_info['Rx_Points']:
                        print(f"  ⚠️  Updating Rx_Points from {predefined_info['Rx_Points']} to {params['Samples']}")
                        seq_params['Rx_Points'] = params['Samples']

                    if params['Bandwidth'] != predefined_info['Rx_SW']:
                        print(f"  ⚠️  Updating Rx_SW from {predefined_info['Rx_SW']} to {params['Bandwidth']}")
                        seq_params['Rx_SW'] = params['Bandwidth']

                else:
                    raise RuntimeError(
                        f"Predefined sequence file not found: {seq_file_path}\n"
                        f"Make sure FSL-MRS submodule is initialized:\n"
                        f"  git submodule update --init --recursive"
                    )
            else:
                # Parameters don't match - warn and use idealized sequence
                print(f"\n{'='*80}")
                print(f"⚠️  WARNING: Parameter mismatch with predefined {params['Sequence']} sequence!")
                print(f"{'='*80}")
                print(f"Predefined file: {predefined_info['description']}")
                print(f"  Fixed parameters: B0={predefined_info['B0']}T, TE={predefined_info['TE']}ms")
                print(f"\nYour parameters:")
                for warning in warnings:
                    print(f"  ⚠️  {warning}")
                print(f"\n⚠️  Using predefined file with different parameters would be INACCURATE!")
                print(f"     (Real pulse shapes are field-strength and TE dependent)")
                print(f"\nGenerating IDEALIZED sequence instead (perfect pulses, no realistic effects)")
                print(f"For accurate simulations at your parameters, provide a custom sequence file.")
                print(f"{'='*80}\n")

                seq_params = self._generate_sequence_json(params)

        else:
            # Generate idealized sequence for 'Custom' or unknown sequences
            print(f"⚠️  No predefined file for '{params['Sequence']}' - generating idealized sequence")
            print("   For accurate simulations, use 'Custom Sequence' parameter with your own JSON file")
            seq_params = self._generate_sequence_json(params)

        # Save sequence file for reference
        seq_file = os.path.join(output_path, f'{params["Sequence"]}_sequence.json')
        with open(seq_file, 'w') as f:
            json.dump(seq_params, f, indent=2)
        print(f"Saved sequence file: {seq_file}")

        # Load spin systems for metabolites
        try:
            spinSystems = simutils.readBuiltInSpins()
            print(f"✓ Loaded {len(spinSystems)} built-in spin systems")
        except Exception as e:
            print(f"⚠️  Could not load spin systems: {e}")
            print("  Returning placeholder basis set")
            # Return placeholder
            basis_set = {}
            for metab in params['Metabolites']:
                basis_set[metab] = np.zeros(int(params['Samples']), dtype=complex)
            return basis_set

        # Run simulation for each metabolite
        basis_set = {}
        total_metabs = len(params['Metabolites'])

        for idx, metab in enumerate(params['Metabolites'], 1):
            if progress_callback:
                progress_callback(idx, total_metabs)

            print(f"\n[{idx}/{total_metabs}] Simulating {metab}...")

            # Get spin system
            sys_name = f'sys{metab}'
            if sys_name not in spinSystems:
                print(f"  ⚠️  Spin system '{sys_name}' not found, skipping")
                continue

            spin_system = spinSystems[sys_name]

            try:
                # denmatsim spin systems are lists of sub-spin-systems
                # (e.g. NAA has acetyl + aspartate groups)
                # simulate each sub-system and sum the FIDs
                if isinstance(spin_system, list):
                    FID = None
                    for sub_idx, sub_sys in enumerate(spin_system):
                        scale = sub_sys.get('scaleFactor', 1.0)
                        sub_fid, ax, pmat = simseq.simseq(sub_sys, seq_params, verbose=False)
                        if FID is None:
                            FID = sub_fid * scale
                        else:
                            FID += sub_fid * scale
                    print(f"  ✓ Simulated {len(spin_system)} sub-systems")
                else:
                    FID, ax, pmat = simseq.simseq(spin_system, seq_params, verbose=False)

                # denmatsim returns FID with negative ppm convention and an
                # arbitrary zero-order phase from the TE evolution.
                # Fix both:
                #   1. Conjugate → flips the frequency axis to the standard
                #      MRS convention (low-field resonances at positive ppm)
                #   2. Zero-order phase correction using FID[0] angle →
                #      puts the absorptive signal in the real channel
                FID = np.conj(FID)
                phi0 = np.angle(FID[0])
                FID = FID * np.exp(-1j * phi0)

                # Store FID
                basis_set[metab] = FID
                print(f"  ✓ Generated FID with {len(FID)} points")

                # Optionally save individual metabolite files
                if params.get('Output Format') == 'JSON':
                    metab_file = os.path.join(output_path, f'{metab}.json')
                    metab_data = {
                        'FID': FID.tolist() if hasattr(FID, 'tolist') else list(FID),
                        'sequence': seq_params,
                        'metabolite': metab,
                    }
                    with open(metab_file, 'w') as f:
                        json.dump(metab_data, f, indent=2)
                    print(f"  Saved: {metab_file}")

                elif params.get('Output Format') == 'LCModel RAW':
                    # Save as LCModel RAW format
                    raw_file = os.path.join(output_path, f'{metab}.RAW')
                    self._save_lcmodel_raw(FID, raw_file, params)
                    print(f"  Saved: {raw_file}")

            except Exception as e:
                print(f"  ✗ Simulation failed: {e}")
                import traceback
                traceback.print_exc()
                continue

        print(f"\n{'='*80}")
        print(f"Simulation complete!")
        print(f"Generated {len(basis_set)}/{total_metabs} metabolite spectra")
        print(f"Output location: {output_path}")
        print(f"{'='*80}\n")

        return basis_set

    def _save_lcmodel_raw(self, fid, filepath, params):
        """Save FID in LCModel RAW format"""
        with open(filepath, 'w') as f:
            # LCModel RAW header
            f.write(f" $NMID\n")
            f.write(f" ID='BasisREMY FSL-MRS simulation'\n")
            f.write(f" FMTDAT='(2E15.6)'\n")
            f.write(f" VOLUME=1.0\n")
            f.write(f" TRAMP=1.0\n")
            f.write(f" $END\n")

            # Write FID data (real, imag pairs)
            for point in fid:
                f.write(f" {point.real:15.6E} {point.imag:15.6E}\n")

