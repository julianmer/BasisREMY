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

from oct2py import Oct2Py

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

        # init fidA
        self.octave = Oct2Py()
        self.octave.eval("warning('off', 'all');")
        self.octave.addpath('./externals/fidA/inputOutput/')
        self.octave.addpath('./externals/fidA/processingTools/')
        self.octave.addpath('./externals/fidA/simulationTools/')

        self.octave.addpath('./externals/jbss/')

        self.octave.addpath(self.octave.genpath('./adapters/'))

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

    def run_simulation(self, params, progress_callback=None):
        # create the output directory if it does not exist
        if not os.path.exists(params['Output Path']):
            os.makedirs(params['Output Path'])

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
            results = self.octave.feval('sLASER_makebasisset_function', curfolder, pathtofida,
                                        system, seq_name, basis_name, B1max, flip_angle, refTp,
                                        Npts, sw, lw, Bfield, thkX, thkY, fovX, fovY, nX, nY, te,
                                        centreFreq, metab_list, tau1, tau2, path_to_pulse,
                                        path_to_save, path_to_spin_system, display,
                                        make_basis, make_raw)
            return metab_list, results

        # prepare tasks for each metabolite
        tasks = [(params['Curfolder'], params['Path to FIA-A'], params['System'],
                  params['Sequence'], params['Basis Name'], params['B1max'],
                  params['Flip Angle'], params['RefTp'], params['Samples'],
                  params['Bandwidth'], params['Linewidth'], params['Bfield'],
                  params['thkX'], params['thkY'], params['fovX'], params['fovY'],
                  params['nX'], params['nY'], params['TE'],
                  params['Center Freq'], params['Metabolites'], params['Tau 1'], params['Tau 2'],
                  params['Path to Pulse'], params['Output Path'],
                  params['Path to Spin System'], params['Display'], params['Make .basis'],
                  params['Make .raw'])]

        basis_set = {}
        total_steps = len(tasks)
        for task_idx, task in enumerate(tasks):
            metab_list, outputs = sLASER_makebasisset_function(*task)
            # outputs is a MATLAB cell array; convert to list and access each cell
            for met_idx, metab_name in enumerate(metab_list):
                # MATLAB cell arrays are 0-indexed when accessed from Python via oct2py
                output_struct = outputs[0, met_idx]
                basis_set[metab_name] = output_struct['fids']
            if progress_callback:
                progress_callback(task_idx + 1, total_steps)
        return basis_set
