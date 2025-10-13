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

        self.octave.addpath('./adapters/')

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
        }

        # add file selection fields
        self.file_selection = ['Path to Pulse']

        # define dictionary of mandatory parameters
        self.mandatory_params = {
            "System": "Philips",
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
        }


        # define dictionary of optional parameters
        self.optional_params = {
            'Nucleus': None,
            'TR': None,
        }

    def parseREMY(self, MRSinMRS):
        # extract as much information as possible from the MRSinMRS dict

        mandatory = {
            'System': MRSinMRS.get('Manufacturer', None),
            'Sequence': self.parseProtocol(MRSinMRS.get('Protocol', None)),
            'Samples': MRSinMRS.get('NumberOfDatapoints', None),
            'Bandwidth': MRSinMRS.get('SpectralWidth', None),
            'Bfield': MRSinMRS.get('B0', None),
            # 'Linewidth': 1,   # TODO: find how to handle or get from REMY
            'TE': MRSinMRS.get('TE', None),
            # 'TE2': 0,   # attention! - only holds for SpinEcho or STEAM
            #             # TODO: find sound solution!
            'Center Freq': MRSinMRS.get('Center Freq', None),
        }

        optional = {
            'Nucleus': MRSinMRS.get('Nucleus', None),
            'TR': MRSinMRS.get('TR', None),
            'Model': MRSinMRS.get('Model', None),
            'SoftwareVersion': MRSinMRS.get('SoftwareVersion', None),
            'BodyPart': MRSinMRS.get('BodyPart', None),
            'VOI': MRSinMRS.get('VOI', None),
            'AnteriorPosteriorSize': MRSinMRS.get('AnteriorPosteriorSize', None),
            'LeftRightSize': MRSinMRS.get('LeftRightSize', None),
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

    def run_simulation(self, params):
        pass

    def run_simulation_with_progress(self, params, progress_callback):
        # create the output directory if it does not exist
        if not os.path.exists(params['Output Path']): os.makedirs(params['Output Path'])

        # fixed parameters
        params.update({
            "Curfolder": os.getcwd() + '/jbss/',
            "Path to FIA-A": os.getcwd() + "/fidA/",
            "Path to Spin System": os.getcwd() + "/jbss/my_mets/",
            "Display": False,
        })

        def sLASER_makebasisset_function(curfolder, pathtofida, system,
                                         seq_name, basis_name, B1max, flip_angle, refTp,
                                         Npts, sw,lw, Bfield, thkX, thkY, fovX, fovY, nX, nY, te,
                                         centreFreq, spinSysList, tau1, tau2, path_to_pulse,
                                         path_to_save, path_to_spin_system, display):
            results = self.octave.feval('sLASER_makebasisset_function', curfolder, pathtofida,
                                        system, seq_name, basis_name, B1max, flip_angle, refTp,
                                        Npts, sw, lw, Bfield, thkX, thkY, fovX, fovY, nX, nY, te,
                                        centreFreq, spinSysList, tau1, tau2, path_to_pulse,
                                        path_to_save, path_to_spin_system, display)
            return metab, results[:, 0] + 1j * results[:, 1]

        # TODO: make sure all parameters are cast properly (int, float, str)
        #       e.g. see LCModelBackend for reference
        tasks = [(params['Curfolder'], params['Path to FIA-A'], params['System'],
                  params['Sequence'], params['Basis Name'], params['B1max'],
                  params['Flip Angle'], params['RefTp'], params['Samples'],
                  params['Bandwidth'], params['Linewidth'], params['Bfield'],
                  params['thkX'], params['thkY'], params['fovX'], params['fovY'],
                  params['nX'], params['nY'], params['TE'],
                  params['Center Freq'], [metab], params['Tau 1'], params['Tau 2'],
                  params['Path to Pulse'], params['Output Path'],
                  params['Path to Spin System'], params['Display'])
                 for metab in params['Metabolites']]

        # initialize the progress bar
        total_steps = len(tasks)
        progress_step = 100 / total_steps

        # run simulations sequentially
        basis_set = {}
        for i, task in enumerate(tasks):
            metab, data = sLASER_makebasisset_function(*task)
            basis_set[metab] = data
            progress_callback(i + 1, total_steps)   # update the progress bar
        return basis_set
