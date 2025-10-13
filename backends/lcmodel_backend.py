####################################################################################################
#                                        lcmodel_backend.py                                        #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 08/10/25                                                                                #
#                                                                                                  #
# Purpose: Defines the LCModelBackend class for simulating MRS basis sets using the                #
#          FID-A sim_lcmrawbasis.m function.                                                       #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import os

from multiprocessing import Pool
from oct2py import Oct2Py

# own
from backends.base import Backend


#**************************************************************************************************#
#                                          LCModelBackend                                          #
#**************************************************************************************************#
#                                                                                                  #
# Implements the basis set simulation backend using the FID-A sim_lcmrawbasis.m function. A very   #
# simplified simulation for SE, PRESS, STEAM, and LASER sequences.                                 #
#                                                                                                  #
#**************************************************************************************************#
class LCModelBackend(Backend):
    def __init__(self):
        super().__init__()
        self.name = 'LCModel'

        # init fidA
        self.octave = Oct2Py()
        self.octave.eval("warning('off', 'all');")
        self.octave.addpath('./externals/fidA/inputOutput/')
        self.octave.addpath('./externals/fidA/processingTools/')
        self.octave.addpath('./externals/fidA/simulationTools/')

        # define possible metabolites
        self.metabs = {
            'Ala': False,
            'Asc': True,
            'Asp': False,
            'Ch': False,
            'Cit': False,
            'Cr': True,
            'EtOH': False,
            'GABA': True,
            'GPC': True,
            'GSH': True,
            'Glc': True,
            'Gln': True,
            'Glu': True,
            'Gly': True,
            'H2O': False,
            'Ins': True,
            'Lac': True,
            'Lip': False,
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
            'Tyros': False,
        }

        # dropdown options
        self.dropdown = {
            'Sequence': ['Spin Echo', 'PRESS', 'STEAM', 'LASER'],
            # 'se' for Spin Echo, 'p' for Press, 'st' for Steam, or 'l' for LASER
            'Add Ref.': ['Yes', 'No'],
            'Make .raw': ['Yes', 'No'],
        }

        # define dictionary of mandatory parameters
        self.mandatory_params = {
            'Sequence': None,
            'Samples': None,
            'Bandwidth': None,
            'Bfield': None,
            'Linewidth': 1,
            'TE': None,
            'TE2': 0,
            'Add Ref.': 'No',  # default to 'No'
            'Make .raw': 'Yes',  # default to 'Yes' (need for .m script to run properly)
            'Output Path': None,
            'Metabolites': [key for key, value in self.metabs.items() if value],
            'Center Freq': None,   # currently used for plotting on the ppm scale
        }

        # define dictionary of optional parameters
        self.optional_params = {
            'Nucleus': None,
            'TR': None,
        }

    def parseREMY(self, MRSinMRS):
        # extract as much information as possible from the MRSinMRS dict

        mandatory = {
            'Sequence': self.parseProtocol(MRSinMRS.get('Protocol', None)),
            'Samples': MRSinMRS.get('NumberOfDatapoints', None),
            'Bandwidth': MRSinMRS.get('SpectralWidth', None),
            'Bfield': MRSinMRS.get('B0', None),
            # 'Linewidth': 1,   # TODO: find how to handle or get from REMY
            'TE': MRSinMRS.get('TE', None),
            # 'TE2': 0,   # attention! - only holds for SpinEcho or STEAM
            #             # TODO: find sound solution!
            'Center Freq': MRSinMRS.get('Center Freq', None),   # currently used for plotting only!
        }

        optional = {
            'System': MRSinMRS.get('Manufacturer', None),
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
        # lcmodel backend supports PRESS, STEAM, Spin Echo, LASER

        if protocol is None:
            return None

        if 'mega' in protocol.lower():  # TODO: better check for editing
            print("Warning: LCModelBackend does not support MEGA sequences. "
                  "Ignoring MEGA part of the protocol.")

        if 'slaser' in protocol.lower():
            print("Warning: LCModelBackend does not support sLASER sequences. "
                  "Recommended to switch backend.")
            return None
        elif 'press' in protocol.lower():
            return 'PRESS'
        elif 'steam' in protocol.lower():
            return 'STEAM'
        elif 'spin' in protocol.lower() or 'se' in protocol.lower():
            return 'Spin Echo'
        elif 'laser' in protocol.lower():
            return 'LASER'
        else:
            return None

    def parse2fidA(self, params):
        # change the parameters to the format used by fidA
        if params['Sequence'] == 'Spin Echo': params['Sequence'] = 'se'
        elif params['Sequence'] == 'PRESS': params['Sequence'] = 'p'
        elif params['Sequence'] == 'STEAM': params['Sequence'] = 'st'
        elif params['Sequence'] == 'LASER': params['Sequence'] = 'l'

        params['Add Ref.'] = params['Add Ref.'][0].lower()
        params['Make .raw'] = params['Make .raw'][0].lower()
        return params

    def run_simulation(self, params):
        # create the output directory if it does not exist
        if not os.path.exists(params['Output Path']): os.makedirs(params['Output Path'])

        def sim_lcmrawbasis(n, sw, Bfield, lb, metab, tau1, tau2, addref, makeraw, seq, out_path):
            results = self.octave.feval('sim_lcmrawbasis', n, sw, Bfield, lb, metab,
                                        tau1, tau2, addref, makeraw, seq, out_path + os.sep)
            return metab, results[:, 0] + 1j * results[:, 1]

        # prepare tasks for each metabolite
        params = self.parse2fidA(params)
        tasks = [(float(params['Samples']), float(params['Bandwidth']), float(params['Bfield']),
                  float(params['Linewidth']), metab, float(params['TE']), float(params['TE2']),
                  params['Add Ref.'], params['Make .raw'], params['Sequence'],
                  params['Output Path']) for metab in params['Metabolites']]

        # run simulations sequentially
        results = [sim_lcmrawbasis(*task) for task in tasks]

        # collect results into a dictionary
        basis_set = {metab: data for metab, data in results}
        return basis_set

    def run_simulation_with_progress(self, params, progress_callback):
        # create the output directory if it does not exist
        if not os.path.exists(params['Output Path']): os.makedirs(params['Output Path'])

        def sim_lcmrawbasis(n, sw, Bfield, lb, metab, tau1, tau2, addref, makeraw, seq, out_path):
            results = self.octave.feval('sim_lcmrawbasis', n, sw, Bfield, lb, metab,
                                        tau1, tau2, addref, makeraw, seq, out_path + os.sep)
            return metab, results[:, 0] + 1j * results[:, 1]

        # prepare tasks for each metabolite
        params = self.parse2fidA(params)
        tasks = [(float(params['Samples']), float(params['Bandwidth']), float(params['Bfield']),
                  float(params['Linewidth']), metab, float(params['TE']), float(params['TE2']),
                  params['Add Ref.'], params['Make .raw'], params['Sequence'],
                  params['Output Path']) for metab in params['Metabolites']]

        # initialize the progress bar
        total_steps = len(tasks)
        progress_step = 100 / total_steps

        # run simulations sequentially
        basis_set = {}
        for i, task in enumerate(tasks):
            metab, data = sim_lcmrawbasis(*task)
            basis_set[metab] = data
            progress_callback(i + 1, total_steps)   # update the progress bar
        return basis_set