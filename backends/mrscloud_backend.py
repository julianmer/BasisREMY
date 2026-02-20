####################################################################################################
#                                        mrscloud_backend.py                                       #
####################################################################################################
#                                                                                                  #
# Authors: G. Simegn (gsimegn1@jh.edu)                                                             #
#          J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 18/02/26                                                                                #
#                                                                                                  #
# Purpose: Defines the MRSCloud backend class for simulating MRS basis sets.                       #
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
#                                         MRSCloudBackend                                          #
#**************************************************************************************************#
#                                                                                                  #
# Implements the basis set simulation backend using the MRSCloud approach.                         #
#                                                                                                  #
#**************************************************************************************************#
class MRSCloudBackend(Backend):
    def __init__(self):
        super().__init__()
        self.name = 'MRSCloud'

        # Mark that this backend requires Octave
        self.requires_octave = True
        # Octave will be initialized lazily when needed

        # define possible metabolites
        self.metabs = {
		    'Ala': True,
		    'Asc': True,
		    'Asp': True,
		    'Ch': False,        # not in provided list
		    'Cit': True,
		    'Cr': True,
		    'EtOH': True,
		    'GABA': True,
		    'GPC': True,
		    'GSH': True,
		    'Glc': False,       # not in provided list
		    'Gln': True,
		    'Glu': True,
		    'Gly': True,
		    'H2O': False,       # not in provided list
		    'Ins': True,
		    'Lac': True,
		    'Lip': False,       # not in provided list
		    'NAA': True,
		    'NAAG': True,
		    'PCh': True,
		    'PCr': True,
		    'PE': False,        # not in provided list
		    'Phenyl': True,
		    'Ref0ppm': False,   # not in provided list
		    'Scyllo': True,
		    'Ser': True,
		    'Tau': True,
		    'Tyros': True,
		    'bHB': True,
		    'bHG': True,
		}


        # dropdown options
        self.dropdown = {
			'System': ['Philips', 'Philips_universal', 'Siemens', 'GE'],
            'Sequence': ['unEdited','MEGA', 'HERMES', 'HERCULES', 'sLASER', 'STEAM' ],
            'Add Ref.': ['Yes', 'No'],
            'Make .raw': ['Yes', 'No'],
        }

        # add file selection fields
        self.file_selection = ['Path to Pulse']
		
        # Default spin system path (MRSCloud style)
		#self.default_spin_system_path = './externals/mrscloud/spinSystems/' TODO-Gize

        # define dictionary of mandatory parameters
        self.mandatory_params = {
            "System": None, # MUST be selected
			'Sequence': None, # MUST be selected
			'Samples': None,
		    'Bandwidth': None,
		    'Bfield': None,
		    'Linewidth': 1,
		    'TE': None,
		    'TE2': 0,
		    'Add Ref.': 'No',
		    'Make .raw': 'Yes',
		    'Output Path': None,
		    'Metabolites': [key for key, value in self.metabs.items() if value],
		    'Center Freq': None,

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
        #  backend supports UnEdited/MEGA/HERMES/HERCULES, slASER

        if protocol is None:
            return None
        protocol = str(protocol).lower()    
        if 'mega' in protocol:
            return 'MEGA'
        elif 'hermes' in protocol:
            return 'HERMES'
        elif 'hercules' in protocol:
            return 'HERCULES'
        elif 'slaser' in protocol:
            return 'sLASER'
        elif 'unedited' in protocol or 'press' in protocol:
            return 'unEdited'
        elif 'steam' in str(protocol).lower():
            return 'STEAM'

        else:
            print("Warning: Unsupported sequence for MRSCloud.")
            return None
        
    def parseSystem(self, system):
        # backend only supports Philips and Siemens systems for now
        if system is None:
            return None

        if 'philips' in system.lower():
            return 'Philips'
        elif 'siemens' in system.lower():
            return 'Siemens'
        elif 'universal_Philips' in system.lower():
            return 'Universal_Philips'
        elif 'ge' in system.lower():
            return 'GE'
        else:
            print("Warning: sLaserBackend only supports Philips, Siemens and GE systems. ")
            return None
    def parse2fidA(self, params):
        # change the parameters to the format used by fidA / MRSCloud
        if params['Sequence'] == 'unEdited': params['Sequence'] = 'p'
        elif params['Sequence'] == 'MEGA': params['Sequence'] = 'mega'
        elif params['Sequence'] == 'HERMES': params['Sequence'] = 'hermes'
        elif params['Sequence'] == 'HERCULES': params['Sequence'] = 'hercules'
        elif params['Sequence'] == 'sLASER': params['Sequence'] = 'sl'
        elif params['Sequence'] == 'STEAM': params['Sequence'] = 'st'
        else: raise ValueError(f"Unsupported sequence: {params['Sequence']}")

        params['Add Ref.'] = params['Add Ref.'][0].lower()
        params['Make .raw'] = params['Make .raw'][0].lower()
        return params

    def setup_octave_paths(self):
        """Setup Octave paths for FID-A toolbox."""
        if self.octave is None:
            raise RuntimeError("Octave not initialized. Call initialize_octave() first.")
         
        # MRSCloud can also use this fidA? or just add all paths in MRSCloud
       # self.octave.eval("warning('off', 'all');")
       # self.octave.addpath('./externals/fidA/inputOutput/')
       # self.octave.addpath('./externals/fidA/processingTools/')
       # self.octave.addpath('./externals/fidA/simulationTools/')
        # the .genpath should add all the subfolders?
        self.octave.addpath(self.octave.genpath('./externals/mrscloud/functions/'))

    

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

        def sim_lcmrawbasis(n, sw, Bfield, lb, metab, tau1, tau2, addref, makeraw, seq, out_path):
            results = self.octave.feval('sim_lcmrawbasis', n, sw, Bfield, lb, metab,
                                        tau1, tau2, addref, makeraw, seq, out_path + os.sep)
            return metab, results[:, 0] + 1j * results[:, 1]

        # prepare tasks for each metabolite
        params = self.parse2fidA(params)
        tasks = [(float(params['Samples']), float(params['Bandwidth']), float(params['Bfield']),
                  float(params['Linewidth']), metab, float(params['TE']), float(params['TE2']),
                  params['Add Ref.'], params['Make .raw'], params['Sequence'],
                  output_path) for metab in params['Metabolites']]

        basis_set = {}
        total_steps = len(tasks)
        for i, task in enumerate(tasks):
            metab, data = sim_lcmrawbasis(*task)

            # Ensure it's a proper numpy array
            if not isinstance(data, np.ndarray):
                data = np.array(data, dtype=complex)

            # Flatten if multidimensional
            if data.ndim > 1:
                data = data.flatten()

            basis_set[metab] = data
            if progress_callback:
                progress_callback(i + 1, total_steps)
        return basis_set
