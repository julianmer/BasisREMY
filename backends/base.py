####################################################################################################
#                                             base.py                                              #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 08/10/25                                                                                #
#                                                                                                  #
# Purpose: Defines the abstract Backend class as a base for all simulation backends.               #
#                                                                                                  #
####################################################################################################


#**************************************************************************************************#
#                                             Backend                                              #
#**************************************************************************************************#
#                                                                                                  #
# Defines the backend structure for the simulations. Inherit from this class to create a new       #
# simulation backend with all mandatory attributes and methods.                                    #
#                                                                                                  #
#**************************************************************************************************#
class Backend:
    def __init__(self):
        self.name = None

        # define possible metabolites
        self.metabs = {}

        # dropdown options
        self.dropdown = {}

        # add file selection fields
        self.file_selection = []

        # define dictionary of mandatory parameters
        self.mandatory_params = {}

        # define dictionary of optional parameters
        self.optional_params = {}

    def update_from_backend(self, backend):
        # update the backend parameters from another backend instance
        # TODO: make sure some parameters are reset (see e.g. Sequence in dropdown)
        self.metabs.update({k: v for k, v in backend.metabs.items() if k in self.metabs})
        self.mandatory_params.update({k: v for k, v in backend.mandatory_params.items()
                                      if k in self.mandatory_params})
        self.optional_params.update({k: v for k, v in backend.optional_params.items()
                                     if k in self.optional_params})
        # note: for now no dropdown options are updated, as they are too specific to the backend
        #       at this point

    def parseREMY(self, MRSinMRS):
        # parse REMY output to parameters and optional parameters for the backend
        raise NotImplementedError("This method should be overridden by subclasses.")

    def parseProtocol(self, protocol):
        # parse the sequence name from the raw MRS data protocol to
        # a backend-specific sequence name if needed
        raise NotImplementedError("This method should be overridden by subclasses.")

    def run_simulation(self, params):
        raise NotImplementedError("This method should be overridden by subclasses.")

    def run_simulation_with_progress(self, params, progress_callback):
        raise NotImplementedError("This method should be overridden by subclasses.")
