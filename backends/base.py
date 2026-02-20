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

        # Octave runtime management
        self.requires_octave = False  # Set to True in subclasses that need Octave
        self.octave = None  # Will be initialized when needed

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

        # Mode support - every backend has at least 'Default'
        # Subclasses override to add more modes
        self.modes = ['Default']
        self.current_mode = 'Default'

    def get_modes(self):
        """Return list of available modes for this backend"""
        return self.modes

    def get_current_mode(self):
        """Return currently selected mode"""
        return self.current_mode

    def set_mode(self, mode):
        """
        Set the current mode and return parameters for that mode.

        Subclasses should override get_params_for_mode() to define
        mode-specific parameter sets.

        Args:
            mode: One of the values from self.modes

        Returns:
            dict: Parameters to display for this mode
        """
        if mode not in self.modes:
            raise ValueError(f"Unknown mode '{mode}'. Available: {self.modes}")
        self.current_mode = mode
        return self.get_params_for_mode(mode)

    def get_params_for_mode(self, mode=None):
        """
        Return parameters to display in GUI for the given mode.

        Default implementation returns all mandatory_params.
        Subclasses override this to return mode-specific parameter sets.

        Args:
            mode: Mode name (uses current_mode if None)

        Returns:
            dict: Parameters to show in GUI
        """
        return dict(self.mandatory_params)

    def initialize_octave(self, prefer_docker=True, verbose=False):
        """
        Initialize Octave runtime if required by this backend.

        Args:
            prefer_docker: If True, try Docker first, otherwise try local first
            verbose: Enable verbose output for debugging

        Returns:
            True if successful, False otherwise

        Raises:
            RuntimeError: If Octave is required but cannot be initialized
        """
        if not self.requires_octave:
            return True

        from core.octave_manager import OctaveManager
        manager = OctaveManager(verbose=verbose)

        try:
            self.octave = manager.initialize_octave(prefer_docker=prefer_docker)

            # If using Docker, check for and clean up old processes
            if hasattr(self.octave, 'check_running_processes'):
                existing = self.octave.check_running_processes()
                if existing:
                    print(f"⚠️  Found {len(existing)} existing Octave process(es)")
                    print(f"   Cleaning up old processes...")
                    self.octave.kill_running_processes()
                    print(f"✓ Ready for new simulation")

            return True
        except RuntimeError as e:
            raise RuntimeError(f"Failed to initialize Octave for {self.name} backend:\n{e}")

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

    def run_simulation(self, params, progress_callback=None):
        raise NotImplementedError("This method should be overridden by subclasses.")
