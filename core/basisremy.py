####################################################################################################
#                                            basisremy.py                                          #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 08/10/25                                                                                #
#                                                                                                  #
# Purpose: Defines the BasisREMY class for extracting REMY parameters from MRS data                #
#          and simulating a basis set using different backends.                                    #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import pathlib
import numpy as np

# own
from backends.lcmodel_backend import LCModelBackend
from backends.slaser_backend import sLaserBackend
from externals.remy.MRSinMRS import DataReaders, Table, setup_log, write_log


#**************************************************************************************************#
#                                            BasisREMY                                             #
#**************************************************************************************************#
#                                                                                                  #
# The BasisREMY class is the main class for the BasisREMY tool. It provides the functionality to   #
# extract REMY parameters from MRS data and simulate a basis.                                      #
#                                                                                                  #
#**************************************************************************************************#
class BasisREMY:
    def __init__(self, backend='sLaserSim'):
        self.DRead = DataReaders()
        self.Table = Table()

        self.backends = {
            'LCModel': LCModelBackend(),
            'sLaserSim': sLaserBackend(),
        }
        self.backend = self.backends[backend]
        self.available_backends = list(self.backends.keys())

    def set_backend(self, backend):
        # set the backend to the selected one
        if backend in self.backends:
            old_backend = self.backend
            self.backend = self.backends[backend]
            self.backend.update_from_backend(old_backend)  # update parameters
            print(f"Backend set to: {self.backend.name}")
        else:
            raise ValueError(f"Unknown backend: {backend}. Available backends: {self.available_backends}")

    def run(self, import_fpath, export_fpath=None, method=None, userParams={}, optionalParams={}):
        # run REMY on the selected file
        MRSinMRS = self.runREMY(import_fpath, method)
        params, opt = self.parseREMY(MRSinMRS)
        params['Output Path'] = export_fpath if export_fpath is not None else './'

        # update the mandatory parameters
        self.backend.mandatory_params.update(params)
        self.backend.mandatory_params.update(userParams)

        # update the optional parameters
        self.backend.optional_params.update(opt)
        self.backend.optional_params.update(optionalParams)

        # run fidA simulation
        basis = self.backend.run_simulation(self.mandatory_params)

        # plot the basis set
        import matplotlib.pyplot as plt
        plt.figure()
        for key, value in basis.items():
            plt.plot(np.fft.fft(value), label=key)
        plt.legend()
        plt.show()

    def runREMY(self, import_fpath, method=None):
        # run REMY datareader on the selected file
        if method is None: suf = pathlib.Path(import_fpath).suffix.lower()
        else: suf = method
        log = None
        if suf == '.dat':   # Siemens Twix file
            write_log(log, 'Data Read: Siemens Twix uses pyMapVBVD ')  # log - pyMapVBVD
            MRSinMRS, log = self.DRead.siemens_twix(import_fpath, log)
            vendor_selection = 'Siemens'
        elif suf == '.ima':  # Siemens Dicom file
            write_log(log, 'Data Read: Siemens Dicom uses pydicom ')  # log - pyDicom
            MRSinMRS, log = self.DRead.siemens_ima(import_fpath, log)
            vendor_selection = 'Siemens'
        elif suf == '.rda':  # Siemens RDA file
            write_log(log, 'Data Read: Siemens RDA directly read with RMY ')  # log - pyDicom
            MRSinMRS, log = self.DRead.siemens_rda(import_fpath,    log)
            vendor_selection = 'Siemens'
        elif suf == '.spar':  # Philips SPAR file
            write_log(log, 'Data Read: Philips SPAR uses spec2nii ')  # log - spec2nii
            MRSinMRS, log = self.DRead.philips_spar(import_fpath, log)
            vendor_selection = 'Philips'
        elif suf == '.7':  # GE Pfile
            write_log(log, 'Data Read: GE Pfile uses spec2nii ')  # log - spec2nii
            MRSinMRS, log = self.DRead.ge_7(import_fpath, log)
            vendor_selection = 'GE'
        elif suf == 'bruker_method':  # Bruker Method file
            write_log(log, 'Data Read: Bruker Method uses spec2nii ')  # log - spec2nii
            MRSinMRS, log = self.DRead.bruker_method(import_fpath, log)
            vendor_selection = 'Bruker'
        elif suf == 'bruker_2dseq':  # Bruker 2dseq file
            write_log(log, 'Data Read: Bruker uses BrukerAPI ' +  # log - BrukerAPI
                      'developed by Tomáš Pšorn\n\t' +
                      'github.com/isi-nmr/brukerapi-python')
            MRSinMRS, log = self.DRead.bruker_2dseq(import_fpath, log)
            vendor_selection = 'Bruker'
        elif suf == '.nii' or suf == '.nii.gz':
            write_log(log, 'Data Read: NIfTI json side car')  # log - NIfTI JSON side car
            MRSinMRS, log = self.DRead.nifti_json(import_fpath, log)
            vendor_selection = 'NIfTI'
        else:
            raise ValueError(f'Unknown file format {suf}! Valid formats are:'
                             f' .dat, .ima, .rda, .spar, .7, bruker_method, bruker_2dseq, .nii, .nii.gz')

        # clean the data
        dtype_selection = suf.replace('bruker_', '').replace('.', '')
        MRSinMRS = self.Table.table_clean(vendor_selection, dtype_selection, MRSinMRS)
        return MRSinMRS

    def parseREMY(self, MRSinMRS):
        # extract as much information as possible from the MRSinMRS dict
        mandatory = {
            # 'Sequence': None,  # TODO: find a way to get this from REMY
            'Samples': MRSinMRS.get('samples', None),
            'Bandwidth': MRSinMRS.get('sample_frequency', None),
            'Bfield': MRSinMRS.get('FieldStrength', None),
            # 'Linewidth': 1,   # TODO: find how to get from REMY
            'TE': MRSinMRS.get('EchoTime', None),
            # 'TE2': 0,   # attention! - only holds for SpinEcho or STEAM
            #             # TODO: find sound solution!
            'Center Freq': MRSinMRS.get('synthesizer_frequency', None),
        }
        optional = {
            'Nucleus': MRSinMRS.get('Nucleus', None),
            'TR': MRSinMRS.get('RepetitionTime', None),
        }
        return mandatory, optional

    def run_gui(self):
        from gui.application import Application
        app = Application()
        app.mainloop()
