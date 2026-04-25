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
import json
import numpy as np
import pathlib

# own
from backends.fslmrs_backend import FSLMRSBackend
from backends.mrscloud_backend import MRSCloudBackend
from backends.custom_backends import CustomSLaser
from backends.fida_backends import FIDA_BACKENDS
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
    # Display order for the top-level Category dropdown.
    CATEGORY_ORDER = ['FID-A', 'Custom', 'MRSCloud', 'FSL-MRS']

    def __init__(self, backend='FidaIdeal'):
        self.DRead = DataReaders()
        self.Table = Table()

        # Build the flat backend registry. The FID-A category contains many
        # entries (FidaIdeal = ex-LCModel, plus PRESS shaped, MEGA-PRESS
        # shaped, …); the other categories currently have a single backend.
        self.backends = {}
        for cls in FIDA_BACKENDS:
            inst = cls()
            self.backends[inst.name] = inst
        # Custom category
        custom = CustomSLaser()
        self.backends[custom.name] = custom
        # MRSCloud + FSL-MRS top-level categories
        mc = MRSCloudBackend(); self.backends[mc.name] = mc
        fm = FSLMRSBackend();   self.backends[fm.name] = fm

        # Group backends by their declared category. Order within a category
        # follows the registration order above.
        self.categories: dict[str, list[str]] = {c: [] for c in self.CATEGORY_ORDER}
        for name, inst in self.backends.items():
            cat = getattr(inst, 'category', 'Other')
            self.categories.setdefault(cat, []).append(name)

        self.set_backend(backend)

    @property
    def available_backends(self):
        return list(self.backends.keys())

    def set_backend(self, backend):
        """Switch to the named backend."""
        if backend not in self.backends:
            raise ValueError(
                f"Unknown backend: {backend}. Available backends: "
                f"{self.available_backends}"
            )
        old_backend = getattr(self, 'backend', None)
        self.backend = self.backends[backend]
        if old_backend is not None and old_backend is not self.backend:
            self.backend.update_from_backend(old_backend)
        print(f"Backend set to: {self.backend.name}")

    def set_category(self, category):
        """Select the first backend in the given category. Convenience for
        the GUI's two-level Category → Backend dropdown."""
        if category not in self.categories or not self.categories[category]:
            raise ValueError(
                f"Unknown / empty category: {category}. "
                f"Available: {list(self.categories.keys())}"
            )
        self.set_backend(self.categories[category][0])

    def get_current_category(self):
        return getattr(self.backend, 'category', 'Other')

    def run(self, import_fpath, export_fpath=None, method=None, userParams={}, optionalParams={}, plot=False):
        # run REMY on the selected file
        MRSinMRS = self.runREMY(import_fpath, method)
        params, opt = self.backend.parseREMY(MRSinMRS)
        params['Output Path'] = export_fpath if export_fpath is not None else './'

        # update the mandatory parameters
        self.backend.mandatory_params.update(params)
        self.backend.mandatory_params.update(userParams)

        # update the optional parameters
        self.backend.optional_params.update(opt)
        self.backend.optional_params.update(optionalParams)

        # run fidA simulation
        basis = self.backend.run_simulation(self.backend.mandatory_params)

        # plot the basis set
        if plot:
            import matplotlib.pyplot as plt
            plt.figure()
            for key, value in basis.items():
                plt.plot(np.fft.fft(value), label=key)
            plt.legend()
            plt.show()

        return basis, params

    def runREMY(self, import_fpath, method=None):
        # run REMY datareader on the selected file
        if method is None: suf = pathlib.Path(import_fpath).suffix.lower()
        else: suf = method

        # check for bruker mehtod or 2dseq (no suffix)
        if suf == '':
            if 'method' in pathlib.Path(import_fpath).name.lower():
                suf = 'method'
            elif '2dseq' in pathlib.Path(import_fpath).name.lower():
                suf = '2dseq'

        if suf == '.gz':  # check for .nii.gz
            if pathlib.Path(import_fpath).name.lower().endswith('.nii.gz'):
                suf = '.nii.gz'

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
        elif suf == 'method':  # Bruker Method file
            write_log(log, 'Data Read: Bruker Method uses spec2nii ')  # log - spec2nii
            MRSinMRS, log = self.DRead.bruker_method(import_fpath, log)
            vendor_selection = 'Bruker'
        elif suf == '2dseq':  # Bruker 2dseq file
            write_log(log, 'Data Read: Bruker uses BrukerAPI ' +  # log - BrukerAPI
                      'developed by Tomáš Pšorn\n\t' +
                      'github.com/isi-nmr/brukerapi-python')
            MRSinMRS, log = self.DRead.bruker_2dseq(import_fpath, log)
            vendor_selection = 'Bruker'
        elif suf == '.nii' or suf == '.nii.gz':   # TODO: might want to adjust to be .nii.gz
            write_log(log, 'Data Read: NIfTI json side car')  # log - NIfTI JSON side car
            # MRSinMRS, log = self.DRead.nifti_json(import_fpath, log)   # TODO: fix for nifti
            try:
                with open(import_fpath.replace(suf, '.json'), 'r') as f:
                    MRSinMRS = json.load(f)
            except:
                from nifti_mrs.nifti_mrs import NIFTI_MRS
                MRSinMRS = NIFTI_MRS(import_fpath).hdr_ext

            # homogenize keys to be strings
            MRSinMRS = {str(k): str(v[0]) if isinstance(v, list) and len(v) == 1 else v for k, v in
                        dict(MRSinMRS).items()}

            vendor_selection = 'NIfTI'
        else:
            raise ValueError(f'Unknown file format {suf}! Valid formats are:'
                             f' .dat, .ima, .rda, .spar, .7, bruker_method, bruker_2dseq, .nii, .nii.gz')

        dtype_selection = suf.replace('.', '')  # remove dot if present
        if suf == '.nii.gz': dtype_selection = 'json'  # special case

        # check for missing MRSinMRS Values that might have different names across versions
        try:
            MRSinMRS = self.Table.table_clean(vendor_selection, dtype_selection, MRSinMRS)
        except Exception as e:
            print(f"Warning: table_clean failed: {e}")

        # populate MRS Table
        try:
            self.Table.populate(vendor_selection, dtype_selection, MRSinMRS)
        except Exception as e:
            print(f"Warning: populate table failed: {e}")

        # get unform MRSinMRS table
        MRSinMRS_unif = self.flatten_mrsinmrs_table(self.Table.MRSinMRS_Table)

        # extend with more info
        MRSinMRS_unif.update(self.extract_more(MRSinMRS, vendor_selection, dtype_selection))
        return MRSinMRS_unif

    def extract_more(self, MRSinMRS, vendor, dtype):
        # extract additional information from the raw MRSinMRS dict if possible
        add_info = {}

        if vendor == 'Philips':
            if dtype == 'spar': # Philips SPAR specific
                if 'synthesizer_frequency' in MRSinMRS:
                    add_info['Center Freq'] = MRSinMRS['synthesizer_frequency']

        elif vendor == 'Siemens':
            # can be 'lFrequency', 'Frequency', 'SpectrometerFrequency', 'MRFrequency'
            if 'lFrequency' in MRSinMRS:
                add_info['Center Freq'] = MRSinMRS['lFrequency']
            elif 'Frequency' in MRSinMRS:
                add_info['Center Freq'] = MRSinMRS['Frequency']
            elif 'SpectrometerFrequency' in MRSinMRS:
                add_info['Center Freq'] = MRSinMRS['SpectrometerFrequency']
            elif 'MRFrequency' in MRSinMRS:
                add_info['Center Freq'] = MRSinMRS['MRFrequency']

        elif vendor == 'GE':
            if dtype == '7': # GE Pfile specific
                if 'synthesizer_frequency' in MRSinMRS:
                    add_info['Center Freq'] = MRSinMRS['rhr_rh_ps_mps_freq']

        elif vendor == 'Bruker':
            pass

        elif vendor == 'NIfTI':
            add_info['Center Freq'] = MRSinMRS['SpectrometerFrequency']
            add_info['ExcitationFlipAngle'] = MRSinMRS['ExcitationFlipAngle']

        return add_info

    def flatten_mrsinmrs_table(self, df):
        flat_dict = {}

        # iterate over all rows
        for idx, row in df.iterrows():
            key = str(row['Generic']).strip()  # lower-level key
            val = row['Values']

            # skip empty keys
            if key != '' and key != 'nan':
                # if value is bytes, decode
                if isinstance(val, bytes):
                    val = val.decode(errors='ignore')
                flat_dict[key] = val

        return flat_dict
