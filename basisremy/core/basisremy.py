####################################################################################################
#                                            basisremy.py                                          #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 08/10/25                                                                                #
#                                                                                                  #
# Purpose: Defines the BasisREMY class for extracting REMY parameters from MRS data                 #
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
from basisremy.backends.fslmrs_backend import FSLMRSBackend
from basisremy.backends.mrscloud_backend import MRSCloudBackend
from basisremy.backends.custom_backends import CustomSLaser
from basisremy.backends.fida_backends import FIDA_BACKENDS
from basisremy.backends.vespa_backend import VespaBackend
from basisremy.remy.MRSinMRS import DataReaders, Table, write_log


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
    CATEGORY_ORDER = ['MRSCloud', 'FID-A', 'FSL-MRS', 'Vespa', 'Custom']

    def __init__(self, backend='MRSCloud'):
        self.DRead = DataReaders()
        self.Table = Table()

        # Cache the last REMY-extracted MRSinMRS dict so that switching
        # backends can re-parse it with the new backend's parseREMY().
        self._last_mrsinmrs = None

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
        # MRSCloud + FSL-MRS + Vespa top-level categories
        mc = MRSCloudBackend(); self.backends[mc.name] = mc
        fm = FSLMRSBackend();   self.backends[fm.name] = fm
        vp = VespaBackend();    self.backends[vp.name] = vp

        # Snapshot each backend's pristine defaults so a new file import can
        # start clean instead of inheriting the previous file's values.
        import copy
        for inst in self.backends.values():
            inst._defaults = (copy.deepcopy(inst.mandatory_params),
                              copy.deepcopy(inst.optional_params))
            inst._default_metabs = dict(inst.metabs)

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

        # Step 1 — if we have cached REMY data, let the NEW backend parse it
        # first so it gets the parameters it specifically needs (e.g. FID-A
        # needs Bfield + Center Freq which MRSCloud doesn't expose, and
        # MRSCloud needs System + Field Strength which FID-A doesn't expose).
        if self._last_mrsinmrs is not None:
            try:
                params, opt = self.backend.parseREMY(self._last_mrsinmrs)
                self.backend.mandatory_params.update(
                    {k: v for k, v in params.items() if v is not None})
                self.backend.optional_params.update(
                    {k: v for k, v in opt.items() if v is not None})
            except Exception as e:
                print(f"Warning: could not re-parse REMY data for {backend}: {e}")

        # Step 2 — overlay user-modified physics values from the old backend
        # (e.g. TE, Samples, Bandwidth the user changed by hand) so manual
        # edits survive the switch.
        #
        # Rules:
        #   * Only copy keys that exist in the NEW backend (don't inject alien keys).
        #   * Skip Sequence and Metabolites — they are backend-specific.
        #   * ONLY copy non-None values — a None on the old backend means
        #     "REMY didn't find this / was never set", and we must not let it
        #     overwrite the freshly re-parsed step-1 value (e.g. Bfield, Center Freq).
        if old_backend is not None and old_backend is not self.backend:
            _skip = {'Sequence', 'Metabolites'}
            self.backend.mandatory_params.update({
                k: v for k, v in old_backend.mandatory_params.items()
                if k in self.backend.mandatory_params
                and k not in _skip
                and v is not None
            })
            self.backend.optional_params.update({
                k: v for k, v in old_backend.optional_params.items()
                if k in self.backend.optional_params
                and v is not None
            })
            # Carry over shared metabolite on/off state (e.g. user turned off GABA).
            self.backend.metabs.update({
                k: v for k, v in old_backend.metabs.items()
                if k in self.backend.metabs
            })
            if hasattr(self.backend, '_refresh_metab_list'):
                self.backend._refresh_metab_list()
            elif 'Metabolites' in self.backend.mandatory_params:
                self.backend.mandatory_params['Metabolites'] = [
                    k for k, v in self.backend.metabs.items() if v]

            # Sequence: translate old backend's sequence into the new backend's
            # vocabulary using map_sequence_in(). This preserves the user's
            # explicit sequence choice (e.g. FID-A 'STEAM' → MRSCloud
            # 'STEAM (7T only)') and takes priority over the REMY-parsed value
            # because it reflects what the user actually *wants* to simulate.
            old_seq = (old_backend.mandatory_params.get('Sequence') or '').strip()
            if old_seq:
                mapped = self.backend.map_sequence_in(old_seq)
                if mapped is not None:
                    self.backend.mandatory_params['Sequence'] = mapped
                # If no mapping (e.g. sLASER → FidaIdeal), leave whatever step 1
                # set so the user is shown the REMY value or None.

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

    def reset_backend_params(self):
        """Restore the active backend's parameter defaults (per new file).

        The metabolite selection is preserved — curating it is user work
        that a new import must not wipe.
        """
        import copy
        defaults = getattr(self.backend, '_defaults', None)
        if defaults is None:
            return
        mandatory, optional = defaults
        keep_metabs = self.backend.mandatory_params.get('Metabolites')
        self.backend.mandatory_params = copy.deepcopy(mandatory)
        self.backend.optional_params = copy.deepcopy(optional)
        if keep_metabs is not None and 'Metabolites' in self.backend.mandatory_params:
            self.backend.mandatory_params['Metabolites'] = list(keep_metabs)

    def run(self, import_fpath, export_fpath=None, method=None, userParams={}, optionalParams={}, plot=False):
        # run REMY on the selected file (starting from clean defaults)
        MRSinMRS = self.runREMY(import_fpath, method)
        self.reset_backend_params()
        params, opt = self.backend.parseREMY(MRSinMRS)
        params['Output Path'] = export_fpath if export_fpath is not None else './'

        # update the mandatory parameters (drop None so REMY gaps don't
        # clobber sensible backend defaults; explicit userParams still win)
        self.backend.mandatory_params.update({k: v for k, v in params.items()
                                              if v is not None})
        self.backend.mandatory_params.update(userParams)

        # update the optional parameters
        self.backend.optional_params.update({k: v for k, v in opt.items()
                                             if v is not None})
        self.backend.optional_params.update(optionalParams)

        # run fidA simulation (mandatory params win on key clashes)
        basis = self.backend.run_simulation({**self.backend.optional_params,
                                             **self.backend.mandatory_params})

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

        # Fresh table per import: clean/populate failures must not leave the
        # previous file's values behind for the next one.
        self.Table = Table()

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
        elif suf == '.nii' or suf == '.nii.gz':
            write_log(log, 'Data Read: NIfTI json side car')  # log - NIfTI JSON side car
            try:
                # replace only the suffix — a '.nii' elsewhere in the path
                # must not be rewritten
                sidecar = import_fpath[:-len(suf)] + '.json'
                with open(sidecar, 'r') as f:
                    MRSinMRS = json.load(f)
            except Exception:
                from nifti_mrs.nifti_mrs import NIFTI_MRS
                hdr = NIFTI_MRS(import_fpath).hdr_ext
                # nifti_mrs >= 1.x returns a Hdr_Ext object, not a dict
                MRSinMRS = hdr.to_dict() if hasattr(hdr, 'to_dict') else hdr

            # homogenize keys to be strings
            MRSinMRS = {str(k): str(v[0]) if isinstance(v, list) and len(v) == 1 else v for k, v in
                        dict(MRSinMRS).items()}

            vendor_selection = 'NIfTI'
        else:
            raise ValueError(f'Unknown file format {suf}! Valid formats are:'
                             f' .dat, .ima, .rda, .spar, .7, bruker_method, bruker_2dseq, .nii, .nii.gz')

        dtype_selection = suf.replace('.', '')  # remove dot if present
        # NIfTI (either suffix) uses the sidecar-JSON labels; a plain '.nii'
        # would otherwise look up a nonexistent 's2nlabel_nifti_nii' column
        # and silently lose every extracted parameter.
        if suf in ('.nii', '.nii.gz'): dtype_selection = 'json'

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

        # Cache for later backend switches
        self._last_mrsinmrs = MRSinMRS_unif

        return MRSinMRS_unif

    @staticmethod
    def _freq_mhz(value):
        # Normalize a spectrometer frequency to MHz ('Center Freq' is MHz
        # everywhere in BasisREMY). Headers store Hz (e.g. 127736713) or
        # MHz (e.g. 127.7); non-numeric values pass through unchanged.
        if isinstance(value, list) and value:
            value = value[0]
        try:
            v = float(value)
        except (TypeError, ValueError):
            return value
        return v / 1e6 if v > 1e6 else v

    def extract_more(self, MRSinMRS, vendor, dtype):
        # extract additional information from the raw MRSinMRS dict if possible
        add_info = {}

        if vendor == 'Philips':
            if dtype == 'spar': # Philips SPAR specific
                if 'synthesizer_frequency' in MRSinMRS:
                    add_info['Center Freq'] = self._freq_mhz(MRSinMRS['synthesizer_frequency'])

        elif vendor == 'Siemens':
            # can be 'lFrequency', 'Frequency', 'SpectrometerFrequency', 'MRFrequency'
            if 'lFrequency' in MRSinMRS:
                add_info['Center Freq'] = self._freq_mhz(MRSinMRS['lFrequency'])
            elif 'Frequency' in MRSinMRS:
                add_info['Center Freq'] = self._freq_mhz(MRSinMRS['Frequency'])
            elif 'SpectrometerFrequency' in MRSinMRS:
                add_info['Center Freq'] = self._freq_mhz(MRSinMRS['SpectrometerFrequency'])
            elif 'MRFrequency' in MRSinMRS:
                add_info['Center Freq'] = self._freq_mhz(MRSinMRS['MRFrequency'])

        elif vendor == 'GE':
            if dtype == '7': # GE Pfile specific
                if 'rhr_rh_ps_mps_freq' in MRSinMRS:
                    # Pfiles store the frequency in 0.1 Hz units (≈1.28e9 at 3T)
                    freq = MRSinMRS['rhr_rh_ps_mps_freq']
                    try:
                        f = float(freq)
                        add_info['Center Freq'] = f / 1e7 if f > 1e8 else self._freq_mhz(f)
                    except (TypeError, ValueError):
                        pass

        elif vendor == 'Bruker':
            pass

        elif vendor == 'NIfTI':
            # SpectrometerFrequency is already MHz per the NIfTI-MRS spec
            # (and may be a list); ExcitationFlipAngle is optional.
            if 'SpectrometerFrequency' in MRSinMRS:
                add_info['Center Freq'] = self._freq_mhz(MRSinMRS['SpectrometerFrequency'])
            if 'ExcitationFlipAngle' in MRSinMRS:
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
