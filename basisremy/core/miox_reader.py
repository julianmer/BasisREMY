####################################################################################################
#                                          miox_reader.py                                          #
####################################################################################################
#                                                                                                  #
# Purpose: MIOX opens the data file; the header comes back in the shape REMY's table mapping       #
#          expects (the keys the REMY readers produced), so extraction stays unchanged while every #
#          format MIOX reads is available. A format REMY had no reader for (Philips DICOM) comes   #
#          back with NIfTI-MRS keys and goes through the NIfTI column of the table. Bruker method  #
#          files stay with REMY's own text parser (MIOX wants the scan folder), NIfTI-MRS files    #
#          are read directly.                                                                      #
#                                                                                                  #
####################################################################################################

#*************#
#   imports   #
#*************#
import re


# MR Spectroscopy Storage: XA DICOM, whose dwell time is not oversampled
_XA_SPECTROSCOPY = '1.2.840.10008.5.1.4.1.1.4.2'


def read(path, suffix):
    """(REMY-style header dict, vendor, file suffix) of ``path`` as MIOX reads it, or None when
    MIOX is not installed, cannot read the file, or the format keeps its own reader."""
    if suffix in ('method', '.nii', '.nii.gz'):
        return None
    try:
        import miox
    except ImportError:
        return None
    try:
        datasets = miox.load(path)
    except Exception:                                      # noqa: BLE001 - REMY's reader is next
        return None
    if not datasets:
        return None
    fmt = str(datasets[0].provenance.format.value)
    if fmt == 'dicom_mrs':
        maker = str(datasets[0].meta.get('Manufacturer', '')).lower()
        return _siemens_dicom(datasets) if 'siemens' in maker else _nifti_style(datasets)
    if fmt == 'bruker' and suffix != '2dseq':
        return None
    builder = _BUILDERS.get(fmt)
    return builder(datasets) if builder else None


def _first(value):
    return value[0] if isinstance(value, (list, tuple)) and value else value


def _zooms(affine):
    """Voxel sizes (mm) along the three image axes."""
    import numpy as np
    return [round(float(np.linalg.norm(affine[:3, i])), 6) for i in range(3)]


def _twix(datasets):
    v = datasets[-1].meta.vendor                           # a multi-measurement file: the scan
    raw = {'Manufacturer': 'Siemens'}
    raw.update(v.get('Dicom', {}))
    raw.update(v.get('Config', {}))
    phoenix = dict(v.get('Phoenix', {}))                   # ASCCONV notation already
    raw['RFPulseNames'] = [str(x).strip('"') for k, x in phoenix.items()
                           if re.fullmatch(r'sTXSPEC\.aRFPULSE\[\d+\]\.tName', k)]
    raw['_ascconv'] = phoenix
    if 'TE_Time' not in raw and 'TE' not in raw and 'alTE[0]' in phoenix:
        raw['TE'] = phoenix['alTE[0]']
    for key in list(v.get('Dicom', {})) + list(v.get('Config', {})):
        if 'dwelltime' in key.lower():
            raw[key] = str(raw[key]).replace(',', '.')
            raw['SpectralWidth'] = 1 / (float(raw[key]) * 1e-9)
    for key in ('TE', 'TR'):
        if key in raw:                                     # XA: an array, first entry
            raw[key] = float(str(raw[key]).replace(',', '.').split()[0]) / 1e3
    return raw, 'Siemens', '.dat'


def _siemens_dicom(datasets):
    d = datasets[0]
    m, v = d.meta, d.meta.vendor
    dcm = v.get('dicom', {})
    raw = {'Manufacturer': 'Siemens', 'Model': m.get('ManufacturersModelName'),
           'SoftwareVersion': m.get('SoftwareVersions'), 'Nucleus': _first(m.get('ResonantNucleus')),
           'Sequence': m.get('ProtocolName', ''), 'dwelltime': d.dwelltime,
           'VectorSize': int(d.shape[3]), 'FieldStrength': _first(m.get('SpectrometerFrequency'))}
    raw['ap_size'], raw['lr_size'], raw['cc_size'] = _zooms(d.affine)
    for key, name in (('TR', 'RepetitionTime'), ('TE', 'EchoTime'), ('TI', 'InversionTime')):
        if m.get(name) is not None:
            raw[key] = float(m[name]) * 1e3
    raw['_ascconv'] = dict(v.get('MrPhoenixProtocol', {}))
    xa = str(dcm.get('SOPClassUID', '')) == _XA_SPECTROSCOPY
    raw['SW'] = 1 / raw['dwelltime'] if xa else 1 / (raw['dwelltime'] / 2)   # VA-VE: oversampled
    averages = v.get('csa_image', {}).get('NumberOfAverages', dcm.get('NumberOfAverages'))
    if averages is None:                                   # XA: the tag sits in a nested sequence
        try:
            import pydicom
            ds = pydicom.dcmread(str(d.provenance.paths[0]), stop_before_pixels=True)
            averages = next(e.value for e in ds.iterall() if e.tag == 0x00180083)
        except Exception:                                  # noqa: BLE001 - stays blank
            pass
    if averages is not None:
        raw['NumberOfAverages'] = float(averages)
    return raw, 'Siemens', '.ima'


def _rda(datasets):
    raw = dict(datasets[0].meta.vendor)
    if 'DwellTime' in raw and 'SpectralWidth' not in raw:
        raw['DwellTime'] = str(raw['DwellTime']).replace(',', '.')
        raw['SW'] = 1 / ((float(raw['DwellTime']) / 2) / 1e6)          # 2x oversampling
    raw['Manufacturer'] = 'Siemens'
    return raw, 'Siemens', '.rda'


def _spar(datasets):
    raw = {'Manufacturer': 'Philips'}
    raw.update(datasets[0].meta.vendor)
    try:                                                   # one transient per row: the rows are
        if (int(float(raw.get('dim2_pnts', 1))) == 1 and float(raw['averages']) == 1   # the averages
                and int(float(raw.get('rows', 1))) > 1):
            raw['averages'] = int(float(raw['rows']))
    except (KeyError, TypeError, ValueError):
        pass
    return raw, 'Philips', '.spar'


def _ge(datasets):
    raw = {'Manufacturer': 'GE'}
    raw.update(datasets[0].meta.vendor)
    for key in ('rhi_tr', 'rhi_te'):
        raw[key] = raw[key] / 1e3                          # ms
    return raw, 'GE', '.7'


def _bruker_2dseq(datasets):
    m = datasets[0].meta
    raw = {'Manufacturer': 'Bruker', 'FieldStrength': m.get('FieldStrength'),
           'Sequence': m.get('SequenceName'), 'Averages': m.get('NumberOfAverages'),
           'SoftwareVersion': m.get('SoftwareVersions')}
    for key, name in (('TR', 'RepetitionTime'), ('TE', 'EchoTime')):
        if m.get(name) is not None:
            raw[key] = float(m[name]) * 1e3
    return raw, 'Bruker', '2dseq'


def _nifti_style(datasets):
    """NIfTI-MRS keys (ms, as the NIfTI branch makes them) for a format REMY has no column for."""
    d = datasets[0]
    raw = {str(k): str(x[0]) if isinstance(x, list) and len(x) == 1 else x for k, x in d.meta.items()}
    for key in ('EchoTime', 'RepetitionTime', 'InversionTime', 'MixingTime'):
        if key in raw:
            try:
                raw[key] = float(raw[key]) * 1e3
            except (TypeError, ValueError):
                pass
    raw.setdefault('NumberOfDatapoints', int(d.shape[3]))
    import nibabel as nib
    axis = {'L': 'lr_size', 'R': 'lr_size', 'A': 'ap_size', 'P': 'ap_size', 'S': 'cc_size', 'I': 'cc_size'}
    for code, size in zip(nib.aff2axcodes(d.affine), _zooms(d.affine)):
        raw[axis[code]] = size
    if 'DIM_DYN' in d.dims:
        raw['averages'] = int(d.shape[d.dims.index('DIM_DYN')])
    return raw, 'NIfTI', '.nii'


_BUILDERS = {'siemens_twix': _twix, 'siemens_rda': _rda, 'philips_spar_sdat': _spar,
             'ge_pfile': _ge, 'bruker': _bruker_2dseq}
