####################################################################################################
#                                         design_fields.py                                         #
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 09/10/26                                                                                #
#                                                                                                  #
# Purpose: The sequence-design values a data header holds beyond the scan parameters: mixing       #
#          time, editing pulses and RF pulses. Only slots whose meaning is confirmed are read;     #
#          the rest stays blank and the sheet recommends a value. The sheet keys (TM, Edit On,     #
#          Edit Off, Edit Tp) go onto the sheet as values from the file; the rest is information.  #
####################################################################################################

from __future__ import annotations

import re

# sheet keys a header can fill (the names sequence_setup.recommended uses)
SHEET_KEYS = ('TM', 'Edit On', 'Edit Off', 'Edit Tp')
# information keys (not sheet values)
EDIT_PPM, RF_PULSES = 'Edit frequencies (header, ppm)', 'RF pulses (header)'
# {name: {'role', 'dur_ms', 'bw_hz', 'flip', 'shape', 'waveform'}}, waveform [(amplitude %, phase deg)]
# or None; for the sequence designer
PULSE_PARAMS = 'RF pulse parameters (header)'
# {'max_mT_m', 'rise_ms'}: the scanner's gradient system (strength, rise time 0 -> full); for the designer
GRADIENTS = 'Gradient system (header)'


def _num(v):
    try:
        return float(str(v).strip().rstrip(';').strip())
    except (TypeError, ValueError):
        return None


#**************************************************************************************************#
#                                     Siemens protocol (ASCCONV)                                   #
#**************************************************************************************************#
# One flat form for twix and DICOM: 'sWipMemBlock.alFree[7]', 'sTXSPEC.aRFPULSE[0].tName'.
def flat_twix(phoenix) -> dict:
    """mapVBVD's Phoenix header (tuple keys) in ASCCONV notation."""
    out = {}
    for key, value in phoenix.items():
        name = ''
        for part in (key if isinstance(key, tuple) else (key,)):
            name += f'[{part}]' if str(part).isdigit() else ('.' if name else '') + str(part)
        out[name] = value
    return out


def ascconv(text: str) -> dict:
    """The ASCCONV block of a Siemens protocol text (DICOM CSA / XA private tag) as a dict."""
    m = re.search(r'### ASCCONV BEGIN.*?###(.*?)### ASCCONV END ###', text, re.S)
    out = {}
    for line in (m.group(1) if m else '').splitlines():
        if '=' not in line:
            continue
        key, value = (s.strip() for s in line.split('=', 1))
        value = value.split('#')[0].strip()
        if value.startswith('"'):
            value = value.strip('"')                       # DICOM doubles the quotes
        elif value.lower().startswith('0x'):
            value = int(value, 16)
        else:
            value = _num(value)
        out[key] = value
    return out


def ascconv_from_dicom(dcm) -> dict:
    """The protocol of a Siemens DICOM: CSA series header (VA-VE) or the XA private tag
    (5200,9229)[0] (0021,10FE)[0] (0021,1019)."""
    try:
        from nibabel.nicom import csareader
        csa = csareader.get_csa_header(dcm, 'series')
        if csa is not None and 'MrPhoenixProtocol' in csa['tags']:
            return ascconv(str(csa['tags']['MrPhoenixProtocol']['items'][0]))
    except Exception:                                      # noqa: BLE001
        pass
    try:
        raw = dcm[0x5200, 0x9229][0][0x0021, 0x10FE][0][0x0021, 0x1019].value
        return ascconv(raw.decode('latin-1') if isinstance(raw, bytes) else str(raw))
    except Exception:                                      # noqa: BLE001
        return {}


def _wip(prot, kind, i):
    for block in ('sWipMemBlock', 'sWiPMemBlock'):         # VE/XA, VB
        if f'{block}.{kind}[{i}]' in prot:
            return _num(prot[f'{block}.{kind}[{i}]'])
    return None


def siemens(prot: dict) -> dict:
    """Editing of the mgs_svs_ed family (svs_edit_mgs_univ, smm_svs_herc), slots as spec2nii
    reads them: alFree[7] mode (absent = 0, MEGA), alFree[12] pulse duration (us),
    adFree[8..11] frequencies (ppm); and the RF pulse names."""
    out = {}
    seq = str(prot.get('tSequenceFileName', ''))
    if re.search(r'mgs_svs_ed|svs_edit_mgs|smm_svs_herc', seq, re.I):
        mode = _wip(prot, 'alFree', 7) or 0.0
        f = [_wip(prot, 'adFree', i) for i in (8, 9, 10, 11)]
        tp = _wip(prot, 'alFree', 12)
        if mode == 0 and f[0] and f[3]:                    # MEGA: ON adFree[8], OFF adFree[11]
            out.update({'Edit On': f[0], 'Edit Off': f[3]})
            if tp and tp > 0:
                out['Edit Tp'] = tp / 1e3
        elif any(f):                                       # HERMES / HERCULES: the frequencies;
            out[EDIT_PPM] = [v for v in f if v]            # duration not confirmed (spec2nii: 20 ms)
    names = [str(v).strip('"') for k, v in prot.items()
             if re.fullmatch(r'sTXSPEC\.aRFPULSE\[\d+\]\.tName', k) and str(v).strip('"')]
    if names:
        out[RF_PULSES] = names
    return out


#**************************************************************************************************#
#                                        other vendors                                             #
#**************************************************************************************************#
def ge(hdr: dict) -> dict:
    """GE MEGA-PRESS psd 'gaba' (Big GABA): rhi_user20 / 21 edit offsets (Hz from the transmitter,
    taken at 4.68 ppm), rhi_user22 pulse width (us); as spec2nii and Gannet read them. Other
    psds use these user CVs for other things, so nothing else is read."""
    psd = hdr.get('rhi_psdname', b'')
    psd = (psd.decode('latin-1') if isinstance(psd, bytes) else str(psd)).strip('\x00').strip().lower()
    on, off, tp = (_num(hdr.get(k)) for k in ('rhi_user20', 'rhi_user21', 'rhi_user22'))
    mhz = (_num(hdr.get('rhr_rh_ps_mps_freq')) or 0) * 1e-7
    if psd != 'gaba' or not on or not off or not mhz:
        return {}
    out = {'Edit On': round(4.68 + on / mhz, 2), 'Edit Off': round(4.68 + off / mhz, 2)}
    if tp and tp > 0:
        out['Edit Tp'] = tp / 1e3
    return out


_BRUKER_PULSES = ('VoxPul1', 'VoxPul2', 'VoxPul3', 'ExcPulse1', 'RefPulse1')
_BRUKER_ROLE = {0: 'exc', 1: 'ref', 2: 'inv'}             # PVM_RF_PULSE Type


def _bruker_waveform(raw):
    """A method's inline shape '( 2N ) a1 p1 a2 p2 ...' as N (amplitude %, phase deg) pairs; checked
    on Dataset_00: the integral factor of the pairs equals the struct's Sint (0.236152) and the
    Bloch bandwidth 8840 Hz the struct's 8400 Hz."""
    m = re.match(r'\s*\(\s*(\d+)\s*\)', str(raw or ''))
    if not m:
        return None
    vals = [float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', str(raw)[m.end():])]
    n = int(m.group(1))
    if n < 4 or n % 2 or len(vals) < n:
        return None
    return [(vals[i], vals[i + 1]) for i in range(0, n, 2)]


def bruker(method: dict) -> dict:
    """Bruker method: StTM (ms), the gradient system (PVM_GradCalConst, PVM_RiseTime) and the PVM
    pulse structs (length ms, bandwidth Hz, flip angle, ..., type, ..., shape), with the shape name
    from <name>Enum and the waveform when the method stores it."""
    out = {}
    tm = _num(str(method.get('$StTM', '')).split(';')[0])
    if tm and tm > 0:
        out['TM'] = tm
    cal, rise = _num(method.get('$PVM_GradCalConst')), _num(method.get('$PVM_RiseTime'))  # Hz/mm, ms
    if cal and rise:
        out[GRADIENTS] = {'max_mT_m': cal * 1e6 / 42.577e6, 'rise_ms': rise}
    pulses, params = [], {}
    for name in _BRUKER_PULSES:
        raw = method.get(f'${name}')
        if not raw:
            continue
        fields = [s.strip() for s in str(raw).replace('; ', ' ').strip('(); ').split(',')]
        length, bw, flip = (_num(v) for v in fields[:3])
        kind = _num(fields[8]) if len(fields) > 8 else None
        shape = str(method.get(f'${name}Enum', '')).split(';')[0].strip('<> ')
        if not (length and bw):
            continue
        pulses.append(f'{name}: {length:g} ms, {bw:g} Hz, {flip:g} deg' + (f', {shape}' if shape else ''))
        params[name] = {'role': _BRUKER_ROLE.get(int(kind)) if kind is not None else None,
                        'dur_ms': length, 'bw_hz': bw, 'flip': flip, 'shape': shape,
                        'waveform': _bruker_waveform(method.get(f'${name}Shape'))}
    if pulses:
        out[RF_PULSES] = pulses
        out[PULSE_PARAMS] = params
    return out


def rda(hdr: dict) -> dict:
    tm = _num(hdr.get('TM'))
    return {'TM': tm} if tm and tm > 0 else {}


def nifti(hdr: dict) -> dict:
    """NIfTI-MRS standard keys: MixingTime (already in ms here) and EditPulse
    ({condition: {PulseOffset: ppm, PulseDuration: s}}; ON / OFF read as MEGA)."""
    out = {}
    tm = _num(hdr.get('MixingTime'))
    if tm and tm > 0:
        out['TM'] = tm
    edit = hdr.get('EditPulse')
    if isinstance(edit, dict) and edit:
        on, off = edit.get('ON'), edit.get('OFF')
        if isinstance(on, dict) and isinstance(off, dict) and _num(on.get('PulseOffset')) is not None \
                and _num(off.get('PulseOffset')) is not None:
            out.update({'Edit On': _num(on['PulseOffset']), 'Edit Off': _num(off['PulseOffset'])})
            if _num(on.get('PulseDuration')):
                out['Edit Tp'] = _num(on['PulseDuration']) * 1e3
        else:
            ppm = []
            for cond in edit.values():
                po = cond.get('PulseOffset') if isinstance(cond, dict) else None
                ppm += [_num(v) for v in (po if isinstance(po, list) else [po]) if _num(v) is not None]
            if ppm:
                out[EDIT_PPM] = sorted(set(ppm), reverse=True)
    return out


def extract(raw: dict, vendor: str, dtype: str) -> dict:
    """Design values of one header (REMY's raw reader output)."""
    vendor = vendor.lower()
    if vendor == 'siemens':
        return rda(raw) if dtype == 'rda' else siemens(raw.get('_ascconv') or {})
    if vendor == 'ge':
        return ge(raw)
    if vendor == 'bruker' and dtype == 'method':
        return bruker(raw)
    if vendor == 'nifti':
        return nifti(raw)
    return {}
