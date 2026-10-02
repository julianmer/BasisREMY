####################################################################################################
#                                          capability.py                                           #
####################################################################################################
#                                                                                                  #
# Purpose: What can BasisREMY be used for? For every local test dataset and every simulation       #
#          entry (engine x pulse model) measure three levels: the file reads, the fields that      #
#          entry marks mandatory are all filled from the header ("complete"), and one metabolite   #
#          simulates - with no typing ("auto") or after the blanks were filled with the            #
#          completion rules below ("completed"). Rows go to a CSV that the poster table and the    #
#          per-field panel are rendered from.                                                      #
#                                                                                                  #
#          python -m tests.validation.capability [--no-sim] [--engines FID-A,spant] [--out x.csv]  #
#                                                                                                  #
####################################################################################################

import argparse
import csv
import glob
import os
import re
import sys
import time
import traceback

from basisremy.core.basisremy import BasisREMY
from basisremy.core.metabolite_identity import translate_metabolites

ROOT = 'example_data/REMY_tests'
BLANK = (None, '', 'missing input', 'select option')

# ---- datasets -------------------------------------------------------------------------------------
FORMATS = [   # (format label, vendor, glob under a dataset folder)
    ('GE P-file',      'GE',      '**/*.7'),
    ('Philips SPAR',   'Philips', '**/*Act.SPAR'),
    ('Siemens twix',   'Siemens', '**/*.dat'),
    ('Siemens DICOM',  'Siemens', '**/*.IMA'),
    ('Siemens RDA',    'Siemens', '**/*.rda'),
    ('Bruker method',  'Bruker',  '**/*method'),
]


def discover():
    """One representative file per dataset folder (+ each NIfTI-MRS file)."""
    rows = []
    for d in sorted(glob.glob(f'{ROOT}/Dataset_*')):
        for fmt, vendor, pat in FORMATS:
            files = [f for f in sorted(glob.glob(os.path.join(d, pat), recursive=True))
                     if '_ecc' not in f and '_quant' not in f and not f.endswith(('.pdf', '.tex', '.log', '.csv'))]
            if files:
                rows.append({'dataset': os.path.basename(d), 'format': fmt, 'vendor': vendor, 'file': files[0]})
                break
    for f in sorted(glob.glob(f'{ROOT}/Datasets_nifti/*.nii.gz')):
        if '_ecc' in f or '_quant' in f:            # companions of the same dataset
            continue
        rows.append({'dataset': os.path.basename(f), 'format': 'NIfTI-MRS', 'vendor': 'NIfTI-MRS', 'file': f})
    return rows


# ---- sequence recogniser (harness only; the product one lives in Gusteau) --------------------------
_SEQ = [
    ('MEGA-sLASER', r'mega.?s?laser|mslaser'),
    ('MEGA-PRESS',  r'mega|mpress|meshcher'),
    ('HERMES',      r'hermes'),
    ('HERCULES',    r'hercules'),
    ('sLASER',      r'slaser|oslaser|semi.?laser'),
    ('LASER',       r'(?<!s)laser'),
    ('SPECIAL',     r'special'),
    ('ISIS',        r'isis'),
    ('CSI',         r'csi'),
    ('STEAM',       r'steam'),
    ('PRESS',       r'press|svs_se|svs_edit'),
]


def recognise(*texts):
    """Canonical sequence from the protocol name (or REMY's own Sequence field)."""
    for text in texts:
        if text:
            t = str(text).lower()
            for name, pat in _SEQ:
                if re.search(pat, t):
                    return name
    return 'unknown'


# ---- engine entries ----------------------------------------------------------------------------------
# label, pulse model, backend, {canonical sequence: params to set on the sheet}
ENTRIES = [
    ('FID-A',    'ideal', 'FidaIdeal',           {'PRESS': {'Sequence': 'PRESS'}, 'STEAM': {'Sequence': 'STEAM'}, 'LASER': {'Sequence': 'LASER'}}),
    ('FID-A',    'ideal', 'FidaMegaPressIdeal',  {'MEGA-PRESS': {}}),
    ('FID-A',    'real',  'FidaPressShaped',     {'PRESS': {}}),
    ('FID-A',    'real',  'FidaSteamShaped',     {'STEAM': {}}),
    ('FID-A',    'real',  'FidaSemiLaserShaped', {'sLASER': {}}),
    ('FID-A',    'real',  'FidaMegaPressShaped', {'MEGA-PRESS': {}}),
    ('FSL-MRS',  'ideal', 'FSL-MRS',             {s: {'Sequence': s} for s in ('PRESS', 'STEAM', 'LASER', 'sLASER', 'MEGA-PRESS', 'HERMES', 'HERCULES', 'MEGA-sLASER')}),
    ('MRSCloud', 'real',  'MRSCloud',            {'PRESS': {'Sequence': 'UnEdited', 'Localization': 'PRESS'},
                                                  'sLASER': {'Sequence': 'UnEdited', 'Localization': 'sLASER'},
                                                  'MEGA-PRESS': {'Sequence': 'MEGA', 'Localization': 'PRESS'},
                                                  'MEGA-sLASER': {'Sequence': 'MEGA', 'Localization': 'sLASER'},
                                                  'HERMES': {'Sequence': 'HERMES', 'Localization': 'PRESS'}}),
    ('Vespa',    'ideal', 'Vespa',               {'PRESS': {'Sequence': 'PRESS'}, 'STEAM': {'Sequence': 'STEAM'}}),
    ('Vespa',    'real',  'Vespa',               {'PRESS': {'Sequence': 'PRESS shaped'}}),
    ('spant',    'ideal', 'Spant',               {'PRESS': {'Sequence': 'PRESS'}, 'STEAM': {'Sequence': 'STEAM'}, 'sLASER': {'Sequence': 'sLASER'}, 'MEGA-PRESS': {'Sequence': 'MEGA-PRESS'}}),
    ('spant',    'real',  'Spant',               {'PRESS': {'Sequence': 'PRESS shaped'}}),
    ('Spinach',  'ideal', 'Spinach',             {'PRESS': {'Sequence': 'PRESS'}, 'STEAM': {'Sequence': 'STEAM'}, 'LASER': {'Sequence': 'LASER'}}),
    ('Spinach',  'real',  'SpinachPressShaped',  {'PRESS': {}}),
    ('Spinach',  'real',  'SpinachSemiLaserShaped', {'sLASER': {}}),
    ('jbss',     'real',  'CustomSLaser',        {'sLASER': {}}),
]

# completion values for blanks (what a user would type / pick); recorded per row
COMPLETION = {
    'TE': 30.0, 'Samples': 2048, 'Bandwidth': 2000.0, 'Bfield': 3.0, 'TM': 10.0, 'Nucleus': '1H',
    'Linewidth': 1.0, 'Edit Pulse Path': 'standard:gauss-edit', 'TR': 2000.0,
}
# speed only (not counted as completion): tiny spatial grids, one metabolite
SPEED = {'nX': 2, 'nY': 2, 'fovX': 1.5, 'fovY': 1.5, 'Spatial Points': 5}


def _blank(v):
    return v is None or (isinstance(v, str) and v.strip().lower() in BLANK)


def complete(params, blanks, sequence, vendor):
    """Fill the listed blanks; return (params, list of filled keys)."""
    filled = []
    params = dict(params)
    for k in blanks:
        v = params.get(k)
        if k == 'Sequence' or not _blank(v):
            if k == 'Sequence':
                filled.append(k)                   # the entry's value is already set
            continue
        if k == 'Path to Pulse':
            v = 'standard:goia-wurst' if sequence in ('sLASER', 'LASER', 'MEGA-sLASER') else 'standard:sinc-ref'
        elif k == 'System':
            v = vendor if vendor in ('GE', 'Philips', 'Siemens') else 'Siemens'
        elif k == 'Center Freq':
            v = float(params.get('Bfield') or COMPLETION['Bfield']) * 42.577
        elif k in ('Tau 1', 'Tau 2'):
            v = float(params.get('TE') or COMPLETION['TE']) / 2.0
        elif k in COMPLETION:
            v = COMPLETION[k]
        else:
            continue                                # stays blank -> the run will raise
        params[k] = v
        filled.append(k)
    return params, filled


def evaluate(ds, entry, simulate=True):
    label, model, backend, seqs = entry
    row = {'dataset': ds['dataset'], 'format': ds['format'], 'vendor': ds['vendor'],
           'engine': label, 'pulse_model': model, 'backend': backend, 'read': False,
           'sequence': '', 'protocol': '', 'supported': False, 'shown': '', 'blanks': '', 'filled': '',
           'level': '', 'error': '', 'seconds': 0.0}
    t0 = time.time()
    try:
        br = BasisREMY(backend)                     # fresh instance: no carry-over
        br.runREMY(import_fpath=ds['file'])
        row['read'] = True
        m = br._last_mrsinmrs
        row['protocol'] = str(m.get('Protocol') or '')
        br.set_backend(backend)
        b = br.backend
        row['sequence'] = recognise(row['protocol'], m.get('Sequence'))
    except Exception as exc:                        # noqa: BLE001
        row.update(level='read-fail', error=f'{type(exc).__name__}: {str(exc)[:160]}', seconds=time.time() - t0)
        return row
    if row['sequence'] not in seqs:
        row.update(level='unsupported', seconds=time.time() - t0)
        return row
    row['supported'] = True
    seq_params = seqs[row['sequence']]
    # the sequence counts as extracted when the backend's mapped value names the same
    # sequence (choosing the pulse model / backend variant is not a blank)
    mapped = f"{b.mandatory_params.get('Sequence') or ''} {b.mandatory_params.get('Localization') or ''}"
    seq_blank = 'Sequence' in seq_params and recognise(mapped) != row['sequence']
    b.mandatory_params.update(seq_params)
    visible = b.get_params_for_mode()               # what the GUI shows and validates
    row['shown'] = ';'.join(sorted(k for k in visible if k != 'Metabolites'))
    blanks = sorted(k for k, v in visible.items() if k != 'Metabolites' and _blank(v))
    if seq_blank:
        blanks = sorted(set(blanks) | {'Sequence'})
    row['blanks'] = ';'.join(blanks)
    params = dict(b.mandatory_params)
    if not simulate:
        row.update(level='complete' if not blanks else 'blank', seconds=time.time() - t0)
        return row
    params, filled = complete(params, blanks, row['sequence'], ds['vendor'])
    row['filled'] = ';'.join(filled)
    params.update({k: v for k, v in SPEED.items() if k in params})
    keys = translate_metabolites(['NAA'], b.metabs.keys())
    params['Metabolites'] = keys or list(b.metabs.keys())[:1]
    try:
        out = b.run_simulation(params)
        ok = bool(out) and all(len(v) > 0 for v in out.values())
        row['level'] = ('auto' if not filled else 'completed') if ok else 'failed'
        if not ok:                                  # MRSCloud records per-metabolite failures
            fails = getattr(b, 'last_failures', None) or {}
            row['error'] = '; '.join(f'{k}: {v}' for k, v in fails.items())[:160] or 'empty result'
    except Exception as exc:                        # noqa: BLE001
        row.update(level='failed', error=f'{type(exc).__name__}: {str(exc)[:160]}')
        traceback.print_exc(limit=1)
    row['seconds'] = time.time() - t0
    return row


FIELDS = ['dataset', 'format', 'vendor', 'engine', 'pulse_model', 'backend', 'read', 'sequence',
          'protocol', 'supported', 'shown', 'blanks', 'filled', 'level', 'error', 'seconds']


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='tests/validation/out/capability.csv')
    ap.add_argument('--no-sim', action='store_true', help='extraction levels only (fast)')
    ap.add_argument('--engines', default='', help='comma list of engine labels to run')
    ap.add_argument('--datasets', default='', help='substring filter on the dataset name')
    ap.add_argument('--fresh', action='store_true', help='ignore rows already in the CSV')
    a = ap.parse_args(argv)
    engines = [e for e in a.engines.split(',') if e]
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    done, fields = set(), FIELDS
    if os.path.exists(a.out) and not a.fresh:
        with open(a.out) as f:
            rd = csv.DictReader(f)
            done = {(r['dataset'], r['backend'], r['pulse_model']) for r in rd}
            fields = rd.fieldnames or FIELDS          # append in the file's own column order
    new = not done
    with open(a.out, 'w' if new else 'a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        if new:
            w.writeheader()
        for ds in discover():
            if a.datasets and a.datasets not in ds['dataset']:
                continue
            for entry in ENTRIES:
                if engines and entry[0] not in engines:
                    continue
                key = (ds['dataset'], entry[2], entry[1])
                if key in done:
                    continue
                row = evaluate(ds, entry, simulate=not a.no_sim)
                w.writerow(row); f.flush()
                print(f"{row['dataset'][:34]:34s} {row['engine']:8s} {row['pulse_model']:5s} "
                      f"{row['sequence']:11s} {row['level']:11s} blanks=[{row['blanks']}] "
                      f"filled=[{row['filled']}] {row['error'][:60]} ({row['seconds']:.0f} s)", flush=True)


if __name__ == '__main__':
    sys.exit(main())
