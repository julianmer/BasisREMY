####################################################################################################
#                                          win_basis.py                                            #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 07/10/26                                                                                #
#                                                                                                  #
# Purpose: Published reference basis sets from the WIN MRS basis-set repository (Will Clarke,      #
#          Oxford; external 'win_mrs_basis'). Every basis JSON there is an FSL-MRS simulation      #
#          output and carries the full sequence description it was simulated from (real RF         #
#          waveforms, delays, B0, readout) next to the FID, so one folder is both a sequence       #
#          input and the reference result for it.                                                  #
#                                                                                                  #
#          basis_sets/ is several hundred MB of Git LFS files, so it is not part of the default    #
#          sparse fetch: folders() lists it from the git tree (nothing downloaded) and fetch()     #
#          checks out one folder and pulls its LFS files on demand (needs git-lfs).                #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import glob
import json
import os
import re
import subprocess

import numpy as np

from basisremy.core import externals

NAME = 'win_mrs_basis'
_LFS_POINTER = b'version https://git-lfs'


def _git(*args) -> str:
    repo = externals.ensure(NAME)
    try:
        return subprocess.run(['git', '-C', repo, *args], check=True, capture_output=True, text=True,
                              env={**os.environ, 'GIT_LFS_SKIP_SMUDGE': '1'}).stdout
    except subprocess.CalledProcessError as exc:
        raise externals.ExternalFetchError(
            f"git {' '.join(args)} failed in {repo}: {exc.stderr.strip()}") from exc


def folders() -> list[str]:
    """Folders under basis_sets/ that hold JSON files, from the git tree (nothing downloaded)."""
    paths = _git('ls-tree', '-r', '--name-only', 'HEAD', 'basis_sets').splitlines()
    return sorted({os.path.dirname(p) for p in paths if p.endswith('.json')})


def fetch(folder: str, files: str = '*.json') -> str:
    """Check out ``folder`` of basis_sets/, download its LFS ``files`` and return the local path."""
    _git('sparse-checkout', 'add', folder)
    repo = externals.ensure(NAME)
    path = os.path.join(repo, folder)
    for f in glob.glob(os.path.join(path, files)):
        with open(f, 'rb') as fh:
            pointer = fh.read()
        if not pointer.startswith(_LFS_POINTER):
            continue
        # One LFS download per file, checked against the pointer's hash. 'git lfs pull'
        # instead scans the whole tree with object sizes, which in this blob-less clone
        # fetches every blob separately (the server answers 429).
        res = subprocess.run(['git', '-C', repo, 'lfs', 'smudge', '--', os.path.relpath(f, repo)],
                             input=pointer, capture_output=True)
        if res.returncode or res.stdout.startswith(_LFS_POINTER):
            raise externals.ExternalFetchError(
                f"{f} is still a Git LFS pointer file (git-lfs needed): {res.stderr.decode().strip()}")
        with open(f, 'wb') as fh:
            fh.write(res.stdout)
    return path


def options(description: str) -> dict:
    """The generate.py options recorded in a sequence description (JSON dict or key=value list)."""
    m = re.search(r'Options: (.*)\.\s*$', description or '', re.S)
    if not m:
        return {}
    try:
        return json.loads(m.group(1))
    except ValueError:
        return dict(kv.split('=', 1) for kv in m.group(1).split(', ') if '=' in kv)


def te_ms(opts: dict) -> float | None:
    """TE from the options: one value, or the sLASER TE1/TE2/TE3 list summed."""
    if 'te' not in opts:
        return None
    te = json.loads(str(opts['te']))
    return float(sum(te)) if isinstance(te, list) else float(te)


def describe(folder: str) -> dict:
    """Sequence, B0, TE, readout and RF count of ``folder``, from its first JSON (one LFS file)."""
    names = [p for p in _git('ls-tree', '--name-only', 'HEAD', f'{folder}/').splitlines()
             if p.endswith('.json')]
    path = fetch(folder, os.path.basename(names[0]))
    with open(os.path.join(path, os.path.basename(names[0]))) as fh:
        seq = json.load(fh).get('seq') or {}
    opts = options(seq.get('description', ''))
    name = seq.get('sequenceName')            # a list in newer files, a string in older ones
    return {'folder': folder, 'sequence': name[0] if isinstance(name, list) else name, 'b0': seq.get('B0'),
            'te_ms': te_ms(opts), 'points': seq.get('Rx_Points'), 'sw': seq.get('Rx_SW'),
            'n_rf': len(seq.get('RF', [])), 'options': opts}


def load(folder: str) -> dict[str, dict]:
    """Metabolite -> {'fid', 'dwell_s', 'centre_mhz', 'seq', 'spins'} for every basis JSON of ``folder``."""
    return read(fetch(folder))


def read(path: str) -> dict[str, dict]:
    """As load(), for a basis folder already on disk."""
    out = {}
    for f in sorted(glob.glob(os.path.join(path, '*.json'))):
        with open(f) as fh:
            d = json.load(fh)
        if 'basis' not in d:
            continue   # e.g. a spin-system or sequence file next to the spectra
        b = d['basis']
        out[b['basis_name']] = {'fid': np.asarray(b['basis_re']) + 1j * np.asarray(b['basis_im']),
                                'dwell_s': b['basis_dwell'], 'centre_mhz': b['basis_centre'],
                                'seq': d.get('seq'), 'spins': d.get('spinSys')}
    return out
