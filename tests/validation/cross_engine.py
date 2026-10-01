####################################################################################################
#                                         cross_engine.py                                          #
####################################################################################################
#                                                                                                  #
# Purpose: Cross-engine validation. Simulate one parameter sheet on several backends and compare   #
#          the metabolites pairwise: Pearson r of the real spectra over a ppm window, the          #
#          least-squares amplitude ratio, the normalised RMS difference after scaling, and the     #
#          position of the largest peak. The MRSCloud paper (Hui 2022) validated its engine the    #
#          same way against FID-A and MARSS.                                                       #
#                                                                                                  #
####################################################################################################

import os
import time

import numpy as np

from basisremy.core.metabolite_identity import translate_metabolites

PANEL = ['NAA', 'Cr', 'PCh', 'Ins', 'Glu', 'Gln', 'GABA', 'GSH', 'Lac']

# Extra sheet values a backend needs beyond the common ones.
_EXTRA = {
    'FidaIdeal': {'TE2': 0},
    'FidaMegaPressIdeal': {'Edit On': 1.9, 'Edit Bandwidth (ppm)': 1.0},
    'FSL-MRS': {'Edit On': 1.9, 'Edit Off': 7.5, 'Edit Tp': 14.0},
    'Spant': {'Edit On': 1.9, 'Edit Off': 7.5, 'Edit Bandwidth (Hz)': 80.0},
}


def run_engine(br, backend, sequence, te, metabs=PANEL, b0=3.0, samples=2048,
               bandwidth=2000.0, linewidth=1.0, tm=10.0, extra=None):
    """Simulate *metabs* on *backend*; return ({metabolite: FID}, seconds).

    Metabolite names are translated into the backend's own keys and the
    result is keyed by the panel name again (sub-spectrum tags kept)."""
    br.set_backend(backend)
    b = br.backend
    keys = translate_metabolites(metabs, b.metabs.keys())
    back = {k: m for k, m in zip(keys, metabs)}
    params = dict(b.mandatory_params)
    params.update({k: v for k, v in b.optional_params.items() if k not in params})
    params.update({'Sequence': sequence, 'TE': te, 'Bfield': b0, 'Samples': samples,
                   'Bandwidth': bandwidth, 'Center Freq': b0 * 42.577,
                   'Linewidth': linewidth, 'TM': tm, 'Metabolites': keys})
    params.update(_EXTRA.get(backend, {}))
    params.update(extra or {})
    t0 = time.time()
    out = b.run_simulation(params)
    dt = time.time() - t0
    result = {}
    for name, fid in out.items():
        base, _, tag = name.partition(' (')
        if base in back:
            result[back[base] + (' (' + tag if tag else '')] = np.asarray(fid, dtype=complex)
    return result, dt


def spectrum(fid, bandwidth, b0, window=(0.2, 4.2)):
    """Real spectrum on the GUI's ppm axis, restricted to *window*."""
    n = fid.size
    spec = np.fft.fftshift(np.fft.fft(fid))
    ppm = np.linspace(-bandwidth / 2, bandwidth / 2, n) / (b0 * 42.577) + 4.65
    sel = (ppm >= window[0]) & (ppm <= window[1])
    return ppm[sel], spec[sel]


def compare(fid_a, fid_b, bandwidth, b0, window=(0.2, 4.2)):
    """Agreement of b with a: r (real parts), scale (least squares b -> a),
    nrmsd after scaling, peak positions [ppm]."""
    ppm, sa = spectrum(fid_a, bandwidth, b0, window)
    _, sb = spectrum(fid_b, bandwidth, b0, window)
    ra, rb = sa.real, sb.real
    r = float(np.corrcoef(ra, rb)[0, 1])
    scale = float(np.dot(rb, ra) / np.dot(rb, rb)) if np.dot(rb, rb) > 0 else np.nan
    nrmsd = float(np.linalg.norm(ra - scale * rb) / np.linalg.norm(ra))
    return {'r': r, 'scale': scale, 'nrmsd': nrmsd,
            'peak_a': float(ppm[np.argmax(np.abs(sa))]),
            'peak_b': float(ppm[np.argmax(np.abs(sb))])}


def agreement_table(results, reference, bandwidth, b0, window=(0.2, 4.2)):
    """results: {engine: {metab: fid}} -> {engine: {metab: compare(...)}} vs reference."""
    table = {}
    for engine, fids in results.items():
        if engine == reference:
            continue
        table[engine] = {m: compare(results[reference][m], fids[m], bandwidth, b0, window)
                         for m in fids if m in results[reference]}
    return table


def format_table(table, timings=None):
    lines = []
    for engine, rows in table.items():
        head = f"## {engine}" + (f"  ({timings[engine]:.1f} s)" if timings and engine in timings else "")
        lines += [head, "| metabolite | r | scale | NRMSD | peak ref | peak |", "|---|---|---|---|---|---|"]
        for m, c in rows.items():
            lines.append(f"| {m} | {c['r']:.4f} | {c['scale']:.3f} | {c['nrmsd']:.3f} | "
                         f"{c['peak_a']:.2f} | {c['peak_b']:.2f} |")
        lines.append("")
    return "\n".join(lines)
