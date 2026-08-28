####################################################################################################
#                                       pygamma_worker.py                                          #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 27/08/26                                                                                #
#                                                                                                  #
# Purpose: Standalone PyGAMMA simulation worker for the Vespa backend. Runs inside the dedicated   #
#          Python 3.9 side-environment (PyGAMMA has no wheels for modern Pythons), so this file    #
#          must stay Python-3.9 compatible and import nothing from basisremy.                      #
#                                                                                                  #
#          Usage: python pygamma_worker.py <job.json> <out.json>                                   #
#                                                                                                  #
#          The job carries acquisition parameters and per-metabolite spin systems (shifts in ppm,  #
#          J-couplings in Hz). Transition-table binning follows Vespa-Simulation's canonical       #
#          binning code (pulse_sequences.xml, "Based on TTable1D::calc_spectra()").                #
#                                                                                                  #
####################################################################################################

import json
import math
import sys

import numpy as np
import pygamma as pg


def build_spin_system(shifts_ppm, j_hz, cf_mhz, centre_ppm):
    """Spin system with shifts relative to the rotating-frame centre."""
    nspins = len(shifts_ppm)
    system = pg.spin_system(nspins)
    system.Omega(cf_mhz)
    for i in range(nspins):
        system.PPM(i, shifts_ppm[i] - centre_ppm)
    for i in range(nspins):
        for k in range(i + 1, nspins):
            jval = j_hz[i][k]
            if jval:
                system.J(i, k, jval)
    return system


def bin_table(mx, system, cf_mhz):
    """Vespa's canonical transition-table binning -> (ppm_rel, area, phase_deg).

    ppm values are relative to the rotating-frame centre (the spin system was
    built with centred shifts).
    """
    nlines = mx.size()
    obs_qn = pg.Isotope('1H').qn()
    qnscale = system.qnStates().size()
    qnscale = qnscale / (2.0 * (2.0 * obs_qn + 1))

    indx = mx.Sort(0, -1, 0)
    ppms, areas, phases = [], [], []
    for i in range(nlines):
        freq_ppm = -1 * mx.Fr(indx[i]) / (2.0 * math.pi * cf_mhz)
        val = mx.I(indx[i])
        area = math.sqrt(val.real() ** 2 + val.imag() ** 2) / qnscale
        phase = -math.degrees(math.atan2(val.imag(), val.real()))
        ppms.append(freq_ppm)
        areas.append(area)
        phases.append(phase)
    return ppms, areas, phases


def synthesize_fid(ppms, areas, phases, cf_mhz, samples, bandwidth, linewidth):
    """Sum of damped complex exponentials (Vespa's calculate_fid formula).

    The sign convention matches FID-A / the BasisREMY plot: a line at
    ppm_rel (relative to centre) oscillates at -ppm_rel * f0 Hz... with the
    binned freq already carrying Vespa's sign, the peak lands at the correct
    ppm on an axis of +freq/f0 + centre.
    """
    t = np.arange(samples) / float(bandwidth)
    fid = np.zeros(samples, dtype=complex)
    for ppm, area, phase in zip(ppms, areas, phases):
        # binned ppm is relative to the centre; the GUI/FID-A axis is
        # ppm = +freq/f0 + centre, so a line at relative ppm oscillates
        # at +ppm*f0 Hz (verified: 2.0 ppm singlet -> -2.65 ppm rel).
        hz = ppm * cf_mhz
        fid += area * np.exp(1j * (2.0 * math.pi * hz * t
                                   + math.radians(phase)))
    fid *= np.exp(-math.pi * float(linewidth) * t)
    return fid


def steam_sigma(system, te, tm, h_op):
    """Ideal STEAM: 90y - TE/2 - 90x - TM - 90x - TE/2.

    Structure follows Vespa-Simulation's 'STEAM Ideal' sequence
    (pulse_sequences.xml, B. Soher): the TE crushers around the 2nd and 3rd
    pulses are emulated by running that part on four copies rotated about z
    by 0/90/180/270 degrees and averaging, which keeps only the pathways with
    p(before pulse 2) = -p(after pulse 3), i.e. FID-A's sim_COF(+1).

    The TM spoiler is emulated with zero_mqc(sys, sigma, 0, 1), which keeps
    only zero-order coherence during TM (longitudinal + zero-quantum), like
    FID-A's sim_COF(0). Vespa's own code calls zero_mqc(sys, sigma, 2, 0)
    there; measured in this PyGAMMA build that leaves magnetisation that stays
    transverse through TM untouched for multi-spin systems, which the z
    rotations cannot remove either (same net phase as the stimulated echo),
    giving e.g. 0.15 instead of 0.5 of the spin-echo amplitude for NAA's CH3.
    With (0, 1) every uncoupled spin gives exactly 0.5, as it must.
    """
    u_half_te = pg.prop(h_op, te / 2.0)
    u_tm = pg.prop(h_op, tm)
    sigma0 = pg.Iypuls(system, pg.sigma_eq(system), 90.0)
    sigma0 = pg.evolve(sigma0, u_half_te)

    dephase_angles = (0.0, 90.0, 180.0, 270.0)
    sigma_res = None
    for angle in dephase_angles:
        riz = pg.gen_op(pg.Rz(system, angle))
        sigma = pg.evolve(pg.gen_op(sigma0), riz)
        sigma = pg.Ixpuls(system, sigma, 90.0)
        pg.zero_mqc(system, sigma, 0, 1)   # keep only p = 0 during TM
        sigma = pg.evolve(sigma, u_tm)
        sigma = pg.Ixpuls(system, sigma, 90.0)
        sigma = pg.evolve(sigma, riz)
        sigma *= 1.0 / float(len(dephase_angles))
        if sigma_res is None:
            sigma_res = pg.gen_op(sigma)
        else:
            sigma_res += sigma
    return pg.evolve(sigma_res, u_half_te)


def shaped_pulse_propagator(system, pulse):
    """Propagator of a shaped RF pulse: the product over waveform steps of
    exp(-i (H0 + Hrf_k) dt_k), with H0 the chemical-shift + J Hamiltonian and
    Hrf_k = amp_k (cos(phi_k) Fx + sin(phi_k) Fy) in Hz, so off-resonance
    spins see the pulse's real bandwidth (Vespa's recipe builds this through
    pg.PulComposite(...).GetUsum(-1); measured in this PyGAMMA build that
    propagator showed no offset dependence at all, so it is built explicitly
    here). Returns (U, pulse_duration_s)."""
    h0 = pg.Hcs(system) + pg.HJ(system)
    fx = pg.gen_op(pg.Fx(system))
    fy = pg.gen_op(pg.Fy(system))
    u = None
    total = 0.0
    for amp, phase_deg, dt in zip(pulse['amp_hz'], pulse['phase_deg'],
                                  pulse['dt_s']):
        phi = math.radians(float(phase_deg))
        hx = pg.gen_op(fx)
        hx *= float(amp) * math.cos(phi)
        hy = pg.gen_op(fy)
        hy *= float(amp) * math.sin(phi)
        hrf = hx + hy
        u_step = pg.prop(h0 + hrf, float(dt))
        u = u_step if u is None else u_step * u
        total += float(dt)
    return u, total


def _crushed_refocus(system, sigma, u180, nsteps=8):
    """Apply a refocusing propagator between emulated crusher gradients.

    Same device as Vespa's STEAM crushers: the block is run on copies of the
    density matrix rotated about z by 0..360 degrees before and (with the
    same rotation) after the pulse, then averaged. Only pathways with
    p(after) = -p(before) survive, i.e. the refocused magnetisation; what the
    pulse failed to refocus is dephased instead of surviving as a phase error.
    Vespa's own 'PRESS with real 180 pulses' has no crushers (single voxel).
    """
    res = None
    for k in range(nsteps):
        riz = pg.gen_op(pg.Rz(system, 360.0 * k / nsteps))
        s = pg.evolve(pg.gen_op(sigma), riz)
        s = pg.evolve(s, u180)
        s = pg.evolve(s, riz)
        s *= 1.0 / float(nsteps)
        res = pg.gen_op(s) if res is None else res + s
    return res


def press_shaped_sigma(system, te, h_op, pulse):
    """PRESS with a shaped refocusing pulse replacing both ideal 180s
    (after Vespa's 'PRESS with real 180 pulses'): each pulse is centred where
    the ideal 180 would be, so the free-evolution delays shrink by its
    length, and is flanked by emulated crusher gradients. Symmetric PRESS
    (TE1 = TE2 = TE/2), on-resonance single voxel position - no spatial
    grid, so the slice profile is not simulated, only the pulse's spectral
    response."""
    u180, tp = shaped_pulse_propagator(system, pulse)
    te1 = te2 = te / 2.0
    if tp >= min(te1, te2):
        raise ValueError("refocusing pulse (%.2f ms) does not fit into TE/2"
                         " (%.2f ms)" % (tp * 1e3, te1 * 1e3))
    u1 = pg.prop(h_op, 0.5 * (te1 - tp))
    u2 = pg.prop(h_op, 0.5 * (te1 - tp + te2 - tp))
    u3 = pg.prop(h_op, 0.5 * (te2 - tp))
    sigma = pg.Iypuls(system, pg.sigma_eq(system), 90.0)
    sigma = pg.evolve(sigma, u1)
    sigma = _crushed_refocus(system, sigma, u180)
    sigma = pg.evolve(sigma, u2)
    sigma = _crushed_refocus(system, sigma, u180)
    return pg.evolve(sigma, u3)


def sequence_sigma(system, sequence, te_ms, h_op, tm_ms=None, pulse=None):
    """Return the density matrix at acquisition start."""
    te = float(te_ms) / 1000.0
    if sequence in ('PRESS shaped',):
        if not pulse:
            raise ValueError("PRESS shaped needs a refocusing pulse waveform")
        return press_shaped_sigma(system, te, h_op, pulse)
    if sequence in ('STEAM',):
        if tm_ms is None:
            raise ValueError("STEAM needs a mixing time (tm_ms)")
        return steam_sigma(system, te, float(tm_ms) / 1000.0, h_op)
    sigma = pg.sigma_eq(system)
    if sequence in ('PRESS',):
        # 90y - TE/4 - 180x - TE/2 - 180x - TE/4 (symmetric PRESS)
        sigma = pg.Iypuls(system, sigma, 90.0)
        u1 = pg.prop(h_op, te / 4.0)
        u2 = pg.prop(h_op, te / 2.0)
        sigma = pg.evolve(sigma, u1)
        sigma = pg.Ixpuls(system, sigma, 180.0)
        sigma = pg.evolve(sigma, u2)
        sigma = pg.Ixpuls(system, sigma, 180.0)
        sigma = pg.evolve(sigma, u1)
        return sigma
    if sequence in ('Spin Echo',):
        sigma = pg.Iypuls(system, sigma, 90.0)
        u = pg.prop(h_op, te / 2.0)
        sigma = pg.evolve(sigma, u)
        sigma = pg.Ixpuls(system, sigma, 180.0)
        sigma = pg.evolve(sigma, u)
        return sigma
    raise ValueError("unsupported sequence: %r" % (sequence,))


def run_job(job):
    cf = float(job['cf_mhz'])
    centre = float(job.get('centre_ppm', 4.65))
    samples = int(job['samples'])
    bandwidth = float(job['bandwidth'])
    linewidth = float(job.get('linewidth', 1.0))
    sequence = job['sequence']
    te = float(job['te_ms'])
    tm = job.get('tm_ms')
    pulse = job.get('pulse')

    basis = {}
    for name, subsystems in job['metabolites'].items():
        # denmatsim convention: a metabolite is a list of sub-spin-systems
        # (e.g. NAA = acetyl + aspartyl groups), summed with scale factors.
        total = np.zeros(samples, dtype=complex)
        for sub in subsystems:
            system = build_spin_system(sub['shifts_ppm'], sub['j_hz'],
                                       cf, centre)
            h_op = pg.Hcs(system) + pg.HJ(system)
            sigma = sequence_sigma(system, sequence, te, h_op, tm_ms=tm,
                                   pulse=pulse)
            detector = pg.gen_op(pg.Fm(system))
            # acq must outlive mx: the table references memory owned by the
            # acquire1D object (a one-liner here corrupts the table/crashes)
            acq = pg.acquire1D(detector, h_op, 0.001)
            mx = acq.table(sigma)
            ppms, areas, phases = bin_table(mx, system, cf)
            del acq
            fid = synthesize_fid(ppms, areas, phases, cf,
                                 samples, bandwidth, linewidth)
            total += float(sub.get('scale', 1.0)) * fid
        basis[name] = {'re': total.real.tolist(), 'im': total.imag.tolist()}
    return basis


def main():
    job_path, out_path = sys.argv[1], sys.argv[2]
    with open(job_path, 'r') as f:
        job = json.load(f)
    try:
        basis = run_job(job)
        result = {'ok': True, 'basis': basis}
    except Exception as exc:  # noqa: BLE001 - the parent surfaces this
        result = {'ok': False, 'error': '%s: %s' % (type(exc).__name__, exc)}
    with open(out_path, 'w') as f:
        json.dump(result, f)


if __name__ == '__main__':
    main()
