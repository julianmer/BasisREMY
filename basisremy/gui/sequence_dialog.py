####################################################################################################
#                                       sequence_dialog.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 07/10/26                                                                                #
#                                                                                                  #
# Purpose: The Sequence panel behind the wand of the parameter step: a timeline of the RF pulses   #
#          the current engine will simulate (their waveforms where known) up to the echo, one      #
#          preview per pulse role (source, duration, bandwidth, waveform) and which engines run    #
#          the sequence with ideal pulses, a waveform per pulse, their own pulse set or a whole    #
#          sequence file. The pulses themselves are chosen in the parameter sheet.                 #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os

import numpy as np
from matplotlib.patches import Rectangle
from nicegui import ui

from basisremy.core import sequence_setup as ss
from basisremy.core import sequence_view as sv
from basisremy.core.pulse_library import bandwidth_hz, is_standard

_ACCENT = "var(--br-primary)"
_AXIS = "#8a95a3"
_ROLE_COLOUR = {'exc': '#15627f', 'ref': '#5b7083', 'edit': '#c2892e'}
_ROLE_LABEL = {'exc': '90°', 'ref': '180°', 'edit': 'edit'}
_ROLE_NAME = {'exc': 'Excitation', 'ref': 'Refocusing', 'edit': 'Editing'}
_TP_KEY = {'exc': 'RefTp', 'ref': 'RefTp', 'edit': 'Edit Tp'}   # duration of a role's waveform
_LEVELS = [('ideal', 'Ideal pulses'), ('shaped', 'Waveform per pulse'), ('own', 'Own pulse set'),
           ('file', 'Whole sequence file')]


def _pulse(spec, kind):
    """A Pulse for a timeline event / sheet value, or None when it cannot be read."""
    try:
        if isinstance(spec, tuple):
            return sv.read_pulse(spec[0], index=spec[1])
        if isinstance(spec, dict):
            return sv.block_pulse(spec, kind)
        if spec and (is_standard(spec) or os.path.exists(str(spec).partition('#')[0])):
            return sv.read_pulse(spec, kind)
    except Exception:                                   # noqa: BLE001 - shown as ideal
        pass
    return None


def _style(ax):
    ax.clear()
    ax.set_facecolor("none")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(_AXIS)
    ax.tick_params(colors=_AXIS, labelsize=8)
    ax.set_yticks([])


def draw_timeline(ax, tl) -> None:
    """RF pulses as their amplitude envelopes (thin bars when ideal), the echo and the readout."""
    _style(ax)
    events, echo = tl['events'], tl['echo_ms']
    if not events or echo is None:
        ax.text(0.5, 0.5, "Set the sequence and TE to see its timeline", ha="center", va="center",
                transform=ax.transAxes, color=_AXIS, fontsize=9)
        return
    span = max(echo, max(e['centre_ms'] + e['dur_ms'] / 2 for e in events)) * 1.25
    height = {'exc': 0.65, 'ref': 1.0, 'edit': 0.5}
    for e in events:
        colour, h = _ROLE_COLOUR[e['role']], height[e['role']]
        c, d = e['centre_ms'], e['dur_ms']
        p = _pulse(e['pulse'], e['role']) if d > 0.05 else None
        if p is not None:
            amp = p.waveform[:, 1] / np.max(p.waveform[:, 1])
            t = c - d / 2 + d * (np.arange(p.n) + 0.5) / p.n
            ax.fill_between(t, 0, h * amp, color=colour, alpha=0.85, lw=0)
        else:
            w = max(d, span * 0.006)
            ax.add_patch(Rectangle((c - w / 2, 0), w, h, color=colour, alpha=0.85, lw=0))
        ax.text(c, h + 0.06, _ROLE_LABEL[e['role']], ha="center", va="bottom", fontsize=7,
                color=colour)
    # readout: a decaying echo from the echo time on
    t = np.linspace(echo, span, 300)
    ax.plot(t, 0.35 * np.exp(-(t - echo) / (0.12 * span)) * np.cos(2 * np.pi * (t - echo) / (0.03 * span)),
            color=_AXIS, lw=0.8)
    ax.axvline(echo, color=_AXIS, lw=0.8, ls="--")
    ax.text(echo + span * 0.01, 1.18, f"echo, TE {echo:.4g} ms", ha="left", fontsize=7, color=_AXIS)
    ax.axhline(0, color=_AXIS, lw=0.6)
    ax.set_xlim(-span * 0.04, span)
    ax.set_ylim(-0.45, 1.35)
    ax.set_xlabel("Time from excitation [ms]", color=_AXIS, fontsize=8)


def draw_pulse(ax, pulse, tp_ms) -> None:
    """Amplitude and phase of one waveform over its duration."""
    _style(ax)
    t = np.linspace(0, tp_ms if tp_ms and np.isfinite(tp_ms) else 1.0, pulse.n)
    amp = pulse.waveform[:, 1] / np.max(pulse.waveform[:, 1])
    ax.fill_between(t, 0, amp, color=_ROLE_COLOUR['ref'], alpha=0.8, lw=0)
    ph = np.unwrap(np.deg2rad(pulse.waveform[:, 0]))
    if np.ptp(ph) > 1e-3:
        ax.plot(t, (ph - ph.min()) / np.ptp(ph), color=_ROLE_COLOUR['edit'], lw=0.9)
    unit = "ms" if tp_ms and np.isfinite(tp_ms) else "samples (duration not set)"
    ax.set_xlabel(f"{unit} · grey amplitude, orange phase", color=_AXIS, fontsize=7)


def open_sequence_dialog(app) -> None:
    """Open the Sequence panel for the app's current engine and sheet."""
    br = app.BasisREMY
    # at the page root: switching engine rebuilds the parameter panel, which must not take the
    # open dialog with it
    with ui.context.client.layout:
        dialog = ui.dialog()
    dialog.on_value_change(lambda e: None if e.value else dialog.delete())
    with dialog, ui.card().classes("w-[900px] max-w-full gap-3"):
        with ui.row().classes("w-full items-center justify-between no-wrap"):
            ui.label("Sequence").classes("text-lg font-bold").style(f"color:{_ACCENT}")
            subtitle = ui.label().classes("text-xs br-muted")
        plot = ui.matplotlib(figsize=(8.4, 2.3)).classes("w-full")
        plot.figure.patch.set_alpha(0.0)
        ax = plot.figure.add_subplot(111)

        def redraw() -> None:
            tl = sv.timeline(br.backend, app.seq_file)
            kind = tl['kind'] or sv.sequence_kind(br.backend) or "sequence not set"
            source = "exact, from the sequence file" if tl['exact'] else "placed from the sheet"
            engine = ss.ENGINE_LABEL.get(br.backend.category, br.backend.category)
            subtitle.text = f"{kind} · {engine} · {source}"
            draw_timeline(ax, tl)
            try:
                plot.figure.tight_layout(pad=0.3)
            except Exception:                           # noqa: BLE001
                pass
            plot.update()

        _pulse_previews(app)
        ui.element("div").classes("br-hairline")
        _runs_on(app)
        with ui.row().classes("w-full justify-end"):
            ui.button("Close", on_click=dialog.close).props("flat color=primary")
    redraw()
    dialog.open()


def _source(spec) -> str:
    if spec is None:
        return "Ideal (instantaneous)"
    if isinstance(spec, tuple):
        return f"{os.path.basename(spec[0])}, RF {spec[1] + 1}"
    if isinstance(spec, dict):
        return "from the sequence description"
    if is_standard(spec):
        return f"Standard: {str(spec).split(':', 1)[1]}"
    return f"File: {os.path.basename(str(spec).partition('#')[0])}"


def _pulse_previews(app) -> None:
    """One read-only row per pulse role: its source, duration, bandwidth and waveform."""
    br = app.BasisREMY
    b = br.backend
    tl = sv.timeline(b, app.seq_file)
    seq, route = ss.current(br)
    first = {}
    for e in tl['events']:
        first.setdefault(e['role'], e)
    if not first:
        ui.label("Pick the sequence and TE in the sheet to see its pulses.").classes("text-sm br-muted")
        return
    for role in ('exc', 'ref', 'edit'):
        if role not in first:
            continue
        e = first[role]
        own = route is not None and route.pulses.get(role) == ss.OWN and not tl['exact']
        p = None if own else _pulse(e['pulse'], role)
        with ui.row().classes("w-full items-center no-wrap gap-4"):
            with ui.column().classes("gap-0 grow min-w-0"):
                ui.label(ss.ROLE_NAME[role]).classes("text-sm font-semibold")
                text = (f"{ss.ENGINE_LABEL.get(b.category, b.category)}'s own pulse" if own
                        else _source(e['pulse']))
                if p is not None and e['dur_ms']:
                    try:
                        bw = bandwidth_hz(p, e['dur_ms'])
                    except Exception:                       # noqa: BLE001
                        bw = None
                    text += f" · {e['dur_ms']:.4g} ms" + (f" · bandwidth {bw:.0f} Hz" if bw else "")
                ui.label(text).classes("text-xs br-muted")
            if p is not None:
                mini = ui.matplotlib(figsize=(2.6, 1.1)).classes("w-56 shrink-0")
                mini.figure.patch.set_alpha(0.0)
                draw_pulse(mini.figure.add_subplot(111), p, e['dur_ms'])
                try:
                    mini.figure.tight_layout(pad=0.2)
                except Exception:                           # noqa: BLE001
                    pass
                mini.update()


def _runs_on(app) -> None:
    br = app.BasisREMY
    seq = sv.timeline(br.backend, app.seq_file)['kind'] or ss.current(br)[0]
    current = ss.ENGINE_LABEL.get(br.backend.category, br.backend.category)
    with ui.column().classes("w-full gap-1"):
        ui.label(f"Runs on{f' ({seq})' if seq else ''}").classes("br-section-title")
        if not seq:
            ui.label("Pick the sequence in the sheet to see which engines run it.").classes(
                "text-xs br-muted")
            return
        levels = ss.levels(seq)
        for level, label in _LEVELS:
            with ui.row().classes("items-center gap-1 w-full"):
                ui.label(label).classes("text-xs br-muted w-40 shrink-0")
                if not levels[level]:
                    ui.label("none yet").classes("text-xs br-muted")
                for engine in levels[level]:
                    on = engine.split(' (')[0] == current.split(' (')[0]
                    ui.badge(engine).props(f"{'' if on else 'outline'} color={'primary' if on else 'grey-7'}")
