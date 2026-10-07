####################################################################################################
#                                       sequence_dialog.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 07/10/26                                                                                #
#                                                                                                  #
# Purpose: The Sequence panel of the parameter step. One view of the sequence the current engine   #
#          will simulate - a timeline of its RF pulses (their waveforms where known) up to the     #
#          echo - with one row per pulse role whose source is ideal, a generated standard pulse    #
#          or any pulse file, and the whole sequence from a Pulseq .seq or an FSL-MRS / WIN .json  #
#          file. Underneath, which engines run the sequence with ideal pulses, a waveform per      #
#          role, or a whole sequence file. The rows edit the same sheet values as the parameter    #
#          list (Path to Pulse, RefTp, Edit Pulse Path, Edit Tp).                                  #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os

import json

import numpy as np
from matplotlib.patches import Rectangle
from nicegui import ui

from basisremy.core import sequence_view as sv
from basisremy.core.pulse_library import STANDARD, bandwidth_hz, is_standard
from basisremy.gui.local_file_picker import LocalFilePicker
from basisremy.gui.ui_state import get_state, set_state

_ACCENT = "var(--br-primary)"
_AXIS = "#8a95a3"
_ROLE_COLOUR = {'exc': '#15627f', 'ref': '#5b7083', 'edit': '#c2892e'}
_ROLE_LABEL = {'exc': '90°', 'ref': '180°', 'edit': 'edit'}
_ROLE_NAME = {'exc': 'Excitation', 'ref': 'Refocusing', 'edit': 'Editing'}
_TP_KEY = {'exc': 'RefTp', 'ref': 'RefTp', 'edit': 'Edit Tp'}   # duration of a role's waveform
_IDEAL = "ideal"
_FILE = "__file__"
# BasisREMY category -> engine name in sequence_view.RUNS_ON
_ENGINE = {'Spant': 'spant', 'Custom': 'jbss'}
_LEVELS = [('ideal', 'Ideal pulses'), ('shaped', 'Pulse per role'), ('file', 'Whole sequence file')]


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
            subtitle.text = f"{kind} · {br.backend.display_name} · {source}"
            draw_timeline(ax, tl)
            try:
                plot.figure.tight_layout(pad=0.3)
            except Exception:                           # noqa: BLE001
                pass
            plot.update()

        def changed(rebuild: bool = False) -> None:
            if rebuild:          # engine, sequence or mode changed: the sheet below changes too
                app._build_tab2()
                rows.refresh()
            redraw()
            runs_on.refresh()

        @ui.refreshable
        def rows() -> None:
            _pulse_rows(app, changed)

        @ui.refreshable
        def runs_on() -> None:
            _runs_on(app)

        rows()
        ui.element("div").classes("br-hairline")
        _file_row(app, dialog)
        runs_on()
        with ui.row().classes("w-full justify-end"):
            ui.button("Close", on_click=dialog.close).props("flat color=primary")
    redraw()
    dialog.open()


def _pulse_rows(app, changed) -> None:
    br = app.BasisREMY
    b = br.backend
    tl = sv.timeline(b, app.seq_file)
    if tl['exact']:
        ui.label("Every pulse comes from the sequence file; remove the file below to set them "
                 "one by one.").classes("text-sm")
        return
    if b.name == 'MRSCloud':
        ui.label("MRSCloud uses its own pulse set: the vendor or universal waveforms chosen by "
                 "System, Localization and Mode in the sheet.").classes("text-sm")
        return
    roles = sv.roles(sv.sequence_kind(b))
    if not roles:
        ui.label("Pick the sequence in the sheet first.").classes("text-sm br-muted")
        return
    for role in roles:
        _pulse_row(app, role, changed)


def _pulse_row(app, role, changed) -> None:
    br = app.BasisREMY
    b = br.backend
    key = sv.pulse_key(b, role)
    tp_key = _TP_KEY[role]
    current = b.mandatory_params.get(key) if key else _IDEAL
    options = {}
    if not key or sv.switch_target(b, role, True):
        options[_IDEAL] = "Ideal (instantaneous)"
    if key or sv.switch_target(b, role, False):
        options.update({f"standard:{n}": f"Standard: {n}" for n in STANDARD})
        if current not in options and current != _IDEAL:
            options[current] = f"File: {os.path.basename(str(current))}"
        options[_FILE] = "Pulse file…"
    with ui.row().classes("w-full items-center no-wrap gap-4"):
        with ui.column().classes("gap-1 grow min-w-0"):
            ui.label(_ROLE_NAME[role]).classes("text-sm font-semibold")
            if len(options) <= 1:
                if role == 'edit' and b.name in ('FSL-MRS', 'Spant'):
                    ui.label(f"{b.display_name}'s own Gaussian editing pulse: Edit On / Off and its "
                             "duration or bandwidth in the sheet").classes("text-xs br-muted")
                elif role == 'exc':
                    ui.label("Ideal (instantaneous). A shaped excitation comes from a whole sequence "
                             "file (FSL-MRS), MRSCloud's pulse set or FID-A's STEAM shaped."
                             ).classes("text-xs br-muted")
                else:
                    ui.label(f"Ideal: {b.display_name} has no waveform for this pulse "
                             f"({_shaped_engines(b, role)})").classes("text-xs br-muted")
                return
            sel = ui.select(options, value=current if current in options else None).props(
                "filled dense").classes("w-full")
            info = ui.label().classes("text-xs br-muted")
            tp = None
            if key:
                tp = ui.number(f"Duration [ms] ({tp_key})", value=sv._num(b.mandatory_params.get(tp_key)),
                               format="%.4g").props("filled dense").classes("w-48")
        mini = ui.matplotlib(figsize=(2.6, 1.1)).classes("w-56 shrink-0")
        mini.figure.patch.set_alpha(0.0)
        mini_ax = mini.figure.add_subplot(111)

    def preview() -> None:
        spec = b.mandatory_params.get(key) if key else None
        p = _pulse(spec, role) if spec else None
        tpv = sv._num(b.mandatory_params.get(tp_key))
        mini.set_visibility(p is not None)
        if p is None:
            _style(mini_ax)
            info.text = "instantaneous rotation" if not spec else "cannot read this pulse"
        else:
            draw_pulse(mini_ax, p, tpv)
            try:
                bw = bandwidth_hz(p, tpv) if tpv else None
            except Exception:                           # noqa: BLE001
                bw = None
            info.text = f"{p.n} points" + (f" · bandwidth {bw:.0f} Hz" if bw else "")
        try:
            mini.figure.tight_layout(pad=0.2)
        except Exception:                               # noqa: BLE001
            pass
        mini.update()

    async def pick(e) -> None:
        value = e.value
        if value == current or value is None:
            return
        if value == _FILE:
            value = await _pick_pulse_file(key or 'Path to Pulse')
            if not value:
                sel.value = current
                return
            if os.path.splitext(value)[1].lower() in ('.seq', '.json'):
                value = f"{value}#{role}"     # whole-sequence file: its pulse of this role
        if value == _IDEAL:
            _apply(app, sv.switch_target(b, role, True))
            changed(rebuild=True)
            return
        target = None if key else sv.switch_target(b, role, False)
        if target:
            _apply(app, target)
        new_key = sv.pulse_key(br.backend, role)
        app._update_param(new_key, value)
        p = _pulse(value, role)
        if p is not None and np.isfinite(p.tp_ms) and not is_standard(value):
            app._update_param(tp_key, round(p.tp_ms, 4))   # the file knows its duration
        changed(rebuild=True)

    def set_tp(e) -> None:
        if e.value is not None:
            app._update_param(tp_key, float(e.value))
            preview()
            changed()

    sel.on_value_change(pick)
    if tp is not None:
        tp.on_value_change(set_tp)
    preview()


def _apply(app, target) -> None:
    """Switch engine / sequence / mode as switch_target() says (sheet values carry over)."""
    br = app.BasisREMY
    if 'backend' in target:
        br.set_backend(target['backend'])
    if 'mode' in target:
        br.backend.set_mode(target['mode'])
    if 'Sequence' in target:
        br.backend.mandatory_params['Sequence'] = target['Sequence']


def _shaped_engines(b, role) -> str:
    kind = sv.sequence_kind(b)
    engines = ", ".join(sv.RUNS_ON.get(kind, {}).get('shaped', [])) or "none yet"
    return f"a {_ROLE_NAME[role].lower()} waveform for {kind}: {engines}"


async def _pick_pulse_file(key) -> str | None:
    if LocalFilePicker.active() is not None:
        return None
    state_key = f"last_dir_{key.lower().replace(' ', '_')}"
    start = get_state(state_key) or "~"
    if not isinstance(start, str) or (start != "~" and not os.path.isdir(start)):
        start = "~"
    path = await LocalFilePicker(start, title="Select pulse file (.pta, .RF, .txt, .seq, .json)")
    if path:
        set_state(state_key, os.path.dirname(path))
    return path


def _runs_on(app) -> None:
    br = app.BasisREMY
    kind = sv.timeline(br.backend, app.seq_file)['kind'] or sv.sequence_kind(br.backend)
    current = _ENGINE.get(br.get_current_category(), br.get_current_category())
    levels = sv.RUNS_ON.get(kind)
    with ui.column().classes("w-full gap-1"):
        ui.label(f"Runs on{f' ({kind})' if kind else ''}").classes("br-section-title")
        if not levels:
            ui.label("Pick the sequence in the sheet to see which engines run it.").classes(
                "text-xs br-muted")
            return
        for level, label in _LEVELS:
            with ui.row().classes("items-center gap-1 w-full"):
                ui.label(label).classes("text-xs br-muted w-40 shrink-0")
                if not levels[level]:
                    ui.label("none yet").classes("text-xs br-muted")
                for engine in levels[level]:
                    on = engine.split(' ')[0] == current
                    ui.badge(engine).props(f"{'' if on else 'outline'} color={'primary' if on else 'grey-7'}")


def _file_row(app, dialog) -> None:
    br = app.BasisREMY
    with ui.row().classes("w-full items-center gap-2"):
        ui.label("Whole sequence from a file").classes("text-sm font-semibold")
        loaded = None
        if app.seq_file:
            loaded = f"Pulseq: {os.path.basename(app.seq_file)}"
        elif br.backend.current_mode == 'Custom' and br.backend.optional_params.get('Custom Sequence'):
            loaded = f"FSL-MRS: {os.path.basename(br.backend.optional_params['Custom Sequence'])}"
        if loaded:
            ui.badge(loaded).props("color=primary")

    async def pick(title, suffix, state_key):
        if LocalFilePicker.active() is not None:
            return None
        start = get_state(state_key) or get_state("last_import_dir") or "~"
        if not isinstance(start, str) or (start != "~" and not os.path.isdir(start)):
            start = "~"
        path = await LocalFilePicker(start, title=title, show_file=lambda p: p.suffix.lower() == suffix)
        if path:
            set_state(state_key, os.path.dirname(path))
        return path

    async def use_pulseq() -> None:
        path = await pick("Select Pulseq sequence file", ".seq", "last_dir_seq")
        if not path:
            return
        app.seq_file = path
        if app._apply_sequence_file():
            dialog.close()
            app._build_tab2()
        else:
            app.seq_file = None

    async def use_fslmrs() -> None:
        path = await pick("Select FSL-MRS / WIN sequence (.json)", ".json", "last_dir_fsl_sequence")
        if not path:
            return
        try:
            seq = sv._fsl_sequence(path)
        except Exception as exc:                        # noqa: BLE001
            ui.notify(str(exc), type="negative")
            return
        app.seq_file = None
        br.set_backend('FSL-MRS')
        br.backend.set_mode('Custom')
        with open(path) as fh:
            nested = 'seq' in json.load(fh)
        if nested:   # a basis JSON: FSL-MRS reads the bare description
            path = os.path.join(br.backend.ensure_workdir(), os.path.basename(path))
            with open(path, 'w') as fh:
                json.dump(seq, fh)
        br.backend.optional_params['Custom Sequence'] = path
        ui.notify(f"FSL-MRS runs {os.path.basename(path)} as given (Custom mode).", type="positive")
        dialog.close()
        app._build_tab2()

    def clear() -> None:
        if app.seq_file:
            app.seq_file = None
        elif br.backend.current_mode == 'Custom':
            br.backend.optional_params['Custom Sequence'] = None
            br.backend.set_mode('Simple')
        dialog.close()
        app._build_tab2()

    with ui.row().classes("gap-2 items-center"):
        ui.button("Pulseq (.seq)", icon="timeline", on_click=use_pulseq).props("outline color=primary dense")
        ui.button("FSL-MRS / WIN (.json)", icon="data_object", on_click=use_fslmrs).props(
            "outline color=primary dense")
        if loaded:
            ui.button("Remove file", icon="close", on_click=clear).props("flat dense color=primary")
    ui.label("A .seq file runs in FID-A's shaped backend for its sequence; an FSL-MRS description "
             "(or a basis JSON, which carries its sequence) runs in FSL-MRS Custom mode."
             ).classes("text-xs br-muted")
