####################################################################################################
#                                       sequence_dialog.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 07/10/26                                                                                #
#                                                                                                  #
# Purpose: The sequence designer behind the wand of the parameter step, the same for every engine: #
#          load any sequence or pulse file, pick the sequence, set its timings, choose each pulse  #
#          (ideal, a standard shape or a pulse file; recommended ones pre-selected, with their     #
#          source), see the timeline and which engines run the design, and save it as a Pulseq     #
#          .seq. The saved file is selected in the sheet's file field, exactly as if the user had  #
#          brought it (core/sequence_design.py).                                                   #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import re

import numpy as np
from matplotlib.patches import Rectangle
from nicegui import ui

from basisremy.core import sequence_design as sd
from basisremy.core import sequence_setup as ss
from basisremy.core import sequence_view as sv
from basisremy.core.pulse_library import bandwidth_hz, is_standard
from basisremy.gui.local_file_picker import LocalFilePicker
from basisremy.gui.ui_state import get_state, set_state

_ACCENT = "var(--br-primary)"
_AXIS = "#8a95a3"
_ROLE_COLOUR = {'exc': '#15627f', 'ref': '#5b7083', 'edit': '#c2892e'}
_ROLE_LABEL = {'exc': '90°', 'ref': '180°', 'edit': 'edit'}
_IDEAL, _PICK = 'ideal', '__file__'
_WHOLE_EXTS = ('.seq', '.json')


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
    te = tl.get('te', echo)
    label = (f"echo, TE {echo:.4g} ms" if abs(echo - te) < 1e-6 else         # STEAM: the echo after TE + TM
             f"echo at {echo:.4g} ms (TE {te:.4g} + TM {echo - te:.4g})")
    ax.text(echo + span * 0.01, 1.18, label, ha="left", fontsize=7, color=_AXIS)
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


#**************************************************************************************************#
#                                            designer                                              #
#**************************************************************************************************#
def _seed(app) -> sd.Design | None:
    """The design the designer opens with: the selected sequence file, else the recommended one
    for the sheet's sequence and TE."""
    br = app.BasisREMY
    b = br.backend
    sheet = {**b.optional_params, **b.mandatory_params}
    if app.seq_file and getattr(app, '_design', None) is not None:
        return sd.Design(**{**app._design.__dict__, 'pulses': {r: dict(p) for r, p in
                                                              app._design.pulses.items()}})
    seq = ss.current(br)[0]
    kind = seq if seq in sd.DESIGNABLE else 'PRESS'
    te = sd._num(sheet.get('TE'))
    return sd.recommend(kind, te, sheet, br._last_mrsinmrs) if te else None


def open_sequence_dialog(app) -> None:
    """Open the sequence designer for the app's current engine and sheet."""
    br = app.BasisREMY
    b = br.backend
    sheet = {**b.optional_params, **b.mandatory_params}
    header = getattr(br, '_last_mrsinmrs', None)
    state = {'d': _seed(app), 'kind': None, 'name': None}
    if state['d'] is None:
        seq = ss.current(br)[0]
        state['kind'] = seq if seq in sd.DESIGNABLE else 'PRESS'
    # at the page root: a rebuild of the parameter panel must not take the open dialog with it
    with ui.context.client.layout:
        dialog = ui.dialog()
    dialog.on_value_change(lambda e: None if e.value else dialog.delete())

    def soon() -> None:
        with dialog:
            ui.timer(0.01, render, once=True)

    def set_te(value) -> None:
        te = sd._num(value)
        d = state['d']
        if te is None or te <= 0:
            return
        if d is None:
            state['d'] = sd.recommend(state['kind'], te, sheet, header)
        else:
            fresh = sd.recommend(d.kind, te, {})
            d.te = te
            for k in list(d.timing):
                if k in d.rec:
                    d.timing[k] = fresh.timing[k]
            if 'dur:edit' in d.rec and 'edit' in fresh.pulses:
                d.pulses['edit']['dur'] = fresh.pulses['edit']['dur']
                d.rec['dur:edit'] = fresh.rec['dur:edit']
        state['name'] = None
        soon()

    def set_kind(kind) -> None:
        d = state['d']
        if d is None:
            state['kind'] = kind
            soon()
            return
        if kind == d.kind:
            return
        new = sd.recommend(kind, d.te, sheet, header)
        new.voxel = d.voxel
        if 'voxel' not in d.rec:
            new.rec.pop('voxel', None)
        for role, p in d.pulses.items():                 # chosen pulses carry over
            if role in new.pulses and f'pulse:{role}' not in d.rec and not \
                    str(p['source']).partition('#')[0].lower().endswith(_WHOLE_EXTS):
                new.pulses[role] = dict(p)
                new.rec.pop(f'pulse:{role}', None)
                new.rec.pop(f'dur:{role}', None)
        state['d'], state['name'] = new, None
        soon()

    def set_timing(key, value) -> None:
        d = state['d']
        v = sd._num(value)
        if v is None or v == d.timing.get(key):
            return
        d.timing[key] = v
        d.rec.pop(key, None)
        keys = sd.TIMING_KEYS[d.kind]
        if keys == ('TE1', 'TE2'):                        # the other half follows
            other = 'TE2' if key == 'TE1' else 'TE1'
            d.timing[other] = round(d.te - v, 6)
            d.rec.pop(other, None)
        soon()

    def set_source(role, source) -> None:
        d = state['d']
        dur, why = sd.default_duration(source, role)
        d.pulses[role] = {'source': source, 'dur': dur}
        d.rec.pop(f'pulse:{role}', None)
        d.rec.pop(f'dur:{role}', None)
        if why and why != "the pulse's own duration":
            d.rec[f'dur:{role}'] = why
        state['name'] = None
        soon()

    def set_dur(role, value) -> None:
        v = sd._num(value)
        d = state['d']
        if v is None or v == d.pulses[role]['dur']:
            return
        d.pulses[role]['dur'] = v
        d.rec.pop(f'dur:{role}', None)
        soon()

    def set_edit(i, value) -> None:
        v = sd._num(value)
        d = state['d']
        if v is None or v == d.edit[i]:
            return
        e = list(d.edit)
        e[i] = v
        d.edit = tuple(e)
        d.rec.pop('edit', None)
        soon()

    def set_voxel(i, value) -> None:
        v = sd._num(value)
        d = state['d']
        if v is None or v <= 0 or v == d.voxel[i]:
            return
        vox = list(d.voxel)
        vox[i] = v
        d.voxel = tuple(vox)
        d.rec.pop('voxel', None)
        state['voxel_user'] = True
        soon()

    def set_targets(label, value) -> None:
        try:
            targets = tuple(float(x) for x in str(value).replace(';', ',').split(',') if x.strip())
        except ValueError:
            ui.notify("Editing targets: ppm values separated by commas.", type="warning")
            return
        d = state['d']
        if not targets or targets == tuple(d.scheme.get(label, ())):
            return
        d.scheme[label] = targets
        d.rec.pop('edit', None)
        soon()

    async def load() -> None:
        path = await _pick("Load a sequence or pulse file", "last_dir_sequence_file")
        if not path:
            return
        if path.lower().endswith(_WHOLE_EXTS):
            try:
                state['d'] = sd.read_design(path, sheet.get('Bfield'))
            except Exception as exc:                         # noqa: BLE001
                ui.notify(f"{os.path.basename(path)}: {exc}", type="negative", multi_line=True)
                return
            state['name'] = os.path.splitext(os.path.basename(path))[0]
            soon()
            return
        d = state['d']
        if d is None:
            ui.notify("Set TE first, then the pulse.", type="warning")
            return
        role = await _ask_role(sd.roles(d.kind), path)
        if role:
            set_source(role, path)

    async def save() -> None:
        d = state['d']
        name = (name_in.value or sd.default_name(d)).strip()
        folder = sd.designs_dir()
        path = os.path.join(folder, f"{name}.seq")
        n = 2
        while os.path.exists(path):                          # never overwrite an earlier design
            path = os.path.join(folder, f"{name}_{n}.seq")
            n += 1
        try:
            sd.write_seq(d, path)
        except Exception as exc:                             # noqa: BLE001
            ui.notify(f"Could not save the design: {exc}", type="negative", multi_line=True)
            return
        with app.panel2:            # the dialog deletes itself on close: work in the page
            app.seq_file = path
            if app._apply_sequence_file():
                ui.notify(f"Saved {os.path.basename(path)} in {folder}", type="positive")
            else:
                app.seq_file = None
                ui.notify(f"Saved {os.path.basename(path)} in {folder}; pick an engine that "
                          "runs it.", type="warning")
            app._rebuild_soon()
        dialog.close()

    with dialog, ui.card().classes("w-[980px] max-w-full gap-3"):
        with ui.row().classes("w-full items-center justify-between no-wrap"):
            ui.label("Sequence designer").classes("text-lg font-bold").style(f"color:{_ACCENT}")
            ui.button("Load file", icon="folder_open", on_click=load).props(
                "flat dense color=primary").tooltip(
                "A whole sequence (Pulseq .seq, sequence .json) or one pulse "
                "(.pta, .RF, .txt, Bruker, .mat)")
        body = ui.column().classes("w-full gap-3")
        ui.element("div").classes("br-hairline")
        with ui.row().classes("w-full items-center justify-end gap-2 no-wrap"):
            name_in = ui.input(placeholder="name").props("filled dense").classes("w-80")
            ui.button("Cancel", on_click=dialog.close).props("flat color=primary")
            save_btn = ui.button("Save", icon="save", on_click=save).props("unelevated color=primary")
            save_btn.mark("designer-save")

    def render() -> None:
        body.clear()
        d = state['d']
        with body:
            _header_row(d, state, set_kind, set_te, set_voxel)
            plot = ui.matplotlib(figsize=(9.0, 2.2)).classes("w-full")
            plot.figure.patch.set_alpha(0.0)
            ax = plot.figure.add_subplot(111)
            if d is None:
                draw_timeline(ax, {'events': [], 'echo_ms': None})
            else:
                ev = [{'role': e['role'], 'centre_ms': e['centre'], 'dur_ms': e['dur'],
                       'pulse': e['source']} for e in sd.events(d)] if not _broken(d) else []
                draw_timeline(ax, {'events': ev, 'echo_ms': sd.echo(d) if ev else None, 'te': d.te})
            try:
                plot.figure.tight_layout(pad=0.3)
            except Exception:                                # noqa: BLE001
                pass
            plot.update()
            if d is None:
                ui.label("TE is not in the sheet: type it above.").classes("text-sm br-muted")
                save_btn.disable()
                return
            _timings(d, set_timing)
            for role in sd.roles(d.kind):
                _pulse_row(d, role, set_source, set_dur, set_edit, set_targets)
            _engines(app, d)
            errs = sd.problems(d)
            for e in errs:
                ui.label(e).classes("text-sm").style("color:#c2453c")
        if state['name'] is None and d is not None:
            state['name'] = sd.default_name(d)
        name_in.value = state['name'] or ""
        if d is not None and not sd.problems(d):
            save_btn.enable()
        else:
            save_btn.disable()

    name_in.on("blur", lambda: state.__setitem__('name', name_in.value))
    render()
    dialog.open()


def _broken(d) -> bool:
    return any(sd._num(d.timing.get(k)) is None for k in sd.TIMING_KEYS[d.kind])


#**************************************************************************************************#
#                                          designer parts                                          #
#**************************************************************************************************#
def _num_field(label, value, on_set, css, tip, width="w-28"):
    with ui.row().classes("items-center gap-2 no-wrap"):
        ui.label(label).classes("text-sm")
        field = ui.input(value="" if value is None else f"{value:.6g}").props("filled dense").classes(
            f"{width} br-v-{css}")
        with field:
            ui.tooltip(tip)
        field.on("blur", lambda: on_set(field.value))
        field.on("keydown.enter", lambda: on_set(field.value))
    return field


def _text_field(label, value, on_set, css, tip, width="w-28"):
    with ui.row().classes("items-center gap-2 no-wrap"):
        ui.label(label).classes("text-sm")
        field = ui.input(value=value).props("filled dense").classes(f"{width} br-v-{css}")
        with field:
            ui.tooltip(tip)
        field.on("blur", lambda: on_set(field.value))
        field.on("keydown.enter", lambda: on_set(field.value))
    return field


def _css(d, key) -> tuple[str, str]:
    if d is not None and key in d.rec:
        return 'rec', f"Recommended: {d.rec[key]}"
    return 'user', "Set by you or from the loaded file"


def _header_row(d, state, set_kind, set_te, set_voxel) -> None:
    with ui.row().classes("w-full items-center gap-6"):
        with ui.row().classes("items-center gap-2 no-wrap"):
            ui.label("Sequence").classes("text-sm font-semibold")
            ui.select(sd.DESIGNABLE, value=d.kind if d else state['kind'],
                      on_change=lambda e: set_kind(e.value)).props("filled dense").classes("w-44")
        te = d.te if d else None
        _num_field("TE [ms]", te, set_te, 'file' if te else 'missing',
                   "The scan's echo time (from the data file)")
        if d is not None and any(p['source'] != 'ideal' for r, p in d.pulses.items() if r != 'edit'):
            css, tip = _css(d, 'voxel')
            if state.get('voxel_user'):
                css, tip = 'user', "Set by you"
            elif 'voxel' not in d.rec:
                css, tip = 'file', "From the data or sequence file"
            ui.label("Voxel [cm]").classes("text-sm")
            for i, axis in enumerate(("L-R", "A-P", "C-C")):
                _num_field(axis, d.voxel[i], lambda v, i=i: set_voxel(i, v), css,
                           tip + " · selected by the excitation, the first and the second "
                           "refocusing pulse", width="w-16")


def _timings(d, set_timing) -> None:
    keys = sd.TIMING_KEYS[d.kind]
    if not keys:
        return
    with ui.row().classes("w-full items-center gap-6"):
        ui.label("Timings").classes("br-section-title w-24")
        for k in keys:
            css, tip = _css(d, k)
            _num_field(f"{sd.TIMING_LABEL[k]} [ms]", d.timing.get(k), lambda v, k=k: set_timing(k, v),
                       css, tip)


def _source_label(src) -> str:
    if src == _IDEAL:
        return "Ideal (instantaneous)"
    if is_standard(src):
        return f"Standard: {src.split(':', 1)[1]}"
    path, _, sel = str(src).partition('#')
    header = re.fullmatch(r'header_(.+)_[0-9a-f]{8}\.exc', os.path.basename(path))   # sd._shape_file
    if header:
        return f"From header: {header.group(1)}"
    return f"From {os.path.basename(path)}" + (f" ({sel})" if sel and not sel.isdigit() else "")


def _pulse_row(d, role, set_source, set_dur, set_edit, set_targets) -> None:
    p = d.pulses[role]
    src = p['source']
    options = {} if role == 'edit' else {_IDEAL: _source_label(_IDEAL)}
    options.update({f"standard:{n}": f"Standard: {n}" for n in sd.STANDARD_FOR[role]})
    if src not in options:
        options[src] = _source_label(src)
    options[_PICK] = "Pulse file…"
    css, tip = _css(d, f'pulse:{role}')
    with ui.row().classes("w-full items-center no-wrap gap-4"):
        with ui.column().classes("gap-1 grow min-w-0"):
            with ui.row().classes("items-center gap-3 no-wrap"):
                ui.label(ss.ROLE_NAME[role]).classes("text-sm font-semibold w-24")
                sel = ui.select(options, value=src).props("filled dense").classes(f"w-72 br-v-{css}")
                with sel:
                    ui.tooltip(tip)

                async def pick(e, current=src) -> None:
                    if e.value in (None, current):
                        return
                    new = e.value
                    if new == _PICK:
                        new = await _pick(f"Select the {ss.ROLE_NAME[role].lower()} pulse",
                                          f"last_dir_pulse_{role}")
                        if not new:
                            e.sender.value = current
                            return
                    set_source(role, new)
                sel.on_value_change(pick)
                if src != _IDEAL:
                    dcss, dtip = _css(d, f'dur:{role}')
                    _num_field("Duration [ms]", p['dur'], lambda v: set_dur(role, v), dcss, dtip,
                               width="w-20")
            if role == 'edit' and d.scheme:
                with ui.row().classes("items-center gap-4 pl-28"):
                    ecss, etip = _css(d, 'edit')
                    for label, targets in d.scheme.items():
                        _text_field(f"{label} [ppm]", ", ".join(f"{p:g}" for p in targets),
                                    lambda v, k=label: set_targets(k, v), ecss,
                                    etip + " · two targets = dual-lobe pulse")
            elif role == 'edit':
                with ui.row().classes("items-center gap-4 pl-28"):
                    ecss, etip = _css(d, 'edit')
                    _num_field("ON [ppm]", d.edit[0], lambda v: set_edit(0, v), ecss, etip, "w-20")
                    _num_field("OFF [ppm]", d.edit[1], lambda v: set_edit(1, v), ecss, etip, "w-20")
            pulse = _pulse(src, role) if src != _IDEAL else None
            info = ""
            if pulse is not None and p['dur']:
                try:
                    info = f"bandwidth {bandwidth_hz(pulse, p['dur']):.0f} Hz"
                except Exception:                            # noqa: BLE001
                    info = ""
            elif src != _IDEAL:
                info = "could not read this pulse"
            if info:
                ui.label(info).classes("text-xs br-muted pl-28")
        if pulse is not None:
            mini = ui.matplotlib(figsize=(2.6, 1.0)).classes("w-56 shrink-0")
            mini.figure.patch.set_alpha(0.0)
            draw_pulse(mini.figure.add_subplot(111), pulse, p['dur'])
            try:
                mini.figure.tight_layout(pad=0.2)
            except Exception:                                # noqa: BLE001
                pass
            mini.update()


def _engines(app, d) -> None:
    br = app.BasisREMY
    with ui.column().classes("w-full gap-1"):
        ui.label("Engines").classes("br-section-title")
        with ui.row().classes("items-center gap-1 w-full"):
            for cat in ss.engines():
                if not br.categories.get(cat):
                    continue
                pl = sd.plan(d, cat, br)
                name = ss.ENGINE_LABEL.get(cat, cat)
                if pl.status == 'ok':
                    badge = ui.badge(name).props("color=primary")
                elif pl.status == 'approx':
                    badge = ui.badge(f"{name} ≈").props("outline color=orange-9")
                else:
                    badge = ui.badge(name).props("outline color=grey-5")
                with badge:
                    ui.tooltip(" ".join(pl.notes) or "Runs this design as it is.").classes("max-w-xs")
                badge.mark(f"engine:{cat}:{pl.status}")
        ui.label("Filled: runs as designed · ≈ runs with differences (hover) · grey: cannot run it "
                 "(hover for why)").classes("text-xs br-muted")


#**************************************************************************************************#
#                                            pickers                                               #
#**************************************************************************************************#
async def _pick(title, state_key) -> str | None:
    if LocalFilePicker.active() is not None:
        return None
    start = get_state(state_key) or get_state("last_dir_sequence_file") or "~"
    if not isinstance(start, str) or (start != "~" and not os.path.isdir(start)):
        start = "~"
    path = await LocalFilePicker(start, title=title)
    if path:
        set_state(state_key, os.path.dirname(path))
    return path


async def _ask_role(roles, path) -> str | None:
    if len(roles) == 1:
        return roles[0]
    with ui.dialog() as dlg, ui.card().classes("gap-3"):
        ui.label(f"Use {os.path.basename(path)} as").classes("text-sm font-semibold")
        with ui.row().classes("gap-2"):
            for r in roles:
                ui.button(ss.ROLE_NAME[r], on_click=lambda _, r=r: dlg.submit(r)).props(
                    "outline color=primary")
    return await dlg
