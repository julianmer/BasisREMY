####################################################################################################
#                                       parameter_sheet.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 07/10/26                                                                                #
#                                                                                                  #
# Purpose: The parameter step, laid out the same for every engine: Engine, Sequence and a          #
#          sequence or pulse file (any format) with the wand opening the Sequence panel; one row   #
#          per pulse role (ideal, a standard pulse or a file); Timings, Acquisition and the        #
#          engine's own settings; the metabolites. Every value shows where it comes from - the     #
#          data file, a recommendation (with its source), the user - or that it is missing.        #
#          What an engine cannot run is greyed out with the reason (core/sequence_setup.py).       #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import json
import os

from nicegui import run, ui

from basisremy.core import sequence_setup as ss
from basisremy.core import sequence_view as sv
from basisremy.core.parameter_registry import get as registry_get
from basisremy.gui.help_widget import label_with_help
from basisremy.gui.local_file_picker import LocalFilePicker
from basisremy.gui.ui_state import get_state, set_state

_UNSET = (None, "", "missing input", "Select option")
_IDEAL, _FILE = "ideal", "__file__"
_PULSE_EXTS = {'.pta', '.rf', '.txt', '.exc', '.rfc', '.inv', '.mat'}

# one name per value, whatever the engine calls it
_LABEL = {
    'TE': 'TE', 'Tau 1': 'TE1 (first echo)', 'Tau 2': 'TE2 (second echo)', 'TE2': 'TE2 (second echo)',
    'sLASER TE1': 'TE1 (first pair)', 'sLASER TE2': 'TE2 (between pairs)', 'sLASER TE3': 'TE3 (second pair)',
    'TM': 'TM (mixing time)', 'Tau': 'Echo spacing', 'Nechoes': 'Echoes', 'Delay': 'ADC delay',
    'Bfield': 'Field strength', 'Center Freq': 'Centre frequency', 'Samples': 'Points',
    'Bandwidth': 'Spectral width', 'Linewidth': 'Linewidth', 'Nucleus': 'Nucleus',
    'System': 'Scanner vendor', 'Edit On': 'ON', 'Edit Off': 'OFF', 'Edit Tp': 'Duration',
    'Edit Bandwidth (Hz)': 'Bandwidth', 'RefTp': 'Duration', 'Flip Angle': 'Flip angle (waveform)',
    'Spatial Points': 'Spatial points', 'thkX': 'Slab thickness x', 'thkY': 'Slab thickness y',
    'fovX': 'Simulated extent x', 'fovY': 'Simulated extent y', 'nX': 'Grid points x', 'nY': 'Grid points y',
    'Sim Centre (ppm)': 'Simulation centre', 'STEAM Variant': 'STEAM coherence selection',
    'Vendor Pulse File': 'Vendor pulse file', 'Pulse Phase': 'Pulse phase',
}
_UNITS = {'Edit On': 'ppm', 'Edit Off': 'ppm', 'Edit Tp': 'ms', 'RefTp': 'ms', 'Edit Bandwidth (Hz)': 'Hz'}
_TIMING = ['TE', 'Tau 1', 'Tau 2', 'TE2', 'sLASER TE1', 'sLASER TE2', 'sLASER TE3', 'TM', 'Tau',
           'Nechoes', 'Delay']
_ACQ = ['Bfield', 'Center Freq', 'Samples', 'Bandwidth', 'Nucleus', 'Linewidth', 'System']
_EDIT_FIELDS = ['Edit On', 'Edit Off', 'Edit Tp', 'Edit Bandwidth (Hz)']
# handled by the selectors, the pulse rows or the metabolite list
_ELSEWHERE = {'Sequence', 'Localization', 'Metabolites', 'Custom Sequence', 'Template File',
              'Path to Pulse', 'Edit Pulse Path', 'RefTp', 'Vendor Pulse File', *_EDIT_FIELDS}
_TP_KEY = {'exc': 'RefTp', 'ref': 'RefTp', 'edit': 'Edit Tp'}
# standard pulses offered per role
_STANDARD_FOR = {'exc': ['sinc-exc'],
                 'ref': ['sinc-ref', 'hs4-ref', 'hs1-inv', 'goia-wurst', 'goia-hs', 'foci'],
                 'edit': ['gauss-edit', 'gauss-edit-20ms']}
# engine options: the backend modes that keep the same pulses
_MODE_LABEL = {
    'MRSCloud': ('Pulse set', {'Universal': 'Universal (open, bundled)',
                               'Non-Universal': "Vendor's product pulses"}),
    'FidaSemiLaserShaped': ('Coherence selection', {'Standard': 'Fast (Zhang 2017)',
                                                    'Phase cycled': '4-step phase cycle'}),
    'FidaOnePulse': ('Variant', {'Ideal': 'Plain', 'Delay': 'ADC delay',
                                 'Arbitrary phase': 'Arbitrary phase'}),
}
_STATE_HELP = {'file': "From the data or sequence file", 'user': "Set by you",
               'missing': "Missing: needed before simulating"}


class _Select(ui.select):
    """A select whose options in ``disabled`` are greyed out (Quasar's option ``disable``)."""

    def __init__(self, options, *, disabled=(), **kwargs):
        self._disabled = set(disabled)
        super().__init__(options, **kwargs)

    def _update_options(self) -> None:
        super()._update_options()
        for option, value in zip(self._props['options'], self._values):
            if value in self._disabled:
                option['disable'] = True


def _display(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6f}".rstrip("0").rstrip(".")
    return str(value)


def _label(key) -> str:
    text = _LABEL.get(key, key)
    units = _UNITS.get(key) or registry_get(key).units
    if units and units not in ('points',) and f'[{units}]' not in text:
        text = f"{text} [{units}]"
    return text


#**************************************************************************************************#
#                                          value states                                            #
#**************************************************************************************************#
def _state(app, key, value) -> tuple[str, str]:
    """(css state, tooltip) of one sheet value."""
    br = app.BasisREMY
    b = br.backend
    st = ss.value_state(b, key, value, app._user_set, br.from_file.get(b.name, {}), app._rec)
    if st == 'recommended':
        return 'rec', f"Recommended: {app._rec[key]}"
    if st == 'default':
        return 'rec', f"{ss.ENGINE_LABEL.get(b.category, b.category)} default"
    return st, _STATE_HELP[st]


def _style(app, el, key, value) -> None:
    css, tip = _state(app, key, value)
    el.classes(remove="br-v-file br-v-rec br-v-user br-v-missing", add=f"br-v-{css}")
    if key in app._tooltips:
        app._tooltips[key].set_text(tip)


def refresh_states(app) -> None:
    """Re-colour every field (after an edit: recommendations follow TE)."""
    b = app.BasisREMY.backend
    seq, _ = ss.current(app.BasisREMY)
    if not _whole_file(app):
        app._rec = ss.apply_recommended(b, seq, _keep(app))
    app._programmatic = True
    try:
        for key, el in app._fields.items():
            value = b.mandatory_params.get(key, b.optional_params.get(key))
            if key in app._rec and getattr(el, 'value', None) != _display(value):
                el.value = _display(value)
            _style(app, el, key, value)
    finally:
        app._programmatic = False


def _keep(app) -> set:
    """Keys a recommendation must not overwrite: set by the user or read from the file."""
    b = app.BasisREMY.backend
    return set(app._user_set) | set(app.BasisREMY.from_file.get(b.name, {}))


#**************************************************************************************************#
#                                             fields                                               #
#**************************************************************************************************#
def _field_row(key, sub):
    row = ui.element("div").classes("br-prow br-prow-sub" if sub else "br-prow")
    with row:
        label_with_help(key, _label(key)).classes("br-prow-label")
    return row


def _register(app, key, el, value) -> None:
    with el:
        app._tooltips[key] = ui.tooltip("")
    app._fields[key] = el
    el.mark(f"param:{key.replace(' ', '_')}")
    _style(app, el, key, value)


def text_field(app, key, value, sub=False) -> None:
    with _field_row(key, sub):
        inp = ui.input(value=_display(value)).props("filled dense").classes("br-pfield")
        inp.on_value_change(lambda e, k=key: app._update_param(k, e.value))
    _register(app, key, inp, value)


def dropdown_field(app, key, value) -> None:
    options = app.BasisREMY.backend.dropdown[key]
    initial = str(value) if value is not None and str(value) in [str(k) for k in options] else None
    with _field_row(key, False):
        sel = ui.select(options, value=initial).props("filled dense").classes("br-pfield")
        sel.on_value_change(lambda e, k=key: app._update_param(k, e.value))
    _register(app, key, sel, value)


def file_field(app, key, value) -> None:
    with _field_row(key, False):
        with ui.row().classes("br-pfield items-center gap-1 no-wrap"):
            inp = ui.input(value="" if value in _UNSET else str(value)).props(
                "filled dense").classes("grow min-w-0")
            inp.on_value_change(lambda e, k=key: app._update_param(k, e.value))

            async def browse(k=key, field=inp) -> None:
                if LocalFilePicker.active() is not None:
                    return
                state_key = f"last_dir_{k.lower().replace(' ', '_')}"
                start = get_state(state_key) or "~"
                if not isinstance(start, str) or (start != "~" and not os.path.isdir(start)):
                    start = "~"
                path = await LocalFilePicker(start, title=f"Select {_LABEL.get(k, k).lower()}")
                if path:
                    set_state(state_key, os.path.dirname(path))
                    field.value = path
            ui.button(icon="folder_open", on_click=browse).props("flat dense round color=primary")
    _register(app, key, inp, value)


#**************************************************************************************************#
#                                            the sheet                                             #
#**************************************************************************************************#
def build(app) -> None:
    """The whole parameter step for the current engine and sheet (rebuilt on every switch)."""
    br = app.BasisREMY
    b = br.backend
    seq, route = ss.current(br)
    whole = _whole_file(app)
    app._fields, app._tooltips = {}, {}
    app._rec = {} if whole else ss.apply_recommended(b, seq, _keep(app))
    shown = b.get_params_for_mode()

    with ui.column().classes("w-full gap-5"):
        _selectors(app, seq, route, whole)
        _legend()
        grid = ui.element("div").classes("br-pgrid")
        with grid:
            with ui.column().classes("gap-4 min-w-0"):
                _pulses_card(app, seq, route, whole, shown)
                timing = [k for k in _TIMING if k in shown]
                if b.name == 'CustomSLaser':          # its Tau 1/2 only set the reference's ppm range
                    timing = [k for k in timing if k == 'TE']
                _card(app, "Timings", timing, note=ss.fixed_timing(b) if seq and not whole else None)
                _card(app, "Acquisition", [k for k in _ACQ if k in shown])
                _card(app, "Engine settings", [k for k in shown if k not in _ELSEWHERE
                                               and k not in timing and k not in _ACQ], modes=route)
            metabs = ui.column().classes("gap-2 min-w-0")
            with metabs:
                app.metabs_col = ui.column().classes("w-full gap-0")
                if 'Metabolites' in shown:
                    app._param_metabolites()
        if not app.metab_checks:
            metabs.set_visibility(False)
            grid.classes(add="br-pgrid--single")


def _legend() -> None:
    with ui.row().classes("items-center gap-4 text-xs br-muted -mt-2"):
        for css, text in (('file', 'from the data file'), ('rec', 'recommended / default'),
                          ('user', 'set by you'), ('missing', 'missing')):
            with ui.row().classes("items-center gap-1.5 no-wrap"):
                ui.element("span").classes(f"br-dot br-dot-{css}")
                ui.label(text)


def _row(label: str, help_key: str | None = None):
    row = ui.element("div").classes("br-prow")
    with row:
        if help_key:
            label_with_help(help_key, label).classes("br-prow-label")
        else:
            ui.label(label).classes("br-prow-label text-sm font-semibold")
    return row


# ---- Engine / Sequence / file ---------------------------------------------------------------------
def _selectors(app, seq, route, whole) -> None:
    br = app.BasisREMY
    cat = br.backend.category
    with ui.column().classes("br-card br-plist w-full"):
        eng_opts = {c: ss.ENGINE_LABEL.get(c, c) for c in ss.engines() if br.categories.get(c)}
        if seq:
            eng_opts = {c: lbl if seq in ss.ROUTES[c] else f"{lbl}  ·  no {seq}"
                        for c, lbl in eng_opts.items()}
        with _row("Engine"):
            eng = ui.select(eng_opts, value=cat).props("filled dense").classes("br-selfield")
        eng.on_value_change(lambda e: _on_engine(app, e.value, eng))

        with _row("Sequence"):
            with ui.icon("info_outline").classes("br-muted cursor-help text-base"):
                ui.tooltip(ss.ENGINE_SCOPE.get(cat, "")).classes("max-w-xs")
            if whole:
                ui.label("from the file").classes("br-selfield text-sm br-muted")
            else:
                options = {s: s if s in ss.ROUTES[cat] else f"{s}  ·  not on {ss.ENGINE_LABEL[cat]}"
                           for s in ss.SEQUENCES}
                sel = _Select(options, disabled=[s for s in ss.SEQUENCES if s not in ss.ROUTES[cat]],
                              value=seq).props("filled dense").classes("br-selfield")
                sel.classes(add="br-v-missing" if seq is None else
                            "br-v-file" if 'Sequence' in br.from_file.get(br.backend.name, {})
                            and 'Sequence' not in app._user_set else "br-v-user")
                sel.on_value_change(lambda e: _on_sequence(app, e.value))
        _file_row(app, whole)
        if seq is None and not whole:
            ui.label("Pick the sequence that acquired the data.").classes(
                "text-xs br-muted px-4 pb-2")


async def _on_engine(app, cat, select) -> None:
    br = app.BasisREMY
    if cat == br.backend.category or getattr(app, "_switching", False):
        return
    seq, _ = ss.current(br)
    route = ss.choose(cat, seq, app._want) if seq else None
    target = br.backends[route.backend if route else br.categories[cat][0]]
    app._switching = True
    try:
        if target.requires_octave and target.octave is None and not await _octave_ready(app):
            select.value = br.backend.category
            return
    finally:
        app._switching = False
    app.seq_file = None
    if route:
        ss.apply(br, cat, seq, app._want)
    else:
        br.set_category(cat)
        if 'Sequence' in br.backend.mandatory_params:
            br.backend.mandatory_params['Sequence'] = None
        if seq:
            ui.notify(f"{ss.ENGINE_LABEL[cat]} does not run {seq}: pick a sequence.", type="warning")
    app._rebuild_soon()


async def _octave_ready(app) -> bool:
    from basisremy.core.octave_manager import OctaveManager
    manager = OctaveManager()
    note = ui.notification("Checking for Docker / Octave…", spinner=True, timeout=None)
    try:
        ok = await run.io_bound(lambda: manager.check_docker_availability()
                                or manager.check_local_octave_availability())
    finally:
        note.dismiss()
    if not ok:
        app._show_octave_instructions(manager)
    return ok


def _on_sequence(app, seq) -> None:
    br = app.BasisREMY
    if seq is None or seq == ss.current(br)[0]:
        return
    app._user_set.add('Sequence')
    app.seq_file = None
    ss.apply(br, br.backend.category, seq, app._want)
    app._rebuild_soon()


def _whole_file(app) -> tuple[str, str] | None:
    """(name, what) of the whole-sequence file the engine runs, None without one."""
    b = app.BasisREMY.backend
    if app.seq_file:
        return os.path.basename(app.seq_file), "Pulseq sequence: timing, pulses and slabs"
    if b.name == 'FSL-MRS' and b.current_mode == 'Custom' and b.optional_params.get('Custom Sequence'):
        return os.path.basename(b.optional_params['Custom Sequence']), "Sequence description, run as given"
    if b.name == 'FSL-MRS' and b.current_mode == 'Template' and b.optional_params.get('Template File'):
        return b.optional_params['Template File'], "FSL-MRS example sequence, run as given"
    return None


def _file_row(app, whole) -> None:
    br = app.BasisREMY
    with _row("Sequence / pulse file"):
        with ui.row().classes("br-selfield items-center gap-1 no-wrap"):
            if whole:
                with ui.column().classes("gap-0 min-w-0 grow"):
                    ui.label(whole[0]).classes("text-sm font-semibold truncate w-full")
                    ui.label(whole[1]).classes("text-xs br-muted truncate w-full")
                ui.button(icon="close", on_click=lambda: _clear_file(app)).props(
                    "flat dense round").classes("br-muted").tooltip("Remove the file")
            else:
                ui.label("drop a file here, or browse").classes("text-xs br-muted grow")
            ui.button(icon="folder_open", on_click=lambda: _browse(app)).props(
                "flat dense round color=primary").tooltip(
                "A whole sequence (Pulseq .seq, sequence .json) or one pulse "
                "(.pta, .RF, .txt, Bruker, .mat)")
            if br.backend.category == 'FSL-MRS':
                with ui.button(icon="library_books").props("flat dense round color=primary") as lib:
                    with ui.menu():
                        for info in br.backend.predefined_sequences.values():
                            ui.menu_item(info['description'],
                                         on_click=lambda _, d=info['description']: _use_template(app, d))
                lib.tooltip("FSL-MRS example sequences")
            ui.button(icon="auto_fix_high", on_click=lambda: _open_panel(app)).props(
                "flat dense round color=primary").mark("sequence-panel").tooltip(
                "Sequence: timeline and pulses")


def _open_panel(app) -> None:
    from basisremy.gui.sequence_dialog import open_sequence_dialog
    open_sequence_dialog(app)


async def _browse(app) -> None:
    if LocalFilePicker.active() is not None:
        return
    start = get_state("last_dir_sequence_file") or get_state("last_import_dir") or "~"
    if not isinstance(start, str) or (start != "~" and not os.path.isdir(start)):
        start = "~"
    path = await LocalFilePicker(start, title="Select a sequence or pulse file")
    if path:
        set_state("last_dir_sequence_file", os.path.dirname(path))
        await use_file(app, path)


async def use_file(app, path: str) -> None:
    """A dropped or picked file: a whole sequence, or one pulse for a role."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.seq':
        app.seq_file = path
        if not app._apply_sequence_file():
            app.seq_file = None
        app._rebuild_soon()
    elif ext == '.json':
        _use_sequence_json(app, path)
    elif ext in _PULSE_EXTS or path.lower().endswith(('.exc', '.rfc', '.inv')):
        await _use_pulse_file(app, path)
    else:
        ui.notify(f"{os.path.basename(path)}: not a sequence or pulse file BasisREMY reads.",
                  type="warning")


def _use_sequence_json(app, path: str) -> None:
    br = app.BasisREMY
    try:
        seq = sv._fsl_sequence(path)
    except Exception as exc:                                # noqa: BLE001
        ui.notify(str(exc), type="negative")
        return
    app.seq_file = None
    if br.backend.name != 'FSL-MRS':
        br.set_backend('FSL-MRS')
    br.backend.set_mode('Custom')
    with open(path) as fh:
        nested = 'seq' in json.load(fh)
    if nested:          # a basis set JSON carries its sequence: FSL-MRS reads the bare description
        path = os.path.join(br.backend.ensure_workdir(), os.path.basename(path))
        with open(path, 'w') as fh:
            json.dump(seq, fh)
    br.backend.optional_params['Custom Sequence'] = path
    ui.notify(f"FSL-MRS runs {os.path.basename(path)} as given.", type="positive")
    app._rebuild_soon()


def _use_template(app, description: str) -> None:
    b = app.BasisREMY.backend
    app.seq_file = None
    b.set_mode('Template')
    b.optional_params['Template File'] = description
    app._rebuild_soon()


def _clear_file(app) -> None:
    b = app.BasisREMY.backend
    if app.seq_file:
        app.seq_file = None
    elif b.name == 'FSL-MRS':
        b.optional_params['Custom Sequence'] = None
        b.optional_params['Template File'] = None
        b.set_mode('Simple')
    app._rebuild_soon()


async def _use_pulse_file(app, path: str) -> None:
    br = app.BasisREMY
    seq, _ = ss.current(br)
    cat = br.backend.category
    if seq is None:
        ui.notify("Pick the sequence first, then the pulse.", type="warning")
        return
    can = [r for r in ss.roles(seq) if ss.SHAPED in ss.pulse_choices(cat, seq, r)]
    if not can:
        role = 'ref' if 'ref' in ss.roles(seq) else 'exc'
        ui.notify(ss.pulse_note(cat, seq, role, ss.SHAPED), type="warning", multi_line=True)
        return
    role = can[0]
    if len(can) > 1:
        with ui.dialog() as dlg, ui.card().classes("gap-3"):
            ui.label(f"Use {os.path.basename(path)} as").classes("text-sm font-semibold")
            with ui.row().classes("gap-2"):
                for r in can:
                    ui.button(ss.ROLE_NAME[r], on_click=lambda _, r=r: dlg.submit(r)).props(
                        "outline color=primary")
        role = await dlg
        if role is None:
            return
    set_pulse(app, seq, role, path)


# ---- pulses -------------------------------------------------------------------------------------
def set_pulse(app, seq, role, value) -> None:
    """Make ``role`` ideal or a waveform (standard name or file) and rebuild the sheet."""
    br = app.BasisREMY
    cat = br.backend.category
    app._want[role] = ss.IDEAL if value == _IDEAL else ss.SHAPED
    app._user_set.add(f"pulse:{role}")
    ss.apply(br, cat, seq, app._want)
    if value != _IDEAL:
        b = br.backend
        key = sv.pulse_key(b, role)
        if key:
            if os.path.splitext(str(value))[1].lower() in ('.seq', '.json'):
                value = f"{value}#{role}"
            app._set_value(key, value)
            try:
                tp = sv.read_pulse(value, sv._KIND[role]).tp_ms
            except Exception:                                # noqa: BLE001 - kept as typed
                tp = float('nan')
            if tp == tp and _TP_KEY[role] in b.get_params_for_mode():
                app._set_value(_TP_KEY[role], round(tp, 4))
    app._rebuild_soon()


def _own_text(b, role) -> str:
    if b.name == 'MRSCloud':
        vendor = b.mandatory_params.get('System') or 'vendor'
        return ("MRSCloud's universal pulse set" if b.current_mode == 'Universal'
                else f"MRSCloud's {vendor} product pulses")
    return f"Gaussian, generated by {ss.ENGINE_LABEL.get(b.category, b.category)}"


def _pulses_card(app, seq, route, whole, shown) -> None:
    b = app.BasisREMY.backend
    with ui.column().classes("gap-2 min-w-0 w-full"):
        ui.label("Pulses").classes("br-section-title")
        with ui.column().classes("br-card br-plist w-full"):
            if whole:
                with _row("All pulses"):
                    ui.label("from the sequence file").classes("text-sm br-muted")
                return
            if route is None:
                with _row("Pulses"):
                    ui.label("follow from the sequence").classes("text-sm br-muted")
                return
            for role in ss.roles(seq):
                _pulse_row(app, seq, route, role, shown)
            if 'Vendor Pulse File' in shown:
                file_field(app, 'Vendor Pulse File', b.mandatory_params.get('Vendor Pulse File'))


def _pulse_row(app, seq, route, role, shown) -> None:
    b = app.BasisREMY.backend
    cat = b.category
    kind = route.pulses.get(role)
    with _row(ss.ROLE_NAME[role]):
        if kind == ss.OWN or kind is None:
            ui.label(_own_text(b, role)).classes("br-selfield text-sm")
        else:
            key = sv.pulse_key(b, role) if kind == ss.SHAPED else None
            value = b.mandatory_params.get(key) if key else _IDEAL
            choices = ss.pulse_choices(cat, seq, role)
            options = {_IDEAL: "Ideal (instantaneous)"}
            options.update({f"standard:{n}": f"Standard: {n}" for n in _STANDARD_FOR[role]})
            if value not in options:
                options[value] = f"File: {os.path.basename(str(value))}"
            options[_FILE] = "Pulse file…"
            disabled = set()
            if ss.IDEAL not in choices:
                disabled.add(_IDEAL)
            if ss.SHAPED not in choices:
                disabled |= set(options) - {_IDEAL}
            sel = _Select(options, disabled=disabled, value=value).props("filled dense").classes(
                "br-selfield")
            user = f"pulse:{role}" in app._user_set or (key in app._user_set if key else False)
            sel.classes(add="br-v-user" if user else "br-v-rec")
            sel.tooltip("Set by you" if user else
                        "Recommended: ideal wherever the engine allows" if value == _IDEAL else
                        f"{ss.ENGINE_LABEL.get(cat, cat)} default / recommended")

            async def pick(e, current=value) -> None:
                if e.value in (None, current):
                    return
                new = e.value
                if new == _FILE:
                    new = await _pick_pulse_file(role)
                    if not new:
                        e.sender.value = current
                        return
                set_pulse(app, seq, role, new)
            sel.on_value_change(pick)
    # duration, editing frequencies: small rows under the role
    extra = [_TP_KEY[role]] if kind == ss.SHAPED and _TP_KEY[role] in shown else []
    if role == 'edit':
        extra = [k for k in _EDIT_FIELDS if k in shown]
    for k in extra:
        text_field(app, k, b.mandatory_params.get(k, b.optional_params.get(k)), sub=True)
    other = ss.SHAPED if kind == ss.IDEAL else ss.IDEAL if kind == ss.SHAPED else None
    note = ss.pulse_note(cat, seq, role, other) if other else None
    if note and not (role == 'edit' and other == ss.IDEAL):
        ui.label(note).classes("text-xs br-muted px-4 pb-2 -mt-1")


async def _pick_pulse_file(role) -> str | None:
    if LocalFilePicker.active() is not None:
        return None
    state_key = f"last_dir_pulse_{role}"
    start = get_state(state_key) or get_state("last_dir_sequence_file") or "~"
    if not isinstance(start, str) or (start != "~" and not os.path.isdir(start)):
        start = "~"
    path = await LocalFilePicker(start, title=f"Select the {ss.ROLE_NAME[role].lower()} pulse")
    if path:
        set_state(state_key, os.path.dirname(path))
    return path


# ---- value cards --------------------------------------------------------------------------------
def _card(app, title, keys, note=None, modes=None) -> None:
    b = app.BasisREMY.backend
    mode_row = modes is not None and len(modes.modes) > 1 and b.name in _MODE_LABEL
    if not keys and not note and not mode_row:
        return
    with ui.column().classes("gap-2 min-w-0 w-full"):
        ui.label(title).classes("br-section-title")
        with ui.column().classes("br-card br-plist w-full"):
            if mode_row:
                label, names = _MODE_LABEL[b.name]
                with _row(label):
                    ui.select({m: names.get(m, m) for m in modes.modes}, value=b.current_mode,
                              on_change=lambda e: app._change_mode(e.value)).props(
                        "filled dense").classes("br-pfield")
            for key in keys:
                value = b.mandatory_params.get(key, b.optional_params.get(key))
                if key in b.file_selection:
                    file_field(app, key, value)
                elif key in b.dropdown:
                    dropdown_field(app, key, value)
                else:
                    text_field(app, key, value)
            if note:
                ui.label(note).classes("text-xs br-muted px-4 py-2")
