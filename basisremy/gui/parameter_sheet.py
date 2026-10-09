####################################################################################################
#                                       parameter_sheet.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 07/10/26                                                                                #
#                                                                                                  #
# Purpose: The parameter step, laid out the same for every engine: Engine, Sequence and the        #
#          sequence file (none = ideal pulses, a dropped or picked file, or a design saved from    #
#          the wand's sequence designer); Timings, Acquisition and the engine's own settings; the  #
#          metabolites. Pulses live in the designer only. Every value shows where it comes from -  #
#          the data file, a recommendation (with its source), the user - or that it is missing.    #
#          What an engine cannot run is greyed out with the reason (core/sequence_setup.py,        #
#          core/sequence_design.py).                                                               #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os

from nicegui import run, ui

from basisremy.core import sequence_design as sd
from basisremy.core import sequence_setup as ss
from basisremy.core.parameter_registry import get as registry_get
from basisremy.gui.help_widget import label_with_help
from basisremy.gui.local_file_picker import LocalFilePicker
from basisremy.gui.ui_state import get_state, set_state

_UNSET = (None, "", "missing input", "Select option")
_PULSE_EXTS = {'.pta', '.rf', '.txt', '.exc', '.rfc', '.inv', '.mat'}
_NONE, _TEMPLATE = "", "template:"

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
# pulses and the slab they select: set in the sequence designer, not in the sheet
_PULSE_FIELDS = {'Path to Pulse', 'Edit Pulse Path', 'RefTp', 'Flip Angle', 'Pulse Phase',
                 'thkX', 'thkY', *_EDIT_FIELDS}
# handled by the selectors, the file field, the designer or the metabolite list
_ELSEWHERE = {'Sequence', 'Localization', 'Metabolites', 'Custom Sequence', 'Template File',
              'Vendor Pulse File', *_PULSE_FIELDS}
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


def text_field(app, key, value, sub=False, readonly=False) -> None:
    with _field_row(key, sub):
        inp = ui.input(value=_display(value)).props("filled dense").classes("br-pfield")
        inp.on_value_change(lambda e, k=key: app._update_param(k, e.value))
        if readonly:
            inp.props("readonly")
    _register(app, key, inp, value)
    if readonly:
        app._tooltips[key].set_text("From the sequence file: change it in the sequence designer (wand)")


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
                timing = [k for k in _TIMING if k in shown]
                if b.name == 'CustomSLaser':          # its Tau 1/2 only set the reference's ppm range
                    timing = [k for k in timing if k == 'TE']
                _card(app, "Timings", timing, note=ss.fixed_timing(b) if seq and not whole else None,
                      readonly=bool(app.seq_file))
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
        greyed = []
        if app.seq_file and app._design is not None:
            for c in eng_opts:
                if c != cat and sd.plan(app._design, c, br).status == 'no':
                    eng_opts[c] += "  ·  cannot run this sequence file"
                    greyed.append(c)
        elif seq:
            eng_opts = {c: lbl if seq in ss.ROUTES[c] else f"{lbl}  ·  no {seq}"
                        for c, lbl in eng_opts.items()}
        with _row("Engine"):
            eng = _Select(eng_opts, disabled=greyed, value=cat).props("filled dense").classes(
                "br-selfield")
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
    app._engine_chosen = True
    seq, _ = ss.current(br)
    design = app._design if app.seq_file else None
    if design is not None:
        pl = sd.plan(design, cat, br)
        route = pl.route
        name = 'FSL-MRS' if cat == 'FSL-MRS' else (route.backend if route else br.categories[cat][0])
    else:
        route = ss.choose(cat, seq, app._want) if seq else None
        name = route.backend if route else br.categories[cat][0]
    target = br.backends[name]
    app._switching = True
    try:
        if target.requires_octave and target.octave is None and not await _octave_ready(app):
            select.value = br.backend.category
            return
    finally:
        app._switching = False
    if design is not None:
        app._plan = sd.apply(br, design, app.seq_file, cat)
        if app._plan.status == 'no':
            _drop_file(app)
            br.set_category(cat)
            ui.notify(" ".join(app._plan.notes), type="warning", multi_line=True)
    elif route:
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
    ss.apply(br, br.backend.category, seq, app._want)
    app._rebuild_soon()


def _whole_file(app) -> tuple[str, str] | None:
    """(name, what) of the whole-sequence file the engine runs, None without one."""
    b = app.BasisREMY.backend
    if app.seq_file:
        return os.path.basename(app.seq_file), "Sequence file: timing, pulses and slabs"
    if b.name == 'FSL-MRS' and b.current_mode == 'Template' and b.optional_params.get('Template File'):
        return b.optional_params['Template File'], "FSL-MRS example sequence, run as given"
    return None


# ---- the sequence file ---------------------------------------------------------------------------
def _file_options(app) -> dict:
    """None (ideal pulses / the engine's own), files used in this session, saved designs, and
    FSL-MRS's example sequences."""
    b = app.BasisREMY.backend
    if b.category == 'MRSCloud':
        opts = {_NONE: "MRSCloud's own pulse set"}
        vendor = b.mandatory_params.get('Vendor Pulse File')
        if vendor not in _UNSET:
            opts[vendor] = f"Vendor pulse: {os.path.basename(str(vendor))}"
        for path in [app.seq_file, *app._seq_files, *sd.saved_designs()]:   # designs: PRESS / MEGA
            if path and path.lower().endswith('.seq'):
                opts.setdefault(path, os.path.basename(path))
        return opts
    opts = {_NONE: "Ideal pulses"}
    for path in [app.seq_file, *app._seq_files, *sd.saved_designs()]:
        if path:
            opts.setdefault(path, os.path.basename(path))
    if b.category == 'FSL-MRS':
        for info in b.predefined_sequences.values():
            opts[_TEMPLATE + info['description']] = f"FSL-MRS example: {info['description']}"
    return opts


def _file_value(app) -> str:
    b = app.BasisREMY.backend
    if app.seq_file:
        return app.seq_file
    if b.name == 'FSL-MRS' and b.current_mode == 'Template' and b.optional_params.get('Template File'):
        return _TEMPLATE + b.optional_params['Template File']
    if b.category == 'MRSCloud' and b.mandatory_params.get('Vendor Pulse File') not in _UNSET:
        return b.mandatory_params['Vendor Pulse File']
    return _NONE


def _file_row(app, whole) -> None:
    br = app.BasisREMY
    b = br.backend
    options = _file_options(app)
    value = _file_value(app)
    with _row("Sequence file"):
        with ui.row().classes("br-selfield items-center gap-1 no-wrap"):
            sel = ui.select(options, value=value if value in options else _NONE).props(
                "filled dense").classes("grow min-w-0")
            needs_vendor = 'Vendor Pulse File' in b.get_params_for_mode() and value == _NONE
            sel.classes(add="br-v-missing" if needs_vendor else
                        "br-v-file" if value != _NONE else "br-v-rec")
            sel.mark("sequence-file")
            with sel:
                ui.tooltip("MRSCloud needs the vendor's refocusing pulse file for this scanner: "
                           "drop or pick it" if needs_vendor else
                           "Drop or pick any sequence (Pulseq .seq, sequence .json) or pulse file "
                           "(.pta, .RF, .txt, Bruker, .mat), or make one with the wand")

            async def chosen(e) -> None:
                if e.value == value:
                    return
                if e.value == _NONE:
                    _clear_file(app)
                elif str(e.value).startswith(_TEMPLATE):
                    _use_template(app, e.value[len(_TEMPLATE):])
                else:
                    await use_file(app, e.value)
            sel.on_value_change(chosen)
            ui.button(icon="folder_open", on_click=lambda: _browse(app)).props(
                "flat dense round color=primary").tooltip("Pick a sequence or pulse file")
            ui.button(icon="auto_fix_high", on_click=lambda: _open_panel(app)).props(
                "flat dense round color=primary").mark("sequence-panel").tooltip(
                "Sequence designer: pulses, timings, save as a sequence file")
    pl = getattr(app, '_plan', None)
    if app.seq_file and pl is not None:
        colour = {'ok': 'br-muted', 'approx': '', 'no': ''}[pl.status]
        style = {'ok': '', 'approx': 'color:#a86d12', 'no': 'color:#c2453c'}[pl.status]
        text = " ".join(pl.notes) if pl.notes else "Runs as given."
        ui.label(text).classes(f"text-xs px-4 pb-2 {colour}").style(style)


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
    """A dropped or picked file: a whole sequence runs as it is; a single pulse becomes a design
    (the recommended one, with this pulse) saved next to the others and selected."""
    ext = os.path.splitext(path)[1].lower()
    b = app.BasisREMY.backend
    if ext in ('.seq', '.json'):
        previous = app.seq_file
        app.seq_file = path
        if not app._apply_sequence_file():
            app.seq_file = previous
        app._rebuild_soon()
    elif ext in _PULSE_EXTS or path.lower().endswith(('.exc', '.rfc', '.inv')):
        if b.category == 'MRSCloud':
            if 'Vendor Pulse File' in b.get_params_for_mode():
                app._set_value('Vendor Pulse File', path)
                app._rebuild_soon()
            else:
                ui.notify("MRSCloud uses its own pulse set; pick another engine for this pulse.",
                          type="warning")
            return
        await _pulse_design(app, path)
    else:
        ui.notify(f"{os.path.basename(path)}: not a sequence or pulse file BasisREMY reads.",
                  type="warning")


async def _pulse_design(app, path: str) -> None:
    from basisremy.gui.sequence_dialog import _ask_role
    br = app.BasisREMY
    b = br.backend
    sheet = {**b.optional_params, **b.mandatory_params}
    seq, _ = ss.current(br)
    te = sd._num(sheet.get('TE'))
    if seq not in sd.DESIGNABLE or te is None:
        ui.notify("Set the sequence and TE first, then the pulse.", type="warning")
        return
    d = app._design if app.seq_file and app._design is not None else sd.recommend(seq, te, sheet, getattr(br, '_last_mrsinmrs', None))
    role = await _ask_role([r for r in sd.roles(d.kind)], path)
    if role is None:
        return
    dur, _why = sd.default_duration(path, role)
    d.pulses[role] = {'source': path, 'dur': dur}
    target = os.path.join(sd.designs_dir(), sd.default_name(d) + ".seq")
    n = 2
    while os.path.exists(target):
        target = os.path.join(sd.designs_dir(), f"{sd.default_name(d)}_{n}.seq")
        n += 1
    try:
        sd.write_seq(d, target)
    except Exception as exc:                                # noqa: BLE001
        ui.notify(f"{os.path.basename(path)}: {exc}", type="negative", multi_line=True)
        return
    ui.notify(f"{os.path.basename(path)} as the {ss.ROLE_NAME[role].lower()} pulse: saved as "
              f"{os.path.basename(target)} (open the wand to change it).", type="positive",
              multi_line=True)
    await use_file(app, target)


def _use_template(app, description: str) -> None:
    _drop_file(app)
    b = app.BasisREMY.backend
    b.set_mode('Template')
    b.optional_params['Template File'] = description
    app._rebuild_soon()


def _drop_file(app) -> None:
    """Forget the sequence file: its values go back to recommendations, ideal pulses."""
    br = app.BasisREMY
    b = br.backend
    pl = getattr(app, '_plan', None)
    if pl is not None:
        given = br.from_file.get(b.name, {})
        for k in pl.values:
            if k in ss.SCAN_KEYS:
                continue
            given.pop(k, None)
            for params in (b.mandatory_params, b.optional_params):
                if k in params and k not in ('Path to Pulse', 'Edit Pulse Path'):
                    params[k] = None
    app.seq_file, app._design, app._plan = None, None, None
    if b.name == 'FSL-MRS':
        b.optional_params['Custom Sequence'] = None
        b.optional_params['Template File'] = None
        b.set_mode('Simple')
    if b.name == 'MRSCloud':
        b.optional_params['Sequence File'] = None
    seq, _ = ss.current(br)
    if seq:
        ss.apply(br, b.category, seq, {})


def _clear_file(app) -> None:
    b = app.BasisREMY.backend
    if b.category == 'MRSCloud' and b.mandatory_params.get('Vendor Pulse File') not in _UNSET:
        b.mandatory_params['Vendor Pulse File'] = None
    _drop_file(app)
    app._rebuild_soon()


# ---- value cards --------------------------------------------------------------------------------
def _card(app, title, keys, note=None, modes=None, readonly=False) -> None:
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
                    text_field(app, key, value, readonly=readonly)
            if note:
                ui.label(note).classes("text-xs br-muted px-4 py-2")
