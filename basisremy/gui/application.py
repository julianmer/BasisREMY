####################################################################################################
#                                           application.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 08/10/25                                                                                #
#                                                                                                  #
# Purpose: Defines the GUI application for the BasisREMY tool. Each tab is a different              #
#          step in the process, starting with the data selection and REMY extraction,              #
#          continuing with the parameter configuration, and ending with the basis set simulation.   #
#                                                                                                  #
#          NiceGUI front-end: a sleek, pure-web UI that runs in a native desktop window            #
#          (via pywebview) and installs cleanly under uvx — no system tcl/tk required. The         #
#          server runs on localhost, so local data files are read in place by path; nothing         #
#          is uploaded or copied.                                                                  #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import os
import threading
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # NiceGUI renders figures to SVG; no interactive backend needed
import matplotlib.pyplot as plt
import numpy as np

from nicegui import app, run, ui

# own
from basisremy.core.basisremy import BasisREMY
from basisremy.core.parameter_registry import get as registry_get, metabolite_full_name
from basisremy.gui.help_widget import label_with_help
from basisremy.gui.local_file_picker import LocalFilePicker
from basisremy.gui.ui_state import get_state, set_state
from basisremy.gui.export_dialog import open_export_dialog


# Brand palette sampled from the BasisREMY mouse logo (deep teal-navy + slate).
PRIMARY = "#15627f"
PRIMARY_DARK = "#0a3a4f"
ASSETS = Path(__file__).resolve().parent.parent / "assets" / "imgs"

# Values that count as "not filled in" for validation / dropdown placeholders.
_UNSET = (None, "", "missing input", "Select option")


# File suffixes / names that BasisREMY's REMY reader can actually parse. Used to
# filter the data-file picker so users can only pick processable files.
_MRS_SUFFIXES = {".dat", ".ima", ".rda", ".spar", ".7", ".nii"}


def _is_mrs_file(p: Path) -> bool:
    name = p.name.lower()
    if name.endswith(".nii.gz"):
        return True
    if p.suffix.lower() in _MRS_SUFFIXES:
        return True
    return "method" in name or "2dseq" in name


# Global stylesheet: minimal, modern, theme-aware (light + system/dark). All
# colours come from CSS variables that flip under Quasar's ``body--dark`` class.
_GLOBAL_CSS = f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap');

:root {{
  /* sampled from the mouse logo: deep teal-navy + slate, on a soft grey page */
  --br-primary: {PRIMARY};
  --br-primary-strong: #0c4257;
  --br-ink: #0a2c3b;
  --br-tint: rgba(21,98,127,0.10);
  --br-tint-strong: rgba(21,98,127,0.16);
  --br-bg: #f5f8fa;
  --br-surface: #ffffff;
  --br-fg: #0e2630;
  --br-muted: #5e7280;
  --br-line: #e7edf2;
  --br-field: #eef3f7;
  --br-hover: rgba(21,98,127,0.06);
  --br-step-bg: #e8eef3;
  --br-step-fg: #93a3af;
  --br-shadow: 0 1px 2px rgba(10,44,59,.05), 0 8px 24px -14px rgba(10,44,59,.20);
  --br-wm-a: #0a3a4f;
  --br-wm-b: #2f93b8;
}}
.body--dark {{
  --br-primary: #5cb6d8;
  --br-primary-strong: #82c9e4;
  --br-ink: #cfe3ee;
  --br-tint: rgba(92,182,216,0.14);
  --br-tint-strong: rgba(92,182,216,0.22);
  --br-bg: #0a1016;
  --br-surface: #121b24;
  --br-fg: #e6eef3;
  --br-muted: #8195a2;
  --br-line: #1e2a34;
  --br-field: #161f29;
  --br-hover: rgba(255,255,255,0.045);
  --br-step-bg: #1b2530;
  --br-step-fg: #6f7d8a;
  --br-shadow: 0 1px 2px rgba(0,0,0,.3), 0 12px 32px -18px rgba(0,0,0,.65);
  --br-wm-a: #8fd2ec;
  --br-wm-b: #5cb6d8;
}}

html, body {{
  font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
  background: var(--br-bg);
  color: var(--br-fg);
  -webkit-font-smoothing: antialiased;
  letter-spacing: -0.01em;
}}
body, .body--dark {{ background: var(--br-bg); color: var(--br-fg); }}
.q-page, .q-layout {{ background: var(--br-bg); }}
.material-icons {{ font-family: 'Material Icons' !important; }}
::selection {{ background: var(--br-tint); }}

/* ---- header ---------------------------------------------------------- */
.br-header {{
  border-bottom: 1px solid var(--br-line);
  background: var(--br-surface);
}}
.br-wordmark {{
  background: linear-gradient(95deg, var(--br-wm-a), var(--br-wm-b));
  -webkit-background-clip: text;
  background-clip: text;
  color: transparent;
  letter-spacing: -0.03em;
}}
/* brand mouse — a dark mouse in light mode, a light mouse in dark mode */
.br-logo {{ height: 40px; width: auto; display: inline-block; }}
.br-logo-dark {{ display: none; }}
.body--dark .br-logo-light {{ display: none; }}
.body--dark .br-logo-dark {{ display: inline-block; }}

/* faint mouse watermark — sits quietly behind the content */
.br-watermark {{
  position: fixed; right: -34px; bottom: -46px;
  width: 320px; height: auto; opacity: .05;
  pointer-events: none; z-index: 0; user-select: none;
}}
.br-wm-dark {{ display: none; }}
.body--dark .br-wm-light {{ display: none; }}
.body--dark .br-wm-dark {{ display: block; }}

/* ---- structure ------------------------------------------------------- */
.br-hairline {{ height: 1px; width: 100%; background: var(--br-line); }}
.br-section-title {{
  font-size: 11px; font-weight: 700; letter-spacing: 0.1em;
  text-transform: uppercase; color: var(--br-muted);
}}
.br-legend {{ border-left: 1px solid var(--br-line); }}
.br-muted {{ color: var(--br-muted); }}

/* keep fixed Quasar grey text readable in dark mode */
.body--dark .text-grey-9, .body--dark .text-grey-8 {{ color: var(--br-fg) !important; }}
.body--dark .text-grey-7, .body--dark .text-grey-6 {{ color: var(--br-muted) !important; }}

/* ---- sleek numbered stepper ----------------------------------------- */
.br-stepper {{ user-select: none; }}
.br-step {{
  display: flex; align-items: center; gap: 10px;
  padding: 7px 12px; border-radius: 999px; cursor: pointer;
  transition: background .15s ease;
}}
.br-step:hover {{ background: var(--br-hover); }}
.br-step-num {{
  width: 26px; height: 26px; border-radius: 50%;
  display: flex; align-items: center; justify-content: center;
  font-size: 12px; font-weight: 700;
  background: var(--br-step-bg); color: var(--br-step-fg);
  transition: all .18s ease;
}}
.br-step-label {{
  font-size: 14px; font-weight: 600; color: var(--br-step-fg);
  transition: color .18s ease;
}}
.br-step.is-active .br-step-num {{
  background: var(--br-primary); color: #fff;
  box-shadow: 0 6px 16px -6px rgba(96,115,137,.85);
}}
.br-step.is-active .br-step-label {{ color: var(--br-fg); }}
.br-step.is-done .br-step-num {{ background: var(--br-primary); color: #fff; }}
.br-step.is-done .br-step-label {{ color: var(--br-fg); opacity: .75; }}
.br-step.is-locked {{ cursor: default; }}
.br-step.is-locked:hover {{ background: transparent; }}
.br-step-line {{
  width: 44px; height: 2px; border-radius: 2px;
  background: var(--br-line); transition: background .18s ease;
}}
.br-step-line.is-done {{ background: var(--br-primary); }}

/* ---- dropzone -------------------------------------------------------- */
.br-drop {{
  border: 1.5px dashed var(--br-line); border-radius: 20px;
  background: var(--br-surface);
  transition: border-color .18s ease, background .18s ease,
              transform .18s ease, box-shadow .18s ease;
}}
.br-drop:hover {{
  border-color: var(--br-primary);
  background: var(--br-tint);
  transform: translateY(-2px);
  box-shadow: var(--br-shadow);
}}
.br-icon-badge {{
  width: 60px; height: 60px; border-radius: 18px;
  display: flex; align-items: center; justify-content: center;
  background: var(--br-tint);
}}

/* ---- selected-file card --------------------------------------------- */
.br-filecard {{
  border: 1px solid var(--br-line); border-radius: 16px;
  background: var(--br-surface); padding: 13px 14px;
  box-shadow: var(--br-shadow);
}}
.br-file-ic {{
  width: 40px; height: 40px; border-radius: 12px; flex-shrink: 0;
  display: flex; align-items: center; justify-content: center;
  background: var(--br-tint);
}}

/* ---- soft filled controls ------------------------------------------- */
.q-field--filled .q-field__control,
.body--dark .q-field--filled .q-field__control {{
  background: var(--br-field);
  border-radius: 10px;
}}
.q-field--filled .q-field__control:before,
.q-field--filled .q-field__control:after {{ display: none; }}
.q-field--filled.q-field--focused .q-field__control {{
  box-shadow: 0 0 0 2px var(--br-tint-strong);
}}
.q-btn {{
  border-radius: 10px; text-transform: none;
  font-weight: 600; letter-spacing: 0;
}}
.q-btn.q-btn--standard {{ box-shadow: 0 12px 24px -16px var(--br-primary); }}
.br-tab-panel-pad {{ padding: 44px 2px 12px; }}

/* ---- settings-style cards & parameter rows -------------------------- */
.br-card {{
  border: 1px solid var(--br-line); border-radius: 16px;
  background: var(--br-surface); overflow: hidden;
  box-shadow: var(--br-shadow);
}}
.br-plist {{ display: flex; flex-direction: column; width: 100%; gap: 0 !important; }}
.br-prow {{
  display: flex; align-items: center; gap: 14px; width: 100%;
  padding: 0 16px; border-top: 1px solid var(--br-line);
  min-height: 38px;
  transition: background .12s ease;
}}
.br-prow:first-child {{ border-top: none; }}
.br-prow:hover {{ background: var(--br-hover); }}
.br-prow-label {{ flex: 1 1 auto; min-width: 0; font-size: 13px; }}
.br-pfield {{ flex: 0 0 auto; width: 184px; }}
/* wider selector fields for the Simulation Software / Mode pickers */
.br-selfield {{ flex: 0 0 auto; width: 280px; }}
/* compact, low-profile fields inside settings rows (overrides Quasar's 40px) */
.br-prow .q-field--filled .q-field__control {{ background: var(--br-field); }}
.br-prow .q-field__control {{
  min-height: 30px !important; height: 30px !important; padding: 0 10px !important;
}}
.br-prow .q-field__control-container {{ padding-top: 0 !important; min-height: 30px; }}
.br-prow .q-field__native, .br-prow .q-field__input {{
  min-height: 30px !important; height: 30px !important;
  padding: 0; font-size: 13px; line-height: 30px;
}}
.br-prow .q-field__marginal {{ height: 30px; }}
/* drop the reserved hint/error row so settings rows stay tight */
.br-prow .q-field__bottom {{ display: none !important; }}
.br-prow .q-field--with-bottom {{ padding-bottom: 0 !important; }}
.br-card-head {{
  display: flex; align-items: center; justify-content: space-between;
  padding: 9px 14px 7px;
}}
/* two-column settings layout: top-aligned so cards keep their natural
   height (no stretched empty space under the shorter column) */
.br-pgrid {{
  display: grid; grid-template-columns: 1fr 1fr;
  gap: 32px; align-items: start; width: 100%;
}}
.br-pgrid--single {{ grid-template-columns: minmax(0, 560px); }}
/* dense metabolite grid */
.br-metab .q-checkbox {{ min-height: 0; }}
.br-metab .q-checkbox__inner {{ font-size: 30px; }}
.br-metab .q-checkbox__label {{ font-size: 13px; padding-left: 2px; }}
</style>
"""


#**************************************************************************************************#
#                                         BasisREMYApp                                              #
#**************************************************************************************************#
#                                                                                                  #
# Builds the whole single-page UI. One instance is created per client connection (in native        #
# mode there is exactly one). All simulation state lives on the instance.                          #
#                                                                                                  #
#**************************************************************************************************#
class BasisREMYApp:

    def __init__(self) -> None:
        # data backend
        self.BasisREMY = BasisREMY()

        # selection / simulation state
        self.selected_file: str | None = None
        self.basis_set: dict | None = None
        self._basis_set_valid = False

        # simulation threading
        self._sim_stop_event = threading.Event()
        self._sim_thread: threading.Thread | None = None
        self._sim_timer = None
        self._sim_step = 0
        self._sim_total = 1
        self._sim_done = False
        self._sim_error: Exception | None = None

        # widget handles rebuilt on every backend/mode change
        self.metab_checks: dict = {}
        self.simulate_button = None

        # step navigation state (custom sleek stepper drives the tab panels)
        self._current_step = "data"
        # Native-window file drops (real paths) arrive as app-level events
        # from NiceGUI's pywebview bridge, not as DOM events.
        app.native.on("drop", self._on_native_drop)
        self._step_unlocked = {"data": True, "params": False, "sim": False}
        self._step_items: dict = {}
        self._step_nums: dict = {}
        self._step_labels: dict = {}
        self._step_lines: list = []

        self._build()

    # ============================================================== top-level UI
    def _build(self) -> None:
        # Theme controller: start on the system preference (auto); the header
        # toggle flips between explicit light / dark.
        self.darkmode = ui.dark_mode(value=None)

        # Decorative brand mouse, faint, fixed in the lower-right corner.
        ui.html(
            '<img src="/assets/basisremy_mouse_all_colors/png/'
            'basisremy_mouse_charcoal.png" class="br-watermark br-wm-light" alt="" />'
            '<img src="/assets/basisremy_mouse_all_colors/png/'
            'basisremy_mouse_light_gray.png" class="br-watermark br-wm-dark" alt="" />'
        )

        self._header()

        with ui.column().classes(
            "w-full max-w-4xl mx-auto px-6 pt-6 pb-10 gap-2 relative z-10"
        ):
            self._stepper()

            # Hidden Quasar tabs act purely as the controller for the panels;
            # navigation happens through the custom stepper above.
            with ui.tabs().classes("hidden") as self.tabs:
                self.tab1 = ui.tab("data")
                self.tab2 = ui.tab("params")
                self.tab3 = ui.tab("sim")

            with ui.tab_panels(self.tabs, value="data").classes(
                "w-full bg-transparent"
            ) as self.panels:
                with ui.tab_panel("data").classes("br-tab-panel-pad"):
                    self._build_tab1()
                self.panel2 = ui.tab_panel("params").classes("br-tab-panel-pad")
                self.panel3 = ui.tab_panel("sim").classes("br-tab-panel-pad")

            self.panels.on_value_change(self._on_tab_changed)

        self._build_tab2()
        self._build_tab3_progress()
        self._refresh_stepper()

    def _header(self) -> None:
        with ui.element("div").classes("br-header w-full relative z-10"):
            with ui.row().classes(
                "w-full max-w-4xl mx-auto items-center no-wrap px-6 py-3 gap-2.5"
            ):
                ui.html(
                    '<img src="/assets/basisremy_mouse_all_colors/png/'
                    'basisremy_mouse_navy_blue.png" class="br-logo br-logo-light" '
                    'alt="BasisREMY" />'
                    '<img src="/assets/basisremy_mouse_all_colors/png/'
                    'basisremy_mouse_sky_blue.png" class="br-logo br-logo-dark" '
                    'alt="BasisREMY" />'
                ).classes("shrink-0 leading-none")
                ui.label("BasisREMY").classes(
                    "br-wordmark text-xl font-extrabold leading-none"
                )
                ui.space()
                toggle = ui.button(
                    icon="dark_mode", on_click=self._toggle_theme
                ).props("flat round dense").classes("br-muted")
                toggle.bind_icon_from(
                    self.darkmode, "value",
                    backward=lambda v: "light_mode" if v else "dark_mode",
                )
                toggle.tooltip("Toggle light / dark theme")

    async def _toggle_theme(self) -> None:
        # ``darkmode`` starts in auto mode (value=None), where Quasar resolves the
        # theme client-side. Read the rendered state and flip to the explicit
        # opposite so the first press always works.
        is_dark = await ui.run_javascript(
            "document.body.classList.contains('body--dark')"
        )
        self.darkmode.value = not is_dark

    # ============================================================== stepper
    def _stepper(self) -> None:
        steps = [("data", "Data"), ("params", "Parameters"), ("sim", "Simulate")]
        with ui.row().classes("br-stepper items-center justify-center w-full py-1"):
            for i, (key, label) in enumerate(steps):
                if i > 0:
                    self._step_lines.append(
                        ui.element("div").classes("br-step-line")
                    )
                item = ui.row().classes("br-step").on(
                    "click", lambda k=key: self._step_click(k)
                )
                with item:
                    self._step_nums[key] = ui.label(str(i + 1)).classes("br-step-num")
                    self._step_labels[key] = ui.label(label).classes("br-step-label")
                self._step_items[key] = item

    def _refresh_stepper(self) -> None:
        order = ["data", "params", "sim"]
        try:
            current = self.panels.value or self._current_step
        except Exception:
            current = self._current_step
        cur_i = order.index(current) if current in order else 0

        for i, key in enumerate(order):
            item = self._step_items.get(key)
            num = self._step_nums.get(key)
            if item is None or num is None:
                continue
            item.classes(remove="is-active is-done is-locked")
            if key == current:
                item.classes(add="is-active")
                num.set_text(str(i + 1))
            elif self._step_unlocked.get(key) and i < cur_i:
                item.classes(add="is-done")
                num.set_text("✓")
            elif self._step_unlocked.get(key):
                num.set_text(str(i + 1))
            else:
                item.classes(add="is-locked")
                num.set_text(str(i + 1))

        for i, line in enumerate(self._step_lines):
            line.classes(remove="is-done")
            if self._step_unlocked.get(order[i + 1]):
                line.classes(add="is-done")

    def _step_click(self, key: str) -> None:
        if self._step_unlocked.get(key):
            self._goto(key)

    def _unlock(self, key: str) -> None:
        self._step_unlocked[key] = True
        self._refresh_stepper()

    def _lock(self, key: str) -> None:
        self._step_unlocked[key] = False
        self._refresh_stepper()

    # ============================================================== tab navigation
    def _on_tab_changed(self, _event=None) -> None:
        """Invalidate stale results when the user steps back to edit params."""
        try:
            current = self.panels.value
        except Exception:
            return
        self._current_step = current

        if current == "params":
            if self._basis_set_valid:
                self._basis_set_valid = False
                self._step_unlocked["sim"] = False
                self.validate_inputs()

            # cancel any in-flight simulation immediately
            if self._sim_thread is not None and self._sim_thread.is_alive():
                print("⏹  Cancelling simulation (user navigated back)…")
                self._sim_stop_event.set()
                # step 3 must be re-earned via Simulate; otherwise it shows a
                # frozen progress page whose poll timer is gone
                self._step_unlocked["sim"] = False
                try:
                    from basisremy.docker.docker_octave import DockerOctave
                    octave = self.BasisREMY.backend.octave
                    # isinstance, not hasattr: oct2py resolves unknown names in
                    # the Octave workspace and raises
                    if isinstance(octave, DockerOctave):
                        octave.kill_running_processes()
                except Exception as e:  # noqa: BLE001
                    print(f"  (Could not kill Docker Octave process: {e})")

        self._refresh_stepper()

    def _goto(self, value: str) -> None:
        self.panels.set_value(value)

    # ============================================================== TAB 1
    def _build_tab1(self) -> None:
        with ui.column().classes("w-full items-center gap-5 pt-6"):
            self._data_body = ui.column().classes(
                "w-full max-w-md items-stretch gap-3"
            )
            self._render_data_body()

            with ui.row().classes("items-center gap-1"):
                self.process_button = ui.button(
                    "Continue", on_click=self._process_file
                ).props("color=primary unelevated icon-right=arrow_forward")
                self.process_button.disable()
                ui.button("Skip", on_click=self._skip_file).props(
                    "flat color=primary"
                )

    def _render_data_body(self) -> None:
        self._data_body.clear()
        with self._data_body:
            if self.selected_file:
                self._file_card()
            else:
                self._dropzone()

    def _dropzone(self) -> None:
        drop = ui.column().classes(
            "br-drop w-full items-center justify-center cursor-pointer py-14 gap-3"
        )
        with drop:
            with ui.element("div").classes("br-icon-badge"):
                ui.icon("upload_file").classes("text-2xl").style("color:var(--br-primary)")
            ui.label("Drop MRS data or click to browse").classes(
                "text-sm font-medium br-muted"
            )
        drop.on("click", self._pick_data_file)
        # In the native window, drops are delivered app-wide by NiceGUI's
        # pywebview bridge with real file paths (see _on_native_drop). These
        # element handlers cover the browser, where real paths are unavailable:
        # prevent the tab from navigating away and explain.
        drop.on("dragover.prevent", js_handler="(e) => e.preventDefault()")
        drop.on("drop.prevent", self._on_browser_drop, [])

    def _on_browser_drop(self) -> None:
        if app.native.main_window is not None:
            return  # the native drop event carries the real path instead
        ui.notify("Drops only work in the desktop window — please click "
                  "the box to browse instead.", type="warning")

    def _on_native_drop(self, e) -> None:
        paths = [p for p in (e.args.get("files") or []) if p and os.path.exists(p)]
        if not paths:
            return
        # An open file picker (data file, pulse waveform, custom sequence,
        # export browse) takes priority: navigate it to the dropped path.
        picker = LocalFilePicker.active()
        if picker is not None:
            picker.show_dropped_path(paths[0])
            return
        if self._current_step != "data":
            return
        mrs = [p for p in paths if _is_mrs_file(Path(p))]
        path = (mrs or paths)[0]
        self.selected_file = path
        set_state("last_import_dir", os.path.dirname(path))
        self.process_button.enable()
        self._render_data_body()

    def _file_card(self) -> None:
        name = Path(self.selected_file).name
        folder = str(Path(self.selected_file).parent)
        with ui.row().classes("br-filecard w-full items-center gap-3 no-wrap"):
            with ui.element("div").classes("br-file-ic"):
                ui.icon("description").classes("text-xl").style("color:var(--br-primary)")
            with ui.column().classes("min-w-0 grow gap-0"):
                ui.label(name).classes("text-sm font-semibold truncate w-full")
                ui.label(folder).classes("text-xs br-muted truncate w-full")
            ui.button(
                icon="close", on_click=self._clear_file
            ).props("flat round dense").classes("br-muted shrink-0")

    def _clear_file(self) -> None:
        self.selected_file = None
        self.process_button.disable()
        self._render_data_body()

    async def _pick_data_file(self) -> None:
        if LocalFilePicker.active() is not None:
            return  # double-click on the dropzone must not stack two dialogs
        start = get_state("last_import_dir") or "~"
        if not isinstance(start, str) or (start != "~" and not os.path.isdir(start)):
            start = "~"
        path = await LocalFilePicker(
            start, title="Select MRS data file", show_file=_is_mrs_file
        )
        if path:
            set_state("last_import_dir", os.path.dirname(path))
            self.selected_file = path
            self.process_button.enable()
            self._render_data_body()

    async def _process_file(self) -> None:
        if not self.selected_file:
            ui.notify("No file selected.", type="warning")
            return
        if getattr(self, "_processing", False):
            return  # a double-click fits inside one client round-trip
        self._processing = True
        picked = self.selected_file
        print(f"Processing file: {picked}")
        self.process_button.props("loading")
        try:
            # Parsing can take seconds for large files — keep the UI alive.
            MRSinMRS = await run.io_bound(self.BasisREMY.runREMY, picked)
            if self.selected_file != picked:
                return  # file cleared or replaced while parsing — discard
            # A new file starts from clean defaults — values from the
            # previous file must not masquerade as this file's metadata.
            self.BasisREMY.reset_backend_params()
            params, opt = self.BasisREMY.backend.parseREMY(MRSinMRS)
            # drop None so REMY gaps don't clobber backend defaults
            self.BasisREMY.backend.mandatory_params.update(
                {k: v for k, v in params.items() if v is not None})
            self.BasisREMY.backend.optional_params.update(
                {k: v for k, v in opt.items() if v is not None})
            self.BasisREMY._last_mrsinmrs = MRSinMRS
        except Exception as exc:  # noqa: BLE001
            ui.notify(f"Could not read file: {exc}", type="negative")
            print(f"REMY error: {exc}")
            return
        finally:
            self._processing = False
            self.process_button.props(remove="loading")

        self._build_tab2()
        self._unlock("params")
        self._goto("params")

    def _skip_file(self) -> None:
        self._unlock("params")
        self._goto("params")

    # ============================================================== TAB 2
    def _build_tab2(self) -> None:
        self.panel2.clear()
        self.metab_checks = {}
        self.simulate_button = None
        backend = self.BasisREMY.backend

        # Settle mode-dependent state (current_mode, dropdown options,
        # file_selection) BEFORE drawing the selectors so the Mode picker and
        # the parameter list below it stay consistent (e.g. a parsed Philips
        # vendor steers MRSCloud to Non-universal up front).
        params_to_show = backend.get_params_for_mode()

        with self.panel2:
            with ui.column().classes("w-full gap-6"):
                self._backend_selectors()

                self._pgrid = ui.element("div").classes("br-pgrid")
                with self._pgrid:
                    with ui.column().classes("gap-2 min-w-0"):
                        ui.label("Parameters").classes("br-section-title")
                        self.params_col = ui.column().classes(
                            "br-card br-plist w-full"
                        )
                    self._metabs_wrap = ui.column().classes("gap-2 min-w-0")
                    with self._metabs_wrap:
                        self.metabs_col = ui.column().classes("w-full gap-0")

                # Per-backend mode selector. For single-backend software
                # (FSL-MRS, MRSCloud) the mode now lives in the selectors card
                # above, so only show it here for multi-backend software
                # (e.g. FID-A's semi-LASER Standard / Phase-cycled).
                cat_backends = self.BasisREMY.categories.get(
                    self.BasisREMY.get_current_category(), []
                )
                if len(cat_backends) > 1 and len(backend.modes) > 1:
                    with self.params_col:
                        with ui.element("div").classes("br-prow"):
                            ui.label("Mode").classes(
                                "br-prow-label text-sm font-semibold"
                            )
                            ui.select(
                                backend.modes,
                                value=backend.current_mode,
                                on_change=lambda e: self._change_mode(e.value),
                            ).props("filled dense").classes("br-pfield")

                for key, value in params_to_show.items():
                    if key in backend.file_selection:
                        self._param_file(key, value)
                    elif key == "Metabolites":
                        self._param_metabolites()
                    elif key in backend.dropdown:
                        self._param_dropdown(key, value)
                    else:
                        self._param_text(key, value)

                # Hide the (empty) metabolite column for backends without a
                # metabolite list so the parameters span a single tidy column.
                if not self.metab_checks:
                    self._metabs_wrap.set_visibility(False)
                    self._pgrid.classes(add="br-pgrid--single")

                ui.element("div").classes("br-hairline")
                with ui.row().classes("w-full justify-between items-center"):
                    ui.button("Back", icon="arrow_back",
                              on_click=lambda: self._goto("data")).props("flat color=primary")
                    self.simulate_button = ui.button(
                        "Simulate basis set", icon="auto_awesome",
                        on_click=self._simulate_basis,
                    ).props("color=primary unelevated")
                    self.simulate_button.disable()

        self.validate_inputs()

    # ---- backend / category selectors -------------------------------------
    def _backend_selectors(self) -> None:
        br = self.BasisREMY
        current_category = br.get_current_category()
        category_options = [c for c in br.CATEGORY_ORDER if br.categories.get(c)]
        for c in br.categories:
            if c not in category_options and br.categories[c]:
                category_options.append(c)

        def backends_for(cat):
            names = br.categories.get(cat, [])
            label_to_name, labels = {}, []
            for n in names:
                b = br.backends[n]
                label = getattr(b, "display_name", None) or b.name
                if getattr(b, "_is_stub", False):
                    label += " (in development)"
                label_to_name[label] = n
                labels.append(label)
            return labels, label_to_name

        async def do_switch(target_name) -> bool:
            if target_name == br.backend.name:
                return True
            new_backend = br.backends[target_name]
            if new_backend.requires_octave and new_backend.octave is None:
                from basisremy.core.octave_manager import OctaveManager
                manager = OctaveManager()
                # probes take seconds — don't freeze the dropdown click
                available = await run.io_bound(
                    lambda: manager.check_docker_availability()
                    or manager.check_local_octave_availability())
                if not available:
                    self._show_octave_instructions(manager)
                    return False
            br.set_backend(target_name)
            return True

        with ui.column().classes("br-card br-plist w-full"):
            with ui.element("div").classes("br-prow"):
                ui.label("Simulation Software").classes(
                    "br-prow-label text-sm font-semibold"
                )
                category_select = ui.select(
                    category_options, value=current_category
                ).props("filled dense").classes("br-selfield")

            labels, self._backend_label_map = backends_for(current_category)
            current_label = next(
                (lbl for lbl, nm in self._backend_label_map.items()
                 if nm == br.backend.name),
                labels[0] if labels else "",
            )
            backend_select = None
            if len(labels) > 1:
                # Multiple backends in this software (FID-A): the backend list
                # IS the "Mode".
                with ui.element("div").classes("br-prow"):
                    ui.label("Mode").classes(
                        "br-prow-label text-sm font-semibold"
                    )
                    backend_select = ui.select(
                        labels, value=current_label
                    ).props("filled dense").classes("br-selfield")
            elif len(br.backend.modes) > 1:
                # Single backend exposing several modes (FSL-MRS, MRSCloud):
                # show the mode selector right under the software picker.
                with ui.element("div").classes("br-prow"):
                    ui.label("Mode").classes(
                        "br-prow-label text-sm font-semibold"
                    )
                    ui.select(
                        br.backend.modes,
                        value=br.backend.current_mode,
                        on_change=lambda e: self._change_mode(e.value),
                    ).props("filled dense").classes("br-selfield")

        async def on_category_change(e) -> None:
            cat = e.value
            new_labels, label_map = backends_for(cat)
            if not new_labels:
                return
            self._backend_label_map = label_map
            target_name = label_map[new_labels[0]]
            if target_name == br.backend.name:
                return  # revert echo — leave the panel (and any open dialog) be
            if getattr(self, "_switching", False):
                return
            self._switching = True
            try:
                ok = await do_switch(target_name)
            finally:
                self._switching = False
            if ok:
                self._build_tab2()
            else:
                category_select.value = br.get_current_category()

        async def on_backend_change(e) -> None:
            target_name = self._backend_label_map.get(e.value)
            if target_name is None or target_name == br.backend.name:
                return  # unknown label or revert echo
            if getattr(br.backends[target_name], "_is_stub", False):
                ui.notify("This shaped sequence is under development — "
                          "simulation support is coming soon.", type="info")
                cur = next((lbl for lbl, nm in self._backend_label_map.items()
                            if nm == br.backend.name), e.value)
                backend_select.value = cur
                return
            if getattr(self, "_switching", False):
                return
            self._switching = True
            try:
                ok = await do_switch(target_name)
            finally:
                self._switching = False
            if ok:
                self._build_tab2()
            else:
                cur = next((lbl for lbl, nm in self._backend_label_map.items()
                            if nm == br.backend.name), e.value)
                backend_select.value = cur

        category_select.on_value_change(on_category_change)
        if backend_select is not None:
            backend_select.on_value_change(on_backend_change)

    def _change_mode(self, mode: str) -> None:
        self.BasisREMY.backend.set_mode(mode)
        self._build_tab2()
        # Some backends veto a mode for the current selection (e.g. MRSCloud
        # forces Non-Universal for GE, which has no Universal waveform set) —
        # say so instead of silently snapping the dropdown back.
        actual = self.BasisREMY.backend.current_mode
        if actual != mode:
            ui.notify(f"'{mode}' is not available for the current System "
                      f"selection — kept '{actual}'.", type="warning")

    # ---- individual parameter widgets -------------------------------------
    def _update_param(self, key: str, value) -> None:
        backend = self.BasisREMY.backend
        if key in backend.mandatory_params:
            backend.mandatory_params[key] = value
        elif key in backend.optional_params:
            backend.optional_params[key] = value
        self.validate_inputs()
        if key in getattr(backend, "schema_affecting_keys", set()):
            self._build_tab2()

    def _param_text(self, key, value) -> None:
        with self.params_col:
            with ui.element("div").classes("br-prow"):
                label_with_help(key).classes("br-prow-label")
                inp = ui.input(
                    value="" if value is None else str(value),
                ).props("filled dense").classes("br-pfield")
                inp.on_value_change(lambda e, k=key: self._update_param(k, e.value))

    def _param_dropdown(self, key, value) -> None:
        backend = self.BasisREMY.backend
        options = backend.dropdown[key]
        # ``options`` may be a list (label == value) or a dict (value -> label).
        keys = list(options)
        initial = str(value) if (
            value is not None and str(value) in [str(k) for k in keys]
        ) else None
        with self.params_col:
            with ui.element("div").classes("br-prow"):
                label_with_help(key).classes("br-prow-label")
                sel = ui.select(options, value=initial).props(
                    "filled dense"
                ).classes("br-pfield")
                sel.on_value_change(lambda e, k=key: self._update_param(k, e.value))

    def _param_file(self, key, value) -> None:
        with self.params_col:
            with ui.element("div").classes("br-prow"):
                label_with_help(key).classes("br-prow-label")
                with ui.row().classes("br-pfield items-center gap-1 no-wrap"):
                    inp = ui.input(
                        value="" if value in _UNSET else str(value),
                    ).props("filled dense").classes("grow min-w-0")
                    inp.on_value_change(lambda e, k=key: self._update_param(k, e.value))

                    async def browse(k=key, field=inp) -> None:
                        if LocalFilePicker.active() is not None:
                            return  # don't stack a second picker
                        # Per-field folder memory: pulse waveforms, sequence
                        # JSONs, etc. live in different places.
                        state_key = f"last_dir_{k.lower().replace(' ', '_')}"
                        cur = field.value or ""
                        start = (os.path.dirname(cur)
                                 if cur and os.path.isdir(os.path.dirname(cur))
                                 else get_state(state_key) or "~")
                        if not isinstance(start, str) or (
                                start != "~" and not os.path.isdir(start)):
                            start = "~"
                        path = await LocalFilePicker(start, title=f"Select {k}")
                        if path:
                            set_state(state_key, os.path.dirname(path))
                            field.value = path
                            self._update_param(k, path)

                    ui.button(icon="folder_open", on_click=browse).props(
                        "flat dense round color=primary"
                    )

    def _param_metabolites(self) -> None:
        backend = self.BasisREMY.backend
        selected = backend.mandatory_params.get("Metabolites", [])
        with self.metabs_col:
            with ui.row().classes("items-center justify-between w-full"):
                ui.label("Metabolites").classes("br-section-title")
                with ui.row().classes("gap-0"):
                    ui.button("Default",
                              on_click=lambda: self._set_metabs(None)).props(
                        "flat dense color=primary").classes("text-xs")
                    ui.button("None",
                              on_click=lambda: self._set_metabs(False)).props(
                        "flat dense color=primary").classes("text-xs")
                    ui.button("All",
                              on_click=lambda: self._set_metabs(True)).props(
                        "flat dense color=primary").classes("text-xs")
            self.metab_checks = {}
            with ui.element("div").classes("br-card br-metab w-full px-3 py-1.5"):
                with ui.grid(columns=2).classes("gap-x-3 gap-y-0 w-full"):
                    for metab in backend.metabs:
                        cb = ui.checkbox(metab, value=metab in selected).props(
                            "dense"
                        ).classes("text-sm")
                        full = metabolite_full_name(metab)
                        if full:
                            cb.tooltip(full)
                        cb.on_value_change(self._update_metabs)
                        self.metab_checks[metab] = cb

    def _set_metabs(self, value: 'bool | None') -> None:
        """Set every metabolite checkbox: True=all, False=none, None=defaults."""
        defaults = getattr(self.BasisREMY.backend, '_default_metabs', {})
        for m, cb in self.metab_checks.items():
            cb.value = defaults.get(m, False) if value is None else value
        self._update_metabs()

    def _update_metabs(self, _event=None) -> None:
        backend = self.BasisREMY.backend
        selected = []
        for m, cb in self.metab_checks.items():
            # keep backend.metabs in sync — set_backend carries the selection
            # over to the next backend from there, not from mandatory_params
            if m in backend.metabs:
                backend.metabs[m] = bool(cb.value)
            if cb.value:
                selected.append(m)
        backend.mandatory_params["Metabolites"] = selected
        self.validate_inputs()

    # ---- validation -------------------------------------------------------
    def validate_inputs(self) -> None:
        if self.simulate_button is None:
            return
        backend = self.BasisREMY.backend

        def _ok(key, value) -> bool:
            if value in _UNSET:
                return False
            # Registry units mark numeric parameters — reject text like "abc"
            # in TE/Bandwidth here instead of crashing mid-simulation.
            if (key not in backend.dropdown
                    and key not in backend.file_selection
                    and registry_get(key).units):
                try:
                    float(value)
                except (TypeError, ValueError):
                    return False
            return True

        # Validate what the current mode actually shows — hidden fields (e.g.
        # Sequence in FSL-MRS Template mode) must not block Simulate.
        visible = backend.get_params_for_mode()
        all_filled = all(
            _ok(key, value)
            for key, value in visible.items()
            if key != "Metabolites"
        )
        if self.metab_checks:
            at_least_one = any(cb.value for cb in self.metab_checks.values())
        else:
            at_least_one = True
        if all_filled and at_least_one:
            self.simulate_button.enable()
        else:
            self.simulate_button.disable()

    # ============================================================== Octave check
    def _show_octave_instructions(self, manager) -> None:
        instructions = manager._get_installation_instructions()
        dialog = ui.dialog()
        with dialog, ui.card().classes(
            "w-[680px] max-w-full gap-3 p-6 rounded-2xl"
        ).style("box-shadow:0 24px 60px -30px rgba(16,24,40,.5)"):
            with ui.row().classes("items-center gap-2 no-wrap"):
                ui.icon("info").classes("text-2xl").style("color:var(--br-primary)")
                ui.label("Octave runtime required").classes(
                    "text-lg font-bold text-grey-9"
                )
            with ui.scroll_area().classes("w-full h-96").style(
                "border:1px solid var(--br-line);border-radius:12px;"
            ):
                ui.label(instructions).classes(
                    "text-xs whitespace-pre font-mono p-3"
                )
            with ui.row().classes("w-full justify-end"):
                ui.button("OK", on_click=dialog.close).props("color=primary unelevated")
        dialog.open()

    # ============================================================== TAB 3
    def _build_tab3_progress(self) -> None:
        self.panel3.clear()
        with self.panel3:
            self.tab3_container = ui.column().classes("w-full items-center gap-6 py-2")
            with self.tab3_container:
                with ui.column().classes(
                    "items-center gap-4 py-12 w-full max-w-sm"
                ) as self._progress_box:
                    ui.spinner("dots", size="lg").style("color:var(--br-primary)")
                    self.sim_status = ui.label("Simulating basis set…").classes(
                        "text-base font-semibold text-grey-9"
                    )
                    self.progress = ui.linear_progress(
                        value=0, show_value=False
                    ).classes("w-full").props("instant-feedback rounded size=8px")
                    self.progress_label = ui.label("0%").classes(
                        "text-sm text-grey-6"
                    )
                self.results_container = ui.column().classes("w-full gap-4")
            ui.element("div").classes("br-hairline")
            with ui.row().classes("w-full justify-start"):
                ui.button("Back", icon="arrow_back",
                          on_click=lambda: self._goto("params")).props("flat color=primary")

    async def _simulate_basis(self) -> None:
        # A cancelled run only stops between metabolites; starting a second
        # thread meanwhile would revive it via _sim_stop_event.clear() and
        # race on the shared progress/result state.
        if self._sim_thread is not None and self._sim_thread.is_alive():
            ui.notify("A simulation is still running — please wait a moment.",
                      type="warning")
            return
        if getattr(self, "_sim_launching", False):
            return  # a second click during the availability probe await
        self._sim_launching = True

        backend = self.BasisREMY.backend
        try:
            if backend.requires_octave and backend.octave is None:
                from basisremy.core.octave_manager import OctaveManager
                manager = OctaveManager()
                # Docker/Octave probes take seconds — don't freeze the click.
                available = await run.io_bound(
                    lambda: manager.check_docker_availability()
                    or manager.check_local_octave_availability())
                if not available:
                    self._show_octave_instructions(manager)
                    return
            if self.BasisREMY.backend is not backend:
                return  # backend switched while the probe was running
        finally:
            self._sim_launching = False

        # reset tab3 to a clean progress state
        self._basis_set_valid = False
        self.basis_set = None
        self._build_tab3_progress()
        self._progress_box.set_visibility(True)

        self._unlock("sim")
        self._goto("sim")

        print("Simulating basis set with the following parameters:")
        for key, value in backend.mandatory_params.items():
            print(f"{key}: {value}")

        metabs = backend.mandatory_params.get("Metabolites", [])
        self._sim_step = 0
        self._sim_total = max(1, len(metabs))
        self._sim_done = False
        self._sim_error = None
        self._sim_failures = {}
        self.progress.set_value(0)
        self.progress_label.set_text("0%")

        self._sim_stop_event.clear()
        self._sim_thread = threading.Thread(target=self._run_simulation, daemon=True)
        self._sim_thread.start()

        self._sim_timer = ui.timer(0.1, self._poll_simulation)

    def _run_simulation(self) -> None:
        backend = self.BasisREMY.backend

        def progress_callback(step, total_steps):
            self._sim_step = step
            self._sim_total = max(1, total_steps)

        try:
            # Mandatory params win on key clashes; optional params (TM,
            # Custom Sequence, Template File, …) must reach the backend too.
            basis = backend.run_simulation(
                {**backend.optional_params, **backend.mandatory_params},
                progress_callback,
                stop_event=self._sim_stop_event,
            )
        except Exception as exc:  # noqa: BLE001
            if self._sim_stop_event.is_set():
                print("⏹  Simulation cancelled.")
            else:
                self._sim_error = exc
            self._sim_done = True
            return

        if self._sim_stop_event.is_set():
            print("⏹  Simulation cancelled.")
            self._sim_done = True
            return

        self._sim_failures = dict(getattr(backend, 'last_failures', None) or {})
        if self._sim_failures and not basis:
            self._sim_error = RuntimeError(
                "all metabolites failed: "
                + ", ".join(sorted(self._sim_failures)))
            self._sim_done = True
            return

        self.basis_set = basis
        self._sim_done = True

    def _poll_simulation(self) -> None:
        # live progress
        frac = self._sim_step / self._sim_total if self._sim_total else 0
        self.progress.set_value(frac)
        self.progress_label.set_text(f"{int(frac * 100)}%")

        if not self._sim_done:
            return

        # finished — stop polling
        if self._sim_timer is not None:
            self._sim_timer.deactivate()
            self._sim_timer = None

        if self._sim_stop_event.is_set():
            return

        if self._sim_error is not None:
            ui.notify(f"Simulation failed: {self._sim_error}", type="negative")
            self._progress_box.set_visibility(False)
            with self.results_container:
                with ui.row().classes("items-center gap-2"):
                    ui.icon("error").classes("text-2xl").style("color:#c2453c")
                    ui.label("Simulation failed — see the console for details.").classes(
                        "text-sm font-semibold text-grey-9"
                    )
            print(f"Simulation error: {self._sim_error}")
            return

        print("Simulation complete.")
        self._basis_set_valid = True
        self._render_results()

    def _render_results(self) -> None:
        # hide the progress widgets
        self._progress_box.set_visibility(False)

        self.checkbox_vars = {}
        self.metab_colors = {}
        default_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

        with self.results_container:
            with ui.row().classes("items-center gap-2 self-start"):
                ui.icon("check_circle").classes("text-xl").style("color:#3f8f5b")
                ui.label("Basis set ready").classes(
                    "text-base font-bold text-grey-9"
                )

            # Edited sequences return '<metab> (ON/OFF/DIFF)' entries — show
            # exactly one sub-spectrum at a time, switched with a toggle.
            self._subspec_toggle = None
            self._legend_rows = {}
            subs = [s for s in ("DIFF", "ON", "OFF")
                    if any(k.endswith(f" ({s})") for k in self.basis_set)]
            if subs:
                with ui.row().classes("items-center gap-4 self-start"):
                    ui.label("Sub-spectrum:").classes(
                        "text-sm font-semibold text-grey-8")
                    self._subspec_toggle = ui.toggle(
                        subs, value=subs[0]).props("dense")
                    self._subspec_toggle.on_value_change(
                        self._on_subspec_switch)

            failures = getattr(self, "_sim_failures", None) or {}
            if failures:
                names = ", ".join(sorted(failures))
                ui.notify(f"{len(failures)} metabolite(s) failed to simulate: "
                          f"{names}", type="warning")
                with ui.row().classes("items-center gap-2 self-start"):
                    ui.icon("warning").classes("text-lg").style("color:#c2892e")
                    ui.label(f"Not simulated (see console): {names}").classes(
                        "text-sm text-grey-8"
                    )

            with ui.row().classes("w-full no-wrap gap-6 items-start"):
                # plot on the left
                with ui.column().classes("grow min-w-0"):
                    self.plot = ui.matplotlib(figsize=(7, 3)).classes("w-full")
                    fig = self.plot.figure
                    fig.patch.set_alpha(0.0)
                    self.ax = fig.add_subplot(111)

                # legend / metabolite toggles on the right
                with ui.column().classes(
                    "br-legend gap-0 max-h-72 overflow-auto pl-4"
                ):
                    for i, metab in enumerate(self.basis_set.keys()):
                        color = default_colors[i % len(default_colors)]
                        self.metab_colors[metab] = color
                        row = ui.row().classes("items-center gap-2 no-wrap")
                        with row:
                            ui.element("div").style(
                                f"width:11px;height:11px;border-radius:3px;"
                                f"background:{color};"
                            )
                            cb = ui.checkbox(metab, value=True).props("dense")
                            cb.on_value_change(self._update_plot)
                            self.checkbox_vars[metab] = cb
                        self._legend_rows[metab] = row
                        row.set_visibility(self._subspec_visible(metab))

            ui.button("Export basis…", icon="download",
                      on_click=self._open_export_dialog).props(
                "color=primary unelevated"
            )

        self._update_plot()

    def _subspec_visible(self, name: str) -> bool:
        """Only the toggled sub-spectrum is shown; plain entries always are."""
        toggle = getattr(self, "_subspec_toggle", None)
        if toggle is None:
            return True
        for sub in ("DIFF", "ON", "OFF"):
            if name.endswith(f" ({sub})"):
                return sub == toggle.value
        return True

    def _on_subspec_switch(self, _event=None) -> None:
        for name, row in getattr(self, "_legend_rows", {}).items():
            row.set_visibility(self._subspec_visible(name))
        self._update_plot()

    def _open_export_dialog(self) -> None:
        if not self.basis_set:
            return
        open_export_dialog(self.basis_set, self.BasisREMY.backend.mandatory_params)

    def _update_plot(self, _event=None) -> None:
        axis_col = "#8a95a3"  # neutral grey that reads on light and dark
        self.ax.clear()
        self.ax.set_facecolor("none")
        self.ax.set_xlabel("Chemical shift [ppm]", color=axis_col)
        self.ax.set_yticks([])
        self.ax.set_yticklabels([])
        for side in ("top", "right", "left"):
            self.ax.spines[side].set_visible(False)
        self.ax.spines["bottom"].set_color(axis_col)
        self.ax.tick_params(colors=axis_col)

        mp = self.BasisREMY.backend.mandatory_params
        cf_raw = mp.get("Center Freq")
        if cf_raw in (None, "", "missing input"):
            b0_raw = mp.get("Bfield") or \
                str(mp.get("Field Strength") or "3T").replace("T", "").strip()
            try:
                b0 = float(b0_raw)
            except (TypeError, ValueError):
                b0 = 3.0
            cf = 42.577e6 * b0
        else:
            cf = float(cf_raw)
            if cf < 20:
                # ppm rotating-frame centre (FID-A shaped backends), not a
                # Larmor frequency — derive the Larmor frequency from B0
                try:
                    cf = 42.577e6 * float(mp.get("Bfield") or 3.0)
                except (TypeError, ValueError):
                    cf = 42.577e6 * 3.0
            elif cf < 1000:
                cf *= 1e6  # MHz → Hz
        bw = float(mp["Bandwidth"])

        for metab, cb in self.checkbox_vars.items():
            if (cb.value and metab in self.basis_set
                    and self._subspec_visible(metab)):
                data = self.basis_set[metab]
                if not isinstance(data, np.ndarray):
                    try:
                        data = np.array(data, dtype=complex)
                    except Exception:  # noqa: BLE001
                        continue
                if data.ndim > 1:
                    data = data.flatten()
                if data.size == 0:
                    continue
                try:
                    ydata = np.real(np.fft.fftshift(np.fft.fft(data)))
                    npts = data.size
                    ppm_axis = np.linspace(-bw / 2, bw / 2, npts) / cf * 1e6 + 4.65
                    self.ax.plot(ppm_axis, ydata, color=self.metab_colors[metab])
                except Exception as e:  # noqa: BLE001
                    print(f"Warning: Could not plot {metab}: {e}")
                    continue

        self.ax.set_xlim(10, 0)
        try:
            self.ax.figure.tight_layout(pad=0.4)
        except Exception:  # noqa: BLE001
            pass
        self.plot.update()


#**************************************************************************************************#
#                                          entry points                                            #
#**************************************************************************************************#
def build_page() -> None:
    """Construct the single-page UI for one client connection."""
    ui.colors(
        primary=PRIMARY,
        secondary=PRIMARY_DARK,
        accent=PRIMARY,
        dark="#0e1822",
        positive="#3f8f5b",
        negative="#c2453c",
        info=PRIMARY,
        warning="#c98a2b",
    )
    ui.add_head_html(_GLOBAL_CSS)
    BasisREMYApp()


def run_app(*, native: bool = True, show: bool = True, port: int = 8080) -> None:
    """Configure routes and start the NiceGUI server.

    In native mode a desktop window is opened via pywebview; otherwise the UI
    is served in the default browser.
    """
    app.add_static_files("/assets", str(ASSETS))
    ui.page("/")(build_page)
    ui.run(
        native=native,
        title="BasisREMY",
        window_size=(1100, 820) if native else None,
        port=port,
        show=show,
        reload=False,
        storage_secret="basisremy",
    )
