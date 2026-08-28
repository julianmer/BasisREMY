####################################################################################################
#                                          export_dialog.py                                        #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 26/06/26                                                                                #
#                                                                                                  #
# Purpose: NiceGUI dialog for exporting a generated basis set in any of the formats                #
#          supported by core.exporters (LCModel .basis, LCModel .RAW per metab,                    #
#          jMRUI .txt, FSL-MRS .json directory, Osprey .mat).                                      #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import re
import traceback

from nicegui import run, ui

from basisremy.core.exporters import (
    export as _export,
    export_subspectra as _export_subspectra,
    subspectrum_path,
    sequence_label,
    SUPPORTED_FORMATS,
    FORMAT_LABELS,
    FORMAT_EXTENSIONS,
)
from basisremy.gui.help_widget import help_icon
from basisremy.gui.local_file_picker import LocalFilePicker
from basisremy.gui.ui_state import get_state, set_state


_ACCENT = "var(--br-primary)"


def _default_basename(params: dict) -> str:
    """Descriptive default file name, e.g. 'PRESS_TE35_3T' (Osprey-style)."""
    parts = []
    seq = sequence_label(params)
    if seq:
        parts.append(seq)
    te = params.get("TE")
    try:
        parts.append(f"TE{float(te):g}")
    except (TypeError, ValueError):
        pass
    b0 = params.get("Bfield") or str(params.get("Field Strength") or "").replace("T", "")
    try:
        parts.append(f"{float(b0):g}T")
    except (TypeError, ValueError):
        pass
    name = "_".join(parts) or "basis"
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def open_export_dialog(basis: dict, params: dict) -> None:
    """Open the unified export dialog (LCModel / jMRUI / FSL-MRS / Osprey)."""
    if not basis:
        ui.notify("Nothing to export yet.", type="warning")
        return

    params = dict(params or {})
    remembered = get_state("last_export_dir")
    default_dir = (remembered
                   if isinstance(remembered, str) and os.path.isdir(remembered)
                   else os.path.abspath("./output"))

    dialog = ui.dialog()
    with dialog, ui.card().classes("w-[600px] max-w-full gap-3"):
        ui.label("Export basis set").classes("text-lg font-bold").style(
            f"color:{_ACCENT}"
        )

        # ---- format ----
        with ui.row().classes("items-center w-full no-wrap gap-2"):
            fmt_select = ui.select(
                {k: FORMAT_LABELS[k] for k in SUPPORTED_FORMATS},
                value="lcmodel_basis",
                label="Format",
            ).classes("grow")
            help_icon("Output Format")

        # ---- sub-spectra (edited sequences only) ----
        has_subspectra = any(
            k.endswith((' (ON)', ' (OFF)', ' (DIFF)')) for k in basis)
        subspec_select = None
        if has_subspectra:
            subspec_select = ui.select(
                ['All', 'DIFF only', 'ON only', 'OFF only'],
                value='All', label='Sub-spectra',
            ).classes("w-full")

        # ---- output path ----
        path_input = ui.input(
            "Output path",
            value=os.path.join(default_dir,
                               _default_basename(params)
                               + FORMAT_EXTENSIONS["lcmodel_basis"]),
        ).classes("w-full")

        async def browse() -> None:
            fmt = fmt_select.value
            ext = FORMAT_EXTENSIONS[fmt]
            cur = path_input.value or default_dir
            if ext:  # single file -> save picker
                start_dir = os.path.dirname(cur) or default_dir
                chosen = await LocalFilePicker(
                    start_dir,
                    save_mode=True,
                    default_name=os.path.basename(cur) or ("basis" + ext),
                    title="Save basis as…",
                )
            else:    # directory output
                chosen = await LocalFilePicker(
                    cur if os.path.isdir(cur) else default_dir,
                    dirs_only=True,
                    title="Select output directory",
                )
            if chosen:
                path_input.value = chosen
                set_state("last_export_dir",
                          chosen if os.path.isdir(chosen)
                          else os.path.dirname(chosen))

        with ui.row().classes("w-full justify-end"):
            ui.button("Browse…", icon="folder_open", on_click=browse).props("flat")

        # keep the suggested filename extension in sync with the format
        def on_fmt_change() -> None:
            fmt = fmt_select.value
            ext = FORMAT_EXTENSIONS[fmt]
            cur = path_input.value or os.path.join(default_dir, "basis")
            base = os.path.splitext(cur)[0] if "." in os.path.basename(cur) else cur
            path_input.value = base + ext if ext else base
        fmt_select.on_value_change(on_fmt_change)

        # ---- summary ----
        n_pts = next(iter(basis.values())).size if basis else 0

        def _summary(selection: str = "All") -> str:
            if subspec_select is None:
                n = f"{len(basis)} metabolites"
            elif selection == "All":
                n = (f"{len(basis)} entries, one basis set per sub-spectrum "
                     f"(_ON / _OFF / _DIFF)")
            else:
                sub = selection.split()[0]
                n = (f"{sum(k.endswith(f'({sub})') for k in basis)} "
                     f"metabolites ({sub} sub-spectrum)")
            return (f"{n} · {n_pts} points · "
                    f"BW = {params.get('Bandwidth', '?')} Hz · "
                    f"TE = {params.get('TE', '?')} ms · "
                    f"Sequence = {sequence_label(params) or '?'}")

        info_label = ui.label(_summary()).classes("text-xs italic text-grey-7")
        if subspec_select is not None:
            subspec_select.on_value_change(
                lambda e: info_label.set_text(_summary(e.value)))

        status = ui.label("").classes("text-sm")

        async def confirm_overwrite(p: str) -> bool:
            with ui.dialog() as confirm, ui.card().classes("gap-3"):
                ui.label(f"“{os.path.basename(p)}” already exists. Overwrite?"
                         ).classes("text-sm")
                with ui.row().classes("w-full justify-end gap-2"):
                    ui.button("Cancel",
                              on_click=lambda: confirm.submit(False)).props("flat")
                    ui.button("Overwrite",
                              on_click=lambda: confirm.submit(True)
                              ).props("color=primary")
            result = bool(await confirm)  # Esc/outside-click → None → False
            confirm.clear()  # closed dialogs are only hidden, not removed
            return result

        async def do_export() -> None:
            if getattr(do_export, "_busy", False):
                return  # second click during a slow write
            fmt = fmt_select.value
            path = (path_input.value or "").strip()
            if not path:
                ui.notify("Please choose an output path.", type="warning")
                return
            # normpath: a trailing separator must not dodge the overwrite
            # check that export() (which abspath-normalizes) would then hit
            path = os.path.normpath(path)
            ext = FORMAT_EXTENSIONS[fmt]
            final = path if (not ext or path.lower().endswith(ext)) else path + ext
            which = None
            targets = [final]
            if subspec_select is not None:
                # edited basis: one output per sub-spectrum, tagged in the name
                which = subspec_select.value.split()[0]   # All / DIFF / ON / OFF
                subs = ("ON", "OFF", "DIFF") if which == "All" else (which,)
                targets = [subspectrum_path(final, fmt, s) for s in subs]
            existing = [p for p in targets if ext and os.path.isfile(p)]
            if existing and not await confirm_overwrite(existing[0]):
                return
            do_export._busy = True
            export_btn.props("loading")
            try:
                status.set_text(f"Writing {FORMAT_LABELS[fmt]}…")
                status.style("color:#607389")
                # slow formats must not freeze the dialog
                if which is not None:
                    outs = await run.io_bound(_export_subspectra, basis, path,
                                              fmt, params, which)
                    out = outs[0]
                    written = ", ".join(os.path.basename(o) for o in outs)
                else:
                    out = await run.io_bound(_export, basis, path, fmt, params)
                    written = out
                set_state("last_export_dir",
                          out if os.path.isdir(out)
                          else os.path.dirname(os.path.abspath(out)))
                ui.notify(f"Exported {written}", type="positive")
                dialog.close()
            except Exception as exc:  # noqa: BLE001
                traceback.print_exc()
                status.set_text(f"✗ {exc}")
                status.style("color:#c25450")
                ui.notify(f"Export failed: {exc}", type="negative")
            finally:
                do_export._busy = False
                export_btn.props(remove="loading")

        with ui.row().classes("w-full justify-end gap-2"):
            ui.button("Cancel", on_click=dialog.close).props("flat")
            export_btn = ui.button("Export", icon="save",
                                   on_click=do_export).props("color=primary")

    dialog.open()
