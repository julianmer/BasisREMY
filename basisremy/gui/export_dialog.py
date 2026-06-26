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
import traceback

from nicegui import ui

from basisremy.core.exporters import (
    export as _export,
    SUPPORTED_FORMATS,
    FORMAT_LABELS,
    FORMAT_EXTENSIONS,
)
from basisremy.gui.help_widget import help_icon
from basisremy.gui.local_file_picker import LocalFilePicker


_ACCENT = "var(--br-primary)"


def open_export_dialog(basis: dict, params: dict) -> None:
    """Open the unified export dialog (LCModel / jMRUI / FSL-MRS / Osprey)."""
    if not basis:
        ui.notify("Nothing to export yet.", type="warning")
        return

    params = dict(params or {})
    label_to_key = {FORMAT_LABELS[k]: k for k in SUPPORTED_FORMATS}
    default_dir = os.path.abspath("./output")

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

        # ---- output path ----
        path_input = ui.input(
            "Output path",
            value=os.path.join(default_dir, "basis" + FORMAT_EXTENSIONS["lcmodel_basis"]),
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
        n_metabs = len(basis)
        n_pts = next(iter(basis.values())).size if basis else 0
        info = (
            f"{n_metabs} metabolites · {n_pts} points · "
            f"BW = {params.get('Bandwidth', '?')} Hz · "
            f"TE = {params.get('TE', '?')} ms · "
            f"Sequence = {params.get('Sequence', '?')}"
        )
        ui.label(info).classes("text-xs italic text-grey-7")

        status = ui.label("").classes("text-sm")

        def do_export() -> None:
            fmt = fmt_select.value
            path = (path_input.value or "").strip()
            if not path:
                ui.notify("Please choose an output path.", type="warning")
                return
            try:
                status.set_text(f"Writing {FORMAT_LABELS[fmt]}…")
                status.style("color:#607389")
                out = _export(basis, path, fmt, params)
                ui.notify(f"Exported to {out}", type="positive")
                dialog.close()
            except Exception as exc:  # noqa: BLE001
                traceback.print_exc()
                status.set_text(f"✗ {exc}")
                status.style("color:#c25450")
                ui.notify(f"Export failed: {exc}", type="negative")

        with ui.row().classes("w-full justify-end gap-2"):
            ui.button("Cancel", on_click=dialog.close).props("flat")
            ui.button("Export", icon="save", on_click=do_export).props("color=primary")

    dialog.open()
