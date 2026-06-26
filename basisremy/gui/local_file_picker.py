####################################################################################################
#                                       local_file_picker.py                                        #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 26/06/26                                                                                #
#                                                                                                  #
# Purpose: A NiceGUI dialog that browses the *local* filesystem of the machine running the app      #
#          and returns an on-disk path. Because the NiceGUI server runs on the user's own          #
#          machine (localhost / native window), this reads files in place — nothing is uploaded     #
#          or copied, preserving the path-based workflow the Octave/Docker backends rely on.        #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

from pathlib import Path
from typing import Callable

from nicegui import ui


_ACCENT = "#607389"


class LocalFilePicker(ui.dialog):
    """Modal local-filesystem browser.

    Usage::

        path = await LocalFilePicker('~')
        if path:
            ...

    Resolves to the selected path string, or ``None`` if cancelled.

    Args:
        directory:    Starting directory.
        dirs_only:    Only allow selecting a directory (hides files).
        save_mode:    Show a filename field; result is ``dir / filename``.
        default_name: Initial filename in save mode.
        title:        Dialog heading.
        show_file:    Optional predicate; files for which it returns ``False``
                      are shown greyed-out and cannot be selected (folders are
                      always selectable).
    """

    def __init__(
        self,
        directory: str = "~",
        *,
        dirs_only: bool = False,
        save_mode: bool = False,
        default_name: str = "",
        title: str = "Select a file",
        show_file: Callable[[Path], bool] | None = None,
    ) -> None:
        super().__init__()

        start = Path(directory).expanduser()
        if not start.is_dir():
            start = start.parent if start.parent.is_dir() else Path.home()
        self._path = start.resolve()
        self._dirs_only = dirs_only
        self._save_mode = save_mode
        self._show_file = show_file

        with self, ui.card().classes("w-[640px] max-w-full gap-2"):
            ui.label(title).classes("text-lg font-bold").style(f"color:{_ACCENT}")

            # current-path bar with an "up" button
            with ui.row().classes("items-center w-full no-wrap gap-2"):
                ui.button(icon="arrow_upward", on_click=self._go_up) \
                    .props("flat dense round").tooltip("Parent folder")
                ui.button(icon="home", on_click=self._go_home) \
                    .props("flat dense round").tooltip("Home")
                self._path_label = ui.label().classes(
                    "text-xs text-grey-7 truncate grow"
                )

            # scrollable directory listing
            self._list = ui.scroll_area().classes(
                "w-full h-72 border rounded"
            )

            if self._save_mode:
                self._name_input = ui.input(
                    "File name", value=default_name
                ).classes("w-full")

            with ui.row().classes("w-full justify-end gap-2"):
                ui.button("Cancel", on_click=lambda: self.submit(None)) \
                    .props("flat")
                if self._dirs_only:
                    ui.button("Select this folder",
                              on_click=lambda: self.submit(str(self._path))) \
                        .props(f"color=primary")
                elif self._save_mode:
                    ui.button("Save", on_click=self._confirm_save) \
                        .props("color=primary")

        self._render()

    # ------------------------------------------------------------------ nav
    def _go_up(self) -> None:
        if self._path.parent != self._path:
            self._path = self._path.parent
            self._render()

    def _go_home(self) -> None:
        self._path = Path.home()
        self._render()

    def _enter(self, p: Path) -> None:
        self._path = p
        self._render()

    def _confirm_save(self) -> None:
        name = (self._name_input.value or "").strip()
        if not name:
            ui.notify("Please enter a file name", type="warning")
            return
        self.submit(str(self._path / name))

    # --------------------------------------------------------------- listing
    def _render(self) -> None:
        self._path_label.set_text(str(self._path))
        self._list.clear()
        try:
            entries = sorted(
                self._path.iterdir(),
                key=lambda e: (not e.is_dir(), e.name.lower()),
            )
        except (PermissionError, OSError) as exc:
            with self._list:
                ui.label(f"Cannot open folder: {exc}").classes("text-red p-2")
            return

        with self._list:
            with ui.list().props("dense").classes("w-full"):
                for entry in entries:
                    if entry.name.startswith("."):
                        continue
                    is_dir = entry.is_dir()
                    if self._dirs_only and not is_dir:
                        continue
                    enabled = (
                        is_dir
                        or self._show_file is None
                        or self._show_file(entry)
                    )
                    self._row(entry, is_dir, enabled)

    def _row(self, entry: Path, is_dir: bool, enabled: bool = True) -> None:
        def on_click() -> None:
            if is_dir:
                self._enter(entry)
            elif not self._save_mode:
                self.submit(str(entry))

        if not enabled:
            with ui.item().classes("opacity-40 cursor-not-allowed"):
                with ui.item_section().props("avatar"):
                    ui.icon("description").style("color:#999")
                with ui.item_section():
                    ui.item_label(entry.name)
                    ui.item_label("Unsupported format").props("caption")
            return

        with ui.item(on_click=on_click).props("clickable"):
            with ui.item_section().props("avatar"):
                ui.icon("folder" if is_dir else "description").style(
                    f"color:{_ACCENT if is_dir else '#999'}"
                )
            with ui.item_section():
                ui.item_label(entry.name)
