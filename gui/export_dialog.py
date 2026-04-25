####################################################################################################
#                                          export_dialog.py                                        #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Purpose: Tk dialog for exporting a generated basis set in any of the formats                     #
#          supported by core.exporters (LCModel .basis, LCModel .RAW per metab,                    #
#          jMRUI .txt, FSL-MRS .json directory, Osprey .mat).                                      #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from core.exporters import (
    export as _export,
    SUPPORTED_FORMATS,
    FORMAT_LABELS,
    FORMAT_EXTENSIONS,
)
from gui.help_widget import HelpIcon


_BG = "#f0f0f0"
_ACCENT = "#607389"
_FG = "#000000"


class ExportDialog(tk.Toplevel):
    def __init__(self, parent, basis: dict, params: dict):
        super().__init__(parent)
        self.basis = basis
        self.params = dict(params or {})

        self.title("Export basis set")
        self.geometry("560x270")
        self.configure(bg=_BG)
        self.transient(parent)
        self.grab_set()

        # ---- format ----
        frm = tk.Frame(self, bg=_BG)
        frm.pack(fill=tk.X, padx=14, pady=(14, 6))
        tk.Label(frm, text="Format:", bg=_BG, fg=_FG, font=("Arial", 12, "bold")).pack(side=tk.LEFT)

        self.fmt_var = tk.StringVar(value=FORMAT_LABELS["lcmodel_basis"])
        labels = [FORMAT_LABELS[k] for k in SUPPORTED_FORMATS]
        self._label_to_key = {FORMAT_LABELS[k]: k for k in SUPPORTED_FORMATS}
        combo = ttk.Combobox(frm, textvariable=self.fmt_var, values=labels,
                             state="readonly", font=("Arial", 12), width=42)
        combo.pack(side=tk.LEFT, padx=(8, 4))
        HelpIcon(frm, "Output Format").pack(side=tk.LEFT)
        combo.bind("<<ComboboxSelected>>", self._on_fmt_change)

        # ---- output path ----
        frm2 = tk.Frame(self, bg=_BG)
        frm2.pack(fill=tk.X, padx=14, pady=(6, 6))
        tk.Label(frm2, text="Output:", bg=_BG, fg=_FG, font=("Arial", 12, "bold")).pack(side=tk.LEFT)

        default_dir = os.path.abspath("./output")
        self.path_var = tk.StringVar(value=os.path.join(default_dir, "basis.basis"))
        entry = tk.Entry(frm2, textvariable=self.path_var, font=("Arial", 11), width=46)
        entry.pack(side=tk.LEFT, padx=(8, 4), fill=tk.X, expand=True)
        tk.Button(frm2, text="Browse…", command=self._browse,
                  bg=_ACCENT, fg=_FG).pack(side=tk.LEFT)
        HelpIcon(frm2, "Output Format").pack(side=tk.LEFT, padx=(4, 0))

        # ---- summary ----
        n_metabs = len(self.basis)
        n_pts = next(iter(self.basis.values())).size if self.basis else 0
        info = (f"{n_metabs} metabolites · {n_pts} points · "
                f"BW = {self.params.get('Bandwidth','?')} Hz · "
                f"TE = {self.params.get('TE','?')} ms · "
                f"Sequence = {self.params.get('Sequence','?')}")
        tk.Label(self, text=info, bg=_BG, fg=_FG, font=("Arial", 10, "italic"),
                 wraplength=520, justify="left").pack(padx=14, pady=(2, 8), anchor="w")

        # ---- buttons ----
        btns = tk.Frame(self, bg=_BG)
        btns.pack(fill=tk.X, padx=14, pady=10)
        tk.Button(btns, text="Cancel", command=self.destroy,
                  bg=_BG, fg=_FG, font=("Arial", 11)).pack(side=tk.RIGHT, padx=(6, 0))
        tk.Button(btns, text="Export", command=self._do_export,
                  bg=_ACCENT, fg=_FG, font=("Arial", 11, "bold")).pack(side=tk.RIGHT)

        self._status = tk.Label(self, text="", bg=_BG, fg=_ACCENT, font=("Arial", 10))
        self._status.pack(padx=14, pady=(0, 4), anchor="w")

    # ---- helpers -----------------------------------------------------------
    def _selected_fmt(self) -> str:
        return self._label_to_key[self.fmt_var.get()]

    def _on_fmt_change(self, _evt=None):
        # Update the suggested filename to match the new format
        fmt = self._selected_fmt()
        ext = FORMAT_EXTENSIONS[fmt]
        cur = self.path_var.get()
        base = os.path.splitext(cur)[0] if "." in os.path.basename(cur) else cur
        if ext:
            self.path_var.set(base + ext)
        else:
            # directory output; strip extension
            self.path_var.set(base)

    def _browse(self):
        fmt = self._selected_fmt()
        ext = FORMAT_EXTENSIONS[fmt]
        if ext:  # single file
            initial = self.path_var.get()
            chosen = filedialog.asksaveasfilename(
                title="Save basis as…",
                defaultextension=ext,
                initialfile=os.path.basename(initial),
                initialdir=os.path.dirname(initial) or ".",
                filetypes=[(FORMAT_LABELS[fmt], f"*{ext}"), ("All files", "*.*")],
            )
        else:    # directory
            chosen = filedialog.askdirectory(title="Select output directory",
                                             initialdir=self.path_var.get() or ".")
        if chosen:
            self.path_var.set(chosen)

    def _do_export(self):
        fmt = self._selected_fmt()
        path = self.path_var.get().strip()
        if not path:
            messagebox.showerror("Export", "Please choose an output path.", parent=self)
            return
        try:
            self._status.config(text=f"Writing {FORMAT_LABELS[fmt]}…", fg=_ACCENT)
            self.update_idletasks()
            out = _export(self.basis, path, fmt, self.params)
            self._status.config(text=f"✓ Exported to {out}", fg="#2c7a3d")
            messagebox.showinfo("Export complete",
                                f"Basis exported as:\n{FORMAT_LABELS[fmt]}\n\n{out}",
                                parent=self)
            self.destroy()
        except Exception as e:
            import traceback
            traceback.print_exc()
            self._status.config(text=f"✗ {e}", fg="#c25450")
            messagebox.showerror("Export failed", str(e), parent=self)

