####################################################################################################
#                                          help_widget.py                                          #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Purpose: A small circled "?" help icon and a `LabelWithHelp` convenience that places the         #
#          icon directly next to a parameter label. Hovering (or clicking) the icon shows a        #
#          richly-styled tooltip pulled from core.parameter_registry.                              #
#                                                                                                  #
#          macOS-safe: Tk's wm_overrideredirect tooltips don't always render on macOS, so we       #
#          additionally call the (unsupported) MacWindowStyle "help / noActivates" hint and lift   #
#          the toplevel after the geometry is realised.                                            #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import platform
import tkinter as tk

from core.parameter_registry import tooltip_text, get as get_param, TODO


# Visual tokens (kept in-sync with gui/application.py palette).
_BG          = "#f0f0f0"
_FG          = "#000000"
_ACCENT      = "#607389"
_TODO_ACCENT = "#c25450"   # subtle red tint when description is still a placeholder
_TIP_BG      = "#ffffe0"
_TIP_FG      = "#1a1a1a"

_IS_MAC = platform.system() == "Darwin"


class HelpIcon(tk.Canvas):
    """A small circled "?" icon.

    Hover (or click) opens a styled tooltip with the parameter's full
    description, units, typical range and source — sourced from
    `core.parameter_registry`.
    """

    def __init__(self, master, param_key: str, *, bg: str = _BG, size: int = 16, **kw):
        info = get_param(param_key)
        is_todo = TODO in info.description
        self._fg = _TODO_ACCENT if is_todo else _ACCENT

        super().__init__(
            master,
            width=size, height=size,
            bg=bg, bd=0, highlightthickness=0,
            cursor="question_arrow",
            **kw,
        )
        self.param_key = param_key
        self._size = size
        self._tip = None
        self._draw(filled=False)

        self.bind("<Enter>", self._on_enter)
        self.bind("<Leave>", self._on_leave)
        self.bind("<Button-1>", self._toggle)

    # ------------------------------------------------------------------ paint
    def _draw(self, filled: bool):
        self.delete("all")
        s = self._size
        pad = 1
        outline = self._fg
        fill = self._fg if filled else ""
        text_color = "#ffffff" if filled else self._fg
        self.create_oval(pad, pad, s - pad, s - pad,
                         outline=outline, fill=fill, width=1)
        # Centre the glyph; macOS metrics differ slightly from X11/Win.
        y_off = 0 if _IS_MAC else -1
        self.create_text(s / 2, s / 2 + y_off, text="?",
                         fill=text_color,
                         font=("Helvetica", max(8, int(s * 0.62)), "bold"))

    # ----------------------------------------------------------------- events
    def _on_enter(self, _evt=None):
        self._draw(filled=True)
        self._show()

    def _on_leave(self, _evt=None):
        self._draw(filled=False)
        self._hide()

    def _toggle(self, _evt=None):
        if self._tip is None:
            self._show()
        else:
            self._hide()

    # -------------------------------------------------------------- lifecycle
    def _show(self):
        if self._tip is not None:
            return
        text = tooltip_text(self.param_key)
        info = get_param(self.param_key)
        is_todo = TODO in info.description

        # Position the tooltip just below-right of the icon.
        x = self.winfo_rootx() + self._size + 6
        y = self.winfo_rooty() + self._size + 4

        tip = tk.Toplevel(self)
        tip.wm_overrideredirect(True)
        try:
            tip.wm_attributes("-topmost", True)
        except tk.TclError:
            pass

        # macOS-specific: ensure the borderless window actually shows.
        if _IS_MAC:
            try:
                tip.tk.call("::tk::unsupported::MacWindowStyle",
                            "style", tip._w, "help", "noActivates")
            except tk.TclError:
                pass

        outer = tk.Frame(tip, bg=_ACCENT, bd=0)
        outer.pack(fill=tk.BOTH, expand=True)

        inner = tk.Frame(outer, bg=_TIP_BG, bd=0)
        inner.pack(fill=tk.BOTH, expand=True, padx=1, pady=1)

        title_color = _TODO_ACCENT if is_todo else _ACCENT
        title = tk.Label(
            inner,
            text=(info.label or self.param_key) + ("  (TODO)" if is_todo else ""),
            bg=_TIP_BG, fg=title_color,
            font=("Arial", 10, "bold"),
            anchor="w", justify="left",
        )
        title.pack(fill=tk.X, padx=8, pady=(6, 0))

        body = tk.Label(
            inner,
            text=text,
            bg=_TIP_BG, fg=_TIP_FG,
            font=("Arial", 9),
            anchor="w", justify="left",
            wraplength=340,
        )
        body.pack(fill=tk.X, padx=8, pady=(2, 6))

        # Position + force visibility (macOS sometimes needs a nudge).
        tip.wm_geometry(f"+{x}+{y}")
        tip.update_idletasks()
        try:
            tip.deiconify()
            tip.lift()
        except tk.TclError:
            pass

        self._tip = tip

    def _hide(self, _event=None):
        if self._tip is not None:
            try:
                self._tip.destroy()
            except tk.TclError:
                pass
            self._tip = None


class LabelWithHelp(tk.Frame):
    """Label + circled "?" packed side-by-side.

    Used in the parameter grid so the help icon sits *next to the parameter
    name*, not at the end of the row.
    """

    def __init__(self, master, param_key: str, *,
                 text: str | None = None,
                 font=("Arial", 12, "bold"),
                 bg: str = _BG,
                 fg: str = _FG,
                 icon_size: int = 15,
                 **kw):
        super().__init__(master, bg=bg, **kw)
        label_text = text if text is not None else f"{param_key}:"
        self.label = tk.Label(self, text=label_text, font=font, bg=bg, fg=fg)
        self.label.pack(side=tk.LEFT)
        self.icon = HelpIcon(self, param_key, bg=bg, size=icon_size)
        self.icon.pack(side=tk.LEFT, padx=(4, 0))


