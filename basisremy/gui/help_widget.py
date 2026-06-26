####################################################################################################
#                                          help_widget.py                                          #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 20/05/26                                                                                #
#                                                                                                  #
# Purpose: A small circled "?" help icon and a `label_with_help` convenience that places the       #
#          icon directly next to a parameter label. Hovering the icon shows a richly-styled        #
#          tooltip pulled from core.parameter_registry.                                            #
#                                                                                                  #
#          NiceGUI re-implementation of the former Tkinter widgets — pure web, works seamlessly    #
#          under uvx (no system tcl/tk).                                                           #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

from nicegui import ui

from basisremy.core.parameter_registry import tooltip_text, get as get_param, TODO


# Visual tokens (kept in-sync with gui/application.py palette).
_TODO_ACCENT = "#c25450"   # subtle red tint when description is still a placeholder


def help_icon(param_key: str):
    """A small circled "?" icon with a rich tooltip for `param_key`.

    The tooltip body comes from :func:`core.parameter_registry.tooltip_text`.
    Parameters that still carry a placeholder description are tinted red.
    """
    info = get_param(param_key)
    is_todo = TODO in info.description
    title_color = _TODO_ACCENT if is_todo else "var(--br-primary)"

    icon = ui.icon("help_outline").classes("cursor-help").style(
        "color:var(--br-muted); font-size:1.05rem;"
    )
    title = (info.label or param_key) + ("  (TODO)" if is_todo else "")
    with icon:
        with ui.tooltip().classes("max-w-sm").style(
            "background: var(--br-surface); color: var(--br-fg);"
            "border: 1px solid var(--br-line);"
            "border-radius: 10px;"
            "box-shadow: 0 10px 30px -12px rgba(0,0,0,.45);"
        ):
            with ui.column().classes("gap-0.5 p-1"):
                ui.label(title).classes("text-sm font-bold").style(
                    f"color:{title_color}"
                )
                ui.label(tooltip_text(param_key)).classes(
                    "text-xs whitespace-pre-line"
                ).style("color: var(--br-fg)")
    return icon


def label_with_help(param_key: str, text: str | None = None):
    """A parameter label with the help icon packed right next to the name."""
    with ui.row().classes("items-center gap-1 no-wrap") as row:
        ui.label(text if text is not None else f"{param_key}:").classes(
            "text-sm font-semibold"
        )
        help_icon(param_key)
    return row
