####################################################################################################
#                                          ui_state.py                                             #
####################################################################################################
#                                                                                                  #
# Purpose: Tiny persistent store for GUI conveniences (e.g. last-used directories). State lives    #
#          in ``$BASISREMY_HOME/gui_state.json`` (default ``~/.basisremy/``) so it survives        #
#          restarts on every OS. All failures are swallowed — a read-only home directory must      #
#          never break the GUI.                                                                    #
#                                                                                                  #
####################################################################################################

import json
import os
from pathlib import Path


def _state_file() -> Path:
    env = os.environ.get("BASISREMY_HOME")
    base = Path(env).expanduser() if env else Path.home() / ".basisremy"
    return base / "gui_state.json"


def get_state(key, default=None):
    """Return the stored value for ``key``, or ``default``."""
    try:
        with open(_state_file(), "r") as f:
            data = json.load(f)
        return data.get(key, default) if isinstance(data, dict) else default
    except Exception:
        return default


def set_state(key, value) -> None:
    """Persist ``key: value``; silently a no-op if the store is unusable."""
    path = _state_file()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        try:
            with open(path, "r") as f:
                data = json.load(f)
        except Exception:
            data = {}
        if not isinstance(data, dict):
            data = {}
        data[key] = value
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
    except Exception:
        pass
