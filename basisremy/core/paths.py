####################################################################################################
#                                            paths.py                                              #
####################################################################################################
#                                                                                                  #
# Purpose: Single source of truth for runtime file-system locations so the backends work both       #
#          from a source checkout and from an installed wheel (pip / uvx), regardless of the       #
#          current working directory.                                                              #
#                                                                                                  #
#          - ADAPTERS_DIR : first-party Octave adapter scripts shipped *inside* the package.        #
#          - runtime_root(): directory that contains ``externals/`` (the third-party submodules).  #
#                            In a source checkout this is the repo root; for an installed copy it  #
#                            is a per-user data directory that ``externals`` are fetched into.     #
#          - externals_root(): ``runtime_root() / "externals"``.                                   #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import shutil
from pathlib import Path

# Directory holding the package's own Octave adapter scripts. These ship inside
# the wheel, so reference them by absolute path rather than relative to cwd.
ADAPTERS_DIR = Path(__file__).resolve().parent.parent / "adapters"

# Package root (the ``basisremy/`` directory) and its parent.
_PACKAGE_DIR = Path(__file__).resolve().parent.parent
_PACKAGE_PARENT = _PACKAGE_DIR.parent


def _looks_like_source_checkout() -> Path | None:
    """Return the repo root if we are running from a source checkout, else None.

    A source checkout has the package's parent directory containing both the
    ``externals`` submodule folder and a ``.git`` directory.
    """
    if (_PACKAGE_PARENT / "externals").is_dir() and (_PACKAGE_PARENT / ".git").exists():
        return _PACKAGE_PARENT
    return None


def runtime_root() -> Path:
    """Return the directory that contains (or will contain) ``externals/``.

    Resolution order:
      1. ``$BASISREMY_HOME`` if set (explicit override).
      2. The repo root when running from a source checkout.
      3. A per-user data directory (``~/.basisremy``) for installed copies.
    """
    env = os.environ.get("BASISREMY_HOME")
    if env:
        return Path(env).expanduser().resolve()

    src = _looks_like_source_checkout()
    if src is not None:
        return src

    return Path.home() / ".basisremy"


def externals_root() -> Path:
    """Return the ``externals`` directory under :func:`runtime_root`."""
    return runtime_root() / "externals"


def ensure_adapters() -> Path:
    """Make the bundled adapter scripts available at ``runtime_root()/adapters``.

    The Octave backends reference ``./adapters/...`` relative to the working
    directory (which is the runtime root and, for Docker, the mounted volume),
    so the package's adapters must be reachable there.

    In a source checkout the adapters already live under the runtime root, so a
    relative symlink is used (single source of truth). For an installed copy the
    adapters sit in ``site-packages`` outside the runtime root, so they are
    copied in.
    """
    root = runtime_root()
    target = root / "adapters"
    if target.exists():
        return target

    root.mkdir(parents=True, exist_ok=True)

    if ADAPTERS_DIR.is_relative_to(root):
        # Source checkout: a relative symlink keeps a single source of truth and
        # works inside the Docker mount (the target is under the same root).
        try:
            target.symlink_to(os.path.relpath(ADAPTERS_DIR, root))
            return target
        except OSError:
            pass  # Fall back to copying if symlinks are not permitted.

    shutil.copytree(ADAPTERS_DIR, target)
    return target
