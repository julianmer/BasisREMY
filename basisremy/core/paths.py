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


def octave_adapters_base(octave) -> str:
    """Return the path prefix the Octave backends should ``addpath`` for adapters.

    The bundled adapter scripts live inside the package (:data:`ADAPTERS_DIR`)
    and never need to be copied or symlinked next to the working directory:

      * Local Octave runs in-process, so it uses the package's absolute path
        directly.
      * The Docker runtime makes the adapters available through its
        ``/workspace`` bind-mount and exposes the matching ``addpath`` prefix
        (relative to the ``/workspace`` working dir) as the ``ADAPTERS_MOUNT``
        instance attribute.

    The Docker runtime sets ``ADAPTERS_MOUNT`` per instance, so read it from the
    instance ``__dict__`` rather than via ``getattr``: the local Octave object
    (``oct2py.Oct2Py``) forwards unknown attribute access into the Octave
    session, which would otherwise mask the absent attribute.
    """
    mount = octave.__dict__.get("ADAPTERS_MOUNT")
    if not mount:
        mount = getattr(type(octave), "ADAPTERS_MOUNT", None)
    if mount:
        return mount
    return str(ADAPTERS_DIR)
