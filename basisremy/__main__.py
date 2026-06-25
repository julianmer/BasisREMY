####################################################################################################
#                                           __main__.py                                            #
####################################################################################################
#                                                                                                  #
# Purpose: Console entry point for the BasisREMY application (the ``basisremy`` command).          #
#                                                                                                  #
#          This thin launcher resolves the project root, makes the vendored ``externals/`` code    #
#          importable, and anchors the working directory at the project root so the Octave/Docker   #
#          backends that use relative ``./externals`` and ``./output`` paths keep working exactly   #
#          as before.                                                                               #
#                                                                                                  #
#          A future NiceGUI front-end can replace ``_start_tk_app`` behind this same entry point   #
#          without changing the user-facing ``basisremy`` command.                                 #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def _find_project_root() -> Path:
    """Return the runtime root (the directory that contains ``externals/``).

    For a source checkout this is the repo root; for an installed copy it is a
    per-user data directory (see :func:`basisremy.core.paths.runtime_root`).
    """
    from basisremy.core.paths import runtime_root

    return runtime_root()


def _prepare_runtime() -> Path:
    """Anchor the working directory and materialize the bundled adapters.

    Returns the resolved runtime root.
    """
    from basisremy.core.paths import ensure_adapters, runtime_root

    root = runtime_root()
    root.mkdir(parents=True, exist_ok=True)

    # The Octave backends reference ``./externals`` and ``./adapters`` relative
    # to the working directory, and the Docker runtime mounts that directory
    # into the container. So make the bundled adapters reachable under the
    # runtime root and anchor the working directory there. (externals/ is
    # fetched lazily per-backend; see basisremy.core.externals.)
    try:
        ensure_adapters()
    except OSError:
        pass

    try:
        os.chdir(root)
    except OSError:
        pass

    return root


def _load_application():
    """Import and return the Tkinter ``Application`` class.

    Kept separate from :func:`main` so smoke tests can verify the import/path
    wiring without opening a window, and so a future NiceGUI front-end can be
    slotted in here behind the same entry point.
    """
    from basisremy.gui.application import Application

    return Application


def _start_tk_app() -> None:
    """Instantiate and run the Tkinter application main loop."""
    Application = _load_application()
    app = Application()
    app.mainloop()


def _run_environment_check() -> int:
    """Print a friendly report of the runtime environment and exit.

    Verifies the package, GUI toolkit and Octave runtimes (Docker / local)
    without launching the GUI. Handy for users and as a non-interactive smoke
    test (e.g. ``uvx --from . basisremy --check``). Always returns 0 so it can
    be wired into CI; problems are reported as human-readable lines.
    """
    from basisremy import __version__

    def _mark(ok: bool) -> str:
        return "OK   " if ok else "MISSING"

    print("BasisREMY environment check")
    print("=" * 40)
    print(f"  package version : {__version__}")
    print(f"  python          : {sys.version.split()[0]} ({sys.executable})")

    # GUI toolkit -----------------------------------------------------------
    try:
        import tkinter  # noqa: F401
        tk_ok = True
    except Exception:
        tk_ok = False
    print(f"  tkinter (GUI)   : {_mark(tk_ok)}")

    # Octave runtimes -------------------------------------------------------
    docker_ok = local_ok = False
    try:
        from basisremy.core.octave_manager import OctaveManager

        manager = OctaveManager()
        docker_ok = bool(manager.check_docker_availability())
        local_ok = bool(manager.check_local_octave_availability())
    except Exception as exc:  # noqa: BLE001
        print(f"  octave check    : could not run ({exc})")
    print(f"  docker octave   : {_mark(docker_ok)}")
    print(f"  local octave    : {_mark(local_ok)}")
    print("=" * 40)

    if not tk_ok:
        print("Note: the GUI needs Tkinter. See the README 'Troubleshooting' section.")
    if not (docker_ok or local_ok):
        print("Note: simulation backends need Docker OR a local Octave install.")
        print("      Data extraction and parameter configuration still work without them.")
    else:
        print("Ready: an Octave runtime is available for simulation backends.")
    return 0


def _build_parser() -> argparse.ArgumentParser:
    from basisremy import __version__

    parser = argparse.ArgumentParser(
        prog="basisremy",
        description="Generate study-specific MR spectroscopy basis sets from raw MRS data.",
    )
    parser.add_argument(
        "--version", action="version", version=f"basisremy {__version__}",
    )
    parser.add_argument(
        "--check", action="store_true",
        help="Report the environment (Python, Tkinter, Docker/Octave) and exit, "
             "without launching the GUI.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Console-script entry point for the ``basisremy`` command."""
    args = _build_parser().parse_args(argv)

    _prepare_runtime()

    if args.check:
        return _run_environment_check()

    try:
        _start_tk_app()
    except ModuleNotFoundError as exc:  # most commonly missing Tkinter
        missing = getattr(exc, "name", "") or str(exc)
        if missing in ("tkinter", "_tkinter") or "tkinter" in str(exc).lower():
            sys.stderr.write(
                "\nERROR: Tkinter is not available in this Python installation.\n"
                "Tkinter ships with CPython but some distributions split it out.\n\n"
                "  - macOS (Homebrew):   brew install python-tk\n"
                "  - Debian/Ubuntu:      sudo apt install python3-tk\n"
                "  - Fedora:             sudo dnf install python3-tkinter\n"
                "  - Windows:            reinstall Python from python.org with the\n"
                "                        'tcl/tk and IDLE' option enabled\n\n"
            )
            return 1
        sys.stderr.write(
            f"\nERROR: A required module could not be imported: {missing}\n"
            "Try reinstalling dependencies (e.g. 'uv sync' in the project).\n\n"
        )
        return 1
    except Exception as exc:  # noqa: BLE001 - surface a friendly startup error
        sys.stderr.write(f"\nERROR: BasisREMY failed to start: {exc}\n\n")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
