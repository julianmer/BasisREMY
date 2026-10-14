####################################################################################################
#                                           __main__.py                                            #
####################################################################################################
#                                                                                                  #
# Purpose: Console entry point for the BasisREMY application (the ``basisremy`` command).          #
#                                                                                                  #
#          This thin launcher resolves the project root, makes the vendored ``externals/`` code    #
#          importable, and anchors the working directory at the project root so the Octave/Docker  #
#          backends that use relative ``./externals`` and ``./output`` paths keep working exactly  #
#          as before.                                                                              #
#                                                                                                  #
#          The front-end is NiceGUI: a pure-web UI shown in a native desktop window (pywebview),   #
#          which installs cleanly under uvx with no system tcl/tk dependency. Without a display    #
#          (headless server, SSH login) or with ``--browser`` the UI is served for a web browser.  #
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
    """Anchor the working directory at the runtime root.

    Returns the resolved runtime root.
    """
    from basisremy.core.paths import runtime_root

    root = runtime_root()
    root.mkdir(parents=True, exist_ok=True)

    # The Octave backends reference ``./externals`` relative to the working
    # directory, and the Docker runtime mounts that directory into the
    # container, so anchor the working directory at the runtime root.
    # (externals/ is fetched lazily per-backend; see basisremy.core.externals.
    # The bundled adapters live inside the package and are referenced directly;
    # see basisremy.core.paths.octave_adapters_base.)
    try:
        os.chdir(root)
    except OSError:
        pass

    return root


def _load_application():
    """Import and return the NiceGUI ``BasisREMYApp`` class.

    Kept separate from :func:`main` so smoke tests can verify the import/path
    wiring (basisremy.gui -> core -> backends -> externals) without opening a
    window.
    """
    from basisremy.gui.application import BasisREMYApp

    return BasisREMYApp


def _has_display() -> bool:
    """Return True when a desktop window can be opened on this machine.

    macOS and Windows always have a display; on Linux/BSD a window needs an
    X11 or Wayland session, which headless servers and SSH logins lack.
    """
    if sys.platform in ("darwin", "win32"):
        return True
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


def _browser_banner(port: int, *, auto: bool) -> str:
    """The notice shown once the browser-mode server is up."""
    import getpass
    import socket

    user = getpass.getuser()
    host = socket.gethostname().split(".")[0]
    why = ("No display found, so BasisREMY runs in your web browser instead of a window."
           if auto else "BasisREMY runs in your web browser (--browser).")
    return (
        "\n"
        f"{why}\n"
        "\n"
        f"  On this machine:    open  http://localhost:{port}\n"
        f"  From your own PC:   ssh -L {port}:localhost:{port} {user}@{host}\n"
        f"                      then open  http://localhost:{port}  there\n"
        "  Stop the server:    Ctrl+C\n"
    )


def _start_gui(*, browser: bool = False, port: int = 8080) -> None:
    """Launch the NiceGUI application.

    By default a native desktop window is opened via pywebview. With
    ``browser=True``, or when no display is available (headless server, plain
    SSH login), the UI is served on ``http://localhost:<port>`` instead and the
    URL is printed so it can be opened locally or through SSH port forwarding.
    """
    from basisremy.gui.application import run_app

    if browser or not _has_display():
        run_app(native=False, show=False, port=port, host="127.0.0.1",
                on_ready=lambda: print(_browser_banner(port, auto=not browser), flush=True))
        return

    run_app(native=True, port=port)


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
        import nicegui  # noqa: F401
        gui_ok = True
    except Exception:
        gui_ok = False
    print(f"  nicegui (GUI)   : {_mark(gui_ok)}")
    try:
        import webview  # noqa: F401  (pywebview, for the native window)
        native_ok = True
    except Exception:
        native_ok = False
    print(f"  pywebview (win) : {_mark(native_ok)}")

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

    if not gui_ok:
        print("Note: the GUI needs NiceGUI (pip install nicegui). "
              "See the README 'Troubleshooting' section.")
    if not native_ok:
        print("Note: the native desktop window needs pywebview (pip install pywebview); "
              "without it the UI opens in your browser.")
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
        help="Report the environment (Python, NiceGUI, Docker/Octave) and exit, "
             "without launching the GUI.",
    )
    parser.add_argument(
        "--browser", action="store_true",
        help="Serve the UI at http://localhost:PORT for a web browser instead of "
             "opening a desktop window (automatic when no display is available, "
             "e.g. on a server; use SSH port forwarding to reach it).",
    )
    parser.add_argument(
        "--port", type=int, default=8080, metavar="PORT",
        help="Port for the UI server (default: 8080).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Console-script entry point for the ``basisremy`` command."""
    args = _build_parser().parse_args(argv)

    _prepare_runtime()

    if args.check:
        return _run_environment_check()

    try:
        _start_gui(browser=args.browser, port=args.port)
    except ModuleNotFoundError as exc:  # most commonly missing nicegui / pywebview
        missing = getattr(exc, "name", "") or str(exc)
        if missing in ("nicegui",) or "nicegui" in str(exc).lower():
            sys.stderr.write(
                "\nERROR: NiceGUI is not installed in this Python environment.\n\n"
                "  pip install nicegui pywebview\n"
                "  (or 'uv sync' / 'uvx --from . basisremy' in the project)\n\n"
            )
            return 1
        if missing in ("webview", "pywebview") or "webview" in str(exc).lower():
            sys.stderr.write(
                "\nERROR: pywebview is required for the native desktop window.\n\n"
                "  pip install pywebview\n\n"
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
