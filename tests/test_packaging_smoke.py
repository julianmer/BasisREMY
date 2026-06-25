####################################################################################################
#                                    test_packaging_smoke.py                                       #
####################################################################################################
#                                                                                                  #
# Purpose: Minimal packaging smoke tests for the ``basisremy`` console entry point.                #
#                                                                                                  #
#          These verify that the package imports, the CLI entry point is wired correctly, and the  #
#          application can initialize far enough to catch import/path errors -- without launching  #
#          the long-running GUI.                                                                   #
#                                                                                                  #
####################################################################################################

import importlib
import importlib.metadata as md

import pytest


def test_launcher_module_imports():
    """The package entry module imports and exposes a callable ``main``."""
    launcher = importlib.import_module("basisremy.__main__")
    assert callable(launcher.main)
    assert callable(launcher._prepare_runtime)
    assert callable(launcher._load_application)


def test_project_root_is_resolvable():
    """The launcher can locate a project root containing the package and externals."""
    launcher = importlib.import_module("basisremy.__main__")
    root = launcher._find_project_root()
    assert (root / "basisremy").is_dir()
    assert (root / "externals").is_dir()


def test_console_entry_point_registered():
    """The ``basisremy`` console script points at ``basisremy.__main__:main``.

    Skipped when the project is not installed (e.g. running pytest against a
    bare source checkout without ``uv sync`` / ``pip install``).
    """
    try:
        eps = md.entry_points(group="console_scripts")
    except TypeError:  # Python < 3.10 selection API
        eps = md.entry_points().get("console_scripts", [])

    match = {ep.name: ep.value for ep in eps}
    if "basisremy" not in match:
        pytest.skip("Project not installed; entry point unavailable.")
    assert match["basisremy"] == "basisremy.__main__:main"


def test_application_class_loads():
    """The GUI application class imports without opening a window.

    This exercises the full import chain (basisremy.gui -> core -> backends ->
    externals) so packaging/path regressions are caught early. It is skipped
    only when the Tk GUI toolkit itself is unavailable in the test environment.
    """
    launcher = importlib.import_module("basisremy.__main__")
    launcher._prepare_runtime()
    try:
        app_cls = launcher._load_application()
    except ModuleNotFoundError as exc:
        if "tkinter" in str(exc).lower():
            pytest.skip(f"Tk GUI toolkit unavailable: {exc}")
        raise
    assert isinstance(app_cls, type)
    assert app_cls.__name__ == "Application"
