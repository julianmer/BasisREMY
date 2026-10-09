"""BasisREMY — study-specific MR spectroscopy basis set generation."""

__version__ = "0.2.0"


def prepare_runtime():
    """Make BasisREMY's runtime folder (where the simulation engines' code is fetched) the working
    directory, as the ``basisremy`` command does, and return it. Call it at the start of a script
    that simulates, and give file paths as absolute paths."""
    from basisremy.__main__ import _prepare_runtime
    return _prepare_runtime()
