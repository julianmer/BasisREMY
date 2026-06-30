####################################################################################################
#                                          externals.py                                            #
####################################################################################################
#                                                                                                  #
# Purpose: Fetch the large third-party Octave/MATLAB libraries that BasisREMY's simulation         #
#          backends rely on. These are NOT shipped in the wheel (they are big and carry their own  #
#          licenses, some of which forbid redistribution), so each is cloned on demand from its    #
#          original source at a pinned commit into ``externals_root()``.                           #
#                                                                                                  #
#          In a source checkout the ``externals/`` git submodules are already present, so          #
#          :func:`ensure` is a no-op there.                                                        #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import shutil
import subprocess
import sys

from basisremy.core.paths import externals_root

# name -> (git url, pinned commit). Commits match the repo's submodule pins.
REGISTRY: dict[str, tuple[str, str]] = {
    "fidA": (
        "https://github.com/CIC-methods/FID-A.git",
        "0c3611c4c5529b8d138317183dc68b4577f7df02",
    ),
    "jbss": (
        "https://github.com/arcj-hub/BasisSetSimulation.git",
        "ad9535eddff9eb2066e1fcdf9e34c5b60d6ae42d",
    ),
    "fsl_mrs": (
        "https://git.fmrib.ox.ac.uk/fsl/fsl_mrs.git",
        "b67e9235c758072c2fbfd36d6dab21d6aac31fc6",
    ),
    "mrscloud": (
        "https://github.com/shui5/MRSCloud.git",
        "8877e3fdc2bf31e85256de25d85ef659d8272cfc",
    ),
    # MRS Basis Set Conversion Toolbox — used by core.exporters to write the
    # various basis-set formats (LCModel / jMRUI / FSL-MRS / Osprey).
    "kbsct": (
        "https://github.com/igweckay/MRS-Basis-Set-Conversion-Toolbox.git",
        "53925137e29fbefd6582171595af59665edd3f9f",
    ),
}


class ExternalFetchError(RuntimeError):
    """Raised when an external library could not be fetched."""


def is_present(name: str) -> bool:
    """Return True if the external ``name`` already exists locally."""
    dest = externals_root() / name
    return dest.is_dir() and any(dest.iterdir())


def ensure(name: str) -> str:
    """Ensure external ``name`` is available locally, cloning it if missing.

    Returns the path to the external as a string. Raises
    :class:`ExternalFetchError` if the name is unknown or the clone fails.
    """
    if name not in REGISTRY:
        raise ExternalFetchError(f"Unknown external '{name}'.")

    dest = externals_root() / name
    if is_present(name):
        return str(dest)

    url, commit = REGISTRY[name]
    dest.parent.mkdir(parents=True, exist_ok=True)

    print(
        f"Fetching '{name}' (one-time download from {url}). This may take a while...",
        file=sys.stderr,
    )
    # Skip Git-LFS smudging during clone/checkout. Some upstreams (notably
    # FID-A) keep large example datasets in LFS and occasionally exceed their
    # LFS budget, which makes the smudge filter fail and aborts the whole
    # clone. We only need the Octave/MATLAB source, not the LFS sample data,
    # so disabling smudge keeps fetching robust regardless of LFS quota.
    env = {**os.environ, "GIT_LFS_SKIP_SMUDGE": "1"}
    try:
        subprocess.run(
            ["git", "clone", "--quiet", url, str(dest)],
            check=True,
            env=env,
        )
        # Pin to the recorded commit. Upstreams occasionally rewrite history and
        # drop the pinned commit (it may also live on a branch the default clone
        # didn't materialise). Try a direct fetch of the commit, and if it is
        # genuinely gone, fall back to the cloned default branch with a warning
        # so the backend still works instead of failing outright.
        if subprocess.run(
            ["git", "-C", str(dest), "checkout", "--quiet", commit], env=env
        ).returncode != 0:
            subprocess.run(
                ["git", "-C", str(dest), "fetch", "--quiet", "origin", commit],
                env=env,
            )
            if subprocess.run(
                ["git", "-C", str(dest), "checkout", "--quiet", commit], env=env
            ).returncode != 0:
                print(
                    f"  ⚠️  pinned commit {commit[:10]} for '{name}' is no longer "
                    "available upstream; using the repository's default branch "
                    "instead.",
                    file=sys.stderr,
                )
    except FileNotFoundError as exc:  # git not installed
        raise ExternalFetchError(
            "git is required to fetch simulation backends but was not found on PATH."
        ) from exc
    except subprocess.CalledProcessError as exc:
        # Remove the partial checkout so the next attempt starts clean instead
        # of tripping the is_present() short-circuit on a half-fetched dir.
        shutil.rmtree(dest, ignore_errors=True)
        raise ExternalFetchError(
            f"Failed to fetch '{name}' from {url} (commit {commit[:10]}): {exc}"
        ) from exc

    return str(dest)
