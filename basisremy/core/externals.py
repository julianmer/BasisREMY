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

import importlib
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
    # Spinach (Kuprov; MIT) — the spin-dynamics kernel behind the Spinach
    # backend. Fetched sparsely (see SPARSE) and patched for Octave (PATCHES).
    "spinach": (
        "https://github.com/IlyaKuprov/Spinach.git",
        "998fbc02777f4f7785757494b43dbb42a0954274",
    ),
}

# Externals that are too big to clone whole: only these top-level directories
# are checked out (git sparse checkout, blobs fetched on demand). Spinach's
# repository is ~520 MB; its kernel is all the backend needs.
SPARSE: dict[str, list[str]] = {
    "spinach": ["kernel"],
}

# Patches applied on top of the pinned commit (idempotent: skipped when the
# tree already carries them). Spinach is written for MATLAB R2024b; the patch
# is the small set of source edits Octave 7 needs (see spinach_octave.patch).
PATCHES: dict[str, str] = {
    "spinach": os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                            "adapters", "backends", "spinach_octave.patch"),
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
    # Skip Git-LFS smudging during clone/checkout. Some upstreams (notably
    # FID-A) keep large example datasets in LFS and occasionally exceed their
    # LFS budget, which makes the smudge filter fail and aborts the whole
    # clone. We only need the Octave/MATLAB source, not the LFS sample data,
    # so disabling smudge keeps fetching robust regardless of LFS quota.
    env = {**os.environ, "GIT_LFS_SKIP_SMUDGE": "1"}

    if is_present(name):
        # Self-heal fetches made before nested-submodule support: e.g.
        # fsl_mrs nests denmatsim as a submodule, which a plain clone leaves
        # as an empty directory. A no-op (~ms) when already complete.
        _init_nested_submodules(dest, env)
        _apply_patch(name, dest, env)
        return str(dest)

    url, commit = REGISTRY[name]
    dest.parent.mkdir(parents=True, exist_ok=True)

    print(
        f"Fetching '{name}' (one-time download from {url}). This may take a while...",
        file=sys.stderr,
    )
    try:
        if name in SPARSE:
            # blob-less sparse clone: only the listed directories are
            # materialised, and only their blobs are downloaded
            subprocess.run(
                ["git", "clone", "--quiet", "--filter=blob:none", "--sparse",
                 url, str(dest)],
                check=True,
                env=env,
            )
            subprocess.run(
                ["git", "-C", str(dest), "sparse-checkout", "set", *SPARSE[name]],
                check=True,
                env=env,
            )
        else:
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

    # Some externals nest submodules of their own (fsl_mrs -> denmatsim);
    # without this their directories stay empty after the clone.
    _init_nested_submodules(dest, env)
    _apply_patch(name, dest, env)

    # A sys.path entry that pointed at this (then-missing) directory has been
    # negatively cached by the import system (sys.path_importer_cache), so
    # Python-package externals (e.g. fsl_mrs's denmatsim) would stay
    # unimportable until the app restarts. Drop the stale caches.
    importlib.invalidate_caches()

    return str(dest)


def _init_nested_submodules(dest, env) -> None:
    """Initialize an external's own nested submodules (best effort)."""
    if not os.path.exists(os.path.join(str(dest), '.git')):
        return
    subprocess.run(
        ['git', '-C', str(dest), 'submodule', 'update', '--init',
         '--recursive', '--quiet'],
        env=env, check=False,
    )


def _apply_patch(name: str, dest, env) -> None:
    """Apply the Octave patch registered for ``name`` (no-op when absent or
    already applied). A checkout that carries an *older* version of the patch
    (a cached CI tree, a user upgrading BasisREMY) is completed file by file:
    every per-file chunk is skipped when already applied and applied when it
    fits. Raises :class:`ExternalFetchError` when a chunk fits neither way —
    the checkout is not the pinned commit."""
    patch = PATCHES.get(name)
    if patch is None or not os.path.exists(os.path.join(str(dest), '.git')):
        return
    patch = os.path.abspath(patch)
    git = ['git', '-C', str(dest)]
    kw = dict(env=env, capture_output=True)

    def fits(path, reverse=False):
        args = [*git, 'apply', '--check'] + (['--reverse'] if reverse else []) + [path]
        return subprocess.run(args, **kw).returncode == 0

    if fits(patch, reverse=True):
        return  # already applied
    if fits(patch):
        subprocess.run([*git, 'apply', patch], env=env, check=True)
        return
    # partially applied: go chunk by chunk
    import tempfile
    with open(patch, encoding='utf-8') as f:
        text = f.read()
    for chunk in _split_patch(text):
        tmp = tempfile.NamedTemporaryFile('w', suffix='.patch', delete=False, encoding='utf-8')
        try:
            tmp.write(chunk)
            tmp.close()
            if fits(tmp.name, reverse=True):
                continue
            if fits(tmp.name):
                subprocess.run([*git, 'apply', tmp.name], env=env, check=True)
                continue
            target = chunk.splitlines()[0]
            raise ExternalFetchError(
                f"The Octave patch for '{name}' does not fit the checkout at {dest} "
                f"({target}; expected commit {REGISTRY[name][1][:10]})."
            )
        finally:
            os.remove(tmp.name)


def _split_patch(text: str) -> list[str]:
    """Split a unified git diff into one chunk per file."""
    chunks, current = [], []
    for line in text.splitlines(keepends=True):
        if line.startswith('diff --git ') and current:
            chunks.append(''.join(current))
            current = []
        current.append(line)
    if current:
        chunks.append(''.join(current))
    return chunks
