####################################################################################################
#                                         spant_manager.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 28/08/26                                                                                #
#                                                                                                  #
# Purpose: Runtime for the spant backend (R). Same shape as the Octave runtime: Docker first (a    #
#          rocker/r-ver image with spant, built on first use), the user's own R installation as    #
#          the fallback, an environment variable to force either. Both runtimes drive the same     #
#          worker script (adapters/backends/spant_worker.R) through Rscript and JSON files, under  #
#          a hard timeout — nothing has to be compiled against the user's R.                       #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import glob
import json
import os
import platform
import shutil
import subprocess
import tempfile
from pathlib import Path

from basisremy.core.paths import ADAPTERS_DIR

_WORKER = ADAPTERS_DIR / 'backends' / 'spant_worker.R'
_DOCKER_IMAGE = 'basisremy-spant:latest'
# rocker/r-ver ships R with its build toolchain (multi-arch); spant's
# compiled dependencies (nloptr, RNifti, ...) need cmake and libcurl on top.
_DOCKERFILE = (
    "FROM rocker/r-ver:4.5.1\n"
    "RUN apt-get update && apt-get install -y --no-install-recommends "
    "cmake libcurl4-openssl-dev libssl-dev && rm -rf /var/lib/apt/lists/*\n"
    "RUN install2.r --error --skipinstalled spant\n"
)


class SpantUnavailable(RuntimeError):
    """spant cannot run on this machine (with guidance in the message)."""


# ------------------------------------------------------------------ discovery
def docker_available() -> bool:
    from basisremy.core import pygamma_manager
    return pygamma_manager.docker_available()


def _windows_rscript_candidates() -> list[str]:
    out = []
    try:
        import winreg  # Windows only
        for root in (winreg.HKEY_LOCAL_MACHINE, winreg.HKEY_CURRENT_USER):
            for key in (r'SOFTWARE\R-core\R', r'SOFTWARE\R-core\R64'):
                try:
                    with winreg.OpenKey(root, key) as k:
                        path, _ = winreg.QueryValueEx(k, 'InstallPath')
                        out.append(os.path.join(path, 'bin', 'Rscript.exe'))
                except OSError:
                    continue
    except ImportError:
        pass
    for base in (os.environ.get('ProgramFiles', r'C:\Program Files'),
                 os.path.join(os.environ.get('LOCALAPPDATA', ''), 'Programs')):
        out.append(os.path.join(base, 'R', 'R-*', 'bin', 'Rscript.exe'))
    return out


def find_rscript() -> str | None:
    """The user's Rscript: $BASISREMY_RSCRIPT, then PATH, then the usual
    install locations of each OS (newest version first)."""
    env = os.environ.get('BASISREMY_RSCRIPT')
    if env and os.path.isfile(env):
        return env
    exe = shutil.which('Rscript')
    if exe:
        return exe
    system = platform.system()
    if system == 'Darwin':
        candidates = ['/Library/Frameworks/R.framework/Resources/bin/Rscript',
                      '/opt/homebrew/bin/Rscript', '/usr/local/bin/Rscript']
    elif system == 'Windows':
        candidates = _windows_rscript_candidates()
    else:
        candidates = ['/usr/bin/Rscript', '/usr/local/bin/Rscript',
                      '/opt/R/*/bin/Rscript']
    for pattern in candidates:
        for path in sorted(glob.glob(pattern), reverse=True):
            if os.path.isfile(path):
                return path
    return None


def spant_version(rscript: str, timeout: float = 120.0) -> str | None:
    """Installed spant version for this Rscript, or None."""
    try:
        proc = subprocess.run(
            [rscript, '-e', 'cat(as.character(packageVersion("spant")))'],
            capture_output=True, text=True, timeout=timeout)
    except Exception:  # noqa: BLE001  (timeout, unrunnable binary, ...)
        return None
    out = proc.stdout.strip()
    return out if proc.returncode == 0 and out else None


def ensure_spant(rscript: str, log=print, timeout: float = 3600.0) -> str:
    """Install spant into the user's R library on first use; return its version."""
    version = spant_version(rscript)
    if version:
        return version
    log("Installing the spant R package into your user library (one-time; "
        "when R has to build it from source this takes 10-40 min)...")
    code = (
        'lib <- Sys.getenv("R_LIBS_USER"); '
        'dir.create(lib, recursive = TRUE, showWarnings = FALSE); '
        '.libPaths(c(lib, .libPaths())); '
        'options(Ncpus = max(1L, parallel::detectCores() - 1L)); '
        'install.packages("spant", lib = lib, repos = "https://cloud.r-project.org")'
    )
    try:
        proc = subprocess.run([rscript, '-e', code], capture_output=True,
                              text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise SpantUnavailable(
            f"Installing spant did not finish within {timeout / 60:.0f} min. "
            f"Install it yourself in R (install.packages('spant')) or use the "
            f"Docker runtime.")
    version = spant_version(rscript)
    if not version:
        raise SpantUnavailable(
            "Installing the spant R package failed:\n"
            + (proc.stderr or proc.stdout or '').strip()[-1500:]
            + "\n\nWhen R builds spant from source it needs CMake (for nloptr): "
            "macOS `brew install cmake`, Debian/Ubuntu `sudo apt install cmake`, "
            "Windows https://cmake.org/download/ — then simulate again.")
    log(f"✓ spant {version} installed")
    return version


def ensure_docker_image(log=print, timeout: float = 3600.0) -> str:
    """Build the spant Docker image on first use."""
    have = subprocess.run(['docker', 'image', 'inspect', _DOCKER_IMAGE],
                          capture_output=True)
    if have.returncode == 0:
        return _DOCKER_IMAGE
    log("Building the spant Docker image (one-time; R compiles spant and its "
        "dependencies, 10-30 min)...")
    try:
        proc = subprocess.run(['docker', 'build', '-t', _DOCKER_IMAGE, '-'],
                              input=_DOCKERFILE, capture_output=True, text=True,
                              timeout=timeout)
    except subprocess.TimeoutExpired:
        raise SpantUnavailable(
            f"Building the spant Docker image did not finish within "
            f"{timeout / 60:.0f} min.")
    if proc.returncode != 0:
        raise SpantUnavailable(
            f"spant Docker image build failed:\n{proc.stderr[-1500:]}")
    log("✓ spant Docker image ready")
    return _DOCKER_IMAGE


# ------------------------------------------------------------------ runtime choice
def installation_instructions() -> str:
    system = platform.system()
    r_hint = {
        'Darwin': "  • macOS: brew install r   (or the CRAN installer, https://cran.r-project.org/bin/macosx/)",
        'Windows': "  • Windows: https://cran.r-project.org/bin/windows/base/",
    }.get(system, "  • Linux: sudo apt-get install r-base   (Debian/Ubuntu) — or https://cran.r-project.org/")
    return "\n".join([
        "spant runtime not available. BasisREMY runs spant (R) either",
        "",
        "OPTION 1 (preferred): Docker — install and start Docker Desktop / Engine;",
        "  the spant image is built automatically on first use.",
        "",
        "OPTION 2: your own R — install R, BasisREMY installs the spant package",
        "  into your user library on first use:",
        r_hint,
        "  (set BASISREMY_RSCRIPT=/path/to/Rscript if R is not on PATH; when R",
        "  builds spant from source it needs CMake: brew/apt install cmake)",
        "",
        "BASISREMY_SPANT_RUNTIME=docker|local forces one of the two.",
    ])


def preferred_runtime(prefer_docker: bool = True) -> str:
    """'docker' or 'local' — Docker first like the Octave runtime, the user's
    own R as the fallback; $BASISREMY_SPANT_RUNTIME forces one."""
    forced = os.environ.get('BASISREMY_SPANT_RUNTIME')
    if forced in ('docker', 'local'):
        return forced
    order = ('docker', 'local') if prefer_docker else ('local', 'docker')
    for runtime in order:
        if runtime == 'docker' and docker_available():
            return 'docker'
        if runtime == 'local' and find_rscript():
            return 'local'
    raise SpantUnavailable(installation_instructions())


def ensure_runtime(runtime: str, log=print) -> str:
    """Make the chosen runtime ready; returns a short description."""
    if runtime == 'docker':
        if not docker_available():
            raise SpantUnavailable("Docker is not running.\n" + installation_instructions())
        ensure_docker_image(log)
        return f"Docker image {_DOCKER_IMAGE}"
    rscript = find_rscript()
    if not rscript:
        raise SpantUnavailable("No R installation found.\n" + installation_instructions())
    version = ensure_spant(rscript, log)
    return f"{rscript} (spant {version})"


# ------------------------------------------------------------------ worker
def run_worker(job: dict, runtime: str | None = None,
               timeout: float = 1800.0) -> dict:
    """Run one job through spant_worker.R; returns the parsed out.json."""
    runtime = runtime or preferred_runtime()
    with tempfile.TemporaryDirectory(prefix='basisremy_spant_') as tmp:
        tmp_path = Path(tmp)
        job = dict(job)
        # stage the waveform next to the job so the container sees it too
        if job.get('pulse_file'):
            src = Path(job['pulse_file'])
            dst = tmp_path / src.name
            shutil.copy2(src, dst)
            job['pulse_file'] = f'/data/{src.name}' if runtime == 'docker' else str(dst)
        (tmp_path / 'job.json').write_text(json.dumps(job))
        shutil.copy2(_WORKER, tmp_path / 'spant_worker.R')
        out_path = tmp_path / 'out.json'
        if runtime == 'docker':
            cmd = ['docker', 'run', '--rm', '-v', f'{tmp}:/data', _DOCKER_IMAGE,
                   'Rscript', '/data/spant_worker.R', '/data/job.json', '/data/out.json']
        else:
            rscript = find_rscript()
            if not rscript:
                raise SpantUnavailable(installation_instructions())
            cmd = [rscript, str(tmp_path / 'spant_worker.R'),
                   str(tmp_path / 'job.json'), str(out_path)]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"spant worker timed out after {timeout:.0f}s.")
        if not out_path.exists():
            raise RuntimeError(
                f"spant worker produced no result (exit {proc.returncode}):\n"
                f"{(proc.stderr or proc.stdout or '').strip()[-1500:]}")
        result = json.loads(out_path.read_text())
    if not result.get('ok'):
        raise RuntimeError(f"spant simulation error: {result.get('error')}")
    return result
