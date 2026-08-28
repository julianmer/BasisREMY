####################################################################################################
#                                       pygamma_manager.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 27/08/26                                                                                #
#                                                                                                  #
# Purpose: Manage the dedicated PyGAMMA side-environment for the Vespa backend. PyGAMMA only       #
#          publishes wheels for Python <= 3.9 (x86_64), so simulations run in a separate           #
#          uv-managed interpreter (x86_64 via Rosetta on Apple Silicon) driven through a           #
#          subprocess worker with a hard timeout.                                                  #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import tempfile
from pathlib import Path

from basisremy.core.paths import runtime_root

_WORKER = (Path(__file__).resolve().parents[1] / 'adapters' / 'backends'
           / 'pygamma_worker.py')

# uv python keys per platform for the newest cpython with a pygamma wheel
_PY_KEYS = {
    ('Darwin', 'arm64'): 'cpython-3.9-macos-x86_64-none',    # Rosetta
    ('Darwin', 'x86_64'): 'cpython-3.9-macos-x86_64-none',
    ('Linux', 'x86_64'): 'cpython-3.9-linux-x86_64-gnu',
    ('Windows', 'AMD64'): 'cpython-3.9-windows-x86_64-none',
}


class PyGammaUnavailable(RuntimeError):
    """PyGAMMA cannot run on this machine (with guidance in the message)."""


def env_dir() -> Path:
    return runtime_root() / 'pygamma-env'


def env_python() -> Path:
    sub = ('Scripts', 'python.exe') if os.name == 'nt' else ('bin', 'python')
    return env_dir().joinpath(*sub)


def is_available() -> bool:
    return env_python().exists()


def _python_key() -> str:
    key = _PY_KEYS.get((platform.system(), platform.machine()))
    if key is None:
        raise PyGammaUnavailable(
            f"PyGAMMA has no wheel for {platform.system()}/{platform.machine()} "
            f"(wheels exist only for Python 3.9 on x86_64/AMD64). A Docker-based "
            f"runtime for this platform is on the roadmap."
        )
    return key


def ensure_env(create: bool = True, log=print) -> Path:
    """Return the side-env's python, creating the env on first use."""
    python = env_python()
    if python.exists():
        return python
    if not create:
        raise PyGammaUnavailable("PyGAMMA environment not set up yet.")
    uv = shutil.which('uv')
    if uv is None:
        raise PyGammaUnavailable(
            "Setting up the PyGAMMA environment needs the 'uv' tool "
            "(https://docs.astral.sh/uv/), which was not found on PATH."
        )
    key = _python_key()
    log(f"Setting up the PyGAMMA environment (one-time, ~1 min): {env_dir()}")
    try:
        subprocess.run([uv, 'python', 'install', key],
                       check=True, capture_output=True, text=True)
        subprocess.run([uv, 'venv', '--python', key, str(env_dir())],
                       check=True, capture_output=True, text=True)
        subprocess.run([uv, 'pip', 'install', '--python', str(python),
                        'pygamma', 'numpy<2'],
                       check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        raise PyGammaUnavailable(
            f"PyGAMMA environment setup failed:\n{exc.stderr or exc}"
        )
    log("✓ PyGAMMA environment ready")
    return python


_DOCKER_IMAGE = 'basisremy-pygamma:latest'


def docker_available() -> bool:
    try:
        return subprocess.run(['docker', 'info'], capture_output=True,
                              timeout=10).returncode == 0
    except Exception:  # noqa: BLE001
        return False


def ensure_docker_image(log=print) -> str:
    """Build the PyGAMMA Docker image on first use (linux/amd64 wheel)."""
    have = subprocess.run(['docker', 'image', 'inspect', _DOCKER_IMAGE],
                          capture_output=True)
    if have.returncode == 0:
        return _DOCKER_IMAGE
    log("Building the PyGAMMA Docker image (one-time, ~1-2 min)...")
    dockerfile = (
        "FROM --platform=linux/amd64 python:3.9-slim\n"
        "RUN pip install --no-cache-dir pygamma 'numpy<2'\n"
    )
    proc = subprocess.run(
        ['docker', 'build', '--platform', 'linux/amd64',
         '-t', _DOCKER_IMAGE, '-'],
        input=dockerfile, capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise PyGammaUnavailable(
            f"PyGAMMA Docker image build failed:\n{proc.stderr[-800:]}")
    log("✓ PyGAMMA Docker image ready")
    return _DOCKER_IMAGE


def _docker_marker() -> Path:
    return env_dir().with_name('pygamma-prefer-docker')


def prefer_docker(reason: str = '') -> None:
    """Remember that the side-env is broken on this machine (e.g. a wedged
    Rosetta translation) so future runs go straight to Docker."""
    marker = _docker_marker()
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(reason or 'side-env failed; using Docker runtime\n')


def _env_healthy(timeout: float = 60.0) -> bool:
    """Does the side-env import pygamma within `timeout` s? A wedged Rosetta
    translation hangs or crashes right here, so this is also the probe that
    lets a stale 'prefer docker' marker heal itself."""
    python = env_python()
    if not python.exists():
        return False
    try:
        return subprocess.run([str(python), '-c', 'import pygamma'],
                              capture_output=True, timeout=timeout).returncode == 0
    except Exception:  # noqa: BLE001  (timeout, unrunnable interpreter, ...)
        return False


def preferred_runtime() -> str:
    """The local side-env is the default; Docker only when the env is known
    broken (marker / env var) or the platform has no PyGAMMA wheel. The
    marker is dropped again as soon as the side-env responds (e.g. after the
    reboot that clears a wedged Rosetta cache)."""
    if os.environ.get('BASISREMY_PYGAMMA_RUNTIME') in ('docker', 'env'):
        return os.environ['BASISREMY_PYGAMMA_RUNTIME']
    marker = _docker_marker()
    if marker.exists() and docker_available():
        if _env_healthy():
            marker.unlink(missing_ok=True)
            print("✓ PyGAMMA side-environment responds again — using it "
                  "(Docker marker removed)")
        else:
            return 'docker'
    try:
        _python_key()
    except PyGammaUnavailable:
        if docker_available():
            return 'docker'
        raise
    return 'env'


def run_worker(job: dict, timeout: float = 600.0,
               python: 'Path | None' = None,
               worker: 'Path | None' = None,
               runtime: 'str | None' = None) -> dict:
    """Run one simulation job; returns the worker's basis dict.

    runtime 'env' uses the uv side-environment; 'docker' runs the worker in
    the linux/amd64 PyGAMMA image (isolated from any host Rosetta state).
    A hard timeout guards against wedged native code — the caller gets a
    clear error, never a hang.
    """
    worker = Path(worker) if worker else _WORKER
    if runtime is None:
        runtime = 'env' if python else preferred_runtime()
    with tempfile.TemporaryDirectory(prefix='basisremy_pygamma_') as tmp:
        tmp_path = Path(tmp)
        job_path = tmp_path / 'job.json'
        out_path = tmp_path / 'out.json'
        job_path.write_text(json.dumps(job))
        if runtime == 'docker':
            image = ensure_docker_image()
            shutil.copy2(worker, tmp_path / 'worker.py')
            cmd = ['docker', 'run', '--rm', '--platform', 'linux/amd64',
                   '-v', f'{tmp}:/data', image,
                   'python', '/data/worker.py', '/data/job.json',
                   '/data/out.json']
        else:
            python = Path(python) if python else ensure_env(create=True)
            cmd = [str(python), str(worker), str(job_path), str(out_path)]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=timeout)
        except subprocess.TimeoutExpired:
            raise RuntimeError(
                f"PyGAMMA worker timed out after {timeout:.0f}s. If this "
                f"happens right after a crash on Apple Silicon, Rosetta's "
                f"translation cache may be wedged — reboot, or use the "
                f"Docker runtime."
            )
        if proc.returncode != 0:
            raise RuntimeError(
                f"PyGAMMA worker failed (exit {proc.returncode}):\n"
                f"{(proc.stderr or proc.stdout or '').strip()[-1000:]}"
            )
        result = json.loads(out_path.read_text())
    if not result.get('ok'):
        raise RuntimeError(f"PyGAMMA simulation error: {result.get('error')}")
    return result['basis']
