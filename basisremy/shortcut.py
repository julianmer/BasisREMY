####################################################################################################
#                                           shortcut.py                                            #
####################################################################################################
#                                                                                                  #
# Purpose: ``basisremy --install-shortcut``: add a clickable BasisREMY launcher with the logo to   #
#          the desktop menus, so the GUI opens like a normal application.                          #
#                                                                                                  #
#          - macOS  : /Applications/BasisREMY.app (Finder, Launchpad, Spotlight, Dock).            #
#          - Linux  : ~/.local/share/applications/basisremy.desktop (app grid / menu).             #
#          - Windows: BasisREMY.lnk in the Start Menu and on the Desktop.                          #
#                                                                                                  #
#          The launcher runs the same Python environment the command was started from, or          #
#          ``uvx basisremy`` when started through uvx (whose cached environment can be removed).   #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import plistlib
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image

from basisremy import __version__

ICON_PNG = Path(__file__).resolve().parent / "assets" / "imgs" / "basisremy_icon.png"


def _launch_command() -> list[str]:
    """Return the command the launcher runs to open the GUI."""
    if "archive-v0" in Path(sys.executable).parts:  # uvx: a throw-away environment in uv's cache
        uvx = shutil.which("uvx")
        if uvx is None:
            raise RuntimeError("started through uvx, but 'uvx' is not on PATH")
        return [uvx, "basisremy"]
    return [sys.executable, "-m", "basisremy"]


def _install_macos(cmd: list[str]) -> list[Path]:
    # /Applications is where Finder and Launchpad look; ~/Applications if not writable.
    folder = Path("/Applications")
    if not os.access(folder, os.W_OK):
        folder = Path.home() / "Applications"
    # Build next to the target and move it in at once: Finder caches a bundle it
    # sees half-written (no executable yet) as broken, with a prohibited sign.
    target = folder / "BasisREMY.app"
    app = folder / ".BasisREMY.app.tmp"
    if app.exists():
        shutil.rmtree(app)
    (app / "Contents" / "MacOS").mkdir(parents=True)
    (app / "Contents" / "Resources").mkdir()

    Image.open(ICON_PNG).save(app / "Contents" / "Resources" / "BasisREMY.icns")
    with open(app / "Contents" / "Info.plist", "wb") as f:
        plistlib.dump({
            "CFBundleName": "BasisREMY",
            "CFBundleDisplayName": "BasisREMY",
            "CFBundleIdentifier": "io.github.julianmer.basisremy",
            "CFBundleExecutable": "BasisREMY",
            "CFBundleIconFile": "BasisREMY",
            "CFBundlePackageType": "APPL",
            "CFBundleShortVersionString": __version__,
            "NSHighResolutionCapable": True,
        }, f)

    # Apps opened from Finder get a bare PATH (/usr/bin:/bin:...), so keep the
    # terminal's PATH to still find docker, octave, Rscript, uvx, ...
    exe = app / "Contents" / "MacOS" / "BasisREMY"
    exe.write_text("#!/bin/sh\n"
                   f"export PATH={shlex.quote(os.environ.get('PATH', ''))}\n"
                   f"exec {shlex.join(cmd)}\n")
    exe.chmod(0o755)

    if target.exists():
        shutil.rmtree(target)
    app.rename(target)
    lsregister = ("/System/Library/Frameworks/CoreServices.framework/Frameworks/"
                  "LaunchServices.framework/Support/lsregister")
    if os.path.exists(lsregister):  # refresh the cached app record and icon
        subprocess.run([lsregister, "-f", str(target)], check=False)
    return [target]


def _install_linux(cmd: list[str]) -> list[Path]:
    data = Path(os.environ.get("XDG_DATA_HOME") or Path.home() / ".local" / "share")
    icon = data / "icons" / "basisremy.png"
    icon.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(ICON_PNG, icon)

    entry = data / "applications" / "basisremy.desktop"
    entry.parent.mkdir(parents=True, exist_ok=True)
    entry.write_text("[Desktop Entry]\n"
                     "Type=Application\n"
                     "Name=BasisREMY\n"
                     "Comment=Generate MRS basis sets from raw MRS data\n"
                     "Exec=" + " ".join(f'"{a}"' for a in cmd) + "\n"
                     f"Icon={icon}\n"
                     "Terminal=false\n"
                     "Categories=Science;Education;\n")
    entry.chmod(0o755)
    return [entry]


def _install_windows(cmd: list[str]) -> list[Path]:
    folder = Path(os.environ.get("LOCALAPPDATA") or Path.home()) / "BasisREMY"
    folder.mkdir(parents=True, exist_ok=True)
    ico = folder / "basisremy.ico"
    Image.open(ICON_PNG).save(ico, sizes=[(16, 16), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)])

    def ps(s: str) -> str:  # PowerShell single-quoted string
        return "'" + s.replace("'", "''") + "'"

    # WindowStyle 7 = minimised: the console window stays out of the way.
    script = (
        "$w = New-Object -ComObject WScript.Shell\n"
        "foreach ($d in @([Environment]::GetFolderPath('Programs'), [Environment]::GetFolderPath('Desktop'))) {\n"
        "  $p = Join-Path $d 'BasisREMY.lnk'\n"
        "  $s = $w.CreateShortcut($p)\n"
        f"  $s.TargetPath = {ps(cmd[0])}\n"
        f"  $s.Arguments = {ps(subprocess.list2cmdline(cmd[1:]))}\n"
        f"  $s.IconLocation = {ps(str(ico))}\n"
        f"  $s.WorkingDirectory = {ps(str(Path.home()))}\n"
        "  $s.WindowStyle = 7\n"
        "  $s.Save()\n"
        "  Write-Output $p\n"
        "}\n"
    )
    out = subprocess.run(["powershell", "-NoProfile", "-NonInteractive", "-Command", script],
                         capture_output=True, text=True, check=True).stdout
    return [Path(line) for line in out.splitlines() if line.strip()]


def install_shortcut() -> int:
    """Install the BasisREMY launcher for this platform and report where it went."""
    try:
        cmd = _launch_command()
        if sys.platform == "darwin":
            made = _install_macos(cmd)
        elif sys.platform == "win32":
            made = _install_windows(cmd)
        else:
            made = _install_linux(cmd)
    except Exception as exc:  # noqa: BLE001 - surface a friendly error
        sys.stderr.write(f"\nERROR: could not install the BasisREMY shortcut: {exc}\n\n")
        return 1

    print("BasisREMY shortcut installed:")
    for p in made:
        print(f"  {p}")
    print(f"It runs: {shlex.join(cmd)}")
    return 0
