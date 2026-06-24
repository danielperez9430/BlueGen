#!/usr/bin/env python3
"""
Centralized tool detection for the PRS Research Pipeline.

Searches for PLINK, bcftools, tabix, and other system tools
with version validation and clear error messages.

Search order:
    1. Project-local tools/ directory
    2. System PATH
    3. Common install locations (/usr/local/bin, /opt/homebrew/bin)

Usage:
    from utils.tool_detection import find_plink, find_bcftools
    plink_path, plink_version = find_plink()  # raises SystemExit if not found
    bcftools_path = find_bcftools()  # returns None if not found
"""

import shutil
import subprocess
import sys
from pathlib import Path


def _check_version(binary: str, min_version: str | None = None) -> str:
    """Get version string from a binary. Returns 'unknown' on failure."""
    try:
        r = subprocess.run(
            [binary, "--version"], capture_output=True, text=True, timeout=5
        )
        first_line = r.stdout.split("\n")[0].strip() if r.stdout else "unknown"
        return first_line
    except Exception:
        return "unknown"


def _find_plink_raw() -> tuple[str | None, str | None]:
    """Find PLINK binary. Returns (path, version) or (None, None) if not found."""
    project_root = Path(__file__).resolve().parent.parent.parent.parent
    candidates = [
        project_root / "tools" / "plink",
        project_root / "tools" / "plink2",
        Path.home() / "bin" / "plink",
        Path("/usr/local/bin/plink"),
        Path("/opt/homebrew/bin/plink"),
    ]

    for c in candidates:
        if c.exists() and c.is_file():
            ver = _check_version(str(c))
            return str(c), ver

    for name in ["plink2", "plink"]:
        path = shutil.which(name)
        if path:
            ver = _check_version(path)
            return path, ver

    return None, None


def find_plink(min_version: str = "1.90") -> tuple[str, str]:
    """
    Find PLINK binary. Returns (path, version_string).

    Search order:
      1. tools/plink (project-local)
      2. tools/plink2 (project-local)
      3. System PATH (plink, plink2)
      4. Common install locations

    Raises SystemExit if not found.
    """
    path, ver = _find_plink_raw()
    if path:
        return path, ver or "unknown"

    print(
        "\n\033[0;31m✗ PLINK not found.\033[0m\n"
        "  Download PLINK 1.9: https://www.cog-genomics.org/plink/\n"
        "  Place the binary in: tools/plink\n"
        "  Or install: brew install plink",
        file=sys.stderr,
    )
    sys.exit(1)


def find_bcftools() -> str | None:
    """Find bcftools. Returns path or None if not found."""
    path = shutil.which("bcftools")
    if path:
        return path

    # Check project tools/
    project_root = Path(__file__).resolve().parent.parent.parent.parent
    local = project_root / "tools" / "bcftools"
    if local.exists():
        return str(local)

    return None


def find_tabix() -> str | None:
    """Find tabix. Returns path or None if not found."""
    path = shutil.which("tabix")
    if path:
        return path

    project_root = Path(__file__).resolve().parent.parent.parent.parent
    local = project_root / "tools" / "tabix"
    if local.exists():
        return str(local)

    return None


def require_bcftools() -> str:
    """Find bcftools or exit with clear error message."""
    path = find_bcftools()
    if not path:
        print(
            "\n\033[0;31m✗ bcftools not found.\033[0m\n"
            "  Install: brew install bcftools\n"
            "  Or: apt install bcftools",
            file=sys.stderr,
        )
        sys.exit(1)
    return path


def require_tabix() -> str:
    """Find tabix or exit with clear error message."""
    path = find_tabix()
    if not path:
        print(
            "\n\033[0;31m✗ tabix not found.\033[0m\n"
            "  Install: brew install htslib\n"
            "  Or: apt install tabix",
            file=sys.stderr,
        )
        sys.exit(1)
    return path


def detect_all() -> dict:
    """Detect all tools and return a summary dict. Graceful when tools missing."""
    plink_path, plink_ver = _find_plink_raw()
    return {
        "plink": {"path": plink_path or "not found",
                  "version": plink_ver or "unknown"},
        "bcftools": {"path": find_bcftools(), "available": find_bcftools() is not None},
        "tabix": {"path": find_tabix(), "available": find_tabix() is not None},
    }


def print_tool_status():
    """Print tool detection status to stdout."""
    tools = detect_all()
    print("\n  Tool Detection:")
    for name, info in tools.items():
        if info.get("path") or (isinstance(info, dict) and info.get("path")):
            path = info["path"]
            ver = info.get("version", "")
            ver_str = f" ({ver})" if ver and ver != "unknown" else ""
            print(f"    ✓ {name}: {path}{ver_str}")
        else:
            print(f"    ✗ {name}: not found")
