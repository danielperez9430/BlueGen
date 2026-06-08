"""Auto-detect and re-exec under project venv if needed.
Import at the top of any script that needs the venv Python:
    import _venv

Place this BEFORE any non-stdlib imports (pandas, numpy, etc.).
Safe to import from any subdirectory depth."""

import sys, os

# Find the venv by walking up from this file's location
_this_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.dirname(_this_dir)  # prs_research_pipeline/
_env_root = os.path.dirname(_project_root)   # project root
_venv_python = os.path.join(_env_root, "venv", "bin", "python3")
_venv_real = os.path.realpath(_venv_python) if os.path.exists(_venv_python) else ""

# Also search parent/../venv for scripts in nested subdirectories
if not _venv_real:
    for _depth in range(1, 5):
        _candidate = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            *([".."] * _depth), "venv", "bin", "python3"
        )
        _candidate = os.path.realpath(_candidate)
        if os.path.exists(_candidate):
            _venv_real = _candidate
            break

# The venv uses the same Python binary as the system (macOS Homebrew).
# Detect by checking if we're actually USING the venv's site-packages.
if _venv_real and os.path.exists(_venv_real):
    _in_venv = any("venv" in p for p in sys.path if "site-packages" in p)
    if not _in_venv:
        _script = sys.argv[0]
        # Only re-exec if running an actual script file (not -c, -m, or interactive)
        if _script and not _script.startswith("-") and os.path.exists(
            _script if os.path.isabs(_script) else os.path.join(os.getcwd(), _script)
        ):
            _script = os.path.abspath(_script) if os.path.isabs(_script) else os.path.join(os.getcwd(), _script)
            os.execv(_venv_real, [_venv_real, _script] + sys.argv[1:])

# Clean up module-level names so they don't leak
del _this_dir, _project_root, _env_root, _venv_python, _venv_real
