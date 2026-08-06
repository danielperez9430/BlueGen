import sys
from pathlib import Path

# scripts/ needs to be on sys.path so submodules here can import sibling
# top-level packages (utils.*, clinical.*) the same way comprehensive_report.py
# itself does - done once here instead of repeating it in every submodule.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
