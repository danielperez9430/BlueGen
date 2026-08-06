"""
BlueGen scoring library (IMPROVEMENT_PLAN.md 2.1).

Importable core of the PRS pipeline's critical path (Stages F-H), extracted
from the subprocess-invoked numbered scripts under scripts/prs/ so the
actual scoring math can be unit-tested without subprocess or (for scoring.py)
a real PLINK binary. The numbered scripts remain the CLI entry points and are
still what prs.py's run_script() invokes - they're now thin wrappers around
this package.
"""
