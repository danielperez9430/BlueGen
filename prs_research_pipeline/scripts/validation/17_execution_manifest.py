#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   PHASE 7 — MODULE 2: PIPELINE EXECUTION MANIFEST                            ║
║   scripts/17_execution_manifest.py                                           ║
║                                                                            ║
║   Records a complete, auditable execution trace of the full pipeline.       ║
║                                                                            ║
║   SCIENTIFIC FREEZE LAYER — Every stage is logged. Every file tracked.      ║
║                                                                            ║
║   Records per stage:                                                        ║
║     • Input and output file paths + hashes                                  ║
║     • Wall-clock execution time                                             ║
║     • Peak memory usage (if available)                                      ║
║     • Sample counts (variants, individuals)                                 ║
║     • SNP retention rates through QC/LD/PCA steps                           ║
║     • Error/warning counts                                                  ║
║                                                                            ║
║   Output:                                                                   ║
║     validation/execution_manifest.json                                      ║
║     validation/execution_timeline.tsv                                       ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
import json
import hashlib
import logging
import time
import functools
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ── Data Structures ────────────────────────────────────────────────────────────

@dataclass
class StageRecord:
    """Record of a single pipeline stage execution."""
    stage_id: str
    stage_name: str
    phase: str = ""
    start_time: str = ""
    end_time: str = ""
    duration_seconds: float = 0.0
    status: str = "pending"  # pending, running, success, failed, skipped
    input_files: Dict[str, str] = field(default_factory=dict)
    output_files: Dict[str, str] = field(default_factory=dict)
    metrics: Dict[str, Any] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


@dataclass
class ExecutionManifest:
    """Complete pipeline execution manifest."""
    run_id: str
    pipeline_version: str = "7.0.0"
    started_at: str = ""
    completed_at: str = ""
    total_duration_seconds: float = 0.0
    total_stages: int = 0
    stages_completed: int = 0
    stages_failed: int = 0
    stages_skipped: int = 0
    stages: List[StageRecord] = field(default_factory=list)
    variant_flow: Dict[str, int] = field(default_factory=dict)
    sample_flow: Dict[str, int] = field(default_factory=dict)
    provenance: Dict[str, str] = field(default_factory=dict)


# ── Execution Tracer ──────────────────────────────────────────────────────────

class ExecutionTracer:
    """
    Traces every pipeline stage execution for audit trail.

    Usage:
        tracer = ExecutionTracer(run_id="run_001")
        tracer.start()

        with tracer.stage("A", "VCF → PLINK"):
            # ... execute stage ...

        with tracer.stage("B", "Quality Control"):
            # ... execute stage ...

        manifest = tracer.finish()
    """

    def __init__(self, run_id: str = "", output_dir: str = "validation"):
        self.run_id = run_id or datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.manifest = ExecutionManifest(run_id=self.run_id)
        self._start_time: Optional[float] = None

    def start(self) -> None:
        """Start the pipeline trace."""
        self._start_time = time.time()
        self.manifest.started_at = datetime.now(timezone.utc).isoformat()
        logger.info("═══ Execution Trace STARTED ═══")

    def stage(self, stage_id: str, stage_name: str, phase: str = "") -> "StageContext":
        """Context manager for a pipeline stage."""
        return StageContext(self, stage_id, stage_name, phase)

    def record_stage(
        self,
        stage_id: str,
        stage_name: str,
        phase: str = "",
        inputs: Optional[Dict[str, str]] = None,
        outputs: Optional[Dict[str, str]] = None,
        metrics: Optional[Dict[str, Any]] = None,
        duration: float = 0.0,
        status: str = "success",
        errors: Optional[List[str]] = None,
        warnings: Optional[List[str]] = None,
    ) -> StageRecord:
        """Record a completed stage."""
        record = StageRecord(
            stage_id=stage_id,
            stage_name=stage_name,
            phase=phase,
            start_time=datetime.now(timezone.utc).isoformat(),
            end_time=datetime.now(timezone.utc).isoformat(),
            duration_seconds=round(duration, 3),
            status=status,
            input_files=inputs or {},
            output_files=outputs or {},
            metrics=metrics or {},
            errors=errors or [],
            warnings=warnings or [],
        )
        self.manifest.stages.append(record)
        self._update_counts(record)
        return record

    def track_variant_flow(self, stage: str, n_variants: int) -> None:
        """Track variant count through pipeline stages."""
        self.manifest.variant_flow[stage] = n_variants

    def track_sample_flow(self, stage: str, n_samples: int) -> None:
        """Track sample count through pipeline stages."""
        self.manifest.sample_flow[stage] = n_samples

    def set_provenance(self, key: str, value: str) -> None:
        """Record provenance metadata."""
        self.manifest.provenance[key] = value

    def finish(self) -> ExecutionManifest:
        """Complete the trace and save manifest."""
        if self._start_time:
            self.manifest.total_duration_seconds = round(time.time() - self._start_time, 1)
        self.manifest.completed_at = datetime.now(timezone.utc).isoformat()

        # Compute summaries
        self.manifest.total_stages = len(self.manifest.stages)
        self.manifest.stages_completed = sum(
            1 for s in self.manifest.stages if s.status == "success"
        )
        self.manifest.stages_failed = sum(
            1 for s in self.manifest.stages if s.status == "failed"
        )
        self.manifest.stages_skipped = sum(
            1 for s in self.manifest.stages if s.status == "skipped"
        )

        self._save()
        self._save_timeline()

        logger.info("═══ Execution Trace COMPLETE ═══")
        logger.info(f"  Stages: {self.manifest.total_stages} total, "
                   f"{self.manifest.stages_completed} completed, "
                   f"{self.manifest.stages_failed} failed")

        return self.manifest

    def _update_counts(self, record: StageRecord) -> None:
        """Update manifest counts from stage record."""
        if record.status == "success":
            self.manifest.stages_completed += 1
        elif record.status == "failed":
            self.manifest.stages_failed += 1
        elif record.status == "skipped":
            self.manifest.stages_skipped += 1

    def _save(self) -> None:
        """Save execution manifest to JSON."""
        path = self.output_dir / "execution_manifest.json"
        with open(path, "w") as fh:
            json.dump({
                "run_id": self.manifest.run_id,
                "pipeline_version": self.manifest.pipeline_version,
                "started_at": self.manifest.started_at,
                "completed_at": self.manifest.completed_at,
                "total_duration_seconds": self.manifest.total_duration_seconds,
                "total_stages": self.manifest.total_stages,
                "stages_completed": self.manifest.stages_completed,
                "stages_failed": self.manifest.stages_failed,
                "stages_skipped": self.manifest.stages_skipped,
                "variant_flow": self.manifest.variant_flow,
                "sample_flow": self.manifest.sample_flow,
                "provenance": self.manifest.provenance,
                "stages": [
                    {
                        "stage_id": s.stage_id,
                        "stage_name": s.stage_name,
                        "phase": s.phase,
                        "duration_seconds": s.duration_seconds,
                        "status": s.status,
                        "input_files": s.input_files,
                        "output_files": s.output_files,
                        "metrics": s.metrics,
                        "errors": s.errors,
                        "warnings": s.warnings,
                    }
                    for s in self.manifest.stages
                ],
            }, fh, indent=2, default=str)
        logger.info(f"  ✅ Execution manifest: {path}")

    def _save_timeline(self) -> None:
        """Save execution timeline as TSV."""
        path = self.output_dir / "execution_timeline.tsv"
        rows = []
        for s in self.manifest.stages:
            rows.append({
                "stage_id": s.stage_id,
                "stage_name": s.stage_name,
                "duration_s": s.duration_seconds,
                "status": s.status,
                "n_input_files": len(s.input_files),
                "n_output_files": len(s.output_files),
                "n_errors": len(s.errors),
                "n_warnings": len(s.warnings),
            })
        df = pd.DataFrame(rows)
        df.to_csv(path, sep="\t", index=False)
        logger.info(f"  ✅ Execution timeline: {path}")


class StageContext:
    """Context manager for tracing a pipeline stage."""

    def __init__(self, tracer: ExecutionTracer, stage_id: str, stage_name: str, phase: str = ""):
        self.tracer = tracer
        self.stage_id = stage_id
        self.stage_name = stage_name
        self.phase = phase
        self._start: float = 0.0
        self.inputs: Dict[str, str] = {}
        self.outputs: Dict[str, str] = {}
        self.metrics: Dict[str, Any] = {}
        self.errors: List[str] = []
        self.warnings: List[str] = []

    def __enter__(self):
        self._start = time.time()
        logger.info(f"── Stage {self.stage_id}: {self.stage_name} ──")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        duration = time.time() - self._start
        status = "failed" if exc_type else "success"

        if exc_type:
            self.errors.append(str(exc_val))

        self.tracer.record_stage(
            stage_id=self.stage_id,
            stage_name=self.stage_name,
            phase=self.phase,
            inputs=self.inputs,
            outputs=self.outputs,
            metrics=self.metrics,
            duration=duration,
            status=status,
            errors=self.errors,
            warnings=self.warnings,
        )

        if exc_type is None:
            logger.info(f"  ✅ Stage {self.stage_id}: {duration:.1f}s")
        else:
            logger.error(f"  ❌ Stage {self.stage_id}: FAILED ({duration:.1f}s)")

        # Don't suppress exceptions
        return False

    def add_input(self, name: str, path: str) -> None:
        """Track an input file."""
        self.inputs[name] = path

    def add_output(self, name: str, path: str) -> None:
        """Track an output file."""
        self.outputs[name] = path
        # Auto-hash output
        if Path(path).exists():
            sha = hashlib.sha256()
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(65536), b""):
                    sha.update(chunk)
            self.outputs[f"{name}_sha256"] = sha.hexdigest()[:16]

    def add_metric(self, key: str, value: Any) -> None:
        """Add a stage metric."""
        self.metrics[key] = value

    def add_warning(self, msg: str) -> None:
        """Add a stage warning."""
        self.warnings.append(msg)

    def add_error(self, msg: str) -> None:
        """Add a stage error."""
        self.errors.append(msg)


# ── CLI Entry Point ──────────────────────────────────────────────────────────

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 7 Module 2: Pipeline Execution Manifest"
    )
    parser.add_argument("--run-id", help="Run identifier")
    parser.add_argument("--output-dir", "-o", default="validation")
    parser.add_argument("--init", action="store_true",
                       help="Initialize a new execution trace")
    parser.add_argument("--record-stage", nargs=4,
                       metavar=("ID", "NAME", "DURATION", "STATUS"),
                       help="Record a completed stage")
    parser.add_argument("--track-variants", nargs=2,
                       metavar=("STAGE", "COUNT"),
                       help="Track variant flow")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    tracer = ExecutionTracer(run_id=args.run_id or "", output_dir=args.output_dir)

    if args.init:
        tracer.start()

    if args.record_stage:
        sid, sname, dur, status = args.record_stage
        tracer.record_stage(
            stage_id=sid, stage_name=sname,
            duration=float(dur), status=status,
        )

    if args.track_variants:
        stage, count = args.track_variants
        tracer.track_variant_flow(stage, int(count))

    manifest = tracer.finish()

    print(f"\n═══ Execution Manifest ═══")
    print(f"  Run ID: {manifest.run_id}")
    print(f"  Stages: {manifest.total_stages}")
    print(f"  Duration: {manifest.total_duration_seconds:.1f}s")
    print(f"\n  Variant flow:")
    for stage, n in manifest.variant_flow.items():
        print(f"    {stage}: {n:,} variants")

    return 0


if __name__ == "__main__":
    sys.exit(main())
