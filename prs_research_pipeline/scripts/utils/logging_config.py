#!/usr/bin/env python3
"""
Unified logging configuration for the PRS Research Pipeline.

Provides consistent logging across all pipeline stages with:
- Colored console output (via ColorFormatter)
- File logging with rotation (5MB default, 3 backups)
- Stage-level progress tracking (via StageProgress context manager)

Usage:
    from utils.logging_config import setup_logging, stage_logger, StageProgress

    # Initialize at pipeline start
    setup_logging(level="INFO", log_file="reports/pipeline.log")

    # Get logger for a specific stage
    logger = stage_logger("Stage A")
    logger.info("Processing VCF...")

    # Track stage progress
    with StageProgress("Stage B", "Quality Control"):
        run_qc()
"""

import logging
import sys
from logging.handlers import RotatingFileHandler
from pathlib import Path


class ColorFormatter(logging.Formatter):
    """Colored console formatter."""

    COLORS = {
        logging.DEBUG: "\033[2m",      # dim
        logging.INFO: "\033[0;32m",    # green
        logging.WARNING: "\033[1;33m", # yellow
        logging.ERROR: "\033[0;31m",   # red
        logging.CRITICAL: "\033[1;31m", # bold red
    }
    RESET = "\033[0m"

    def format(self, record):
        color = self.COLORS.get(record.levelno, "")
        record.msg = f"{color}{record.msg}{self.RESET}"
        return super().format(record)


def setup_logging(
    level: str = "INFO",
    log_file: str | None = None,
    max_bytes: int = 5 * 1024 * 1024,
    backup_count: int = 3,
) -> logging.Logger:
    """
    Configure pipeline-wide logging.

    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_file: Optional file path for log output
        max_bytes: Max log file size before rotation
        backup_count: Number of rotated log files to keep

    Returns:
        Root pipeline logger
    """
    root = logging.getLogger("bluegen")
    root.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Clear existing handlers
    root.handlers.clear()

    # Console handler with colors
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(ColorFormatter("%(message)s"))
    root.addHandler(console)

    # File handler with rotation
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = RotatingFileHandler(
            str(log_path), maxBytes=max_bytes, backupCount=backup_count
        )
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")
        )
        root.addHandler(file_handler)

    return root


def stage_logger(name: str) -> logging.Logger:
    """Get a logger for a specific pipeline stage."""
    return logging.getLogger(f"bluegen.{name}")


class StageProgress:
    """Context manager for tracking stage progress."""

    def __init__(self, name: str, description: str = ""):
        self.name = name
        self.description = description
        self.logger = stage_logger(name)

    def __enter__(self):
        self.logger.info(f"\n\033[1;36m═══ {self.name}: {self.description} ═══\033[0m")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            self.logger.info(f"\033[0;32m✓ {self.name} complete\033[0m")
        else:
            self.logger.error(f"\033[0;31m✗ {self.name} failed: {exc_val}\033[0m")
        return False  # Don't suppress exceptions
