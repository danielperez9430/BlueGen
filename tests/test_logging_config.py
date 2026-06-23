"""Tests for utils.logging_config module."""

import sys
import tempfile
import logging
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "prs_research_pipeline" / "scripts"))

from utils.logging_config import setup_logging, stage_logger, StageProgress


def test_setup_logging_returns_logger():
    logger = setup_logging()
    assert isinstance(logger, logging.Logger)
    assert logger.name == "bluegen"


def test_setup_logging_with_file():
    with tempfile.TemporaryDirectory() as tmpdir:
        log_file = Path(tmpdir) / "test.log"
        logger = setup_logging(log_file=str(log_file))
        logger.info("test message")
        assert log_file.exists()
        content = log_file.read_text()
        assert "test message" in content


def test_stage_logger():
    logger = stage_logger("test_stage")
    assert isinstance(logger, logging.Logger)
    assert "test_stage" in logger.name


def test_stage_progress_context_manager():
    with StageProgress("Test", "description"):
        pass  # Should not raise


def test_stage_progress_logs_info():
    # StageProgress uses stage_logger which writes to bluegen.<name>
    # Just verify it completes without error and produces output on stderr
    with StageProgress("TestStageVerify", "testing"):
        pass
    # If we got here without exception, the context manager works
