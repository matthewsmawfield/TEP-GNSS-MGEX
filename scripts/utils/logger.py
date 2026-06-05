#!/usr/bin/env python3
"""TEP-GNSS-MGEX logging utilities.

Aligned with TEP-GNSS-II (Cairo) logger style:
- Color-coded console output
- Custom levels: PROCESS, SUCCESS, TEST, TITLE
- Real-time flushing
- Memory usage tracking
"""

import sys
import os
import logging
from pathlib import Path
from typing import Optional

# Custom log levels (must be registered before first use)
logging.addLevelName(25, "PROCESS")
logging.addLevelName(26, "SUCCESS")
logging.addLevelName(27, "TEST")
logging.addLevelName(28, "TITLE")

STEP_LOGGER = None


class TEPFormatter(logging.Formatter):
    """Formatter with color support."""

    COLORS = {
        'SUCCESS': '\033[1;32m',
        'WARNING': '\033[1;33m',
        'ERROR': '\033[1;31m',
        'INFO': '\033[0;37m',
        'DEBUG': '\033[0;90m',
        'PROCESS': '\033[0;34m',
        'TEST': '\033[1;35m',
        'TITLE': '\033[1;36m',
        'CRITICAL': '\033[1;41m',
    }
    RESET = '\033[0m'

    def __init__(self, fmt=None, datefmt=None, use_colors=True):
        super().__init__(fmt, datefmt='%H:%M:%S')
        self.use_colors = use_colors

    def format(self, record):
        message = record.getMessage()
        level_mapping = {
            25: ('PROCESS', self.COLORS['PROCESS']),
            26: ('SUCCESS', self.COLORS['SUCCESS']),
            27: ('TEST', self.COLORS['TEST']),
            28: ('TITLE', self.COLORS['TITLE']),
            logging.INFO: ('INFO', self.COLORS['INFO']),
            logging.WARNING: ('WARNING', self.COLORS['WARNING']),
            logging.ERROR: ('ERROR', self.COLORS['ERROR']),
            logging.DEBUG: ('DEBUG', self.COLORS['DEBUG']),
            logging.CRITICAL: ('CRITICAL', self.COLORS['CRITICAL']),
        }
        level_name, color = level_mapping.get(record.levelno, ('INFO', self.COLORS['INFO']))
        timestamp = self.formatTime(record, self.datefmt)
        if self.use_colors:
            return f"{color}[{timestamp}] [{level_name}] {message}{self.RESET}"
        return f"[{timestamp}] [{level_name}] {message}"


class TEPFileFormatter(logging.Formatter):
    """Clean formatter for file output without ANSI color codes."""

    def __init__(self, fmt=None, datefmt=None):
        super().__init__(fmt, datefmt='%H:%M:%S')

    def format(self, record):
        message = record.getMessage()
        level_mapping = {
            25: 'PROCESS', 26: 'SUCCESS', 27: 'TEST', 28: 'TITLE',
            logging.INFO: 'INFO', logging.WARNING: 'WARNING',
            logging.ERROR: 'ERROR', logging.DEBUG: 'DEBUG',
            logging.CRITICAL: 'CRITICAL',
        }
        level_name = level_mapping.get(record.levelno, 'INFO')
        timestamp = self.formatTime(record, self.datefmt)
        return f"[{timestamp}] [{level_name}] {message}"


class TEPLogger:
    """Structured logger with file and console output."""

    LEVEL_MAP = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
        "CRITICAL": logging.CRITICAL,
    }

    def __init__(self, name: str, log_file_path: Optional[Path] = None, level: str = "INFO", reset_log: bool = False):
        # Use a unique internal name to avoid handler sharing between instances
        self._internal_name = f"tep_{name}_{id(self)}"
        self.logger = logging.getLogger(self._internal_name)
        self.logger.setLevel(self.LEVEL_MAP.get(level, logging.INFO))
        # Remove any stale handlers from previous instances with same id
        self.logger.handlers.clear()

        # Console handler with colors
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        try:
            ch.stream.reconfigure(line_buffering=True)
        except AttributeError:
            pass  # Python <3.9 fallback
        ch.setFormatter(TEPFormatter())
        self.logger.addHandler(ch)
        self.logger.propagate = False

        # File handler
        if log_file_path is None:
            log_file_path = Path(__file__).resolve().parents[2] / "logs" / "general_tep_gnss.log"

        self.log_file_path = log_file_path
        log_file_path.parent.mkdir(parents=True, exist_ok=True)
        if reset_log:
            try:
                with open(log_file_path, 'w') as f:
                    f.write("")
            except Exception:
                pass

        fh = logging.FileHandler(log_file_path, mode='a', encoding='utf-8')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(TEPFileFormatter())
        self.logger.addHandler(fh)

    def clear_log(self):
        """Close file handlers, truncate the log file, and reopen fresh handlers."""
        # Close and remove all file handlers
        for handler in list(self.logger.handlers):
            if isinstance(handler, logging.FileHandler):
                handler.close()
                self.logger.removeHandler(handler)
        # Truncate the file
        try:
            with open(self.log_file_path, 'w') as f:
                f.write("")
        except Exception:
            pass
        # Reopen a fresh file handler
        fh = logging.FileHandler(self.log_file_path, mode='a', encoding='utf-8')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(TEPFileFormatter())
        self.logger.addHandler(fh)

    def _flush_all(self):
        for handler in self.logger.handlers:
            if hasattr(handler, 'stream') and hasattr(handler.stream, 'flush'):
                handler.stream.flush()

    def debug(self, msg: str):
        self.logger.debug(msg)

    def info(self, msg: str):
        self.logger.info(msg)

    def warning(self, msg: str):
        self.logger.warning(msg)

    def error(self, msg: str):
        self.logger.error(msg)

    def critical(self, msg: str):
        self.logger.critical(msg)

    def process(self, msg: str):
        self.logger.log(25, msg)
        self._flush_all()

    def success(self, msg: str):
        self.logger.log(26, msg)
        self._flush_all()

    def test(self, msg: str):
        self.logger.log(27, msg)
        self._flush_all()

    def title(self, msg: str):
        self.logger.log(28, msg)
        self._flush_all()


def set_step_logger(logger: TEPLogger):
    global STEP_LOGGER
    STEP_LOGGER = logger


def print_status(msg: str, level: str = "INFO"):
    """Print and optionally log a status message."""
    if STEP_LOGGER:
        if level == "SUCCESS":
            STEP_LOGGER.success(msg)
        elif level == "ERROR":
            STEP_LOGGER.error(msg)
        elif level == "WARNING":
            STEP_LOGGER.warning(msg)
        elif level == "PROCESS":
            STEP_LOGGER.process(msg)
        elif level == "DEBUG":
            STEP_LOGGER.debug(msg)
        elif level == "TITLE":
            STEP_LOGGER.title("=" * 80)
            STEP_LOGGER.title(msg)
            STEP_LOGGER.title("=" * 80)
        else:
            STEP_LOGGER.info(msg)
        # Ensure file handler flushes so logs are not lost on crash
        if STEP_LOGGER and STEP_LOGGER.logger.handlers:
            for handler in STEP_LOGGER.logger.handlers:
                if isinstance(handler, logging.FileHandler):
                    handler.flush()
    else:
        import datetime
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")
        print(f"[{timestamp}] [{level}] {msg}")
    sys.stdout.flush()


def print_table(headers, rows, title=None):
    """Print a simple ASCII table."""
    if title:
        print_status(f"\n{title}")
        print_status("=" * len(title))
    col_widths = [max(len(str(h)), max(len(str(r[i])) for r in rows)) for i, h in enumerate(headers)]
    header_line = " | ".join(h.ljust(col_widths[i]) for i, h in enumerate(headers))
    print_status(header_line)
    print_status("-" * len(header_line))
    for row in rows:
        print_status(" | ".join(str(row[i]).ljust(col_widths[i]) for i in range(len(row))))
    print_status("")


def check_memory_usage(context: str = "Unknown") -> None:
    """Log current memory usage."""
    try:
        import psutil
        import gc
        process = psutil.Process(os.getpid())
        mem_info = process.memory_info()
        rss_mb = mem_info.rss / (1024 * 1024)
        vms_mb = mem_info.vms / (1024 * 1024)
        msg = f"Memory usage in {context}: RSS={rss_mb:.2f} MB, VMS={vms_mb:.2f} MB"
        if STEP_LOGGER:
            STEP_LOGGER.debug(msg)
        else:
            print(f"[DEBUG] {msg}")
        gc.collect()
    except ImportError:
        pass
