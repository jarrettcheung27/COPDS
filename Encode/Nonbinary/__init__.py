from __future__ import annotations

from pathlib import Path
import importlib
import sys

_THIS_DIR = Path(__file__).resolve().parent

# Prefer local prebuilt extension, fall back to build output if present.
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

_build_dir = _THIS_DIR / "build" / f"lib.win-amd64-cpython-{sys.version_info.major}{sys.version_info.minor}"
if _build_dir.exists() and str(_build_dir) not in sys.path:
    sys.path.insert(0, str(_build_dir))

fftqspa = importlib.import_module("fftqspa")

__all__ = ["fftqspa"]
