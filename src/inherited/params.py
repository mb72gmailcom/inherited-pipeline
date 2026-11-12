"""Write ``params.json`` files that record how an output directory was produced."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from inherited import __version__


def resolved(path: Path | None) -> str | None:
    """Return an absolute path string, or ``None`` when no path was given."""
    return str(path.resolve()) if path is not None else None


def write_params(path: Path, payload: dict[str, Any]) -> Path:
    """Write ``payload`` plus ``package_version`` and ``run_at`` to ``path``."""
    path.parent.mkdir(parents=True, exist_ok=True)
    record = {
        **payload,
        "package_version": __version__,
        "run_at": datetime.now(timezone.utc).isoformat(),
    }
    with path.open("w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
    return path
