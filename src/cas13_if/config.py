"""Configuration loading, resolution, and canonical hashing."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, cast

import yaml


class ConfigError(ValueError):
    """Raised when a configuration cannot be loaded or validated."""


ConfigDict = dict[str, Any]


def load_config(path: Path) -> ConfigDict:
    """Load a YAML mapping and return a recursively resolved plain dictionary."""
    if not path.is_file():
        raise ConfigError(f"configuration does not exist: {path}")
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigError(f"invalid YAML in {path}: {exc}") from exc
    if raw is None:
        return {}
    if not isinstance(raw, dict):
        raise ConfigError(f"configuration root must be a mapping: {path}")
    return _resolve_interpolations(cast(ConfigDict, raw))


def _resolve_interpolations(config: ConfigDict) -> ConfigDict:
    """Resolve the limited ${a.b} interpolation used by project path configs."""
    resolved = cast(ConfigDict, json.loads(json.dumps(config)))
    for _ in range(20):
        changed = _walk_and_resolve(resolved, resolved)
        if not changed:
            break
    unresolved = _find_unresolved(resolved)
    if unresolved:
        names = ", ".join(sorted(unresolved))
        raise ConfigError(f"unresolved configuration interpolation(s): {names}")
    return resolved


def _walk_and_resolve(value: Any, root: ConfigDict) -> bool:
    changed = False
    if isinstance(value, dict):
        for key, item in value.items():
            if isinstance(item, str):
                replacement = _resolve_string(item, root)
                if replacement != item:
                    value[key] = replacement
                    changed = True
            else:
                changed = _walk_and_resolve(item, root) or changed
    elif isinstance(value, list):
        for index, item in enumerate(value):
            if isinstance(item, str):
                replacement = _resolve_string(item, root)
                if replacement != item:
                    value[index] = replacement
                    changed = True
            else:
                changed = _walk_and_resolve(item, root) or changed
    return changed


def _resolve_string(value: str, root: ConfigDict) -> Any:
    start = value.find("${")
    if start < 0:
        return value
    end = value.find("}", start)
    if end < 0:
        raise ConfigError(f"unterminated interpolation: {value}")
    key = value[start + 2 : end]
    replacement: Any = root
    for component in key.split("."):
        if not isinstance(replacement, dict) or component not in replacement:
            return value
        replacement = replacement[component]
    token = value[start : end + 1]
    if value == token:
        return replacement
    if not isinstance(replacement, (str, int, float, bool)):
        raise ConfigError(f"cannot embed non-scalar interpolation {token}")
    return value.replace(token, str(replacement))


def _find_unresolved(value: Any) -> set[str]:
    found: set[str] = set()
    if isinstance(value, str) and "${" in value:
        found.add(value)
    elif isinstance(value, dict):
        for item in value.values():
            found.update(_find_unresolved(item))
    elif isinstance(value, list):
        for item in value:
            found.update(_find_unresolved(item))
    return found


def canonical_json(config: ConfigDict) -> str:
    """Return stable compact JSON for hashing and manifests."""
    return json.dumps(config, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def config_hash(config: ConfigDict, length: int = 10) -> str:
    """Hash resolved configuration deterministically."""
    if length < 7 or length > 64:
        raise ValueError("config hash length must be between 7 and 64")
    return hashlib.sha256(canonical_json(config).encode("utf-8")).hexdigest()[:length]
