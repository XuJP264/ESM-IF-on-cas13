from pathlib import Path

import pytest

from cas13_if.config import ConfigError, config_hash, load_config


def test_config_resolution_and_hash(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text("root: /x\nnested:\n  path: ${root}/data\n", encoding="utf-8")
    first = load_config(path)
    second = load_config(path)
    assert first["nested"]["path"] == "/x/data"
    assert config_hash(first) == config_hash(second)


def test_config_rejects_unresolved_and_non_mapping(tmp_path: Path) -> None:
    unresolved = tmp_path / "unresolved.yaml"
    unresolved.write_text("path: ${missing.value}\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="unresolved"):
        load_config(unresolved)
    sequence = tmp_path / "sequence.yaml"
    sequence.write_text("- one\n- two\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="mapping"):
        load_config(sequence)


def test_config_hash_length_validation() -> None:
    with pytest.raises(ValueError, match="between"):
        config_hash({}, length=2)
