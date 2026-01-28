"""Configuration management using OmegaConf."""

from pathlib import Path
from typing import Optional
from dataclasses import dataclass
from omegaconf import OmegaConf, MISSING, DictConfig


@dataclass
class OffTargetsConfig:
    """Schema for off-targets command configuration."""

    risearch_file: str = MISSING  # Required
    transcriptome: Optional[str] = None
    feature: str = "exon"
    expression_metric: str = "RPKM"
    accessibility_dir: Optional[str] = None
    output: Optional[str] = None
    fasta: Optional[str] = None
    window: int = 80
    span: int = 40
    unpaired: int = 30
    on_target: Optional[str] = None
    on_target_risearch_file: Optional[str] = None
    query: Optional[str] = None
    on_target_expression: float = 1000.0
    on_target_accessibility: Optional[str] = None
    legacy_format: bool = False
    detailed_report: bool = False
    sense_only: bool = False
    type: str = "gw"
    verbose: bool = False


@dataclass
class AccessibilityConfig:
    """Schema for accessibility command configuration."""

    fasta: str = MISSING  # Required
    output: str = MISSING  # Required
    window: int = 80
    span: int = 40
    unpaired: int = 30
    verbose: bool = False


@dataclass
class PipelineConfig:
    """Top-level config schema."""

    command: str = MISSING  # "off-targets" or "accessibility"
    off_targets: Optional[OffTargetsConfig] = None
    accessibility: Optional[AccessibilityConfig] = None


def load_config(config_path: Path) -> DictConfig:
    """
    Load and validate a YAML config file.

    Resolves relative paths relative to the config file's directory.

    Args:
        config_path: Path to YAML configuration file.

    Returns:
        Validated OmegaConf DictConfig merged with schema defaults.

    Raises:
        omegaconf.MissingMandatoryValue: If required fields are missing.
        ValueError: If command is invalid.
    """
    # Load schema as structured config
    schema = OmegaConf.structured(PipelineConfig)

    # Load user config
    user_cfg = OmegaConf.load(config_path)

    # Merge (user values override schema defaults)
    cfg = OmegaConf.merge(schema, user_cfg)

    # Validate command
    if cfg.command not in ("off-targets", "accessibility"):
        raise ValueError(
            f"Unknown command: {cfg.command}. Must be 'off-targets' or 'accessibility'."
        )

    # Resolve paths relative to config file directory
    config_dir = config_path.parent.resolve()
    cfg = _resolve_paths(cfg, config_dir)

    return cfg


def _resolve_paths(cfg: DictConfig, base_dir: Path) -> DictConfig:
    """
    Resolve relative paths in config relative to base_dir.

    Only resolves string fields that look like paths (contain / or end with common extensions).
    """
    path_fields = {
        "off_targets": [
            "risearch_file",
            "transcriptome",
            "accessibility_dir",
            "output",
            "fasta",
            "on_target",
            "on_target_risearch_file",
            "query",
            "on_target_accessibility",
        ],
        "accessibility": ["fasta", "output"],
    }

    for section, fields in path_fields.items():
        section_cfg = getattr(cfg, section, None)
        if section_cfg is None:
            continue
        for field_name in fields:
            value = getattr(section_cfg, field_name, None)
            if value is not None and isinstance(value, str):
                p = Path(value)
                if not p.is_absolute():
                    resolved = (base_dir / p).resolve()
                    OmegaConf.update(cfg, f"{section}.{field_name}", str(resolved))

    return cfg


def config_to_kwargs(cfg: DictConfig, command: str) -> dict:
    """
    Convert OmegaConf section to kwargs dict for command function.

    Converts string paths to Path objects for Typer compatibility.
    """
    section = getattr(cfg, command.replace("-", "_"))
    kwargs = OmegaConf.to_container(section, resolve=True)

    # Convert path strings to Path objects
    path_fields = [
        "risearch_file",
        "transcriptome",
        "accessibility_dir",
        "output",
        "fasta",
        "on_target",
        "on_target_risearch_file",
        "query",
        "on_target_accessibility",
    ]
    for key in path_fields:
        if key in kwargs and kwargs[key] is not None:
            kwargs[key] = Path(kwargs[key])

    # Rename 'type' to 'predictions_type' to match function signature
    if "type" in kwargs:
        kwargs["predictions_type"] = kwargs.pop("type")

    # Argument Mapping Logic
    # ----------------------

    key_mapping = {
        "window": "window_size",
        "span": "max_span",
        "unpaired": "unpaired_prob",
    }

    if command == "off-targets":
        # off_targets.run specific mappings
        key_mapping.update(
            {
                "transcriptome": "gtf_file",
                "feature": "feature_type",
                "output": "output_file",
                "fasta": "genome_file",
                "on_target": "on_target_file",
                "query": "query_file",
            }
        )
    elif command == "accessibility":
        # accessibility.run specific mappings
        key_mapping.update(
            {
                "fasta": "genome",
                "output": "output_dir",
            }
        )

    for old_key, new_key in key_mapping.items():
        if old_key in kwargs:
            kwargs[new_key] = kwargs.pop(old_key)

    return kwargs
