"""Configuration management using OmegaConf."""

from pathlib import Path
from typing import Optional
from dataclasses import dataclass
from omegaconf import OmegaConf, MISSING, DictConfig


@dataclass
class OffTargetsConfig:
    risearch_file: str = MISSING
    input_dir: Optional[str] = None
    transcriptome: Optional[str] = None
    transcriptome_format: str = "auto"
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
    alpha: str = "1.0"
    gamma: str = "1.0"
    theta: str = ""
    legacy_format: bool = False
    detailed_report: bool = False
    sense_only: bool = False
    type: str = "gw"
    verbose: bool = False
    chunk_mode: bool = False
    batch_size: int = 50


@dataclass
class AccessibilityConfig:
    fasta: str = MISSING
    output: str = MISSING
    window: int = 80
    span: int = 40
    unpaired: int = 30
    verbose: bool = False


@dataclass
class OrthologsConfig:
    target_gene: str = MISSING
    species_list: str = MISSING
    output: str = MISSING
    email: str = MISSING
    insecta_level: int = 50557
    plot_lengths: bool = False
    run_msa: bool = False


@dataclass
class PipelineConfig:
    command: str = MISSING  # "off-targets", "accessibility", or "orthologs"
    verbose: bool = False
    off_targets: Optional[OffTargetsConfig] = None
    accessibility: Optional[AccessibilityConfig] = None
    orthologs: Optional[OrthologsConfig] = None


def load_config(config_path: Path) -> DictConfig:
    """Load and validate a YAML config file.

    Resolves relative paths relative to the config file's directory.
    """
    schema = OmegaConf.structured(PipelineConfig)
    cfg = OmegaConf.merge(schema, OmegaConf.load(config_path))

    if cfg.command not in ("off-targets", "accessibility", "orthologs"):
        raise ValueError(
            f"Unknown command: {cfg.command}. Must be 'off-targets', 'accessibility', or 'orthologs'."
        )

    return _resolve_paths(cfg, config_path.parent.resolve())


def _resolve_paths(cfg: DictConfig, base_dir: Path) -> DictConfig:
    path_fields = {
        "off_targets": [
            "risearch_file", "transcriptome", "accessibility_dir",
            "output", "fasta", "on_target", "on_target_risearch_file",
            "query", "on_target_accessibility",
        ],
        "accessibility": ["fasta", "output"],
        "orthologs": ["target_gene", "species_list", "output"],
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
                    OmegaConf.update(cfg, f"{section}.{field_name}", str((base_dir / p).resolve()))

    return cfg


def config_to_kwargs(cfg: DictConfig, command: str) -> dict:
    """Convert OmegaConf section to kwargs dict for command function."""
    section = getattr(cfg, command.replace("-", "_"))
    kwargs = OmegaConf.to_container(section, resolve=True)

    path_fields = [
        "risearch_file", "transcriptome", "accessibility_dir", "output",
        "fasta", "on_target", "on_target_risearch_file", "query",
        "on_target_accessibility", "target_gene", "species_list",
    ]
    for key in path_fields:
        if key in kwargs and kwargs[key] is not None:
            kwargs[key] = Path(kwargs[key])

    if "type" in kwargs:
        kwargs["predictions_type"] = kwargs.pop("type")

    key_mapping = {
        "window": "window_size",
        "span": "max_span",
        "unpaired": "unpaired_prob",
    }

    if command == "off-targets":
        key_mapping.update({
            "transcriptome": "gtf_file",
            "feature": "feature_type",
            "output": "output_file",
            "fasta": "genome_file",
            "on_target": "on_target_file",
            "query": "query_file",
        })
        for arg in ["sirna_fasta", "target_fasta", "target_index", "workers", "on_target_ids_file"]:
            if arg not in kwargs:
                kwargs[arg] = None
    elif command == "accessibility":
        key_mapping.update({"fasta": "genome", "output": "output_dir"})
    elif command == "orthologs":
        key_mapping.update({
            "target_gene": "target_gene_file",
            "species_list": "species_list_file",
            "output": "output_dir",
        })

    for old_key, new_key in key_mapping.items():
        if old_key in kwargs:
            kwargs[new_key] = kwargs.pop(old_key)

    return kwargs
