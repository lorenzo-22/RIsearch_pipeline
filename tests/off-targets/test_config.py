"""Tests for config module."""

import pytest
from omegaconf import MissingMandatoryValue
from riot.config import load_config, config_to_kwargs


class TestLoadConfig:
    def test_load_off_targets_config(self, tmp_path):
        config_file = tmp_path / "config.yaml"
        # Create dummy path for risearch_file
        (tmp_path / "test.tsv").touch()

        config_file.write_text("""
command: off-targets
off_targets:
  risearch_file: test.tsv
  output: results.tsv
""")
        cfg = load_config(config_file)
        assert cfg.command == "off-targets"
        # Path resolution check
        assert cfg.off_targets.risearch_file == str(tmp_path / "test.tsv")
        assert cfg.off_targets.output == str(tmp_path / "results.tsv")

    def test_load_accessibility_config(self, tmp_path):
        config_file = tmp_path / "config.yaml"
        config_file.write_text("""
command: accessibility
accessibility:
  fasta: genome.fa
  output: output_dir
""")
        cfg = load_config(config_file)
        assert cfg.command == "accessibility"
        assert cfg.accessibility.fasta == str(tmp_path / "genome.fa")

    def test_missing_required_field(self, tmp_path):
        config_file = tmp_path / "config.yaml"
        config_file.write_text("""
command: off-targets
off_targets: {}
""")
        # load_config triggers validation when resolving paths, so error is raised there
        with pytest.raises(MissingMandatoryValue):
            load_config(config_file)

    def test_invalid_command(self, tmp_path):
        config_file = tmp_path / "config.yaml"
        config_file.write_text("command: invalid")
        with pytest.raises(ValueError, match="Unknown command"):
            load_config(config_file)

    def test_config_to_kwargs_mappings(self, tmp_path):
        config_file = tmp_path / "config.yaml"
        config_file.write_text("""
command: off-targets
off_targets:
  risearch_file: test.tsv
  transcriptome: ann.gtf
  feature: CDS
""")
        cfg = load_config(config_file)
        kwargs = config_to_kwargs(cfg, "off-targets")

        assert "gtf_file" in kwargs
        assert kwargs["gtf_file"] == tmp_path / "ann.gtf"
        assert kwargs["feature_type"] == "CDS"
        assert "transcriptome" not in kwargs

    def test_config_to_kwargs_accessibility(self, tmp_path):
        config_file = tmp_path / "config.yaml"
        config_file.write_text("""
command: accessibility
accessibility:
  fasta: genome.fa
  output: out
""")
        cfg = load_config(config_file)
        kwargs = config_to_kwargs(cfg, "accessibility")

        assert "genome" in kwargs
        assert kwargs["genome"] == tmp_path / "genome.fa"
        assert "output_dir" in kwargs
        assert "fasta" not in kwargs
