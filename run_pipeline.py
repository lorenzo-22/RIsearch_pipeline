#!/usr/bin/env python3
"""
Orchestrator for the RIsearch siRNA off-target discovery pipeline.

Runs pipeline stages in dependency order:
  1. (optional) index       -- risearch-pipeline index
  2. (optional) convert     -- convert_risearch_to_parquet.py
  3.            accessibility -- risearch-pipeline accessibility
  4.            off-targets  -- risearch-pipeline off-targets

Local mode (default): steps run sequentially; logs go to logs/<timestamp>/<step>.log.
Slurm mode (--slurm): each step submitted with sbatch, chained via --dependency=afterok.
Dry-run (--dry-run): prints commands without executing; Slurm dry-run shows full sbatch
  commands with <JOBID_stepname> placeholders so the dependency chain is visible.

Config format (YAML):
  steps: [accessibility, off-targets]   # default; add index/convert to enable them

  # Per-step Slurm resource overrides (merged with built-in defaults)
  slurm:
    partition: batch
    account: myaccount
    accessibility:
      time: "08:00:00"
      mem: 32G
      cpus_per_task: 8
    off_targets:
      time: "04:00:00"
      mem: 64G
      cpus_per_task: 16

  index:
    target: data/genome.fa          # required
    output: data/genome.idx         # optional

  convert:
    input_dir: data/raw_risearch/   # required
    out_dir: data/parquet/          # optional (default: same as input_dir)
    ids_file: data/sirna_ids.txt    # optional
    workers: 32

  accessibility:
    fasta: data/genome.fa           # required
    output: data/accessibility/     # required
    window: 80
    span: 40
    unpaired: 30
    temperature: 37.0
    workers: 8

  off_targets:
    risearch_file: data/parquet/    # required (file or parquet directory)
    transcriptome: data/anno.gtf    # optional
    accessibility_dir: data/accessibility/
    output: results/off_targets.tsv
    alpha: "0.8;1.0"
    gamma: "1.0"
    type: gw
"""

import argparse
import shlex
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import yaml


STEP_ORDER = ["index", "convert", "accessibility", "off-targets"]
DEFAULT_STEPS = ["accessibility", "off-targets"]

# Built-in per-step Slurm resource defaults.
_DEFAULT_RESOURCES: dict[str, dict] = {
    "index":         {"time": "01:00:00", "mem": "8G",  "cpus_per_task": 4},
    "convert":       {"time": "02:00:00", "mem": "16G", "cpus_per_task": 16},
    "accessibility": {"time": "08:00:00", "mem": "32G", "cpus_per_task": 8},
    "off-targets":   {"time": "04:00:00", "mem": "64G", "cpus_per_task": 16},
}


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

def _load_config(path: Path) -> dict:
    return yaml.safe_load(path.read_text()) or {}


def _resolve(value: str | None, base: Path) -> str | None:
    if value is None:
        return None
    p = Path(value)
    return str((base / p).resolve() if not p.is_absolute() else p)


# ---------------------------------------------------------------------------
# Command builders: config section dict → CLI argv list
# ---------------------------------------------------------------------------

def _build_index(cfg: dict, base: Path) -> list[str]:
    target = _resolve(cfg.get("target"), base)
    if not target:
        raise ValueError("'target' is required")
    cmd = ["risearch-pipeline", "index", target]
    output = _resolve(cfg.get("output"), base)
    if output:
        cmd += ["--output", output]
    if cfg.get("verbose"):
        cmd.append("--verbose")
    return cmd


def _build_convert(cfg: dict, base: Path) -> list[str]:
    script = Path(__file__).parent / "convert_risearch_to_parquet.py"
    input_dir = _resolve(cfg.get("input_dir"), base)
    if not input_dir:
        raise ValueError("'input_dir' is required")
    cmd = [sys.executable, str(script), input_dir]
    out_dir = _resolve(cfg.get("out_dir"), base)
    if out_dir:
        cmd += ["--out-dir", out_dir]
    ids_file = _resolve(cfg.get("ids_file"), base)
    if ids_file:
        cmd += ["--ids-file", ids_file]
    if "workers" in cfg:
        cmd += ["--workers", str(cfg["workers"])]
    if cfg.get("skip_existing"):
        cmd.append("--skip-existing")
    return cmd


def _build_accessibility(cfg: dict, base: Path) -> list[str]:
    fasta = _resolve(cfg.get("fasta"), base)
    if not fasta:
        raise ValueError("'fasta' is required")
    output = _resolve(cfg.get("output"), base)
    if not output:
        raise ValueError("'output' is required")
    cmd = ["risearch-pipeline", "accessibility", "--fasta", fasta, "--output", output]
    for key, flag in [("window", "--window"), ("span", "--span"), ("unpaired", "--unpaired"),
                      ("temperature", "--temperature"), ("workers", "--workers")]:
        if key in cfg and cfg[key] is not None:
            cmd += [flag, str(cfg[key])]
    if cfg.get("verbose"):
        cmd.append("--verbose")
    return cmd


def _build_off_targets(cfg: dict, base: Path) -> list[str]:
    cmd = ["risearch-pipeline", "off-targets"]

    # Path options
    for key, flag in [
        ("risearch_file",          "--risearch-file"),
        ("sirna_fasta",            "--sirna-fasta"),
        ("target_fasta",           "--target-fasta"),
        ("target_index",           "--target-index"),
        ("transcriptome",          "--transcriptome"),
        ("accessibility_dir",      "--accessibility-dir"),
        ("output",                 "--output"),
        ("fasta",                  "--fasta"),
        ("on_target",              "--on-target"),
        ("on_target_risearch_file","--on-target-risearch-file"),
        ("query",                  "--query"),
        ("on_target_accessibility","--on-target-accessibility"),
        ("on_target_ids",          "--on-target-ids"),
    ]:
        val = _resolve(cfg.get(key), base)
        if val is not None:
            cmd += [flag, val]

    # Scalar options
    for key, flag in [
        ("workers",              "--workers"),
        ("feature",              "--feature"),
        ("expression_metric",    "--expression-metric"),
        ("transcriptome_format", "--transcriptome-format"),
        ("window",               "--window"),
        ("span",                 "--span"),
        ("unpaired",             "--unpaired"),
        ("temperature",          "--temperature"),
        ("on_target_expression", "--on-target-expression"),
        ("alpha",                "--alpha"),
        ("gamma",                "--gamma"),
        ("theta",                "--theta"),
        ("type",                 "--type"),
        ("output_format",        "--output-format"),
    ]:
        if key in cfg and cfg[key] is not None:
            cmd += [flag, str(cfg[key])]

    # Boolean flags
    for key, flag in [
        ("legacy_format",   "--legacy-format"),
        ("detailed_report", "--detailed-report"),
        ("sense_only",      "--sense-only"),
        ("verbose",         "--verbose"),
        ("profile",         "--profile"),
        ("summary_only",    "--summary-only"),
    ]:
        if cfg.get(key):
            cmd.append(flag)

    return cmd


_BUILDERS = {
    "index":         _build_index,
    "convert":       _build_convert,
    "accessibility": _build_accessibility,
    "off-targets":   _build_off_targets,
}

_CFG_KEYS = {
    "index":         "index",
    "convert":       "convert",
    "accessibility": "accessibility",
    "off-targets":   "off_targets",
}

# Dependency graph: index/convert/accessibility are independent of each other.
# off-targets waits for accessibility (always) and convert (if convert is in the run).
_DEPS: dict[str, list[str]] = {
    "index":         [],
    "convert":       [],
    "accessibility": [],
    "off-targets":   ["accessibility", "convert"],
}


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------

def _run_local(steps: list[str], commands: list[list[str]], log_dir: Path) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    for step, cmd in zip(steps, commands):
        log_file = log_dir / f"{step}.log"
        print(f"\n[{step}] {shlex.join(cmd)}")
        print(f"[{step}] log → {log_file}")
        with log_file.open("w") as fh:
            fh.write(f"CMD: {shlex.join(cmd)}\n\n")
            fh.flush()
            rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
        if rc != 0:
            print(f"[{step}] FAILED (exit {rc}) — see {log_file}", file=sys.stderr)
            sys.exit(rc)
        print(f"[{step}] OK")


def _run_slurm(
    steps: list[str],
    commands: list[list[str]],
    slurm_cfg: dict,
    resources: dict[str, dict],
    log_dir: Path,
    dry_run: bool,
) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    job_ids: dict[str, str] = {}  # step → job ID (real or placeholder)

    for step, cmd in zip(steps, commands):
        res = resources[step]
        sbatch = [
            "sbatch", "--parsable",
            "--time",          res.get("time", "04:00:00"),
            "--mem",           res.get("mem", "16G"),
            "--cpus-per-task", str(res.get("cpus_per_task", 4)),
            "--job-name",      f"rip_{step.replace('-', '_')}",
            "--output",        str(log_dir / f"{step}_%j.out"),
            "--error",         str(log_dir / f"{step}_%j.err"),
        ]
        partition = slurm_cfg.get("partition")
        if partition:
            sbatch += ["--partition", partition]
        account = slurm_cfg.get("account")
        if account:
            sbatch += ["--account", account]

        active_deps = [s for s in _DEPS.get(step, []) if s in steps]
        if active_deps:
            dep_str = ":".join(job_ids[d] for d in active_deps)
            sbatch += ["--dependency", f"afterok:{dep_str}"]

        sbatch += ["--wrap", shlex.join(cmd)]

        if dry_run:
            print(shlex.join(sbatch))
            print()
            job_ids[step] = f"<JOBID_{step.replace('-', '_')}>"
        else:
            result = subprocess.run(sbatch, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"[{step}] sbatch failed: {result.stderr.strip()}", file=sys.stderr)
                sys.exit(result.returncode)
            job_ids[step] = result.stdout.strip().split(";")[0]
            print(f"[{step}] submitted → job {job_ids[step]}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", "-c", required=True, type=Path,
                        metavar="YAML",
                        help="Pipeline orchestrator config file.")
    parser.add_argument("--slurm", action="store_true",
                        help="Submit steps as Slurm jobs.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print commands without executing.")
    parser.add_argument("--steps", metavar="STEP[,STEP...]",
                        help=f"Steps to run (default: {','.join(DEFAULT_STEPS)}). "
                             f"Available: {','.join(STEP_ORDER)}.")
    # Global Slurm resource overrides (applied to every step)
    parser.add_argument("--partition", metavar="NAME",
                        help="Slurm partition (overrides config).")
    parser.add_argument("--account", metavar="NAME",
                        help="Slurm account (overrides config).")
    parser.add_argument("--time", metavar="HH:MM:SS",
                        help="Wall time for every step (overrides per-step defaults).")
    parser.add_argument("--mem", metavar="SIZE",
                        help="Memory for every step, e.g. 32G (overrides per-step defaults).")
    parser.add_argument("--cpus-per-task", type=int, metavar="N",
                        help="CPUs per step (overrides per-step defaults).")
    args = parser.parse_args()

    config_path = args.config.resolve()
    if not config_path.exists():
        sys.exit(f"Config not found: {config_path}")

    cfg = _load_config(config_path)
    base = config_path.parent

    # Resolve steps
    if args.steps:
        steps = [s.strip() for s in args.steps.split(",")]
    else:
        steps = cfg.get("steps", DEFAULT_STEPS)

    unknown = [s for s in steps if s not in STEP_ORDER]
    if unknown:
        sys.exit(f"Unknown step(s): {unknown}. Available: {STEP_ORDER}")

    # Preserve canonical dependency order
    steps = [s for s in STEP_ORDER if s in steps]

    # Build argv lists
    commands: list[list[str]] = []
    for step in steps:
        step_cfg = cfg.get(_CFG_KEYS[step]) or {}
        try:
            commands.append(_BUILDERS[step](step_cfg, base))
        except ValueError as exc:
            sys.exit(f"Config error for step '{step}': {exc}")

    # Log directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = Path(cfg.get("log_dir", "logs")) / timestamp

    # Resolve Slurm config: file → CLI overrides
    slurm_cfg: dict = dict(cfg.get("slurm") or {})
    if args.partition:
        slurm_cfg["partition"] = args.partition
    if args.account:
        slurm_cfg["account"] = args.account

    # Per-step resources: built-in defaults → per-step config block → global CLI overrides
    resources: dict[str, dict] = {}
    for step in steps:
        res = dict(_DEFAULT_RESOURCES.get(step, {}))
        step_key = step.replace("-", "_")
        step_slurm = (slurm_cfg.get(step_key) or {})
        res.update(step_slurm)
        if args.time:
            res["time"] = args.time
        if args.mem:
            res["mem"] = args.mem
        if args.cpus_per_task:
            res["cpus_per_task"] = args.cpus_per_task
        resources[step] = res

    # Execute
    if args.dry_run and not args.slurm:
        for step, cmd in zip(steps, commands):
            print(f"# {step}")
            print(shlex.join(cmd))
            print()
        return

    if args.slurm:
        _run_slurm(steps, commands, slurm_cfg, resources, log_dir, args.dry_run)
    else:
        _run_local(steps, commands, log_dir)


if __name__ == "__main__":
    main()
