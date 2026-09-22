"""Tests for the Slurm submission path of scripts/run_pipeline.py.

Slurm cannot run in CI, so the thing worth pinning is the *sbatch command that
would be submitted*. `--slurm --dry-run` prints exactly that (README documents it
as the way to verify a pipeline before submitting), which makes the whole path
testable without a scheduler: assert on the generated argv, and use dry-run's
`<JOBID_*>` placeholders to check that `--dependency=afterok:` chaining is wired
between the right jobs.

`scripts/` is not an installed package (no console script, not in the wheel), so
the module is loaded by path.
"""

import importlib.util
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
RUN_PIPELINE = REPO_ROOT / "scripts" / "run_pipeline.py"


def _load_run_pipeline():
    spec = importlib.util.spec_from_file_location("run_pipeline", RUN_PIPELINE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["run_pipeline"] = module
    spec.loader.exec_module(module)
    return module


rp = _load_run_pipeline()


def _sbatch_lines(capsys) -> list[list[str]]:
    """Parse captured dry-run stdout back into one argv list per sbatch call."""
    out = capsys.readouterr().out
    return [shlex.split(line) for line in out.splitlines() if line.startswith("sbatch")]


def _flag(argv: list[str], name: str) -> str | None:
    """Value following `name` in an argv list, or None if the flag is absent."""
    return argv[argv.index(name) + 1] if name in argv else None


class TestRunSlurmDryRun:
    def test_emits_one_sbatch_per_job(self, tmp_path, capsys):
        jobs = [
            rp.Job(
                key="accessibility",
                step="accessibility",
                cmd=["sioff", "accessibility"],
            ),
            rp.Job(key="off-targets", step="off-targets", cmd=["sioff", "off-targets"]),
        ]
        resources = {
            s: dict(rp._DEFAULT_RESOURCES[s]) for s in ("accessibility", "off-targets")
        }

        rp._run_slurm(jobs, {}, resources, tmp_path / "logs", dry_run=True)

        assert len(_sbatch_lines(capsys)) == 2

    def test_resource_flags_come_from_the_step_defaults(self, tmp_path, capsys):
        job = rp.Job(
            key="off-targets", step="off-targets", cmd=["sioff", "off-targets"]
        )

        rp._run_slurm(
            [job],
            {},
            {"off-targets": dict(rp._DEFAULT_RESOURCES["off-targets"])},
            tmp_path / "logs",
            dry_run=True,
        )

        (argv,) = _sbatch_lines(capsys)
        assert argv[:2] == ["sbatch", "--parsable"]
        assert _flag(argv, "--time") == "04:00:00"
        assert _flag(argv, "--mem") == "64G"
        assert _flag(argv, "--cpus-per-task") == "16"
        assert _flag(argv, "--job-name") == "rip_off_targets"

    def test_falls_back_to_builtin_defaults_when_resources_are_empty(
        self, tmp_path, capsys
    ):
        """`res.get(..., default)` fallbacks in _run_slurm, not the _DEFAULT_RESOURCES table.

        Every step in _DEFAULT_RESOURCES sets all three keys, so these fallbacks are
        unreachable via the normal path and would otherwise go untested.
        """
        job = rp.Job(key="custom", step="custom", cmd=["sioff", "off-targets"])

        rp._run_slurm([job], {}, {"custom": {}}, tmp_path / "logs", dry_run=True)

        (argv,) = _sbatch_lines(capsys)
        assert _flag(argv, "--time") == "04:00:00"
        assert _flag(argv, "--mem") == "16G"
        assert _flag(argv, "--cpus-per-task") == "4"

    def test_per_step_resource_overrides_win(self, tmp_path, capsys):
        job = rp.Job(
            key="accessibility", step="accessibility", cmd=["sioff", "accessibility"]
        )

        rp._run_slurm(
            [job],
            {},
            {"accessibility": {"time": "12:00:00", "mem": "128G", "cpus_per_task": 32}},
            tmp_path / "logs",
            dry_run=True,
        )

        (argv,) = _sbatch_lines(capsys)
        assert _flag(argv, "--time") == "12:00:00"
        assert _flag(argv, "--mem") == "128G"
        assert _flag(argv, "--cpus-per-task") == "32"

    def test_log_paths_carry_the_slurm_jobid_pattern(self, tmp_path, capsys):
        log_dir = tmp_path / "logs"
        job = rp.Job(
            key="off-targets:human", step="off-targets", cmd=["sioff", "off-targets"]
        )

        rp._run_slurm(
            [job],
            {},
            {"off-targets": dict(rp._DEFAULT_RESOURCES["off-targets"])},
            log_dir,
            dry_run=True,
        )

        (argv,) = _sbatch_lines(capsys)
        # %j is expanded by Slurm, so it must survive shlex.join unmangled.
        assert _flag(argv, "--output") == str(log_dir / "off_targets_human_%j.out")
        assert _flag(argv, "--error") == str(log_dir / "off_targets_human_%j.err")

    def test_partition_and_account_are_passed_through_when_set(self, tmp_path, capsys):
        job = rp.Job(
            key="accessibility", step="accessibility", cmd=["sioff", "accessibility"]
        )
        resources = {"accessibility": dict(rp._DEFAULT_RESOURCES["accessibility"])}

        rp._run_slurm(
            [job],
            {"partition": "batch", "account": "myaccount"},
            resources,
            tmp_path / "logs",
            dry_run=True,
        )

        (argv,) = _sbatch_lines(capsys)
        assert _flag(argv, "--partition") == "batch"
        assert _flag(argv, "--account") == "myaccount"

    @pytest.mark.parametrize("slurm_cfg", [{}, {"partition": "", "account": ""}])
    def test_partition_and_account_are_omitted_when_unset_or_empty(
        self, tmp_path, capsys, slurm_cfg
    ):
        job = rp.Job(
            key="accessibility", step="accessibility", cmd=["sioff", "accessibility"]
        )
        resources = {"accessibility": dict(rp._DEFAULT_RESOURCES["accessibility"])}

        rp._run_slurm([job], slurm_cfg, resources, tmp_path / "logs", dry_run=True)

        (argv,) = _sbatch_lines(capsys)
        assert "--partition" not in argv
        assert "--account" not in argv

    def test_command_is_wrapped_and_quoted(self, tmp_path, capsys):
        job = rp.Job(
            key="off-targets",
            step="off-targets",
            cmd=["sioff", "off-targets", "-o", "/out/has space.tsv"],
        )

        rp._run_slurm(
            [job],
            {},
            {"off-targets": dict(rp._DEFAULT_RESOURCES["off-targets"])},
            tmp_path / "logs",
            dry_run=True,
        )

        (argv,) = _sbatch_lines(capsys)
        # --wrap must be last, and the inner command must survive a shell round-trip
        # with its spaces intact.
        assert argv[-2] == "--wrap"
        assert shlex.split(argv[-1]) == [
            "sioff",
            "off-targets",
            "-o",
            "/out/has space.tsv",
        ]


class TestDependencyChaining:
    def test_first_job_has_no_dependency(self, tmp_path, capsys):
        jobs = [
            rp.Job(
                key="accessibility",
                step="accessibility",
                cmd=["sioff", "accessibility"],
            ),
            rp.Job(
                key="off-targets",
                step="off-targets",
                cmd=["sioff", "off-targets"],
                dep_keys=["accessibility"],
            ),
        ]
        resources = {
            s: dict(rp._DEFAULT_RESOURCES[s]) for s in ("accessibility", "off-targets")
        }

        rp._run_slurm(jobs, {}, resources, tmp_path / "logs", dry_run=True)

        first, second = _sbatch_lines(capsys)
        assert "--dependency" not in first
        assert _flag(second, "--dependency") == "afterok:<JOBID_accessibility>"

    def test_multiple_dependencies_are_colon_joined(self, tmp_path, capsys):
        jobs = [
            rp.Job(key="index", step="index", cmd=["sioff", "index"]),
            rp.Job(
                key="accessibility",
                step="accessibility",
                cmd=["sioff", "accessibility"],
            ),
            rp.Job(
                key="off-targets",
                step="off-targets",
                cmd=["sioff", "off-targets"],
                dep_keys=["index", "accessibility"],
            ),
        ]
        resources = {
            s: dict(rp._DEFAULT_RESOURCES[s])
            for s in ("index", "accessibility", "off-targets")
        }

        rp._run_slurm(jobs, {}, resources, tmp_path / "logs", dry_run=True)

        *_, last = _sbatch_lines(capsys)
        assert (
            _flag(last, "--dependency") == "afterok:<JOBID_index>:<JOBID_accessibility>"
        )

    def test_fanout_groups_do_not_cross_depend(self, tmp_path, capsys):
        """Each transcriptome group chains only within itself."""
        jobs = [
            rp.Job(
                key="accessibility:human",
                step="accessibility",
                cmd=["sioff", "accessibility"],
            ),
            rp.Job(
                key="accessibility:mouse",
                step="accessibility",
                cmd=["sioff", "accessibility"],
            ),
            rp.Job(
                key="off-targets:human",
                step="off-targets",
                cmd=["sioff", "off-targets"],
                dep_keys=["accessibility:human"],
            ),
            rp.Job(
                key="off-targets:mouse",
                step="off-targets",
                cmd=["sioff", "off-targets"],
                dep_keys=["accessibility:mouse"],
            ),
        ]
        resources = {
            s: dict(rp._DEFAULT_RESOURCES[s]) for s in ("accessibility", "off-targets")
        }

        rp._run_slurm(jobs, {}, resources, tmp_path / "logs", dry_run=True)

        _, _, human, mouse = _sbatch_lines(capsys)
        assert _flag(human, "--dependency") == "afterok:<JOBID_accessibility_human>"
        assert _flag(mouse, "--dependency") == "afterok:<JOBID_accessibility_mouse>"


class TestDryRunSubmitsNothing:
    def test_no_subprocess_is_spawned(self, tmp_path, monkeypatch, capsys):
        """The whole point of --dry-run: print, never submit."""

        def _fail(*args, **kwargs):  # pragma: no cover - only runs on regression
            raise AssertionError(f"dry-run must not execute anything, got: {args!r}")

        monkeypatch.setattr(rp.subprocess, "run", _fail)
        job = rp.Job(
            key="accessibility", step="accessibility", cmd=["sioff", "accessibility"]
        )

        rp._run_slurm(
            [job],
            {},
            {"accessibility": dict(rp._DEFAULT_RESOURCES["accessibility"])},
            tmp_path / "logs",
            dry_run=True,
        )

        assert _sbatch_lines(capsys)  # it did print one


class TestEndToEndDryRun:
    """The documented route: `run_pipeline.py --config … --slurm --dry-run`."""

    def test_shipped_example_config_generates_a_chained_two_job_pipeline(
        self, tmp_path
    ):
        result = subprocess.run(
            [
                sys.executable,
                str(RUN_PIPELINE),
                "--config",
                str(REPO_ROOT / "example_yaml" / "run-pipeline.example.yaml"),
                "--slurm",
                "--dry-run",
            ],
            capture_output=True,
            text=True,
            cwd=tmp_path,
        )

        assert result.returncode == 0, result.stderr
        lines = [
            shlex.split(x) for x in result.stdout.splitlines() if x.startswith("sbatch")
        ]
        assert len(lines) == 2

        accessibility, off_targets = lines
        assert _flag(accessibility, "--job-name") == "rip_accessibility"
        assert "--dependency" not in accessibility
        assert _flag(accessibility, "--partition") == "batch"
        assert _flag(accessibility, "--account") == "myaccount"

        assert _flag(off_targets, "--job-name") == "rip_off_targets"
        assert _flag(off_targets, "--dependency") == "afterok:<JOBID_accessibility>"

    def test_dry_run_writes_no_log_files(self, tmp_path):
        """Dry-run must not produce output; it does still mkdir the log directory.

        `_run_slurm` calls `log_dir.mkdir(parents=True, exist_ok=True)` before the
        dry_run branch, so an empty `logs/<timestamp>/` is left behind. That is
        arguably a wart — a dry run touching the filesystem at all — but it is the
        current contract, so it is pinned here rather than silently tolerated.
        """
        subprocess.run(
            [
                sys.executable,
                str(RUN_PIPELINE),
                "--config",
                str(REPO_ROOT / "example_yaml" / "run-pipeline.example.yaml"),
                "--slurm",
                "--dry-run",
            ],
            capture_output=True,
            text=True,
            cwd=tmp_path,
            check=True,
        )

        assert (tmp_path / "logs").is_dir()
        assert list((tmp_path / "logs").iterdir()), (
            "expected a timestamped subdirectory"
        )
        assert not list((tmp_path / "logs").rglob("*.out"))
        assert not list((tmp_path / "logs").rglob("*.err"))
        assert not any(p.is_file() for p in (tmp_path / "logs").rglob("*"))
