# CI for RIOT — GitHub Actions (lint + test)

**Date:** 2026-07-27
**Status:** approved design

## Goal

Add Continuous Integration: on every pull request and every push to `main`, a
fresh cloud runner installs the project and runs the linter and the test suite.
Failures show as a red check on the PR, catching problems before merge.

CD (publishing to PyPI) is **out of scope** — it is blocked anyway because the
`risearch` dependency is a `git+ssh` direct URL that PyPI rejects. Revisit once
`risearch` is published to PyPI.

## The core constraint: the `risearch` dependency

`risearch` is a **private** GitHub repo, fetched over **SSH**, built from
**Rust/PyO3**. A CI runner cannot install it without an SSH deploy key *and* a
Rust toolchain. We chose **not** to give CI that access.

Only 1 test file (`tests/off-targets/test_risearch_service.py`) and 2 methods of
`src/riot/services/risearch_service.py` (`index_target`, `run_search`) actually
use it. Everything else — off-targets, accessibility, parsing, probability,
intersection — does not.

Today, however, `import riot` transitively imports `risearch` at module top
(`risearch_service.py:13` ← `core/risearch.py` ← `api.py`), so the package is a
hard requirement even for code paths that never call RIsearch. The README claims
the opposite ("imported lazily"). We fix that.

## Design

### Enabling code changes

1. **Lazy import** — in `src/riot/services/risearch_service.py`, remove the
   top-level `import risearch` and import it locally inside `index_target()` and
   `run_search()` (the only two users). Result: `import riot` and all non-RIsearch
   code work with `risearch` absent.

2. **Optional extra** — in `pyproject.toml`, move `risearch` out of
   `[project].dependencies` into
   `[project.optional-dependencies]` under a `risearch` extra. Relock `uv.lock`.
   - `uv sync` → core install, no SSH/Rust needed (this is what CI runs).
   - `uv sync --extra risearch` → full install for `index`/`search` users.

3. **Test guard** — top of `tests/off-targets/test_risearch_service.py`:
   `pytest.importorskip("risearch")`. The file's tests skip cleanly when the
   package is absent; they still run locally where it is installed.

4. **README** — update the Installation section to document `uv sync` (core) vs
   `uv sync --extra risearch` (full), and correct the lazy-import description so it
   is now accurate.

### CI workflow — `.github/workflows/ci.yml`

Triggers: `pull_request` (any branch) and `push` to `main`.

Two independent jobs, run in parallel on `ubuntu-latest`:

| Job | Steps |
|-----|-------|
| `lint` | checkout → `astral-sh/setup-uv` → `uv sync --frozen` → `uv run ruff check src/ tests/` → `uv run ruff format --check src/ tests/` |
| `test` | checkout → `astral-sh/setup-uv` (Python 3.14) → `uv sync --frozen` → `uv run pytest -q` |

- `--frozen` uses the committed `uv.lock` exactly (reproducible, no re-resolution;
  never touches the private git source since `risearch` is now an unrequested extra).
- risearch-dependent tests auto-skip via the `importorskip` guard.
- Python version: **3.14 only** (the project's required minimum).

### Deferred (YAGNI)

- **Type-check in CI** — repo has `ty` + `pyrefly` as dev deps. Deferred: type
  checkers on Python 3.14 + these libraries may flag many pre-existing issues and
  paint CI red on day one. Add in a follow-up after confirming type-cleanliness.
- **CD / PyPI publish** — blocked by the git+ssh dep; add when `risearch` is on PyPI.
- **Python version matrix** — single version (3.14) for now.

## Verification

- Locally simulate CI's environment: install without the `risearch` extra and
  confirm `import riot` works and `pytest` passes with `test_risearch_service.py`
  skipped.
- Confirm `uv run ruff check` / `ruff format --check` pass on the current tree
  (fix formatting if not).
- After merge, confirm the Actions run shows green lint + test on the PR.

## Success criteria

- A PR shows two passing checks (lint, test) with no manual steps.
- `uv sync` (no extra) succeeds on a machine with no SSH access to the private
  `risearch` repo and no Rust toolchain.
- `import riot` succeeds without `risearch` installed.
