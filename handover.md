# Refactoring Handover Plan

- [x] **Step 1: Rename `TranscriptomeParser` → `AnnotationParser`**
  - Rename file: `services/transcriptome_parser.py` → `services/annotation_parser.py`
  - Rename class: `TranscriptomeParser` → `AnnotationParser`
  - Rename schema: `TRANSCRIPTOME_SCHEMA` → `GTF_SCHEMA` in `models.py`
  - Rename method: keep `load_gtf` (it parses both GTF and BED — name is accurate)
  - Update import in `commands/off_targets.py` (3 occurrences)
  - Rename test file: `tests/off-targets/test_transcriptome_parser.py` → `test_annotation_parser.py`
  - Update import in `tests/off-targets/test_transcriptome_format.py`
  - Update `CLAUDE.md` architecture section
  - Update `README.md` mermaid diagram

- [x] **Step 2: Unify Execution Modes — remove `--input-dir` and `--chunk-mode`**
  - Remove `--input-dir` / `-d` CLI flag from `off_targets.py`. The existing redirect
    (lines 446-448) already handles directories passed to `--risearch-file`; this just
    deletes the redundant separate flag and its code branch.
  - Remove `--chunk-mode` / `--batch-size` CLI flags and the corresponding chunk-mode
    code block (~200 lines). Large single-file inputs should use
    `convert_risearch_to_parquet.py` first; chunk-mode is a workaround, not a feature.
  - Remove `input_dir` and `chunk_mode` / `batch_size` fields from `OffTargetsConfig`
    in `config.py` and the `config_to_kwargs` mapping.
  - Run `pytest tests/` to verify nothing is broken after removal.

- [ ] **Step 3: Remove Useless Flags**
  - After Step 2, audit remaining CLI flags for anything now dead:
    - `--profile` — check if used / tested; remove if only wiring is the flag itself.
    - `--output-format` — still live (tsv/csv/parquet in dir mode), keep.
    - `--summary-only` — still live in dir mode, keep.
    - Any `input_dir` YAML config keys in `config/` example files — clean those up.
  - Clean up associated dead code paths in `config.py` and `off_targets.py`.

- [ ] **Step 4: General Code Polish**
  - Simplify the `_init_worker` / `_process_one_sirna` module-level global pattern:
    use a proper dataclass or namedtuple for worker state instead of bare globals.
  - Deduplicate the alpha-gamma pair construction logic (copy-pasted in 2 places now).
  - Remove dead `is_multi_sirna` branch in standard-mode path if chunk-mode removal
    makes it always True; otherwise consolidate.
  - Remove unused imports revealed by `ruff check`.
  - Ensure `ruff format` and `ty check`/`pyrefly check` pass clean.

- [ ] **Step 5: Update the README**
  - Update mermaid diagram: `TranscriptomeParser` → `AnnotationParser`.
  - Update key-options table: remove `--input-dir`, `--chunk-mode`, `--batch-size`.
  - Update any prose referencing "directory mode" or "chunk mode".
  - Ensure usage examples use `--risearch-file` (not `--input-dir`) for directory input.
