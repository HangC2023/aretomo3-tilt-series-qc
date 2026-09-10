# Known issues / findings log

Code-level bugs and gaps found while using this pipeline, logged here as
they're found (separate from any specific project's own processing notes).
Not a substitute for GitHub issues on anything nontrivial — small/trivial
fixes are just made directly and logged here for the record.

## 2026-09-09 — first pre-gain-corrected MRC dataset (project bi30960_6)

Context: bi30960_6 is the first dataset run through this pipeline that is
(a) raw MRC movies rather than TIFF/EER, and (b) already gain-corrected.
Neither condition had been exercised before, so this was a natural place
to find gaps. Processing notes for that project itself live in its own
`NOTES.md`, not here.

### Fixed

1. **Missing `tqdm` dependency.** `commands/run_aretomo3.py` does
   `from tqdm import tqdm` at module level; `cli.py` imports every command
   eagerly, so the *entire* CLI failed to import in an env without `tqdm`
   installed (not declared in `environment.yml` or `pyproject.toml`).
   Fixed: added to both files, `pip install tqdm`'d into the existing
   `aretomo3-preprocess` env directly.

2. **`--gain` unconditionally required for `--cmd 0`, no escape hatch.**
   `run_aretomo3.py::_preflight()` (~line 1322, pre-fix):
   ```python
   if args.gain is not None:
       if not Path(args.gain).exists():
           errors.append(...)
   elif args.cmd == 0:
       errors.append('--gain is required for --cmd 0 (motion correction mode)')
   ```
   This lands in `errors`, not `warnings` -> unconditional `sys.exit(1)`
   at the end of `run()`. Confirmed `--force` only bypasses `warnings`,
   never `errors` -- there was no existing way to run `--cmd 0` without a
   gain reference, even though the actual command-building code already
   handled `args.gain is None` correctly (just omits `-Gain` from the
   AreTomo3 command line). The bug was purely in this preflight check,
   which assumed gain correction is always needed for cmd=0 -- never
   exercised before since every prior dataset needed gain correction.

   Fixed: added a `--no-gain` boolean flag (mutually exclusive with
   `--gain`). When passed, the preflight skips the "gain required" error
   for `--cmd 0`. Also fixed a related cosmetic bug this exposed: the
   "Gain check: skipped" print unconditionally said `(no gain — cmd != 0)`
   even when the real reason was `--no-gain` with `cmd==0`; now reports
   the actual reason. Verified against real project data via `--dry-run`
   -- preflight now passes cleanly (Gain check, PixelSpacing, Voltage,
   FmInt all green).

### Checked, not bugs

3. `--fm-int` vs. movie-type preflight check (~line 1529) only
   special-cases `tif`/`tiff` and `eer` explicitly; any other extension
   (including `mrc`) falls through to a plain "FmInt: OK" print with no
   dedicated validation. Harmless in practice since the default
   (`--fm-int 1`) is the correct value for already-rendered-frame movies
   (MRC, same as TIFF) -- right answer by default, not because of an
   actual MRC-aware check. Worth tightening later (explicit `mrc` case
   in the same style as `tif`/`tiff`) but not urgent.

### Found, not yet fixed

4. **`shared/discovery.py::mrc_pixel_size()` gives silently wrong results
   on this dataset.** It computes pixel size as `cella.x / nx`, documented
   as relying on "nx and mx agree in every file this codebase produces,
   confirmed against real data" -- that assumption does NOT hold for
   bi30960_6's raw movie MRC headers: `mrc_pixel_size()` returns 50000 Å
   for files whose real pixel size is 2.67 Å (confirmed by IMOD's own
   `header` tool separately, which reports a placeholder 1.0 Å for the
   same files -- neither header-derived path is usable for this dataset).
   Not yet clear whether this is specific to raw camera MRC movies (vs.
   AreTomo3's own MRC *output*, e.g. tomogram volumes, which is what most
   current callers of `mrc_pixel_size()` actually use it for --
   `ctf_handedness.py`, `imod_mtffilter.py`, `pytom_ribo_auto.py`,
   `gapstop_match.py`). Needs checking against a real AreTomo3 output
   volume before deciding whether this needs a fix or just a documented
   caveat for raw-movie use. Not blocking bi30960_6 processing since
   pixel size there comes from `--calibrated-apix`/user input, not this
   helper.

5. **`shared/discovery.py::parse_fraction_filename()` is TIFF-only.**
   Expected pattern `..._NNN_TILT_YYYYMMDD_HHMMSS_fractions.tif[f]` doesn't
   match this (or presumably any future) MRC/bracket-tilt-angle naming
   convention. Not fixed -- bi30960_6 used a standalone script
   (`build_mdoc.py`, in that project's own directory, not this repo) to
   build mdocs for this dataset instead of extending this pipeline's own
   filename parsing. Worth reconsidering if a second MRC-based dataset
   shows up -- at that point a real second regex/parser in
   `shared/discovery.py` (alongside the existing one, per this repo's own
   "three mdoc parsers exist on purpose" precedent) would be justified;
   one-off scripts per dataset naming convention won't scale.

## Not committed/pushed yet

Fixes 1 and 2 above are local edits in this working tree, not committed.
