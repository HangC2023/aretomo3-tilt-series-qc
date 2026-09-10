"""
Tests for validate-mdoc's SubFramePath-slash check.

Root cause (confirmed via strace + reading AreTomo3's C++ source directly,
project bi30960_6, 2026-09-09): CReadMdoc::mExtractFramePath() has a
genuine off-by-one bug in its "strip leading whitespace" loop -- for
exactly one leading space (the normal "Field = value" mdoc format), the
computed start index is the space's own position, not the position after
it, so the strip is a silent no-op and the extracted SubFramePath value
keeps ONE stray leading space. This is invisible with real SerialEM
mdocs, where SubFramePath is a full Windows path and
CReadMdoc::GetFrameFileName() strips everything up to the LAST
backslash/slash -- discarding the stray leading space as a side effect.
A bare basename (no slash at all) has nothing for that stripping logic
to find, so GetFrameFileName() returns the whole string, leading space
included, and the resulting open() call fails with ENOENT. See
CLAUDE.md and build_mdoc.py (project bi30960_6) for the full story --
this was the actual root cause of a "parses fine, then intermittently
fails to open files partway through a run" failure mode that looked
like a GPU memory or AreTomo3-version regression before it was traced
all the way down.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from aretomo3_preprocess.commands.validate_mdoc import validate_file


def _section(zval, tilt, subframe_path, dose=4.16, exptime=1.0):
    # Trailing DateTime line matters -- see test_validate_mdoc_frames.py's
    # _section docstring: _simulate_aretomo3 commits a section only on the
    # line AFTER ExposureTime, so without a filler line here the next
    # section's own "[ZValue = N]" marker gets silently consumed instead.
    return (
        f'[ZValue = {zval}]\n'
        f'TiltAngle = {tilt}\n'
        f'ExposureDose = {dose}\n'
        f'SubFramePath = {subframe_path}\n'
        f'ExposureTime = {exptime}\n'
        f'DateTime = 01-Jan-2026  00:00:00\n'
    )


def _write_mdoc(path: Path, subframe_paths):
    n = len(subframe_paths)
    text = ''.join(_section(i, float(i), p) for i, p in enumerate(subframe_paths))
    path.write_text(text)


# _MIN_TILTS is 7 -- use 7 sections so these mdocs are otherwise fully valid.
_BASENAMES = [f'Position_1_{i:03d}_0.0_fractions.mrc' for i in range(7)]


def test_bare_basename_subframepath_fails():
    """A SubFramePath with no directory separator at all is flagged and
    fails validation -- it's a guaranteed AreTomo3 open() failure, not
    just a style nit."""
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        mdoc_path = Path(tmp) / 'Position_1.mdoc'
        _write_mdoc(mdoc_path, _BASENAMES)
        r = validate_file(str(mdoc_path))
        assert r['success'] is False
        assert r['no_slash_paths'] == _BASENAMES
        assert any('directory separator' in issue for issue in r['issues'])


def test_unix_path_subframepath_passes():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        mdoc_path = Path(tmp) / 'Position_1.mdoc'
        paths = [f'/net/data/acquisition/{f}' for f in _BASENAMES]
        _write_mdoc(mdoc_path, paths)
        r = validate_file(str(mdoc_path))
        assert r['success'] is True
        assert r['no_slash_paths'] == []


def test_windows_unc_path_subframepath_passes():
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        mdoc_path = Path(tmp) / 'Position_1.mdoc'
        paths = [rf'\\GATANCUSTOMER\DoseFractions\{f}' for f in _BASENAMES]
        _write_mdoc(mdoc_path, paths)
        r = validate_file(str(mdoc_path))
        assert r['success'] is True
        assert r['no_slash_paths'] == []


def test_dummy_prefix_subframepath_passes():
    """The documented fix for a synthetic mdoc: any real path segment
    before the basename is enough, it doesn't need to be a real directory
    (AreTomo3 ignores SubFramePath's directory component entirely)."""
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        mdoc_path = Path(tmp) / 'Position_1.mdoc'
        paths = [f'src/{f}' for f in _BASENAMES]
        _write_mdoc(mdoc_path, paths)
        r = validate_file(str(mdoc_path))
        assert r['success'] is True
        assert r['no_slash_paths'] == []
