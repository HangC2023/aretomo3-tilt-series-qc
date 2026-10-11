#!/usr/bin/env python3
"""auto_threshold_mask.py -- derive a mask threshold from the map itself
instead of guessing a fixed constant.

Built 2026-09-15, Sean's proposal: "find a threshold that is the lowest
you can go and still have a single density (no dust or isolated blobs)...
or using the voxel histogram go as close as you get until you see the
curve start to increase and get the largest connected blob from that."

These are the same idea: as threshold T decreases from the map's max, the
count of voxels above T rises slowly while T is still in the structure's
own density range (sparse, compact), then rises sharply once T crosses
into the background/noise floor (background voxels vastly outnumber
structure voxels) -- and that's exactly the point where a second,
disconnected blob (noise/dust) first appears alongside the main
structure. So: sweep T downward, track connected-component count
(scipy.ndimage.label, full 26-connectivity), stop at the last T where the
whole map still has exactly one component. Take the largest connected
component there as the binary core (safety net in case dust and
structure appear in the same step).

Rationale for using the built TEMPLATE (not the raw map) as input: masks
are applied in the template's own frame (same box/centering/apix as the
picking template), so thresholding there keeps everything in one
consistent frame, and matches what build_reference_set.sh already
produces.

This is a genuine alternative to scripts/pdb_to_reference_mask.sh's
molmap-based masks -- no fitted-PDB dependency at all, derives entirely
from the map's own noise floor. Cross-check the two methods against each
other (Dice/Jaccard overlap) rather than trusting either exclusively.

Usage:
  auto_threshold_mask.py <template.mrc> <output_dir> <prefix> \
      [n_steps=300] [extend_inimask=2] [width_soft_edge=4]
"""
import sys

import mrcfile
import numpy as np
from scipy.ndimage import label


def find_auto_threshold(data, n_steps=1000, start_percentile=None,
                         end_percentile=50.0, max_components=20,
                         max_volume_fraction=0.25, min_blob_voxels=50):
    """start_percentile=None (default, recommended) auto-derives the
    sweep's starting point instead of using a fixed percentile: the
    midpoint, in raw intensity, between the histogram's solvent peak
    (its mode -- the dominant background bin) and a robust stand-in for
    the map's brightest real signal (the 99.5th percentile, NOT
    data.max() itself).

    Replaces an earlier fixed start_percentile=99.5 default (2026-09-15),
    which was not universally safe. Failure mode #1 below (isolated
    atomic-scale peaks not yet merged) can extend further into the
    percentile range for some maps than others -- hit this for real on
    the 80S reference (EMD-4263, bi30960_6, 2026-10-02): at the very
    first sweep point (the 99.5th-percentile threshold itself) there were
    already 29 disconnected components, exceeding max_components(20) on
    the FIRST iteration, so the old code's hard break fired immediately
    with last_good_core still None -- RuntimeError, no mask produced.
    40S/60S (EMD-14317/EMD-18765) worked fine at the fixed 99.5 default;
    80S's safe window was only [98.0, 99.0].

    The fix (Sean, same session): rather than hand-pick a percentile per
    map, derive the starting VALUE from the histogram's own bimodal shape
    -- halfway between the solvent peak and the signal's own brightest
    range. Tried halfway-to-data.max() first; failed on 80S (106
    components) because a handful of extreme single-voxel outliers (e.g.
    an ordered ion) pull the literal max up, biasing "halfway" back into
    the same fragmented-peak zone. Using the 99.5th percentile instead of
    the literal max as the "brightest signal" anchor fixed it: validated
    against all three of this project's real reference maps --
    80S: 2 components, 40S: 1, 60S: 1, at the computed halfway point. A
    full sweep from that auto-derived starting point reproduced the
    80S fix's manually-found result almost exactly (5414 vs 5412 core
    voxels from a one-off manual start_percentile=99.0 override).

    An explicit start_percentile can still be passed to bypass
    auto-derivation (e.g. to reproduce a specific earlier run).

    Sweep threshold down through the value range between
    `start_percentile` and `end_percentile`; return the lowest threshold
    for which the volume is still a single, reasonably-sized connected
    component (full 26-connectivity), plus a binary mask of that
    component.

    Design, tuned against real data (caught these failure modes
    empirically, not just reasoned about -- see git history of this file
    for the two earlier broken versions):

    1. Don't start the sweep at data.max(). Right at the very top of the
       range, isolated atomic-scale density peaks (individual
       backbone/sidechain features) are disconnected from EACH OTHER
       before merging into one solid blob further down -- real, benign,
       and irrelevant to mask-building. Starting at `start_percentile`
       (default 99.5) skips straight past this region. (Sean's
       suggestion: "you might have to start closer to the solvent peak
       as at high thresholds it might start with disconnected density.")

    2. Don't sweep across the FULL value range with percentile-spaced
       steps -- percentile spacing gives roughly equal voxel-count
       increments per step BY CONSTRUCTION, which washes out exactly the
       signal being looked for (the rate voxel count grows as threshold
       drops). Instead take LINEAR steps in raw intensity between
       `start_percentile` and `end_percentile` (default: down to the
       median, comfortably past the solvent/background peak) -- narrow
       enough that linear spacing stays informative.

    3. Stop the sweep -- and use the immediately preceding threshold --
       as soon as EITHER: the connected-component count exceeds
       `max_components` (dust appearing; Sean's suggestion: "sweep down
       until you get like 20 blobs"), OR the foreground volume fraction
       exceeds `max_volume_fraction` (guards against the background
       noise floor percolating into one giant connected blob, which is
       trivially "1 component" but meaningless -- this exact failure
       mode broke an earlier version of this function: it found T deep
       in the background where 99.99% of the box was "connected").
       `min_blob_voxels` is a matching lower-bound safety net (Sean:
       "set a size for the biggest blob that it must pass before
       stopping") so a spuriously-early stop can't return a near-empty
       mask.
    """
    if start_percentile is None:
        counts, edges = np.histogram(data, bins=500)
        centers = (edges[:-1] + edges[1:]) / 2
        solvent_peak = float(centers[np.argmax(counts)])
        robust_max = float(np.percentile(data, 99.5))
        start_T = (solvent_peak + robust_max) / 2
    else:
        start_T = float(np.percentile(data, start_percentile))
    end_T = float(np.percentile(data, end_percentile))
    thresholds = np.linspace(start_T, end_T, n_steps)
    structure = np.ones((3, 3, 3), dtype=int)  # full connectivity

    last_good_T = None
    last_good_core = None
    for T in thresholds:
        binary = data > T
        n_vox = int(binary.sum())
        if n_vox == 0:
            continue
        labeled, n_components = label(binary, structure=structure)
        volume_fraction = n_vox / data.size
        sizes = np.bincount(labeled.ravel())
        sizes[0] = 0
        largest_size = int(sizes.max())

        if n_components > max_components or volume_fraction > max_volume_fraction:
            break
        if largest_size < min_blob_voxels:
            continue  # too small/fragmented still, keep going

        last_good_T = float(T)
        last_good_core = labeled == int(np.argmax(sizes))

    if last_good_core is None:
        raise RuntimeError("never found a stable single-blob threshold in "
                            f"[{end_percentile}, {start_percentile}] percentile range "
                            "-- map may need different start/end percentiles, inspect manually")

    return last_good_T, last_good_core


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    template_path, out_dir, prefix = sys.argv[1:4]
    n_steps = int(sys.argv[4]) if len(sys.argv) > 4 else 1000

    with mrcfile.open(template_path, permissive=True) as m:
        data = m.data.astype(np.float64)
        apix = float(m.voxel_size.x)

    T, core = find_auto_threshold(data, n_steps=n_steps)
    print(f"auto threshold: {T:.6f}  (core voxels: {int(core.sum())} = "
          f"{core.sum() / core.size:.4f} of box)")

    import os
    os.makedirs(out_dir, exist_ok=True)
    core_path = os.path.join(out_dir, f"{prefix}_AUTOTHRESH_core.mrc")
    with mrcfile.new(core_path, core.astype(np.float32), overwrite=True) as m:
        m.voxel_size = apix
    print(f"wrote binary core: {core_path}")
    print(f"threshold: {T:.6f}")


if __name__ == "__main__":
    main()
