#!/usr/bin/env python3
"""check_mask_coverage.py -- does a mask actually cover the signal in a map?

Built 2026-09-15 alongside build_reference_set.sh's mask-creation option and
scripts/pdb_to_reference_mask.sh (model-derived masks), per Sean's request:
"We could just test that the resulting mask covers the signal in the map."

Two independent coverage metrics are reported, since either alone can be
misleading:

  power_fraction: sum(v^2 * mask) / sum(v^2) over the whole box. This is the
    physically relevant one for template matching -- the mask multiplies the
    density before the FFT cross-correlation, so this is literally the
    fraction of correlation signal the mask keeps. But it's dominated by the
    brightest few voxels (v^2 weighting), so a mask could score well here
    while still clipping a real but dim lobe of the structure.

  top_percentile_coverage: of the voxels in the top `--top-percentile` % by
    |density| (default top 1%, i.e. the clearest structural signal), what
    fraction fall inside mask > 0.5. Catches exactly the clipping case above.

Usage:
  check_mask_coverage.py --map <map_or_template.mrc> --mask <mask.mrc>
      [--power-fraction-threshold 0.95] [--top-percentile 1.0]
      [--top-percentile-threshold 0.98]

Map and mask must have identical box shape (same convention as the rest of
this project's reference-building scripts -- build them from the same
box-size/apix first, e.g. via build_reference_set.sh / pdb_to_reference_mask.sh).
"""
import argparse
import sys

import mrcfile
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--map", required=True, help="reference map or built template (.mrc)")
    p.add_argument("--mask", required=True, help="mask to evaluate (.mrc)")
    p.add_argument("--power-fraction-threshold", type=float, default=0.95)
    p.add_argument("--top-percentile", type=float, default=1.0,
                    help="top X%% of voxels by |density| used for the hard-coverage check")
    p.add_argument("--top-percentile-threshold", type=float, default=0.98)
    args = p.parse_args()

    with mrcfile.open(args.map, permissive=True) as f:
        v = f.data.astype(np.float64)
    with mrcfile.open(args.mask, permissive=True) as f:
        m = f.data.astype(np.float64)

    if v.shape != m.shape:
        sys.exit(f"ERROR: shape mismatch -- map {v.shape} vs mask {m.shape}. "
                  f"Build both at the same box size/apix before comparing.")
    if m.min() < -1e-4 or m.max() > 1.0001:
        print(f"WARNING: mask values outside [0,1] (min={m.min():.4f}, max={m.max():.4f}) "
              f"-- is this really a mask file?", file=sys.stderr)

    total_power = float(np.sum(v ** 2))
    captured_power = float(np.sum((v ** 2) * m))
    power_fraction = captured_power / total_power if total_power > 0 else float("nan")

    cutoff = np.percentile(np.abs(v), 100.0 - args.top_percentile)
    top_mask = np.abs(v) >= cutoff
    n_top = int(top_mask.sum())
    n_top_covered = int(np.sum(top_mask & (m > 0.5)))
    top_coverage = n_top_covered / n_top if n_top > 0 else float("nan")

    mask_volume_voxels = int(np.sum(m > 0.5))
    mask_volume_frac_of_box = mask_volume_voxels / m.size

    print(f"map:  {args.map}  shape={v.shape}")
    print(f"mask: {args.mask}  shape={m.shape}")
    print()
    print(f"power_fraction (sum v^2*mask / sum v^2):        {power_fraction:.4f}  "
          f"(threshold {args.power_fraction_threshold})")
    print(f"top-{args.top_percentile:g}%-intensity coverage (frac inside mask>0.5): "
          f"{top_coverage:.4f}  (threshold {args.top_percentile_threshold}, "
          f"{n_top_covered}/{n_top} voxels)")
    print(f"mask volume: {mask_volume_voxels} voxels = {mask_volume_frac_of_box:.4f} of box "
          f"(sanity check -- should be a small fraction of the box, not most of it)")

    ok = (power_fraction >= args.power_fraction_threshold
          and top_coverage >= args.top_percentile_threshold)
    print()
    if ok:
        print("PASS: mask captures the map's signal within the given thresholds.")
    else:
        print("WARN: mask may be clipping real signal -- inspect visually before trusting it.")
        sys.exit(1)


if __name__ == "__main__":
    main()
