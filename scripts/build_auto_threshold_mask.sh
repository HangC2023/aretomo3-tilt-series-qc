#!/bin/bash
# build_auto_threshold_mask.sh -- mask from the map's own noise floor
# (scripts/auto_threshold_mask.py: connected-component threshold sweep),
# extended+soft-edged via relion_mask_create for consistency with the
# rest of this project's mask outputs, mirror built by flipping (never
# independently re-thresholded, same discipline as build_reference_set.sh
# and pdb_to_reference_mask.sh).
#
# Built 2026-09-15, Sean's proposal (see auto_threshold_mask.py docstring
# for the full rationale) as an independent alternative to the PDB/molmap-
# derived masks -- no fitted-PDB dependency, derives entirely from the
# map's own signal-vs-noise transition.
#
# Usage:
#   build_auto_threshold_mask.sh <template.mrc> <output_dir> <prefix> \
#       [extend_inimask=2] [width_soft_edge=4] [n_steps=300] \
#       [compare_mask.mrc] [--mirror]
#
# compare_mask.mrc: if given (e.g. the existing PDB-derived mask), reports
# Dice/Jaccard overlap between the two independently-derived masks.
#
# --mirror: also write a _MASK_mirror.mrc (exact axis flip of the mask,
# never independently re-thresholded). Off by default -- pass this flag
# (anywhere in the arguments) to opt in.

set -eo pipefail

MAKE_MIRROR=false
ARGS=()
for arg in "$@"; do
    if [ "$arg" = "--mirror" ]; then
        MAKE_MIRROR=true
    else
        ARGS+=("$arg")
    fi
done
set -- "${ARGS[@]}"

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <template.mrc> <output_dir> <prefix> [extend_inimask] [width_soft_edge] [n_steps] [compare_mask.mrc] [--mirror]" >&2
    exit 1
fi
TEMPLATE="$1"; OUT_DIR="$2"; PREFIX="$3"
EXTEND_INIMASK="${4:-2}"
WIDTH_SOFT_EDGE="${5:-4}"
N_STEPS="${6:-1000}"
COMPARE_MASK="${7:-}"

mkdir -p "$OUT_DIR"
cd "$(dirname "$0")/.."
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate pytom_tm

echo "=== Step 1: connected-component threshold sweep ==="
python3 scripts/auto_threshold_mask.py "$TEMPLATE" "$OUT_DIR" "$PREFIX" "$N_STEPS"
CORE="${OUT_DIR}/${PREFIX}_AUTOTHRESH_core.mrc"

echo "=== Step 2: relion_mask_create (extend+soft-edge on the binary core) ==="
MASK="${OUT_DIR}/${PREFIX}_AUTOTHRESH_MASK.mrc"
relion_mask_create --i "$CORE" --o "$MASK" \
    --ini_threshold 0.5 --extend_inimask "$EXTEND_INIMASK" \
    --width_soft_edge "$WIDTH_SOFT_EDGE" --j 1
echo "wrote: $MASK"

if [ "$MAKE_MIRROR" = true ]; then
    echo "=== Step 3: mirror mask (flip the normal mask, axis=2 -- never re-thresholded) ==="
    python3 - "$MASK" "${OUT_DIR}/${PREFIX}_AUTOTHRESH_MASK_mirror.mrc" <<'PYEOF'
import sys
import mrcfile
import numpy as np
normal_path, mirror_path = sys.argv[1], sys.argv[2]
with mrcfile.open(normal_path, permissive=True) as m:
    data = m.data.copy()
    apix = float(m.voxel_size.x)
mirrored = np.flip(data, axis=2).copy()
with mrcfile.new(mirror_path, mirrored, overwrite=True) as m:
    m.voxel_size = apix
print(f"wrote {mirror_path} (exact flip of {normal_path})")
PYEOF
else
    echo "=== Step 3: skipped (pass --mirror to also write a mirrored mask) ==="
fi

echo "=== Step 4: coverage check against ${TEMPLATE} ==="
python3 scripts/check_mask_coverage.py --map "$TEMPLATE" --mask "$MASK" || true

if [ -n "$COMPARE_MASK" ]; then
    echo "=== Step 5: cross-check against ${COMPARE_MASK} (independent method) ==="
    python3 - "$MASK" "$COMPARE_MASK" <<'PYEOF'
import sys
import mrcfile
import numpy as np
a_path, b_path = sys.argv[1], sys.argv[2]
with mrcfile.open(a_path, permissive=True) as m:
    a = m.data.astype(np.float64) > 0.5
with mrcfile.open(b_path, permissive=True) as m:
    b = m.data.astype(np.float64) > 0.5
intersection = np.logical_and(a, b).sum()
union = np.logical_or(a, b).sum()
dice = 2 * intersection / (a.sum() + b.sum()) if (a.sum() + b.sum()) > 0 else float('nan')
jaccard = intersection / union if union > 0 else float('nan')
print(f"auto-threshold mask volume: {a.sum()} voxels")
print(f"compare mask volume:        {b.sum()} voxels")
print(f"Dice coefficient:    {dice:.4f}")
print(f"Jaccard index:       {jaccard:.4f}")
if dice < 0.7:
    print("NOTE: two independent masking methods disagree substantially (Dice<0.7) -- inspect visually.")
else:
    print("Two independent masking methods agree well (Dice>=0.7).")
PYEOF
fi

echo "=== Done. ==="
