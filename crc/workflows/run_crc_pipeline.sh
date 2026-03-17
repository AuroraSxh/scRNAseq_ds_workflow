#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

source .venv/bin/activate

python src/01_prepare_qc.py
python src/02_integrate_cluster.py
python src/03_annotate_major.py
python src/04_cd8_subtype.py
python src/07_state_umap_and_bubbleplots.py
python src/09_refine_subgroups_and_interactions.py
