#!/bin/bash

set -e  # zakończ przy błędzie
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"  # zakładam że skrypt jest np. w preprocessing/pilotDexCortCtrl

cd "$PROJECT_DIR"

echo "run analysis" 


# BG60BLUE_L4
split-pipe \
  --mode all \
  --kit WT_mini \
  --chemistry v3 \
  -f raw/pilotDexCortCtrl/BG60BLUE_L4/BG60BLUE_MKDL250008382-1A_22VVMKLT4_L4_1.fq.gz \
  --fq2 raw/pilotDexCortCtrl/BG60BLUE_L4/BG60BLUE_MKDL250008382-1A_22VVMKLT4_L4_2.fq.gz \
  --nthreads 16 \
  --output_dir data/pilotDexCortCtrl/BG60BLUE_L4-out \
  --genome_dir data/ref_hg38/GRCh38_ref \
  --sample all-well A1-A12 \
  --chem_score_skip \
  --kit_score_skip

# BG60BLUE_L6
split-pipe \
  --mode all \
  --kit WT_mini \
  --chemistry v3 \
  -f raw/pilotDexCortCtrl/BG60BLUE_L6/BG60BLUE_MKDL250008382-1A_22VVLVLT4_L6_1.fq.gz \
  --fq2 raw/pilotDexCortCtrl/BG60BLUE_L6/BG60BLUE_MKDL250008382-1A_22VVLVLT4_L6_2.fq.gz \
  --nthreads 16 \
  --output_dir data/pilotDexCortCtrl/BG60BLUE_L6-out \
  --genome_dir data/ref_hg38/GRCh38_ref \
  --sample all-well A1-A12 \
  --chem_score_skip \
  --kit_score_skip

# BG60GREEN_L4
split-pipe \
  --mode all \
  --kit WT_mini \
  --chemistry v3 \
  --fq1 raw/pilotDexCortCtrl/BG60GREEN_L4/BG60GREEN_MKDL250008381-1A_22VVMKLT4_L4_1.fq.gz \
  --fq2 raw/pilotDexCortCtrl/BG60GREEN_L4/BG60GREEN_MKDL250008381-1A_22VVMKLT4_L4_2.fq.gz \
  --nthreads 16 \
  --output_dir data/pilotDexCortCtrl/BG60GREEN_L4-out \
  --genome_dir data/ref_hg38/GRCh38_ref \
  --sample all-well A1-A12 \
  --chem_score_skip \
  --kit_score_skip

# BG60GREEN_L6
split-pipe \
  --mode all \
  --kit WT_mini \
  --chemistry v3 \
  -f raw/pilotDexCortCtrl/BG60GREEN_L6/BG60GREEN_MKDL250008381-1A_22VVLVLT4_L6_1.fq.gz \
  --fq2 raw/pilotDexCortCtrl/BG60GREEN_L6/BG60GREEN_MKDL250008381-1A_22VVLVLT4_L6_2.fq.gz \
  --nthreads 16 \
  --output_dir data/pilotDexCortCtrl/BG60GREEN_L6-out \
  --genome_dir data/ref_hg38/GRCh38_ref \
  --sample all-well A1-A12 \
  --chem_score_skip \
  --kit_score_skip

# Połączenie lane'ów BLUE
split-pipe \
  --mode comb \
  --output_dir data/pilotDexCortCtrl/combined_BG60BLUE \
  --sublibraries \
    data/pilotDexCortCtrl/BG60BLUE_L4-out \
    data/pilotDexCortCtrl/BG60BLUE_L6-out

# Połączenie lane'ów GREEN
split-pipe \
  --mode comb \
  --output_dir data/pilotDexCortCtrl/combined_BG60GREEN \
  --sublibraries \
    data/pilotDexCortCtrl/BG60GREEN_L4-out \
    data/pilotDexCortCtrl/BG60GREEN_L6-out

# Połączenie BLUE i GREEN (całościowe połączenie)
split-pipe \
  --mode comb \
  --output_dir data/pilotDexCortCtrl/combined_all \
  --sublibraries \
    data/pilotDexCortCtrl/combined_BG60BLUE \
    data/pilotDexCortCtrl/combined_BG60GREEN

