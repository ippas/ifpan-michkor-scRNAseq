# Project title (institute-name-subject)

#### Project logline (technique, organism, tissue type)
Short description of treatment groups/subjects


## Methods
This sections should be a description of preprocessin and analysis ready to be included in the publication


## Preprocessing
Details of file preprocessing

## Analysis
Details of analysis

*notes: all files included in the repo need to be referenced, either in README or other .md files. The analysis has to be fully reproducible, in principle the repo should contain code + description of how to run it while data and results kept outside*

## About this template
Directories:
- _root_ - README.md, *.Rproj, general configuration files, etc.
- raw - raw data
- preprocessing - scripts
- data - useful data, created by scripts/tools/preprocessing
- analysis - analysis source code
- results - output ready to present

Comand to prepare ref genome
```
split-pipe   \
  --mode mkref   \
--genome_name GRCh38   \
--fasta data/ref_hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz   \
--genes data/ref_hg38/Homo_sapiens.GRCh38.113.chr.gtf.gz   \
--nthreads 16   \
--output_dir data/ref_hg38/GRCh38_ref

```

First trial
```
split-pipe \
  --mode all \
  --kit WT_mega \
  --chemistry v2 \
  -f raw/pilotDexCortCtrl/BG60BLUE_L6/BG60BLUE_MKDL250008382-1A_22VVLVLT4_L6_1.fq.gz \
  --fq2 raw/pilotDexCortCtrl/BG60BLUE_L6/BG60BLUE_MKDL250008382-1A_22VVLVLT4_L6_2.fq.gz \
  --nthreads 16 \
  --output_dir data/pilotDexCortCtrl/BG60BLUE_L6-out \
  --genome_dir data/ref_hg38/GRCh38_ref \
  --sample all-well A1-H12 \
  --chem_score_skip \
  --kit_score_skip
```
