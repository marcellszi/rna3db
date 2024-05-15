#!/bin/bash
#SBATCH -c 64
#SBATCH -t 0
#SBATCH -p <insert partition here>
#SBATCH --mem=64000
#SBATCH -o logs/rna3db_full_release_%j.out
#SBATCH -e logs/rna3db_full_release_%j.err
#SBATCH --mail-user=<insert email address here>
#SBATCH --mail-type=ALL

# the output dir, along with the date of the last release
OUTPUT_DIR=""
PREVIOUS_RELEASE_DATE="2024-04-26"

# you set these once and forget
RNA3DB_ROOT_DIR=""
MMCIF_DIR=""
CMSCAN=""
CMDB=""

NEW_RELEASE_DATE=$(date +"%Y-%m-%d")
mkdir -p $OUTPUT_DIR/$NEW_RELEASE_DATE

# prepare the env
mamba activate rna3db

# download latest mmCIF files
bash $RNA3DB_ROOT_DIR/scripts/download_pdb_mmcif.sh $MMCIF_DIR

# run parse
python -m rna3db parse \
    $MMCIF_DIR \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-jsons/parse.json

# run filter
python -m rna3db filter $MMCIF_DIR $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-jsons/filter.json

# write only the new sequences to a FASTA file
python $RNA3DB_ROOT_DIR/scripts/build_incremental_release_fasta.py \
    $OUTPUT_DIR/$OLD_RELEASE_DATE/rna3db-jsons/parse.json \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-jsons/parse.json \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/$NEW_RELEASE_DATE.fasta

# do cmscan on all new sequences
CMSCAN --cpu 64 \
    -o $OUTPUT_DIR/$NEW_RELEASE_DATE/$NEW_RELEASE_DATE.o \
    -tbl $OUTPUT_DIR/$NEW_RELEASE_DATE/$NEW_RELEASE_DATE.tbl \
    $CMDB \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/$NEW_RELEASE_DATE.fasta

# find new sequences that did not get a hit
python $RNA3DB_ROOT_DIR/scripts/get_nohits.py \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/$NEW_RELEASE_DATE.fasta \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/$NEW_RELEASE_DATE-nohits.fasta \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/

# re-scan new sequences with --max that did not get a hit
CMSCAN --max --cpu 64 \
    -o $OUTPUT_DIR/$NEW_RELEASE_DATE/$NEW_RELEASE_DATE-nohits.o \
    -tbl $OUTPUT_DIR/$NEW_RELEASE_DATE/$NEW_RELEASE_DATE-nohits.tbl \
    $CMDB \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/$NEW_RELEASE_DATE.fasta

# run cluster
python -m rna3db cluster \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-jsons/filter.json \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-jsons/cluster.json \
    --tbl_dir $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-cmscans

# run split
python -m rna3db split \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-jsons/cluster.json \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-jsons/split.json

# make mmCIFs
python scripts/json_to_mmcif.py \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-jsons/split.json \
    /Users/marcell/Documents/rna3db/data/pdb_mmcif/ \
    $MMCIF_DIR \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-mmcifs

# compress files ready for release
tar -czvf \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-cmscans.tar.gz \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-cmscans
tar -czvf \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-cmscans.tar.gz \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-cmscans
tar -cfvJ \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-mmcifs.tar.xz \
    $OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-mmcifs
