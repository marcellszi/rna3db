#!/bin/bash
#SBATCH -c 64
#SBATCH -t 0
#SBATCH -p <insert partition here>
#SBATCH --mem=64000
#SBATCH -o logs/rna3db_incremental_release_%j.out
#SBATCH -e logs/rna3db_incremental_release_%j.err
#SBATCH --mail-user=<insert email address here>
#SBATCH --mail-type=ALL

# where you want the release to be output to
OUTPUT_DIR=""
# where the latest release is located
OLD_RELEASE=""

# you set these once and forget
RNA3DB_ROOT_DIR=""
PDB_MMCIF_DIR=""
CMSCAN=""
CMDB=""

# set useful constants
NEW_RELEASE_DATE=$(date +"%Y-%m-%d")
JSON_PATH=$OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-jsons
MMCIF_PATH=$OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-mmcifs
CMSCAN_PATH=$OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-cmscans
FASTA_PATH=$OUTPUT_DIR/$NEW_RELEASE_DATE/rna3db-fastas
UPLOAD_PATH=$OUTPUT_DIR/$NEW_RELEASE_DATE/uploads

# make directories
mkdir -p $JSON_PATH
mkdir -p $MMCIF_PATH
mkdir -p $CMSCAN_PATH
mkdir -p $FASTA_PATH
mkdir -p $UPLOAD_PATH

# prepare the env
mamba activate rna3db

# download latest mmCIF files
echo "Downloading..."
#bash $RNA3DB_ROOT_DIR/scripts/download_pdb_mmcif.sh $PDB_MMCIF_DIR
echo "Finished downloading."

# run parse
echo "Parsing..."
python -m rna3db parse \
    $PDB_MMCIF_DIR \
    $JSON_PATH/parse.json
echo "Finished parsing."

# run filter
echo "Filtering..."
python -m rna3db filter \
    $JSON_PATH/parse.json \
    $JSON_PATH/filter.json
echo "Finished filtering."

# write only the new sequences to a FASTA file
echo "Building FASTA..."
python $RNA3DB_ROOT_DIR/scripts/build_incremental_release_fasta.py \
    $OLD_RELEASE/rna3db-jsons/parse.json \
    $JSON_PATH/parse.json \
    $FASTA_PATH/$NEW_RELEASE_DATE.fasta
echo "Finished building FASTA."

echo "Copying old cmscans..."
cp $OLD_RELEASE/rna3db-cmscans/*.tbl $CMSCAN_PATH/
cp $OLD_RELEASE/rna3db-cmscans/*.o $CMSCAN_PATH/
echo "Copied old cmscans."

# do cmscan on new sequences
echo "Running cmscan..."
$CMSCAN --cpu 64 \
    -o $CMSCAN_PATH/$NEW_RELEASE_DATE.o \
    --tbl $CMSCAN_PATH/$NEW_RELEASE_DATE.tbl \
    $CMDB \
    $FASTA_PATH/$NEW_RELEASE_DATE.fasta
echo "Finished cmscan."

# find sequences that did not get a hit
echo "Finding sequences without hits..."
python $RNA3DB_ROOT_DIR/scripts/get_nohits.py \
    $FASTA_PATH/$NEW_RELEASE_DATE.fasta \
    $FASTA_PATH/$NEW_RELEASE_DATE-nohits.fasta \
    $CMSCAN_PATH
echo "Finished finding sequences."

# re-scan sequences with --max that did not get a hit
echo "Running cmscan (no hits)..."
$CMSCAN --max --cpu 64 \
    -o $CMSCAN_PATH/$NEW_RELEASE_DATE-nohits.o \
    --tbl $CMSCAN_PATH/$NEW_RELEASE_DATE-nohits.tbl \
    $CMDB \
    $FASTA_PATH/$NEW_RELEASE_DATE-nohits.fasta
echo "Finished cmscan (no hits)."

# run cluster
echo "Clustering..."
python -m rna3db cluster \
    $JSON_PATH/filter.json \
    $JSON_PATH/cluster.json \
    --tbl_dir $CMSCAN_PATH
echo "Cleaning up..."
rm -r $JSON_PATH/mmseqs2_*
echo "Finished clustering."

# run split
echo "Splitting..."
python -m rna3db split \
    $JSON_PATH/cluster.json \
    $JSON_PATH/split.json
echo "Finished splitting."

# make mmCIFs
echo "Making mmCIFs..."
python scripts/json_to_mmcif.py \
    $JSON_PATH/split.json \
    $PDB_MMCIF_DIR \
    $MMCIF_PATH
echo "Finished writing mmCIFs."

# compress files ready for release
echo "Compressing..."
tar -czf \
    $UPLOAD_PATH/rna3db-cmscans.tar.gz \
    $CMSCAN_PATH
tar -czf \
    $UPLOAD_PATH/rna3db-jsons.tar.gz \
    $JSON_PATH
tar -cJf \
    $UPLOAD_PATH/rna3db-mmcifs.tar.xz \
    $MMCIF_PATH

echo "Finished."
echo $NEW_RELEASE_DATE "release ready for upload."

