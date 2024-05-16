# RNA3DB scripts
Below are brief descriptions of the scripts in this folder.

- `scripts/slurm` a directory containing useful SLURM scripts.
- `scripts/build_incremental_release_fasta.py` can be used to extract the different chains from two `parse.json` files. Useful for incramental releases.
- `scripts/download_pdb_mmcif.sh` a script for downloading the latest version of the PDB.
- `scripts/fasta_to_json.py` take a [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) file and create a [JSON](https://en.wikipedia.org/wiki/JSON) usable by RNA3DB.
    - **Note:** that since FASTA files don't contain this information, the `release_date` is set to 1970-01-01, `structure_method` to "", and `resolution` to 0.0.
- `scripts/generate_modifications_cache.py` used to generate a modifications cache. See [Downloading required data](https://github.com/marcellszi/rna3db/wiki/Building-RNA3DB-from-scratch#downloading-required-data) on the RNA3DB Wiki.
- `scripts/get_nohits.py` looks at a FASTA file and `.tbl` file(s) and identifies entries in the FASTA file that get no hits in any of the `.tbl` file(s). Useful for the second `cmscan`. See [Homology Search](https://github.com/marcellszi/rna3db/wiki/Building-RNA3DB-from-scratch#homology-search) on the RNA3DB Wiki.
- `scripts/json_to_fasta.py` converts an RNA3DB [JSON](https://en.wikipedia.org/wiki/JSON) to a [FASTA file](https://en.wikipedia.org/wiki/FASTA_format).
- `scripts/json_to_mmcif.py` is used to build the single-chain [mmCIFs](https://en.wikipedia.org/wiki/Macromolecular_Crystallographic_Information_File). This script re-reads the chains from a `split.json` and writes them to a hierarchial folder, with each [mmCIF](https://en.wikipedia.org/wiki/Macromolecular_Crystallographic_Information_File) file containing a single chain.
