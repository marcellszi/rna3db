# RNA3DB SLURM scripts

These [SLURM](https://slurm.schedmd.com/documentation.html) scripts will eventually be used to build releases automatically.

> **Note:** The scripts are experimental as they haven't been rigorously tested.


## Getting started
The first of these script, `build_full_release.slurm`, builds an entire release from the start. This script does a homology search on all chains found in the PDB, so it takes a long time to run.

The second script, `build_incremental_release.slurm` adds new chains (added to the PDB since last release) to an existing release.

Both files start with a number of [sbatch](https://slurm.schedmd.com/sbatch.html) SLURM commands:
```sh
#SBATCH -c 64
#SBATCH -t 0
#SBATCH -p <insert partition here>
#SBATCH --mem=64000
#SBATCH -o logs/rna3db_full_release_%j.out
#SBATCH -e logs/rna3db_full_release_%j.err
#SBATCH --mail-user=<insert email address here>
#SBATCH --mail-type=ALL
```
You will likely need to edit some of these options if you want to use these scripts. Please see the [SLURM documentation for sbatch](https://slurm.schedmd.com/sbatch.html) on what each line means. At least you will need to either enter a partition, or remove the `-p` option. Similarly, you will need to edit the `--mail-user` option.

Next, there are a number of paths you need to set in both files:
```sh
# where you want the release to be output to
OUTPUT_DIR=""
# where the latest release is located
OLD_RELEASE=""

# you set these once and forget
RNA3DB_ROOT_DIR=""
PDB_MMCIF_DIR=""
CMSCAN=""
CMDB=""
```
- `OUTPUT_DIR` specifies the root directory where the release will be placed
- `OLD_RELEASE` is the path to the directory of the release you want to add the new PDB chains to. This is only needed when you are trying to build an incremental release.
- `RNA3DB_ROOT_DIR` path to the rna3db repository. Scripts are called from `$RNA3DB_ROOT_DIR/scripts/`.
- `CMSCAN` is the path to the `cmscan` executable.
- `CMDB` is the path to the covariance models you want to use for the homology search (`cmscan`). Usually this would come from [Rfam](https://rfam.org/) in the form of `Rfam.cm`.

Once you have set the required paths and edited the sbatch commands as required, you can simply run the jobs via:
```sh
$ sbatch build_full_release.slurm
```
Or:
```sh
$ sbatch build_incremental_release.slurm
```
