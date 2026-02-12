# RNA3DB Test Data

## `mmcifs`
These files are used in `test_pipeline.py` to make sure everything works correctly end-to-end.
- `1ehz.cif`: A single-chain tRNA, good easy example.
- `3tup.cif`: A tRNA + protein. 
- `3cgs.cif`: Two short chains, should probably be filtered.
- `5m0h.cif`: Short fragment in component #0.
- `1y27.cif`: A riboswitch in component #2. 

## `tbls`
- Rows from 2026-01-05 release extracted for the structures found in `mmcifs`. 

## Miscellaneous
- `old_1ehz_A.cif`: For ensuring backwards compatibility with older versions of RNA3DB.