# Internal QC and mzML precursor refining pipeline

## Usage
Mainly internal, but if you must:

```
nextflow run -resume -profile docker qc.nf \
    --db your-peptides.fa \
    --instrument [timstof, qe] \ # pick one of these
    --raw your-QC-run.raw 
```

When testing, instead of `--raw`, we use `--mzml` to keep files small.


## Development

For local work:
```
bash run_tests.sh
```

On github actions we run on PR/push to master.


## Releasing
Create a PR, make changes, update the version etc. Make sure you publish a release to update the container.


## Todo:
- Block PRs where there is no version change
- Fix local testing script
