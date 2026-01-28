# Internal QC and mzML precursor refining pipeline

A nextflow pipeline to run QC for proteomics instruments. Handles Thermo and Bruker DDA and DIA runs,
the output is a `qc.json` file with some statistics, such as nr of scans, peptides, precursor error
(quartiles), retention times, etc. 

## Usage
Mainly internal to our lab, but an example is:
```
nextflow run -resume -profile docker qc.nf \
    --dia \
    --db your-peptides.fa \
    --library your-peptides.speclib \
    --instrument [timstof, qe] \ # pick one of these
    --raw your-QC-run.raw 
```

Use `--dda` instead of `--dia` and remove the `--library` param when running DDA data. When testing on github, instead of `--raw`, we use `--mzml` to keep files small.


## Development

Install Nextflow and Docker. For local work, test with:
```
# This also runs on github actions
bash run_tests.sh

# Locally, supply your own raw files which are too big for GHA
bash run_tests_local.sh
```

On github actions we run on PR/push to master.


## Releasing
Create a PR, make changes, update the version etc. Make sure you publish a release to update the container (or push that by hand).


## Todo:
- Block PRs where there is no version change
