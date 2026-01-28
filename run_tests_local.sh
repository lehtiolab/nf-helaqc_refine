#!/usr/bin/env bash
set -euo pipefail

# Local test script that contains raw files, since we dont have small raw files for on github

rundir=$(pwd)
export repodir=$(dirname "$(realpath -s "$0")")
export testdir="${repodir}/tests/"
testdata_base="${rundir}/static-resources"
export testdata="${testdata_base}/test-data/ddamsproteomics"

## Fill these in yourself:
BRUKER_DIA="/path/to/file.d"
BRUKER_DDA="..."
THERMO_DIA="/path/to/file.raw"
THERMO_DDA="..."
THERMO_DIA_MZML="/path/to/file.mzML"

# Thermo mzML DIA (dont have a small one)
NXF_VER=24.04.3 nextflow run -resume -profile bigtest "${repodir}/qc.nf" -with-trace \
        --dia \
	--mzml "$THERMO_DIA_MZML" \
	--db "${testdata}/lf.fa" \
	--library "lf.fa.predicted.speclib" \
	--instrument qe 


### TIMS raw DIA (raw directly analyzed)
NXF_VER=24.10.3 nextflow run -resume -profile bigtest "${repodir}/qc.nf" -with-trace \
        --dia \
	--raw "$BRUKER_DIA" \
	--instrument timstof \
	--library "lf.fa.predicted.speclib" \
	--trackedpeptides "LGGNEQVTR_2;YILAGVENSK_2" \
	--db "${testdata}/lf.fa" \

# TIMS raw DDA (runs tdf2mzml)
nextflow run -resume -profile bigtest "${repodir}/qc.nf" -with-trace \
        --dda \
	--raw "$BRUKER_DDA" \
	--db "${testdata}/lf.fa" \
	--instrument timstof

# Thermo raw DIA
NXF_VER=24.04.3 nextflow run -resume -profile bigtest "${repodir}/qc.nf" -with-trace \
        --dia \
	--raw "$THERMO_DIA" \
	--db "${testdata}/lf.fa" \
	--library "lf.fa.predicted.speclib" \
	--instrument qe 

# Thermo raw DDA
# How can it be such a bad score?
nextflow run -resume -profile bigtest "${repodir}/qc.nf" -with-trace \
        --dda \
	--raw "$THERMO_DDA" \
	--db "${testdata}/lf.fa" \
	--pepconf 0.5 --psmconf 0.5 \
	--instrument qe 

bash "${repodir}/run_tests.sh"
