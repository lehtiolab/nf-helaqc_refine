#!/usr/bin/env bash
set -euo pipefail

rundir=$(pwd)
export repodir=$(dirname "$(realpath -s "$0")")
export testdir="${repodir}/tests/"
testdata_base="${rundir}/static-resources"
export testdata="${testdata_base}/test-data/ddamsproteomics"

if [ -e "${testdata_base}" ]
then
    cd "${testdata_base}" && git checkout ddamsproteomics_test_data
else
    git clone --single-branch --branch ddamsproteomics_test_data https://github.com/lehtiolab/static-resources
fi
cd "${rundir}"


docker buildx build -t nfhelaqc_test \
	-f ${repodir}/Dockerfile \
	--load \
       	--cache-to type=gha \
       	--cache-from type=gha \
       	${repodir} 

nextflow run -resume -profile test "${repodir}/qc.nf" \
        --dda \
	--mzml  "${testdata}/lf_phos_fr11_500.mzML" \
	--noquant \
	--pepconf 0.05 --psmconf 0.05 \
	--db "${testdata}/lf.fa" \
	--instrument qe

nextflow run -resume -profile test "${repodir}/refine_mzml.nf" \
	--input <(cat "${repodir}/test/refine_mzml.txt" | envsubst) \
	--db "${testdata}/lf.fa" \
	--instrument qe
