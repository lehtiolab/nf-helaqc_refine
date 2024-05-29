#!/usr/bin/env bash

rundir=$(pwd)
export repodir=$(dirname "$(realpath -s "$0")")
export testdir="${repodir}/tests/"
export testdata="${rundir}/static-resources/test-data/ddamsproteomics"

if [ -e "${testdata}" ]
then
    cd "${testdata}" && git checkout ddamsproteomics_test_data
else
    git clone --single-branch --branch ddamsproteomics_test_data https://github.com/lehtiolab/static-resources
fi
cd "${rundir}"


export NXF_VER=23.10.1

curl -s https://get.nextflow.io | bash
docker buildx -t nfhelaqc_test \
	-f ${repodir}/Dockerfile \
       	--cache-to type=gha \
       	--cache-from type=gha ${repodir} 

./nextflow run -resume -profile test "${repodir}/qc.nf" \
	--mzml  "${testdata}/lf_phos_fr11_500.mzML" \
	--noquant \
	--pepconf 0.05 --psmconf 0.05 \
	--db "${testdata}/lf.fa" \
	--instrument qe

./nextflow run -resume -profile test "${repodir}/refine_mzml.nf" \
	--input <(cat "${repodir}/test/refine_mzml.txt" | envsubst) \
	--db "${testdata}/lf.fa" \
	--instrument qe
