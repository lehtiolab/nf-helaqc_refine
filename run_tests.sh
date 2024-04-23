#!/usr/bin/env bash

rundir=$(pwd)
export repodir=$(dirname "$(realpath -s "$0")")
export testdir="${repodir}/tests/"
export testdata="${rundir}/static-resources/test-data/ddamsproteomics"

if [ -e "${testdata}" ]
then
    cd "${testdata}" && git checkout ddamsproteomics_test_data && cd "${rundir}"
else
    git clone --single-branch --branch ddamsproteomics_test_data https://github.com/lehtiolab/static-resources
fi


export NXF_VER=23.10.1

curl -s https://get.nextflow.io | bash
nextflow run -resume -profile docker "${repodir}/refine_mzml.nf" \
	--input <(cat "${repodir}/test/refine_mzml.txt" | envsubst) \
	--db "${testdata}/lf.fa" \
	--instrument qe
