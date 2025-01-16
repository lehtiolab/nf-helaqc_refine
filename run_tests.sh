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


#docker buildx build -t nfhelaqc_test \
#	-f ${repodir}/Dockerfile \
#	--load \
#       	--cache-to type=gha \
#       	--cache-from type=gha \
#       	${repodir} 

nextflow run -resume -profile test "${repodir}/qc.nf" \
        --dda \
	--mzml  "${testdata}/lf_phos_fr11_500.mzML" \
	--pepconf 0.05 --psmconf 0.05 \
	--db "${testdata}/lf.fa" \
	--instrument qe
	#--noquant \

nextflow run -resume -profile test "${repodir}/qc.nf" -with-trace \
        --dia \
	--raw "raws/241216_KN_MetEval_DIA_44min_200ng_HeLaiRT_QC_2024-12-16_6596.d" \
	--library "newlib-mods.speclib" \
	--db "up_tdb.fa"
	#--db "${testdata}/lf.fa" \
	#--pepconf 0.05 --psmconf 0.05 \

	#--raw ddraws/MZ_HeLa200ng_44min_dda_QC_2024-08-19_5739.d" \

# In testing mzML, explicitly pass bruker param
nextflow run -resume -profile test "${repodir}/qc.nf" -with-trace \
        --dda --bruker \
	--raw "raws/MZ_HeLa200ng_44min_dda_QC_2024-08-19_5739.d" \
	--db "up_tdb.fa" \
	--instrument timstof

nextflow run -resume -profile test "${repodir}/qc.nf" -with-trace \
        --dia \
	--raw "raws/20241128_QC_Figaro_HeLa_200ng_DIA_41p5min_2Th_3ms_r2.raw" \
	--library "newlib-mods.speclib" \
	--db "up_tdb.fa" \
	--instrument qe 

#nextflow run -resume -profile test "${repodir}/refine_mzml.nf" \
#	--input <(cat "${repodir}/test/refine_mzml.txt" | envsubst) \
#	--db "${testdata}/lf.fa" \
#	--instrument qe
