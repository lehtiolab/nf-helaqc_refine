mkdir testrun
for fn in a b c d; do [ -e testrun/${fn}.mzML ] || cp testdata/fr08_veryminimal.mzML testrun/${fn}.mzML; done
echo `pwd`/testrun/a.mzML'\t'setA'\t'37-49'\t'08 > testrun/mzmldef.txt
echo `pwd`/testrun/b.mzML'\t'setA'\t'37-49'\t'09 >> testrun/mzmldef.txt
echo `pwd`/testrun/c.mzML'\t'setB'\t'37-49'\t'08 >> testrun/mzmldef.txt
echo `pwd`/testrun/d.mzML'\t'setB'\t'37-49'\t'09 >> testrun/mzmldef.txt

nextflow run quant_proteomics.nf --mzmldef testrun/mzmldef.txt --mods ../Mods.txt --ddb ../decoy_Homo_sapiens.GRCh38.ENS90.pep.all.fa --instrument qe  --tdb ../Homo_sapiens.GRCh38.ENS90.pep.all.fa  -profile testing --normalize --denoms 'setA:126 setB:131' --martmap ../Ensembl90_info.txt --outdir test_setmerge --genes --symbols  --pipep ~/formatted_known_peptides_ENSUniRefseq_TMT_predpi_20150825.txt --isobaric tmt10plex -resume 

