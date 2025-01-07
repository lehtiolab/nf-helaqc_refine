/*
===================================
Longitudinal instrument QC pipeline
===================================
@Authors
Jorrit Boekel @glormph

Usage:
nextflow run qc.nf 
*/

nextflow.enable.dsl = 2

include { DDAQC } from './workflows/ddaqc.nf'
include { DIAQC } from './workflows/diaqc.nf'


process peptidesProteinsReport {
  container 'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple path('tpsms'), path('dpsms'), path('peptable.txt'), path('dpeptides'), path('db.sqlite'), val(noquant)

  output:
  tuple path('peptable.txt'), path('prottable.txt'), path('qc.json')

  script:
  scorecolpat = '^q-value'
  """
  sed -i 's/PepQValue/q-value/;s/QValue/PSM q-value/' "peptable.txt"
  sed -i 's/PepQValue/q-value/;s/QValue/PSM q-value/' dpeptides

  # score col is linearmodel_qval or q-value, but if the column only contains 0.0 or NA (no linear modeling possible due to only q<10e-04), we use Sage score instead
  tscol=\$(head -1 peptable.txt | tr '\\t' '\\n' | grep -n "${scorecolpat}" | cut -f 1 -d':')
  dscol=\$(head -1 dpeptides | tr '\\t' '\\n' | grep -n "${scorecolpat}" | cut -f 1 -d':')
  if [ -n "\$(cut -f \$tscol peptable.txt | tail -n+2 | egrep -v '(NA\$|0\\.0\$)')" ] && [ -n "\$(cut -f \$dscol dpeptides | tail -n+2 | egrep -v '(NA\$|0\\.0\$)')" ]
    then
      scpat="${scorecolpat}"
      logflag="--logscore"
    else
      scpat="sage_discriminant_score"
      logflag=""
      echo 'Not enough q-values or linear-model q-values for peptides to calculate FDR, using Sage disciminant score instead.' >> warnings
  fi
  msstitch proteins -i peptable.txt --decoyfn dpeptides -o tprots --scorecolpattern "\$scpat" \$logflag ${noquant ? '' : '--ms1quant'} --psmtable tpsms
  msstitch conffilt -i tprots -o prottable.txt --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'q-value'

  protcol=\$(head -1 peptable.txt | tr '\\t' '\\n' | grep -n "Master" | cut -f 1 -d':')
  parse_output.py db.sqlite "\$(wc -l tpsms)" "\$(wc -l peptable.txt)" "\$(cut -f \$protcol peptable.txt | grep -v ';' | wc -l)" "\$(wc -l tprots)"
  """
}


workflow {

  DDAQC(params.raw, params.mzml, params.instrument, params.db, params.mods, params.prectol, params.filters,
    params.options, params.noquant, params.psmconf, params.pepconf)

}
