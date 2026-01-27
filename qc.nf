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


process reportingQC {
  container "ghcr.io/lehtiolab/nfhelaqc:${workflow.manifest.version}"

  input:
  tuple val(acq_method), path('tpsms'), path('peptable.txt'), val(nrprots), val(nrpsms), val(nrpeps), val(nrunipeps), val(fwhmscans), path(scan_db), val(trackedpeptides)

  output:
  path('qc.json')

  script:
  """
  parse_output.py --acquisition $acq_method \
    --scandb $scan_db --nrpsms $nrpsms --nrpeps $nrpeps --peaks_on_lc $fwhmscans \
    --nruni $nrunipeps \
    --nrprot $nrprots \
    ${trackedpeptides ? "--trackedpeptides ${trackedpeptides.join(' ')}" :''}
  """
}


workflow {

  trackedpeptides = params.trackedpeptides ? params.trackedpeptides.tokenize(';') : false
  if (params.dda) {
    DDAQC(params.raw, params.mzml, params.instrument, params.db, params.mods, params.prectol, params.filters,
      params.options, params.psmconf, params.pepconf)
    | map { ['dda', it].flatten() + [trackedpeptides] }
    | reportingQC
  
  } else if (params.dia) {
  
    DIAQC(params.raw,  params.mzml, params.library, params.db, params.instrument, trackedpeptides)
    | map { ['dia', it].flatten() + [trackedpeptides] }
    | reportingQC
  }

  reportingQC.out
  | subscribe { it.copyTo("${params.outdir}/${it.baseName}.${it.extension}") }  
}
