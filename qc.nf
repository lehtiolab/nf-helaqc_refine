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
  container 'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'

  input:
  tuple val(acq_method), path('tpsms'), path('peptable.txt'), val(nrprots), val(nrpsms), val(nrpeps), path(scan_db)

  output:
  path('qc.json')

  script:
  protfield = acq_method == 'dia' ? 'Genes' : 'Master protein(s)'
  """
  parse_output.py --acquisition $acq_method \
    --scandb $scan_db --nrpsms $nrpsms --nrpeps $nrpeps \
    --nruni "\$(cut -f${Utils.get_field_nr('peptable.txt', protfield)} peptable.txt | grep -v ';' | wc -l)" \
    --nrprot $nrprots
  """
}


workflow {

  if (params.dda) {
    DDAQC(params.raw, params.mzml, params.instrument, params.db, params.mods, params.prectol, params.filters,
      params.options, params.noquant, params.psmconf, params.pepconf)
    | map { ['dda', it].flatten() }
    | reportingQC
  
  } else if (params.dia) {
  
    DIAQC(params.raw, params.library, params.db)
    | map { ['dia', it].flatten() }
    | reportingQC
  }

}
