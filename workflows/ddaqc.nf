include { msconvert; createNewSpectraLookup } from '../modules.nf' 


process dinosaur {
  container params.test ? 'nfhelaqc_test' : \
    "ghcr.io/lehtiolab/nfhelaqc:${workflow.manifest.version}"

  input:
  path(mzml)

  output:
  path("dinosaur.features.tsv")

  script:
  """
  dinosaur -Xmx${task.memory.toMega()}M --concurrency=${task.cpus} --outName=dinosaur ${mzml}
  """
}


process makeDDB {
  container 'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'
 
  input:
  path(db)

  output:
  path('ddb')

  """
  msstitch makedecoy -i "$db" -o ddb --scramble prot_rev --ignore-target-hits
  """
}


process sagePrepare {
  container params.test ? 'nfhelaqc_test' : \
    "ghcr.io/lehtiolab/nfhelaqc:${workflow.manifest.version}"

  input:
  tuple path(tdb), path(ddb), val(prectol), val(fragtol), path('sage.json')
  output:
  tuple path('td_db'), path('config.json')
  
  script:
  """
  cat $tdb $ddb > td_db
  export PRECTOL=${prectol}
  export FRAGTOL=${fragtol}
  cat sage.json | envsubst > config.json
  """
}


process sage {
  container 'ghcr.io/lazear/sage:v0.14.7'

  input:
  tuple path(db), path('config.json'), path(mzml), path(mods)

  output:
  path('results.sage.pin'), emit: perco
  path('results.sage.tsv'), emit: tsv

  script:
  """
  export RAYON_NUM_THREADS=${task.cpus}
  export SAGE_LOG=trace
  sage --disable-telemetry-i-dont-want-to-improve-sage --write-pin -f $db config.json $mzml
  """
}


process percolator {
  container workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/percolator:3.5--hfd1433f_1' :
    'quay.io/biocontainers/percolator:3.6.5--h6351f2a_0'

  input:
  path('percoin.tsv')

  output:
  path('perco.xml')

  script:
  """
  percolator -j percoin.tsv -X perco.xml -N 500000 --decoy-xml-output -Y --num-threads ${task.cpus}
  """
}


process createPSMTable {
  container 'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "tpsmtable" ? "psmtable.txt" : null }

  input:
  tuple path('perco'), path('psms'), path(lookup), path(db), path(ddb), val(noquant), val(psmconf), val(pepconf)

  output:
  tuple path('tpsmtable'), path('dpsmtable'), path('tpeptides'), path('dpeptides'), emit: tables
  path('tpsmlookup'), emit: lookup

  script:
  """
  msstitch perco2psm --perco perco -i psms -o psms_perco --filtpsm ${psmconf} --filtpep ${pepconf}
  msstitch split -i psms_perco --splitcol TD
  cp $lookup tpsmlookup
  cp $lookup dpsmlookup
  msstitch psmtable -i target.tsv --dbfile tpsmlookup -o tpsmtable --addmiscleav --fasta "$db" ${noquant ? '' : '--ms1quant'} --proteingroup --spectracol ${Utils.get_field_nr('target.tsv', 'filename')}
  msstitch psmtable -i decoy.tsv --dbfile dpsmlookup -o dpsmtable --fasta "$ddb" --proteingroup --spectracol ${Utils.get_field_nr('decoy.tsv', 'filename')}
  msstitch peptides -i tpsmtable -o tpeptides --scorecolpattern sage_discriminant --spectracol ${Utils.get_field_nr('tpsmtable', 'filename')} ${noquant ? '' : '--ms1quantcolpattern area'}
  msstitch peptides -i dpsmtable -o dpeptides --scorecolpattern sage_discriminant --spectracol ${Utils.get_field_nr('dpsmtable', 'filename')}
  """
}


process proteinTables {
  container 'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple path('tpsms'), path('dpsms'), path('peptable.txt'), path('dpeptides'), val(noquant)

  output:
  tuple path('tpsms'), path('peptable.txt'), eval('wc -l < prottable.txt'),
    eval('wc -l < tpsms'), eval('wc -l < peptable.txt')

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
  """
}


workflow DDAQC {

  take:
  raw
  mzml
  instrument
  db
  mods
  prectol
  pwiz_filters
  pwiz_options
  noquant
  psmconf
  pepconf

  main:
  // Set fragment and precursor tolerance, if not passed - defaults:
  def tolerances = [precursor: [timstof: 30, qe: 10], fragment: [timstof: 30, qe: 20]]
  prectol = prectol ?: tolerances.precursor[instrument]
  fragtol = prectol ?: tolerances.fragment[instrument]

  if (raw) {
    channel.fromPath(raw)
    | map { [it, pwiz_filters, pwiz_options] }
    | msconvert
    | set { mzml }

  } else if (mzml) {
    channel.fromPath(mzml)
    | set { mzml }

  } else {
    exit 1, 'Must either input a --raw file.raw, or an --mzml file.mzML'
  }

  channel
    .fromPath(db)
    .set { tdb }

  makeDDB(tdb)

  if (!noquant) {
    mzml
    | dinosaur
    | set { dino }
  } else {
    channel.fromPath('NO__FILE')
    | set { dino }
  }

  mzml
  | combine(dino)
  | createNewSpectraLookup
    

  tdb
  | combine(makeDDB.out)
  | map { it + [prectol, fragtol, file("${baseDir}/assets/sage.json")] }
  | sagePrepare
  | combine(mzml)
  | combine(channel.fromPath(mods))
  | sage
  
  sage.out.perco
  | percolator
  | combine(sage.out.tsv)
  | combine(createNewSpectraLookup.out)
  | combine(tdb)
  | combine(makeDDB.out)
  | map { it + [noquant, psmconf, pepconf] }
  | createPSMTable

  createPSMTable.out.tables
  | map { it + noquant }
  | proteinTables 

  emit:
  proteinTables.out
  | combine(createPSMTable.out.lookup)
}
