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



process msconvert {

  cpus = 4 // FIXME 4 for TIMSTOF, XX for normal?

  input:
  tuple path(raw), val(filters), val(options)

  output:
  path(outfile)

  script:
  // the string "infile" does not have NF escaping characters like & (e.g. in FAIMS 35&65),
  // which it does to "raw". That would work fine but not if the files are quoted in the 
  // script, then they cant be found when there is \&.
  infile = "${raw.baseName}.${raw.extension}"
  outfile = "${raw.baseName}.mzML"
  """
  # Resolve directory if necessary, pwiz tries to read NF soft links as if they are files, which
  # does not work in case of directory
  ${raw.isDirectory() ?  "mv '${infile}' tmpdir && cp -rL tmpdir '${infile}'" : ''}
  wine msconvert "${infile}" ${filters} ${options}
  """
}


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
  container 'quay.io/biocontainers/msstitch:3.15--pyhdfd78af_0'
 
  input:
  path(db)

  output:
  path('ddb')

  """
  msstitch makedecoy -i "$db" -o ddb --scramble prot_rev --ignore-target-hits
  """
}


process createSpectraLookup {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple path(mzml), path(dino)
  
  output:
  path('mslookup_db.sqlite')

  script:
  """
  msstitch storespectra --spectra "${mzml}" --setnames 'QC'
  ${dino.baseName != 'NO__FILE' ? "msstitch storequant --dbfile mslookup_db.sqlite --dinosaur \"${dino}\" --spectra \"${mzml}\" --mztol 20.0 --mztoltype ppm --rttol 5.0" : ''}
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


def get_field_nr(fn, pattern) {
  return "\$(head -n1 ${fn} | tr '\t' '\n' | grep -wn '^${pattern}\$' | cut -f1 -d':')"
}

process createPSMTable {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "tpsmtable" ? "psmtable.txt" : null }

  input:
  //tuple path('tpsms'), path('dpsms'), path(lookup), path(db), path(ddb), val(noquant), val(psmconf), val(pepconf)
  tuple path('perco'), path('psms'), path(lookup), path(db), path(ddb), val(noquant), val(psmconf), val(pepconf)

  output:
  tuple path('tpsmtable'), path('dpsmtable'), path('tpeptides'), path('dpeptides'), path('tpsmlookup')

  script:
  """
  msstitch perco2psm --perco perco -i psms -o psms_perco --filtpsm ${psmconf} --filtpep ${pepconf}
  msstitch split -i psms_perco --splitcol TD
  cp $lookup tpsmlookup
  cp $lookup dpsmlookup
  msstitch psmtable -i target.tsv --dbfile tpsmlookup -o tpsmtable --addmiscleav --fasta "$db" ${noquant ? '' : '--ms1quant'} --proteingroup --spectracol ${get_field_nr('target.tsv', 'filename')}
  msstitch psmtable -i decoy.tsv --dbfile dpsmlookup -o dpsmtable --fasta "$ddb" --proteingroup --spectracol ${get_field_nr('decoy.tsv', 'filename')}
  msstitch peptides -i tpsmtable -o tpeptides --scorecolpattern sage_discriminant --spectracol ${get_field_nr('tpsmtable', 'filename')} ${noquant ? '' : '--ms1quantcolpattern area'}
  msstitch peptides -i dpsmtable -o dpeptides --scorecolpattern sage_discriminant --spectracol ${get_field_nr('dpsmtable', 'filename')}
  """
}


process peptidesProteinsReport {

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

  // Set fragment and precursor tolerance, if not passed - defaults:
  def tolerances = [precursor: [timstof: 30, qe: 10], fragment: [timstof: 30, qe: 20]]
  prectol = params.prectol ?: tolerances.precursor[params.instrument]
  fragtol = params.prectol ?: tolerances.fragment[params.instrument]

  if (params.raw) {
    filters = params.filters.tokenize(';').collect() { x -> "--filter ${x}" }.join(' ')
    options = params.options.tokenize(';').collect() {x -> "--${x}"}.join(' ')
    channel.fromPath(params.raw)
    | map { [it, filters, options] }
    | msconvert
    | set { mzml }

  } else if (params.mzml) {
    channel.fromPath(params.mzml)
    | set { mzml }

  } else {
    exit 1, 'Must either input a --raw file.raw, or an --mzml file.mzML'
  }

  channel
    .fromPath(params.db)
    .set { tdb }

  makeDDB(tdb)

  if (!params.noquant) {
    mzml
    | dinosaur
    | set { dino }
  } else {
    channel.fromPath('NO__FILE')
    | set { dino }
  }

  mzml
  | combine(dino)
  | createSpectraLookup
    

  tdb
  | combine(makeDDB.out)
  | map { it + [prectol, fragtol, file("${baseDir}/assets/sage.json")] }
  | sagePrepare
  | combine(mzml)
  | combine(channel.fromPath(params.mods))
  | sage
  
  sage.out.perco
  | percolator
  | combine(sage.out.tsv)
  | combine(createSpectraLookup.out)
  | combine(tdb)
  | combine(makeDDB.out)
  | map { it + [params.noquant, params.psmconf, params.pepconf] }
  | createPSMTable
  | map { it + params.noquant }
  | peptidesProteinsReport
}
