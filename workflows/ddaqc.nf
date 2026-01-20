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


process tdf2Mzml {
  container 'mfreitas/tdf2mzml'
  
  input:
  path(rawfile)

  output:
  path(mzmlfile)

  script:
  mzmlfile = "${rawfile.baseName}.mzml"
  """
  tdf2mzml.py -i $rawfile
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
  tuple path(db), path('config.json'), path(specfile), val(instrumenttype), path(mods)

  output:
  path('results.sage.pin'), emit: perco
  tuple path('results.sage.tsv'), val(instrumenttype), emit: tsv

  script:
  remove_scan_index_str = instrumenttype == 'bruker'
  """
  export RAYON_NUM_THREADS=${task.cpus}
  export SAGE_LOG=trace
  sage --disable-telemetry-i-dont-want-to-improve-sage --write-pin -f $db config.json $specfile
  ${remove_scan_index_str ? "sed -i 's/index=//' results.sage.tsv" : ''}
  ${remove_scan_index_str ? "sed -i 's/index=//' results.sage.pin" : ''}
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

  input:
  tuple path('perco'), path('psms'), val(instrumenttype), path(lookup), path(db), path(ddb), val(psmconf), val(pepconf)

  output:
  tuple path('tpsmtable'), path('dpsmtable'), path('tpeptides'), path('dpeptides'), emit: tables
  path('tpsmlookup'), emit: lookup

  script:
  add_scan_index_str = instrumenttype == 'bruker'
  """
  msstitch perco2psm --perco perco -i psms -o psms_perco --filtpsm ${psmconf} --filtpep ${pepconf}
  ${add_scan_index_str ? "mv psms_perco ppnoi && head -n1 ppnoi > psms_perco" : ''}
  ${add_scan_index_str ? "paste <(awk -F'\t' -v OFS='\t' '{print \$1,\$2,\$3,\$4,\$5,\"index=\"\$6}' ppnoi) <(cut -f7- ppnoi) | tail -n+2 >> psms_perco" : ''}

  msstitch split -i psms_perco --splitcol TD
  cp $lookup tpsmlookup
  cp $lookup dpsmlookup
  msstitch psmtable -i target.tsv --dbfile tpsmlookup -o tpsmtable --fasta "$db" --ms1quant --proteingroup --spectracol ${Utils.get_field_nr('target.tsv', 'filename')}
  msstitch psmtable -i decoy.tsv --dbfile dpsmlookup -o dpsmtable --fasta "$ddb" --proteingroup --spectracol ${Utils.get_field_nr('decoy.tsv', 'filename')}
  msstitch peptides -i tpsmtable -o tpeptides --scorecolpattern sage_discriminant --spectracol ${Utils.get_field_nr('tpsmtable', 'filename')} --ms1quantcolpattern area
  msstitch peptides -i dpsmtable -o dpeptides --scorecolpattern sage_discriminant --spectracol ${Utils.get_field_nr('dpsmtable', 'filename')}
  """
}


process proteinTables {
  container 'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'

  input:
  tuple path('tpsms'), path('dpsms'), path('peptable.txt'), path('dpeptides')

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
  msstitch proteins -i peptable.txt --decoyfn dpeptides -o tprots --scorecolpattern "\$scpat" \$logflag --ms1quant --psmtable tpsms
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
  psmconf
  pepconf

  main:
  // Set fragment and precursor tolerance, if not passed - defaults:
  def tolerances = [precursor: [timstof: 30, qe: 10, velos: 10, astral: 10], fragment: [timstof: 30, qe: 20, velos: 20, astral: 20]]
  prectol = prectol ?: tolerances.precursor[instrument]
  fragtol = prectol ?: tolerances.fragment[instrument]

  if (raw) {
    channel.fromPath(raw)
|view()
    | branch { 
      thermo: it.extension == 'raw' 
      bruker: it.extension == 'd'
      }
    | set { raw_c }

    raw_c.thermo
    | map { [it, instrument, pwiz_filters, pwiz_options] }
    | msconvert
    | map { [it, 'thermo'] }
    | set { thermo_mzml }

    raw_c.bruker
    | tdf2Mzml
    | map { [it, 'bruker'] }
    | concat(thermo_mzml)
    | set { spectrafn }

  } else if (mzml) {
    // Thermo mzML is passed, testing only
    channel.fromPath(mzml)
    | map { [it, instrument == 'timstof' ? 'bruker': 'thermo' ] }
    | set { spectrafn }

  } else {
    exit 1, 'Must either input a --raw file.raw, a --raw file.d, or an --mzml file.mzML'
  }

  channel
    .fromPath(db)
    .set { tdb }

  makeDDB(tdb)

  spectrafn
  | map { it[0] }
  | dinosaur
  | set { dino }

  spectrafn
  | map { it[0] }
  | combine(dino)
  | createNewSpectraLookup
    

  tdb
  | combine(makeDDB.out)
  | map { it + [prectol, fragtol, file("${baseDir}/assets/sage.json")] }
  | sagePrepare
  | combine(spectrafn)
  | combine(channel.fromPath(mods))
  | sage
  
  sage.out.perco
  | percolator
  | combine(sage.out.tsv)
  | combine(createNewSpectraLookup.out)
  | combine(tdb)
  | combine(makeDDB.out)
  | map { it + [psmconf, pepconf] }
  | createPSMTable

  createPSMTable.out.tables
  | proteinTables 

  emit:
  proteinTables.out
  | map { [it, 0].flatten() }
  | combine(createPSMTable.out.lookup)
}
