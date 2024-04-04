/*
===================================
Longitudinal instrument QC pipeline
===================================
@Authors
Jorrit Boekel @glormph

Usage:
nextflow run qc.nf 
*/

nextflow.enable.dsl = 1

params.mzml = false
params.db = false
params.mods = 'data/labelfreemods.txt'
params.instrument = false
params.noquant = false
params.qval_modelthreshold = false
params.outdir = 'results'
params.prectol = '10.0ppm'
params.threadspercore = 1
params.filters = ''
params.options = ''

filters = params.filters.tokenize(';').collect() { x -> "--filter ${x}" }.join(' ')
options = params.options.tokenize(';').collect() {x -> "--${x}"}.join(' ')

instrument = [qe: 3, timstof: 2, velos:1][params.instrument]

process msconvert {

  cpus = 4 // FIXME 4 for TIMSTOF, XX for normal?

  input:
  path(raw) from Channel.fromPath(params.raw)

  output:
  file(outfile) into (mzml_msgf, mzml_mss, mzml_dino)

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

  when: !params.noquant

  input:
  file mzml from mzml_dino

  output:
  file "dinosaur.features.tsv" into dino_out

  script:
  """
  dinosaur --concurrency=${task.cpus * params.threadspercore} --outName="dinosaur" "${mzml}"
  """
}


process makeDDB {
 
  input:
  path(db) from Channel.of(params.db)
  output:
  file 'ddb' into ddb
  file(db) into (targetdb, psm_tdb)
  """
  msstitch makedecoy -i "$db" -o ddb --scramble prot_rev --ignore-target-hits
  """
}


process createSpectraLookup {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  file mzml from mzml_mss
  file dino from dino_out.ifEmpty('NA')
  
  output:
  file 'mslookup_db.sqlite' into spec_lookup

  script:
  """
  msstitch storespectra --spectra "${mzml}" --setnames 'QC'
  ${!params.noquant ? "msstitch storequant --dbfile mslookup_db.sqlite --dinosaur \"${dino}\" --spectra \"${mzml}\" --mztol 20.0 --mztoltype ppm --rttol 5.0" : ''}
  """
}


process msgfPlus {

  input:
  file mzml from mzml_msgf
  file('db.fa') from targetdb
  file mods  from Channel.fromPath(params.mods)

  output:
  tuple file('tpsms'), file('dpsms') into tdpsms
  
  """
  msgf_plus -Xmx16G -d "db.fa" -s "$mzml" -o "${mzml}.mzid" -thread ${task.cpus * params.threadspercore} -mod $mods -tda 1 -decoy decoy -t ${params.prectol} -ti -1,2 -m 0 -inst $instrument -e 1 -protocol 0 -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
  msgf_plus -Xmx8G edu.ucsd.msjava.ui.MzIDToTsv -i "${mzml}.mzid" -o out.tsv -showDecoy 1
  head -n 1 out.tsv > dpsms
  grep -v decoy_ out.tsv >> tpsms
  grep decoy_ out.tsv |grep -v \$'\\t'ENSP| grep -v '\\;ENSP' | grep -v \$'\\t''sp|' | grep -v '\\;sp|' >> dpsms
  """
}


process createPSMTable {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "tpsmtable" ? "psmtable.txt" : null }

  input:
  tuple file('tpsms'), file('dpsms') from tdpsms
  file 'lookup' from spec_lookup
  file db from psm_tdb 
  file ddb from ddb

  output:
  tuple file('tpsmtable'), file('dpsmtable'), file('tpeptides'), file('dpeptides'), file('tpsmlookup') into peptides_report

  """
  msstitch conffilt -i tpsms -o tfiltpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern '^QVal'
  msstitch conffilt -i dpsms -o dfiltpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern '^QVal'
  msstitch conffilt -i tfiltpsm -o tfiltpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PepQValue'
  msstitch conffilt -i dfiltpsm -o dfiltpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PepQValue'
  cp lookup tpsmlookup
  cp lookup dpsmlookup

  msstitch psmtable -i tfiltpep --dbfile tpsmlookup -o tpsmtable --addmiscleav --fasta "$db" ${params.noquant ? '' : '--ms1quant'} --proteingroup
  msstitch psmtable -i dfiltpep --dbfile dpsmlookup -o dpsmtable --fasta "$ddb" --proteingroup
  msstitch peptides -i tpsmtable -o tpeptides --scorecolpattern MSGFScore --spectracol 1 ${params.noquant ? '' : '--ms1quantcolpattern area'}
  msstitch peptides -i dpsmtable -o dpeptides --scorecolpattern MSGFScore --spectracol 1
  """
}


process peptidesProteinsReport {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple file('tpsms'), file('dpsms'), file('peptable.txt'), file('dpeptides'), file('db.sqlite') from peptides_report

  output:
  tuple file('peptable.txt'), file('prottable.txt'), file('qc.json') into outfiles

  script:
  scorecolpat = '^q-value'
  """
  sed -i 's/PepQValue/q-value/;s/QValue/PSM q-value/' "peptable.txt"
  sed -i 's/PepQValue/q-value/;s/QValue/PSM q-value/' dpeptides

  # score col is linearmodel_qval or q-value, but if the column only contains 0.0 or NA (no linear modeling possible due to only q<10e-04), we use MSGFScore instead
  tscol=\$(head -1 peptable.txt | tr '\\t' '\\n' | grep -n "${scorecolpat}" | cut -f 1 -d':')
  dscol=\$(head -1 dpeptides | tr '\\t' '\\n' | grep -n "${scorecolpat}" | cut -f 1 -d':')
  if [ -n "\$(cut -f \$tscol peptable.txt | tail -n+2 | egrep -v '(NA\$|0\\.0\$)')" ] && [ -n "\$(cut -f \$dscol dpeptides | tail -n+2 | egrep -v '(NA\$|0\\.0\$)')" ]
    then
      scpat="${scorecolpat}"
      logflag="--logscore"
    else
      scpat="MSGFScore"
      logflag=""
      echo 'Not enough q-values or linear-model q-values for peptides to calculate FDR, using MSGFScore instead.' >> warnings
  fi
  msstitch proteins -i peptable.txt --decoyfn dpeptides -o tprots --scorecolpattern "\$scpat" \$logflag ${params.noquant ? '' : '--ms1quant'} --psmtable tpsms
  msstitch conffilt -i tprots -o prottable.txt --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'q-value'

  protcol=\$(head -1 peptable.txt | tr '\\t' '\\n' | grep -n "Master" | cut -f 1 -d':')
  parse_output.py db.sqlite "\$(wc -l tpsms)" "\$(wc -l peptable.txt)" "\$(cut -f \$protcol peptable.txt | grep -v ';' | wc -l)" "\$(wc -l tprots)"
  """
}
