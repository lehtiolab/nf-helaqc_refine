/*
===================================
Longitudinal instrument QC pipeline
===================================
@Authors
Jorrit Boekel @glormph

Usage:
nextflow run longqc.nf 

TODO list
- add dinosaur to get MS1 and LC peaks
*/

params.mzml = false
params.db = false
params.mods = false
params.instrument = false
params.qval_modelthreshold = false
params.outdir = 'results'
params.filters = ''
params.options = ''

filters = params.filters.tokenize(';').collect() { x -> "--filter ${x}" }.join(' ')
options = params.options.tokenize(';').collect() {x -> "--${x}"}.join(' ')

raw = file(params.raw)
//db = file(params.db)
mods = file(params.mods)
instrument = [qe: 3, velos:1][params.instrument]


process msconvert {
  container 'chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:3.0.20066-729ef9c41'

  cpus = 4 // FIXME 4 for TIMSTOF, XX for normal?

  input:
  file raw

  output:
  file(outfile) into (mzml_hk, mzml_msgf, mzml_mss)

  script:
  outfile = "${raw.baseName}.mzML"
  """
  # Resolve directory if necessary, pwiz tries to read NF soft links as if they are files, which
  # does not work in case of directory
  ${raw.isDirectory() ?  "mv ${raw} tmpdir && cp -rL tmpdir ${raw}" : ''}
  wine msconvert ${raw} ${filters} ${options}
  """
}


process hardklor {

  input:
  file mzml from mzml_hk
  file(hkconf) from Channel.fromPath("$baseDir/data/hardklor.conf").first()
  
  output:
  file 'hardklor.out' into hk_out

  """
  hardklor <(cat $hkconf <(echo "$mzml" hardklor.out))
  """
}

process kronik {

  input:
  file 'hardklor.out' from hk_out 
  output:
  file 'kronik.out' into kronik_out
  """
  kronik -c 5 -d 3 -g 1 -m 8000 -n 600 -p 10 hardklor.out kronik.out
  """
}

process makeDDB {
 
  input:
  path(db) from Channel.of(params.db)
  output:
  file 'ddb' into ddb
  file(db) into (targetdb, psm_tdb)
  """
  msslookup makedecoy -i "$db" -o ddb --scramble prot_rev --ignore-target-hits
  """
}


process createSpectraLookup {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  file mzml from mzml_mss
  file kronik from kronik_out
  
  output:
  file 'mslookup_db.sqlite' into spec_lookup

  script:
  """
  msslookup spectra -i "${mzml}" --setnames 'QC'
  msslookup ms1quant --dbfile mslookup_db.sqlite -i ${kronik} --spectra "${mzml}" --quanttype kronik --mztol 20.0 --mztoltype ppm --rttol 5.0 
  """
}


process msgfPlus {

  input:
  file mzml from mzml_msgf
  file('db.fa') from targetdb
  file mods

  output:
  tuple file('tpsms'), file('dpsms') into tdpsms
  
  """
  msgf_plus -Xmx16G -d "db.fa" -s "$mzml" -o "${mzml}.mzid" -thread ${task.cpus * 2} -mod $mods -tda 1 -decoy decoy -t 10.0ppm -ti -1,2 -m 0 -inst $instrument -e 1 -protocol 0 -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
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
  tuple file('tpsmtable'), file('dpsmtable') into psmtable
  tuple file('tpeptides'), file('dpeptides') into prepeptable

  """
  msspsmtable conffilt -i tpsms -o tfiltpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern '^QVal'
  msspsmtable conffilt -i dpsms -o dfiltpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern '^QVal'
  msspsmtable conffilt -i tfiltpsm -o tfiltpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PepQValue'
  msspsmtable conffilt -i dfiltpsm -o dfiltpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PepQValue'
  cp lookup tpsmlookup
  cp lookup dpsmlookup
  msslookup psms -i tfiltpep --dbfile tpsmlookup --spectracol 1 --fasta "$db"
  msslookup psms -i dfiltpep --dbfile dpsmlookup --spectracol 1 --fasta $ddb
  msspsmtable specdata -i tfiltpep --dbfile tpsmlookup -o trtpsms --addmiscleav
  msspsmtable quant -i trtpsms -o tquant --dbfile tpsmlookup --precursor
# FIXME only one command for pgrouping
  msslookup proteingroup -i tquant --dbfile tpsmlookup || exit 3
  msslookup proteingroup -i dfiltpep --dbfile dpsmlookup || exit 3
  msspsmtable proteingroup -i tquant -o tpsmtable --dbfile tpsmlookup
  msspsmtable proteingroup -i dfiltpep -o dpsmtable --dbfile dpsmlookup
  msspeptable psm2pep -i tpsmtable -o tpeptides --scorecolpattern MSGFScore --spectracol 1 --ms1quantcolpattern area
  msspeptable psm2pep -i dpsmtable -o dpeptides --scorecolpattern MSGFScore --spectracol 1
  """
}


process createPeptideProteinTable{

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple file('peptable.txt'), file('dpeptides') from prepeptable
  tuple file('tpsms'), file('dpsms') from psmtable

  output:
  tuple file('peptable.txt'), file('prottable.txt') into outfiles

  script:
  scorecolpat = 'linear model'
  """
  sed -i 's/PepQValue/q-value/;s/QValue/PSM q-value/' "peptable.txt"
  sed -i 's/PepQValue/q-value/;s/QValue/PSM q-value/' dpeptides

  echo Protein ID > tproteins
  echo Protein ID > dproteins
  tail -n+2 tpsms|cut -f 15 |grep -v '\\;'|grep -v "^\$"|sort|uniq >> tproteins
  tail -n+2 dpsms|cut -f 11 |grep -v '\\;'|grep -v "^\$"|sort|uniq >> dproteins
  mssprottable ms1quant -i tproteins -o ms1quant --psmtable tpsms --protcol 13 
  msspeptable modelqvals -i "peptable.txt" -o tlinmodpep --scorecolpattern 'MSGFScore' --fdrcolpattern '^q-value' ${params.qval_modelthreshold ? "--qvalthreshold ${params.qval_modelthreshold}" : ''}
  msspeptable modelqvals -i dpeptides -o dlinmodpep --scorecolpattern 'MSGFScore' --fdrcolpattern '^q-value' ${params.qval_modelthreshold ? "--qvalthreshold ${params.qval_modelthreshold}" : ''}

  # score col is linearmodel_qval or q-value, but if the column only contains 0.0 or NA (no linear modeling possible due to only q<10e-04), we use MSGFScore instead
  tscol=\$(head -1 tlinmodpep | tr '\\t' '\\n' | grep -n "${scorecolpat}" | cut -f 1 -d':')
  dscol=\$(head -1 tlinmodpep | tr '\\t' '\\n' | grep -n "${scorecolpat}" | cut -f 1 -d':')
  if [ -n "\$(cut -f \$tscol tlinmodpep | tail -n+2 | egrep -v '(NA\$|0\\.0\$)')" ] && [ -n "\$(cut -f \$dscol dlinmodpep | tail -n+2 | egrep -v '(NA\$|0\\.0\$)')" ]
    then
      scpat="${scorecolpat}"
      logflag="--logscore"
    else
      scpat="MSGFScore"
      logflag=""
      echo 'Not enough q-values or linear-model q-values for peptides to calculate FDR, using MSGFScore instead.' >> warnings
  fi
  mssprottable bestpeptide -i ms1quant -o tbestpep --peptable tlinmodpep --scorecolpattern "\$scpat" \$logflag --protcolpattern 'Master' 
  mssprottable bestpeptide -i dproteins -o dbestpep --peptable dlinmodpep --scorecolpattern "\$scpat" \$logflag --protcolpattern 'Master'
  mssprottable pickedfdr --picktype result -i tbestpep --decoyfn dbestpep -o fdrprots
  msspsmtable conffilt -i fdrprots -o prottable.txt --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'q-value'
  """
}
