/*
===================================
Longitudinal instrument QC pipeline
===================================
@Authors
Jorrit Boekel @glormph

Usage:
nextflow run longqc.nf 

TODO list
- output missed cleavage?
- add dinosaur to get MS1 and LC peaks
*/

params.mzml = false
params.db = false
params.mods = false
params.instrument = false
params.qval_modelthreshold = false
params.outdir = 'results'

mzml = file(params.mzml)
db = file(params.db)
mods = file(params.mods)
instrument = [qe: 3, velos:1][params.instrument]


process hardklor {
  container 'quay.io/biocontainers/hardklor:2.3.0--0'
  input:
  file mzml
  file(hkconf) from Channel.fromPath("$baseDir/data/hardklor.conf").first()
  
  output:
  file 'hardklor.out' into hk_out

  """
  hardklor <(cat $hkconf <(echo "$mzml" hardklor.out))
  """
}

process kronik {

  container 'quay.io/biocontainers/kronik:2.20--0'
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
  file db
  output:
  file 'ddb' into ddb
  file(db) into targetdb
  """
  msslookup makedecoy -i "$db" -o ddb --scramble prot_rev --ignore-target-hits
  """
}

process createSpectraLookup {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  file mzml
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
  file mzml
  file('db.fa') from targetdb
  file mods

  output:
  set file('tpsms'), file('dpsms') into tdpsms
  
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
  set file('tpsms'), file('dpsms') from tdpsms
  file 'lookup' from spec_lookup
  file db
  file ddb from ddb

  output:
  set file('tpsmtable'), file('dpsmtable') into psmtable
  set file('tpeptides'), file('dpeptides') into prepeptable

  """
  msspsmtable conffilt -i tpsms -o tfiltpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern '^QVal'
  msspsmtable conffilt -i dpsms -o dfiltpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern '^QVal'
  msspsmtable conffilt -i tfiltpsm -o tfiltpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PepQValue'
  msspsmtable conffilt -i dfiltpsm -o dfiltpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PepQValue'
  cp lookup tpsmlookup
  cp lookup dpsmlookup
  msslookup psms -i tfiltpep --dbfile tpsmlookup --spectracol 1 --fasta "$db"
  msslookup psms -i dfiltpep --dbfile dpsmlookup --spectracol 1 --fasta $ddb
  msspsmtable specdata -i tfiltpep --dbfile tpsmlookup -o trtpsms
  msspsmtable quant -i trtpsms -o tquant --dbfile tpsmlookup --precursor
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
  set file('tprepep'), file('dprepep') from prepeptable
  set file('tpsms'), file('dpsms') from psmtable

  output:
  set file('peptable.txt'), file('prottable.txt') into outfiles

  script:
  scorecolpat = 'linear model'
  """
  paste <( cut -f 12 tprepep) <( cut -f 1-11,13-23 tprepep) > "peptable.txt"
  paste <( cut -f 10 dprepep) <( cut -f 1-9,11-20 dprepep) > dpeptides
  sed -i 's/PepQValue/q-value/;s/QValue/PSM q-value/' "peptable.txt"
  sed -i 's/PepQValue/q-value/;s/QValue/PSM q-value/' dpeptides
  echo Protein ID > tproteins
  echo Protein ID > dproteins
  tail -n+2 tpsms|cut -f 13 |grep -v '\\;'|grep -v "^\$"|sort|uniq >> tproteins
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
