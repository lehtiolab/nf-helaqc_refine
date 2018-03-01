/*
===================================
Longitudinal instrument QC pipeline
===================================
@Authors
Jorrit Boekel @glormph

Usage:
nextflow run longqc.nf 

FIXME list
= PSM table, one file, make set column
- mslookup spectra, ms1 area, one set, add to psms
- create peptide table from psms (one only)
- decoy psm/peptides/proteins needed, protein groups
- FIXME quant, then column recheck
  FIXME peptide table cut -f is wrong bc will get extra col with area quant
- create protein table with fdr
- run QC on missed cleavage
- output how to get it? SCP *but that has keys* SMB but then need analysis partition*
  - in celery though
*/

mzml = file(params.mzml)
db = file(params.db)
mods = file(params.mods)
hkconf = file('data/hardklor.conf')
instrument = [qe: 3, velos:1][params.instrument]


process hardklor {
  container 'quay.io/biocontainers/hardklor:2.3.0--0'
  input:
  file mzml
  file hkconf
  
  output:
  file 'hardklor.out' into hk_out

  """
  cp $hkconf config
  echo "$mzml" hardklor.out >> config
  hardklor config
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
  container 'biopython/biopython'
  input:
  file db
  output:
  file 'ddb' into ddb
  """
  #!/usr/bin/env python3
  from Bio import SeqIO
  with open('$db') as fp, open('ddb', 'w') as wfp:
    for decoy in (x[::-1] for x in SeqIO.parse(fp, 'fasta')):
      decoy.description = decoy.description.replace('ENS', 'XXX_ENS')
      decoy.id = 'XXX_{}'.format(decoy.id)
      SeqIO.write(decoy, wfp, 'fasta')
  
  """
}

process createSpectraLookup {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'
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

  container 'quay.io/biocontainers/msgf_plus:2017.07.21--py27_0'

  input:
  file mzml
  file db
  file mods

  output:
  file 'out.mzid.tsv' into mzidtsv
  
  """
  msgf_plus -Xmx16G -d $db -s "$mzml" -o "${mzml}.mzid" -thread ${task.cpus * 2} -mod $mods -tda 1 -t 10.0ppm -ti -1,2 -m 0 -inst $instrument -e 1 -protocol 5 -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
  msgf_plus -Xmx3500M edu.ucsd.msjava.ui.MzIDToTsv -i "${mzml}.mzid" -o out.mzid.tsv -showDecoy 1
  """
}


process determineTDPSMS {
  container 'ubuntu:latest'

  input:
  file 'psms' from mzidtsv

  output:
  set file('tpsms'), file('dpsms') into tdpsms
  """
  head -n 1 psms > dpsms
  grep -v XXX_ psms >> tpsms
  grep XXX_ psms|grep -v \$'\\t'ENSP| grep -v '\\;ENSP' >> dpsms
  """
}


process createPSMPeptideTable {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'
  
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
  msslookup psms -i tfiltpep --dbfile tpsmlookup --spectracol 1 --fasta $db
  msslookup psms -i dfiltpep --dbfile dpsmlookup --spectracol 1 --fasta $ddb
  msspsmtable specdata -i tfiltpep --dbfile tpsmlookup -o trtpsms
  msspsmtable quant -i trtpsms -o tquant --dbfile tpsmlookup --precursor
  msslookup proteingroup -i tquant --dbfile tpsmlookup
  msslookup proteingroup -i dfiltpep --dbfile dpsmlookup
  msspsmtable proteingroup -i tquant -o tpsmtable --dbfile tpsmlookup
  msspsmtable proteingroup -i dfiltpep -o dpsmtable --dbfile dpsmlookup
  msspeptable psm2pep -i tpsmtable -o tpeptides --scorecolpattern MSGFScore --spectracol 1 --ms1quantcolpattern area
  msspeptable psm2pep -i dpsmtable -o dpeptides --scorecolpattern MSGFScore --spectracol 1
  """
}

psmtable
  .into { psm2prottable; psmprotquant }

process cutPasteReplacePeptideProteinTable{

  container 'ubuntu:latest'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "tpep" ? "peptable.txt" : null }

  input:
  set file('tprepep'), file('dprepep') from prepeptable
  set file('tpsm'), file('dpsm') from psm2prottable 

  output:
  set file('tpep'), file('dpep') into peptable
  set file('tprot'), file('dprot') into proteinlist

  """
  paste <( cut -f 12 tprepep) <( cut -f 1-11,13-23 tprepep) > tpep
  paste <( cut -f 10 dprepep) <( cut -f 1-9,11-20 dprepep) > dpep
  sed -i 's/PepQValue/q-value/;s/QValue/PSM q-value/' tpep
  sed -i 's/PepQValue/q-value/;s/QValue/PSM q-value/' dpep
  echo Protein accession > tprot
  echo Protein accession > dprot
  tail -n+2 tpsm|cut -f 13 |grep -v '\\;'|grep -v "^\$"|sort|uniq >> tprot
  tail -n+2 dpsm|cut -f 11 |grep -v '\\;'|grep -v "^\$"|sort|uniq >> dprot
  """
}

process createProteinTable {
  
  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set file('tproteins'), file('dproteins') from proteinlist
  set file('tpsms'), file('dpsms') from psmprotquant
  set file('tpeptides'), file('dpeptides') from peptable

  output:
  file 'prottable.txt' into proteins

  """
  mssprottable ms1quant -i tproteins -o ms1quant --psmtable tpsms --protcol 13 
  msspeptable modelqvals -i tpeptides -o tlinmodpep --scorecolpattern 'MSGFScore' --fdrcolpattern '^q-value'
  msspeptable modelqvals -i dpeptides -o dlinmodpep --scorecolpattern 'MSGFScore' --fdrcolpattern '^q-value'
  mssprottable bestpeptide -i ms1quant -o tbestpep --peptable tlinmodpep --scorecolpattern 'linear model' --logscore --protcolpattern 'Master' 
  mssprottable bestpeptide -i dproteins -o dbestpep --peptable dlinmodpep --scorecolpattern 'linear model' --logscore --protcolpattern 'Master'
  mssprottable pickedfdr --picktype result -i tbestpep --decoyfn dbestpep -o fdrprots
  msspsmtable conffilt -i fdrprots -o prottable.txt --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'q-value'
  """
}
