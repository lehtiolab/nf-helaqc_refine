/*
==============================
QUANTITATIVE PROTEOMICS PIPELINE
==============================
@Authors
Jorrit Boekel @glormph

Usage:

textfile:
/abspath/to/filefr01.mzML\tsetname\tplateid\tfrnr
/abspath/to/filefr02.mzML\tsetname\tplateid\tfrnr

*/

nf_required_version = '0.26.0'
if( ! nextflow.version.matches(">= ${nf_required_version}") ){
  println("Nextflow version too old, ${nf_required_version} required")
  exit(1)
}

/* SET DEFAULT PARAMS */
params.isobaric = false
params.activation = 'hcd'
params.mods = 'Mods.txt'
params.outdir = 'results'
params.normalize = false
params.genes = false
params.symbols = false
params.fastadelim = false
params.genefield = false
params.speclookup = false
params.quantlookup = false
params.hirief = false
params.onlypeptides = false
params.noquant = false

mods = file(params.mods)
tdb = file(params.tdb)
ddb = file(params.ddb)
martmap = file(params.martmap)
qcknitrpsms = file('qc/knitr_psms.Rhtml')
qcknitrplatepsms = file('qc/knitr_psms_perplate.Rhtml')
qcknitrnofrpsms = file('qc/knitr_psms_nofr.Rhtml')
qcknitrprot = file('qc/knitr_prot.Rhtml')
qcknitrnormfac = file('qc/knitr_iso_norm.Rhtml')
qctemplater = file('qc/collect.py')

piannotscript = file('scripts/peptide_pi_annotator.py')
trainingpep = file(params.pipep)


accolmap = [peptides: 12, proteins: 14, genes: 17, assoc: 18]

setdenoms = [:]
if (params.isobaric) {
  params.denoms.tokenize(' ').each{ it -> x=it.tokenize(':'); setdenoms.put(x[0], x[1..-1])}
}

normalize = (params.normalize && params.isobaric) ? true: false

hkconf = file('data/hardklor.conf')
activations = [hcd:'High-energy collision-induced dissociation', cid:'Collision-induced dissociation', etd:'Electron transfer dissociation']
activationtype = activations[params.activation]
plextype = params.isobaric ? params.isobaric.replaceFirst(/[0-9]+plex/, "") : 'false'
massshift = [tmt:0.0013, itraq:0.00125, false:0][plextype]
msgfprotocol = [tmt:4, itraq:2, false:0][plextype]
instrument = params.instrument ? params.instrument : false
msgfinstrument = [velos:1, qe:3, false:0][instrument]

/* PIPELINE START */

// Either feed an mzmldef file (tab separated lines with filepath\tsetname), or /path/to/\*.mzML
if (!params.mzmldef) {
Channel
  .fromPath(params.mzmls)
  .map { it -> [it, 'NA'] }
  .set { mzml_in }
} else {
Channel
  .from(file("${params.mzmldef}").readLines())
  .map { it -> it.tokenize('\t') }
  .set { mzml_in }
}

mzml_in
  .tap { sets }
  .map { it -> [file(it[0]), it[1], it[2] ? it[2] : it[1], it[3] ? it[3] : 'NA' ]} // create file, set plate to setname, and fraction to NA if there is none
  .tap { strips }
  .map { it -> [it[1], it[0].baseName.replaceFirst(/.*\/(\S+)\.mzML/, "\$1"), it[0], it[2], it[3]] }
  .tap{ mzmlfiles; mzml_isobaric; mzml_hklor; mzml_msgf }
  .count()
  .set{ amount_mzml }

sets
  .map{ it -> it[1] }
  .unique()
  .tap { setnames_psm } 
  .collect()
  .map { it -> [it] }
  .into { setnames_featqc; setnames_psmqc }

strips
  .map { it -> it[2] }
  .unique()
  .toList()
  .set { strips_for_deltapi }


process IsobaricQuant {

  container 'quay.io/biocontainers/openms:2.2.0--py27_boost1.64_0'

  when: params.isobaric && !params.quantlookup

  input:
  set val(setname), val(sample), file(infile), val(platename), val(fraction)from mzml_isobaric

  output:
  set val(sample), file("${infile}.consensusXML") into isobaricxml

  """
  IsobaricAnalyzer  -type $params.isobaric -in $infile -out "${infile}.consensusXML" -extraction:select_activation "$activationtype" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true 
  """
}


process hardklor {
  container 'quay.io/biocontainers/hardklor:2.3.0--0'
  when: !params.quantlookup && !params.noquant

  input:
  set val(setname), val(sample), file(infile), val(platename), val(fraction) from mzml_hklor
  file hkconf

  output:
  set val(sample), file('hardklor.out'), file(infile) into hk_out

  """
  cp $hkconf config
  echo "$infile" hardklor.out >> config
  hardklor config
  """
}


process kronik {

  container 'quay.io/biocontainers/kronik:2.20--0'
  when: !params.quantlookup && !params.noquant

  input:
  set val(sample), file('hardklor.out'), file(mzml) from hk_out 

  output:
  set val(sample), file("${sample}.kr"), file(mzml) into kronik_out

  """
  kronik -c 5 -d 3 -g 1 -m 8000 -n 600 -p 10 hardklor.out ${sample}.kr
  """
}

mzmlfiles
  .buffer(size: amount_mzml.value)
  .map { it.sort( {a, b -> a[1] <=> b[1]}) } // sort on sample for consistent .sh script in -resume
  .map { it -> [it.collect() { it[0] }, it.collect() { it[2] }, it.collect() { it[3] } ] } // lists: [sets], [mzmlfiles], [plates]
  .into { mzmlfiles_all; mzmlfiles_all_count }


if (params.speclookup && !params.quantlookup) {
  Channel
    .fromPath(params.speclookup)
    .into{ spec_lookup; countlookup }
} 
if (!params.speclookup && params.quantlookup) {
  Channel
    .fromPath(params.quantlookup)
    .into { spec_lookup; quant_lookup; countlookup }
} 


process createSpectraLookup {

  container 'quay.io/biocontainers/msstitch:2.7--py36_0'

  when: !(params.speclookup || params.quantlookup)

  input:
  set val(setnames), file(mzmlfiles), val(platenames) from mzmlfiles_all

  output:
  file 'mslookup_db.sqlite' into newspeclookup 

  script:
  """
  msslookup spectra -i ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
  """
}

isoquant_amount = params.isobaric ? amount_mzml.value : 1
isobaricxml
  .ifEmpty(['NA', 'NA', 'NA'])
  .buffer(size: isoquant_amount)
  .map { it.sort({a, b -> a[0] <=> b[0]}) }
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }] } // samples, isoxml
  .set { isofiles_sets }

kronik_out
  .ifEmpty(['NA', 'NA'])
  .buffer(size: amount_mzml.value)
  .map { it.sort({a, b -> a[0] <=> b[0]}) }
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }, it.collect() { it[2] }] } // samples, kronikout, mzml
  .set { krfiles_sets }


if (params.noquant && !(params.speclookup || params.quantlookup)) {
  newspeclookup
    .into { quant_lookup; spec_lookup; countlookup }
} else if (!(params.speclookup || params.quantlookup)) {
  newspeclookup
    .into { spec_lookup; countlookup }
}

process quantLookup {

  container 'quay.io/biocontainers/msstitch:2.7--py36_0'
  
  when: !params.quantlookup && !params.noquant

  input:
  file lookup from spec_lookup
  set val(isosamples), file(isofns) from isofiles_sets
  set val(krsamples), file(krfns), file(mzmls) from krfiles_sets

  output:
  file 'db.sqlite' into newquantlookup

  script:
  if (params.isobaric)
  """
  cp $lookup db.sqlite
  msslookup ms1quant --dbfile db.sqlite -i ${krfns.join(' ')} --spectra ${mzmls.join(' ')} --quanttype kronik --mztol 20.0 --mztoltype ppm --rttol 5.0 
  msslookup isoquant --dbfile db.sqlite -i ${isofns.join(' ')} --spectra ${isosamples.collect{ x -> x + '.mzML' }.join(' ')}
  """
  else
  """
  cp $lookup db.sqlite
  msslookup ms1quant --dbfile db.sqlite -i ${krfns.join(' ')} --spectra ${mzmls.join(' ')} --quanttype kronik --mztol 20.0 --mztoltype ppm --rttol 5.0 
  """
}


if (!params.quantlookup && !params.noquant) {
  newquantlookup
    .into { quant_lookup }
} 

mzmlfiles_all_count
  .merge(countlookup)
  .set { specfilein }


process countMS2perFile {

  container 'quay.io/biocontainers/msstitch:2.7--py36_0'

  input:
  set val(setnames), file(mzmlfiles), val(platenames), file(speclookup) from specfilein

  output:
  set val(setnames), file(mzmlfiles), val(platenames), file('amount_spectra_files') into specfilems2

  script:
  """
  sqlite3 $speclookup "SELECT mzmlfilename, COUNT(*) FROM mzml JOIN mzmlfiles USING(mzmlfile_id) JOIN biosets USING(set_id) GROUP BY mzmlfilename" > amount_spectra_files
  """
}


if (params.hirief) {
  specfilems2.set { scans_platecount }
} else {
  specfilems2
    .map { it -> [it[3], ['noplates']] }
    .into { scans_platecount; scans_result }
}


process countMS2sPerPlate {

  container 'biopython/biopython:latest'
  
  publishDir "${params.outdir}", mode: 'copy', overwrite: true 
  when: params.hirief

  input:
  set val(setnames), file(mzmlfiles), val(platenames), file('nr_spec_per_file') from scans_platecount

  output:
  set file('scans_per_plate'), val(splates) into scans_perplate

  script:
  splates = [setnames, platenames].transpose().collect() { "${it[0]}_${it[1]}" }
  """
  #!/usr/bin/env python
  platesets = [\"${splates.join('", "')}\"]
  platescans = {p: 0 for p in platesets}
  fileplates = {fn: p for fn, p in zip([\"${mzmlfiles.join('", "')}\"], platesets)}
  with open('nr_spec_per_file') as fp:
      for line in fp:
          fn, scans = line.strip('\\n').split('|')
          platescans[fileplates[fn]] += int(scans)
  with open('scans_per_plate', 'w') as fp:
      for plate, scans in platescans.items():
          fp.write('{}\\t{}\\n'.format(plate, scans))
  """
}

if (params.hirief) {
  scans_perplate.set { scans_result }
}

process concatTDFasta {
 
  container 'ubuntu:latest'

  input:
  file(tdb)
  file(ddb)

  output:
  file('db.fa') into concatdb

  script:
  """
  cat $tdb $ddb > db.fa
  """
}


process msgfPlus {
  container 'quay.io/biocontainers/msgf_plus:2017.07.21--py27_0'

  input:
  set val(setname), val(sample), file(x), val(platename), val(fraction) from mzml_msgf
  file(db) from concatdb
  file mods

  output:
  set val(setname), val(sample), file("${sample}.mzid") into mzids
  set val(setname), file("${sample}.mzid"), file('out.mzid.tsv'), val(platename), val(fraction) into mzidtsvs
  
  """
  msgf_plus -Xmx16G -d $db -s $x -o "${sample}.mzid" -thread 12 -mod $mods -tda 0 -t 10.0ppm -ti -1,2 -m 0 -inst ${msgfinstrument} -e 1 -protocol ${msgfprotocol} -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
  msgf_plus -Xmx3500M edu.ucsd.msjava.ui.MzIDToTsv -i "${sample}.mzid" -o out.mzid.tsv
  """
}

mzids
  .groupTuple()
  .set { mzids_2pin }


process percolator {

  container 'quay.io/biocontainers/percolator:3.1--boost_1.623'

  input:
  set val(setname), val(samples), file('mzid?') from mzids_2pin

  output:
  set val(setname), file('perco.xml') into percolated

  """
  echo $samples
  mkdir mzids
  count=1;for sam in ${samples.join(' ')}; do ln -s `pwd`/mzid\$count mzids/\${sam}.mzid; echo mzids/\${sam}.mzid >> metafile; ((count++));done
  msgf2pin -o percoin.xml -e trypsin -P "decoy_" metafile
  percolator -j percoin.xml -X perco.xml -N 500000 --decoy-xml-output -y
  """
}


mzidtsvs
  .groupTuple()
  .join(percolated)
  .set { mzperco }


process svmToTSV {

  container 'quay.io/biocontainers/msstitch:2.7--py36_0'

  input:
  set val(setname), file('mzident????'), file('mzidtsv????'), val(platenames), val(fractions), file(perco) from mzperco 

  output:
  set val(setname), val('target'), file('tmzidperco') into tmzidtsv_perco
  set val(setname), val('decoy'), file('dmzidperco') into dmzidtsv_perco

  script:
  """
#!/usr/bin/env python
import re
def count_missed_cleavage(full_pepseq, count=0):
    '''Regex .*[KR][^P] matches until the end and checks if there is a final
    charachter so this will not match the tryptic residue'''
    pepseq = re.sub('[\\+\\-]\\d*.\\d*', '', full_pepseq)
    match = re.match('.*[KR][^P]', pepseq)
    if match:
        count += 1
        return count_missed_cleavage(match.group()[:-1], count)
    else:
        return count

from glob import glob
mzidtsvfns = sorted(glob('mzidtsv*'))
mzidfns = sorted(glob('mzident*'))
fractions = [\"${fractions.join('", "')}\"]
plates = [\"${platenames.join('", "')}\"]
from app.readers import pycolator, xml, tsv, mzidplus
import os
ns = xml.get_namespace_from_top('$perco', None) 
psms = {p.attrib['{%s}psm_id' % ns['xmlns']]: p for p in pycolator.generate_psms('$perco', ns)}
decoys = {True: 0, False: 0}
for psm in sorted([(pid, float(p.find('{%s}svm_score' % ns['xmlns']).text), p) for pid, p in psms.items()], reverse=True, key=lambda x:x[1]):
    pdecoy = psm[2].attrib['{%s}decoy' % ns['xmlns']] == 'true'
    decoys[pdecoy] += 1
    psms[psm[0]] = {'decoy': pdecoy, 'svm': psm[1], 'qval': decoys[True]/decoys[False]}  # T-TDC
decoys = {'true': 0, 'false': 0}
for svm, pep in sorted([(float(x.find('{%s}svm_score' % ns['xmlns']).text), x) for x in pycolator.generate_peptides('$perco', ns)], reverse=True, key=lambda x:x[0]):
    decoys[pep.attrib['{%s}decoy' % ns['xmlns']]] += 1
    [psms[pid.text].update({'pepqval': decoys['true']/decoys['false']}) for pid in pep.find('{%s}psm_ids' % ns['xmlns'])]
oldheader = tsv.get_tsv_header(mzidtsvfns[0])
header = oldheader + ['percolator svm-score', 'PSM q-value', 'peptide q-value', 'Strip', 'Fraction', 'missed_cleavage']
with open('tmzidperco', 'w') as tfp, open('dmzidperco', 'w') as dfp:
    tfp.write('\\t'.join(header))
    dfp.write('\\t'.join(header))
    for fnix, mzidfn in enumerate(mzidfns):
        mzns = mzidplus.get_mzid_namespace(mzidfn)
        siis = (sii for sir in mzidplus.mzid_spec_result_generator(mzidfn, mzns) for sii in sir.findall('{%s}SpectrumIdentificationItem' % mzns['xmlns']))
        for specidi, psm in zip(siis, tsv.generate_tsv_psms(mzidtsvfns[fnix], oldheader)):
            # percolator psm ID is: samplename_SII_scannr_rank_scannr_charge_rank
            scan, rank = specidi.attrib['id'].replace('SII_', '').split('_')
            outpsm = {k: v for k,v in psm.items()}
            spfile = os.path.splitext(psm['#SpecFile'])[0]
            try:
                percopsm = psms['{fn}_SII_{sc}_{rk}_{sc}_{ch}_{rk}'.format(fn=spfile, sc=scan, rk=rank, ch=psm['Charge'])]
            except KeyError:
                continue
            outpsm.update({'percolator svm-score': percopsm['svm'], 'PSM q-value': percopsm['qval'], 'peptide q-value': percopsm['pepqval'], 'Strip': plates[fnix], 'Fraction': fractions[fnix], 'missed_cleavage': count_missed_cleavage(outpsm['Peptide'])})
            if percopsm['decoy']:
                dfp.write('\\n')
                dfp.write('\\t'.join([str(outpsm[k]) for k in header]))
            else:
                outpsm['Protein'] = ';'.join([x for x in outpsm['Protein'].split(';') if 'decoy_' not in x])
                tfp.write('\\n')
                tfp.write('\\t'.join([str(outpsm[k]) for k in header]))
  """
}

tmzidtsv_perco
  .concat(dmzidtsv_perco)
  .groupTuple(by: 1)
  .combine(quant_lookup)
  .set { prepsm }

if (params.hirief) {
  strips_for_deltapi
    .map { it -> [it, trainingpep, piannotscript]}
    .set { stripannot }
} else {
  strips_for_deltapi
    .map { it -> [it, trainingpep, piannotscript]}
    .set { stripannot }
}

process createPSMTable {

  container 'quay.io/biocontainers/msstitch:2.7--py36_0'
  
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {["target_psmlookup.sql", "target_psmtable.txt", "decoy_psmtable.txt"].contains(it) ? it : null}

  input:
  set val(setnames), val(td), file('psms?'), file('lookup') from prepsm
  set file(tdb), file(ddb), file(mmap) from Channel.value([tdb, ddb, martmap])
  set val(allstrips), file(trainingpep), file(piannotscript) from stripannot

  output:
  set val(td), file("${td}_psmtable.txt") into psm_result
  set val(td), file({setnames.collect() { it + '.tsv' }}) into setpsmtables
  set val(td), file("${td}_psmlookup.sql") into psmlookup

  script:
  """
  msspsmtable merge -i psms* -o psms.txt
  msspsmtable conffilt -i psms.txt -o filtpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PSM q-value'
  msspsmtable conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'peptide q-value'
  cp lookup psmlookup
  msslookup psms -i filtpep --dbfile psmlookup ${params.onlypeptides ? '' : "--fasta ${td == 'target' ? tdb : "${ddb} --decoy"}"} ${params.martmap ? "--map ${mmap}" : ''}
  msspsmtable specdata -i filtpep --dbfile psmlookup -o prepsms.txt
  ${!params.noquant ? "msspsmtable quant -i prepsms.txt -o qpsms.txt --dbfile psmlookup --precursor ${params.isobaric && td=='target' ? '--isobaric' : ''}" : 'mv prepsms.txt qpsms.txt'}
  sed 's/\\#SpecFile/SpectraFile/' -i qpsms.txt
  ${!params.onlypeptides ? "msspsmtable genes -i qpsms.txt -o gpsms --dbfile psmlookup" : ''}
  ${!params.onlypeptides ? "msslookup proteingroup -i qpsms.txt --dbfile psmlookup" : ''}
  ${!params.onlypeptides ? "msspsmtable proteingroup -i gpsms -o pgpsms --dbfile psmlookup" : 'mv qpsms.txt pgpsms'}
  ${params.hirief ? "python $piannotscript -i $trainingpep -p pgpsms --o dppsms --stripcolpattern Strip --pepcolpattern Peptide --fraccolpattern Fraction --strippatterns ${allstrips.join(' ')} --intercepts ${allstrips.collect() { params.strips[it].intercept}.join(' ')} --widths ${allstrips.collect() { params.strips[it].fr_width}.join(' ')} --ignoremods \'*\'" : ''}
  msspsmtable split -i ${params.hirief ? 'dppsms' : 'pgpsms'} --bioset
  mv ${params.hirief ? 'dppsms' : 'pgpsms'} ${td}_psmtable.txt
  mv psmlookup ${td}_psmlookup.sql
  """
}

setnames_psm
  .toList()
  .map { it -> [it.sort()]}
  .set { setlist_psm }

setpsmtables
  .map { it -> [it[0], it[1] instanceof java.util.List ? it[1] : [it[1]] ] }
  .map{ it -> [it[0], it[1].sort { a, b -> a.baseName.tokenize('.')[0] <=> b.baseName.tokenize('.')[0] }] } // names are setnames, sort on them then merge with sorted setnames
  .merge(setlist_psm)
  .transpose()
  .set { psm_pep }


process psm2Peptides {
  container 'quay.io/biocontainers/msstitch:2.7--py36_0'

  input:
  set val(td), file('psms'), val(setname) from psm_pep
  
  output:
  set val(td), val(setname), file('psms'), file('peptides') into prepep
  """
  msspeptable psm2pep -i psms -o peptides --scorecolpattern svm --spectracol 1 ${params.isobaric && td == 'target' ? "--isobquantcolpattern plex" : "" } ${!params.noquant ? "--ms1quantcolpattern area" : ""}
  """
}


process shuffleHeaderPeptidesMakeProttables {
 
  container 'ubuntu:latest'

  input:
  set val(td), val(setname), file(psms), file('peptides') from prepep

  output:
  set val(setname), val(td), file(psms), file('peptide_table.txt') into peptides
  set val(setname), val(td), file(psms), file('proteins'), val('proteins') into proteins
  set val(setname), val(td), file(psms), file('genes'), val('genes') into genes
  set val(setname), val(td), file(psms), file('symbols'), val('assoc') into symbols

  // protein, gene and symbol table are outputted regardless of necessity
  script:
  col = accolmap.peptides + 1  // psm2pep adds a column
  """
  paste <( cut -f ${col} peptides) <( cut -f 1-${col-1},${col+1}-500 peptides) > peptide_table.txt
  echo Protein accession |tee proteins genes symbols
  tail -n+2 psms|cut -f ${accolmap.proteins}|grep -v '\\;'|grep -v "^\$"|sort|uniq >> proteins
  tail -n+2 psms|cut -f ${accolmap.genes}|grep -v '\\;'|grep -v "^\$"|sort|uniq >> genes
  tail -n+2 psms|cut -f ${accolmap.assoc}|grep -v '\\;'|grep -v "^\$"|sort|uniq >> symbols
  """
}

proteins
  .tap { proteins_pre }
  .filter { it[1] == 'target' }
  .set { tprot_norm }

process ratioNormalizeProteins {
  container 'quay.io/biocontainers/msstitch:2.7--py36_0'

  when: normalize

  input:
  set val(setname), val(td), file('psms'), file('proteins'), val(acctype) from tprot_norm
  output:
  set val(setname), file('proteinratios') into proteinratios
  """
  msspsmtable isoratio -i psms -o proteinratios --protcol ${accolmap.proteins} --targettable proteins --isobquantcolpattern plex --minint 0.1 --denompatterns ${setdenoms[setname].join(' ')}
  """
}

tpep = Channel.create()
dpep = Channel.create()
peptides.choice(tpep, dpep) { it -> it[1] == 'target' ? 0 : 1}

if (normalize) {
  tpep
    .join(proteinratios)
    .concat(dpep)
    .set { peptable_quant }
} else {
  tpep
    .concat(dpep)
    .set { peptable_quant }
}

process postprodPeptideTable {

  container 'quay.io/biocontainers/msstitch:2.7--py36_0'
  
  input:
  set val(setname), val(td), file('psms'), file('peptides'), file(pratios) from peptable_quant

  output:
  set val(setname), val(td), file("${setname}_linmod"), file(pratios) into pepslinmod
  set val(setname), val('peptides'), val(td), file("${setname}_linmod") into peptides_out
  set val(setname), file('normratiosused') optional true into normratios

  script:
  if (params.isobaric && td=='target')
  """
  msspsmtable isoratio -i psms -o pepisoquant --targettable peptides --protcol ${accolmap.peptides} --isobquantcolpattern plex --minint 0.1 --denompatterns ${setdenoms[setname].join(' ')} ${normalize ? "--normalize median --norm-ratios $pratios" : ''} > normratiosused
  mv pepisoquant peptide_table.txt
  msspeptable modelqvals -i peptide_table.txt -o ${setname}_linmod --scorecolpattern svm --fdrcolpattern '^q-value'
  """
  else
  """
  mv peptides peptide_table.txt
  msspeptable modelqvals -i peptide_table.txt -o ${setname}_linmod --scorecolpattern svm --fdrcolpattern '^q-value'
  """
}


if (params.genes && params.symbols) { 
  pepslinmod.tap { pepsg; pepss }.concat(pepsg, pepss).set { pepslinmod_prot }
  proteins_pre.concat(genes).concat(symbols).set { pgstables } 
} else if (params.genes) { 
  pepslinmod.tap { pepsg }.concat(pepsg).set { pepslinmod_prot }
  proteins_pre.concat(genes).set { pgstables } 
} else { 
  pepslinmod.set { pepslinmod_prot }
  proteins_pre.set { pgstables }
}


pgstables
  .join(pepslinmod_prot, by: [0,1])
  .set { prepgs_in }

process prepProteinGeneSymbolTable {

  container 'quay.io/biocontainers/msstitch:2.7--py36_0'
 
  when: !params.onlypeptides

  input:
  set val(setname), val(td), file('psms'), file('proteins'), val(acctype), file('peplinmod'), file('pratios') from prepgs_in

  output:
  set val(setname), val(acctype), val(td), file('bestpeptides') into bestpep

  script:
  if (params.isobaric && td == 'target')
  """
  mssprottable ms1quant -i proteins -o protms1 --psmtable psms --protcol ${accolmap[acctype]}
  msspsmtable isoratio -i psms -o proteintable --protcol ${accolmap[acctype]} --targettable protms1 --isobquantcolpattern plex --minint 0.1 --denompatterns ${setdenoms[setname].join(' ')} ${normalize ? '--norm-ratios pratios --normalize median': ''}
  mssprottable bestpeptide -i proteintable -o bestpeptides --peptable peplinmod --scorecolpattern ${acctype == 'proteins' ? '\'^q-value\'' : '\'linear model\''} --logscore --protcol ${accolmap[acctype] + 1}
  """
  else
  """
  ${td == 'target' && !params.noquant ? "mssprottable ms1quant -i proteins -o proteintable --psmtable psms --protcol ${accolmap[acctype]}" : 'mv proteins proteintable'}
  mssprottable bestpeptide -i proteintable -o bestpeptides --peptable peplinmod --scorecolpattern ${acctype == 'proteins' ? '\'^q-value\'' : '\'linear model\''} --logscore --protcol ${accolmap[acctype] + 1}
  """
}

tbestpep = Channel.create()
dbestpep = Channel.create()
bestpep
  .groupTuple(by: [0,1])
  .transpose()
  .choice(tbestpep, dbestpep) { it[2] == 'target' ? 0 : 1 }


process proteinFDR {
  container 'quay.io/biocontainers/msstitch:2.7--py36_0'
  when: !params.onlypeptides
  input:
  set val(setname), val(acctype), val(td), file('tbestpep') from tbestpep
  set val(setname), val(acctype), val(td), file('dbestpep') from dbestpep
  set file(tfasta), file(dfasta) from Channel.value([tdb, ddb])

  output:
  set val(setname), val(acctype), file("${setname}_protfdr") into protfdrout
  script:
  if (acctype == 'genes')
  """
  mssprottable pickedfdr --picktype fasta --targetfasta $tfasta --decoyfasta $dfasta ${params.fastadelim ? "--fastadelim \'${params.fastadelim}\' --genefield ${params.genefield}" : ''} -i tbestpep --decoyfn dbestpep -o ${setname}_protfdr
  """
  else
  """
  mssprottable ${acctype == 'proteins' ? 'protfdr' : 'pickedfdr --picktype result'} -i tbestpep --decoyfn dbestpep -o ${setname}_protfdr
  """
}

peptides_out
  .filter { it[2] == 'target' }
  .map { it -> [it[0], it[1], it[3]] }
  .set { peptides_to_merge }

if (!params.onlypeptides) {
  peptides_to_merge
    .concat(protfdrout)
    .groupTuple(by: 1)
    .set { ptables_to_merge }
} else {
  peptides_to_merge
    .groupTuple(by: 1)
    .set { ptables_to_merge }
}

psmlookup
  .filter { it[0] == 'target' }
  .collect()
  .map { it[1] }
  .set { tlookup }

process proteinPeptideSetMerge {

  container 'quay.io/biocontainers/msstitch:2.7--py36_0'
  
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "proteintable" ? "${outname}_table.txt": null}

  input:
  set val(setnames), val(acctype), file(tables) from ptables_to_merge
  file(lookup) from tlookup
  
  output:
  set val(acctype), file('proteintable') into featuretables

  script:
  outname = (acctype == 'assoc') ? 'symbols' : acctype
  """
  cp $lookup db.sqlite
  msslookup ${acctype == 'peptides' ? 'peptides --fdrcolpattern \'^q-value\' --peptidecol' : 'proteins --fdrcolpattern \'q-value\' --protcol'} 1 --dbfile db.sqlite -i ${tables.join(' ')} --setnames ${setnames.join(' ')} ${!params.noquant ? "--ms1quantcolpattern area" : ""}  ${params.isobaric ? '--psmnrcolpattern quanted --isobquantcolpattern plex' : ''} ${acctype in ['genes', 'assoc'] ? "--genecentric ${acctype}" : ''}
  ${acctype == 'peptides' ? 'msspeptable build' : 'mssprottable build --mergecutoff 0.01'} --dbfile db.sqlite -o proteintable ${params.isobaric ? '--isobaric' : ''} ${!params.noquant ? "--precursor": ""} --fdr ${acctype in ['genes', 'assoc'] ? "--genecentric ${acctype}" : ''} ${params.onlypeptides ? "--noncentric" : ''}
  sed -i 's/\\#/Amount/g' proteintable
  """
}

psm_result
  .filter { it[0] == 'target' }
  .merge(scans_result)
  .map { it -> [it[0], it[1], it[2], it[3].unique()] }
  .set { targetpsm_result }


process psmQC {
  container 'r_qc_ggplot'
  input:
  set val(td), file('psms'), file('scans'), val(plates) from targetpsm_result
  val(setnames) from setnames_psmqc
  file(qcknitrplatepsms)
  file(qcknitrnofrpsms)
  file(qcknitrpsms)
  output:
  set val('psms'), file('knitr.html') into psmqccollect
  file('*_psms.html') optional true into platepsmscoll
  // TODO no proteins == no coverage for pep centric
  script:
  if (params.hirief)
  """
  setcol=`python -c 'with open("psms") as fp: h=next(fp).strip().split("\\t");print(h.index("Biological set")+1)'`
  stripcol=`python -c 'with open("psms") as fp: h=next(fp).strip().split("\\t");print(h.index("Strip")+1)'`
  paste -d _ <(cut -f \$setcol psms) <( cut -f \$stripcol psms) | sed 's/Biological set_Strip/plateID/' > platecol
  paste psms platecol > feats
  Rscript -e 'library(ggplot2); library(reshape2); library(knitr); nrsets=${setnames[0].size()}; feats = read.table("feats", header=T, sep="\\t", comment.char = "", quote = ""); amount_ms2 = read.table("scans"); knitr::knit2html("$qcknitrpsms", output="knitr.html"); ${plates.collect() { "plateid=\"${it}\"; knitr::knit2html(\"$qcknitrplatepsms\", output=\"${it}_psms.html\")"}.join('; ')}'
  rm knitr_psms.html
  """
  else
  """
  Rscript -e 'library(ggplot2); library(reshape2); library(knitr); nrsets=${setnames[0].size()}; feats = read.table("psms", header=T, sep="\\t", comment.char = "", quote = ""); amount_ms2 = read.table("scans", sep="|", header=F); knitr::knit2html(\"$qcknitrnofrpsms\", output=\"knitr.html\")'
  """
}

featuretables
  .merge(setnames_featqc)
  .set { featqcinput }

normratios
  .toList()
  .map { it -> [it.collect() { it[0] }.sort(), it.collect() { it[1] }.sort()] }
  .set{ allsetnormratios }


process normRatioParse {
  container 'ubuntu:latest'

  input:
  set val(setnames), file('norm?') from allsetnormratios

  output:
  file('normtable_sets') into normtable
  """
  count=1;for setn in ${setnames.join(' ')}; do echo "" >> norm"\$count" ; tail -n+2 norm"\$count" | sed \$'s/ - /\t'"\$setn"\$'\t/'; done >> normtable_sets
  """
}


process featQC {
  container 'r_qc_ggplot'
  input:
  set val(acctype), file('feats'), val(setnames) from featqcinput
  file('normtable') from normtable
  file(qcknitrprot)
  file(qcknitrnormfac)
  output:
  set val(acctype), file('knitr.html') into qccollect
  file('normalizefactors.html') optional true into normhtml

  script:
  """
  touch normalizefactors.html
  Rscript -e 'library(ggplot2); library(forcats); library(grid); library(reshape2); library(knitr); nrsets=${setnames.size()}; feats = read.table("feats", header=T, sep="\\t", comment.char = "", quote = ""); feattype="$acctype"; knitr::knit2html("$qcknitrprot", output="knitr.html"); ${normalize ? "normtable=\"normtable\"; knitr::knit2html(\"$qcknitrnormfac\", output=\"normalizefactors.html\");": ''}'
  """
}

qccollect
  .concat(psmqccollect)
  .toList()
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }] }
  .set { collected_feats_qc }

if (!params.hirief) {
  Channel.from([1]).set { platepsmscoll }
}
process collectQC {

  container 'r_qc_ggplot'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(acctypes), file('feat?') from collected_feats_qc
  file('norm.html') from normhtml
  file(ppsms) from platepsmscoll
  file(qctemplater)

  output:
  file('qc.html')

  script:
  if (params.hirief)
  """
  count=1; for ac in ${acctypes.join(' ')}; do mv feat\$count \$ac.html; ((count++)); done
  python $qctemplater $params.searchname hirief ${ppsms.join(' ')}
  """
  else
  """
  count=1; for ac in ${acctypes.join(' ')}; do mv feat\$count \$ac.html; ((count++)); done
  python $qctemplater $params.searchname nofrac
  """
}
