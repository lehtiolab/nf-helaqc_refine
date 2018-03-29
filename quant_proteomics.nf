/*
==============================
IPAW: HiRIEF II varDB pipeline
==============================
@Authors
Jorrit Boekel @glormph

Usage:

textfile:
/abspath/to/filefr01.mzML\tsetname\tfrnr
/abspath/to/filefr02.mzML\tsetname\tfrnr

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

mods = file(params.mods)
tdb = file(params.tdb)
ddb = file(params.ddb)
martmap = file(params.martmap)

accolmap = [peptides: 11, proteins: 13, genes: 16, assoc: 17]

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
  .map { it -> [file(it[0]), it[1]]}
  .map { it -> [it[1], it[0].baseName.replaceFirst(/.*\/(\S+)\.mzML/, "\$1"), it[0]] }
  .tap{ mzmlfiles; mzml_isobaric; mzml_hklor; mzml_msgf }
  .count()
  .set{ amount_mzml }

sets
  .map{ it -> it[1] }
  .unique()
  .tap { setnames_psm }
  .collect()
  .subscribe { println "Detected setnames: ${it.join(', ')}" }



process IsobaricQuant {

  container 'quay.io/biocontainers/openms:2.2.0--py27_boost1.64_0'

  when: params.isobaric

  input:
  set val(setname), val(sample), file(infile) from mzml_isobaric

  output:
  set val(sample), file("${infile}.consensusXML") into isobaricxml

  """
  IsobaricAnalyzer  -type $params.isobaric -in $infile -out "${infile}.consensusXML" -extraction:select_activation "$activationtype" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true 
  """
}


process hardklor {
  container 'quay.io/biocontainers/hardklor:2.3.0--0'

  input:
  set val(setname), val(sample), file(infile) from mzml_hklor
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
  .map { it -> [it.collect() { it[0] }, it.collect() { it[2] }] } // lists: [sets], [mzmlfiles]
  .set{ mzmlfiles_all }

/*
// This is not sorted which makes -resume problematic due to chagning bash command.sh script
mzml_lookup
  .reduce([sets:[], files:[]]) { a, b -> a.sets.add(b[3]); a.files.add(b[2]); return a }
  .map { it -> [it.files, it.sets]}
  .set { mzmlfns_sets }
*/

process createSpectraLookup {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  set val(setnames), file(mzmlfiles) from mzmlfiles_all

  output:
  file 'mslookup_db.sqlite' into spec_lookup

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


process quantLookup {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  file lookup from spec_lookup
  set val(isosamples), file(isofns) from isofiles_sets
  set val(krsamples), file(krfns), file(mzmls) from krfiles_sets

  output:
  file 'db.sqlite' into quant_lookup

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
  // FIXME add plate and fr to mzidtsv in this process, get them from input channel!
  container 'quay.io/biocontainers/msgf_plus:2017.07.21--py27_0'

  input:
  set val(setname), val(sample), file(x) from mzml_msgf
  file(db) from concatdb
  file mods

  output:
  set val(setname), val(sample), file("${sample}.mzid") into mzids
  set val(setname), file("${sample}.mzid"), file('out.mzid.tsv') into mzidtsvs
  //set val(setname), val(sample), file("${sample}.mzid"), val(td) into mzids
  //set val("${setname}_${td}"), val(sample), file('out.mzid.tsv') into mzidtsvs
  
  """
  msgf_plus -Xmx16G -d $db -s $x -o "${sample}.mzid" -thread 12 -mod $mods -tda 0 -t 10.0ppm -ti -1,2 -m 0 -inst 3 -e 1 -protocol ${msgfprotocol} -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
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

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  set val(setname), file('mzident?'), file('mzidtsv?'), file(perco) from mzperco 

  output:
  set val(setname), val('target'), file('tmzidperco') into tmzidtsv_perco
  set val(setname), val('decoy'), file('dmzidperco') into dmzidtsv_perco

  script:
  """
#!/usr/bin/env python
from glob import glob
mzidtsvfns = sorted(glob('mzidtsv*'))
mzidfns = sorted(glob('mzident*'))
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
header = oldheader + ['percolator svm-score', 'PSM q-value', 'peptide q-value']
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
            outpsm.update({'percolator svm-score': percopsm['svm'], 'PSM q-value': percopsm['qval'], 'peptide q-value': percopsm['pepqval']})
            if percopsm['decoy']:
                dfp.write('\\n')
                dfp.write('\\t'.join([str(outpsm[k]) for k in header]))
            else:
                tfp.write('\\n')
                tfp.write('\\t'.join([str(outpsm[k]) for k in header]))
  """
}

tmzidtsv_perco
  .concat(dmzidtsv_perco)
  .groupTuple(by: 1)
  .combine(quant_lookup)
  .set { prepsm }

process createPSMTable {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  set val(setnames), val(td), file('psms?'), file('lookup') from prepsm
  set file(tdb), file(ddb), file(mmap) from Channel.value([tdb, ddb, martmap])

  output:
  set val(td), file("${td}_psmtable.txt") into psmtables
  set val(td), file({setnames.collect() { it + '.tsv' }}) into setpsmtables
  set val(td), file('psmlookup') into psmlookup

  script:
  """
  msspsmtable merge -i psms* -o psms.txt
  msspsmtable conffilt -i psms.txt -o filtpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PSM q-value'
  msspsmtable conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'peptide q-value'
  cp lookup psmlookup
  msslookup psms -i filtpep --dbfile psmlookup --fasta ${td == 'target' ? tdb : "${ddb} --decoy"} --map ${mmap} 
  msspsmtable specdata -i filtpep --dbfile psmlookup -o prepsms.txt
  msspsmtable quant -i prepsms.txt -o qpsms.txt --dbfile psmlookup --precursor ${params.isobaric && td=='target' ? "--isobaric" : ""}
  sed 's/\\#SpecFile/SpectraFile/' -i qpsms.txt
  msspsmtable genes -i qpsms.txt -o gpsms --dbfile psmlookup
  msslookup proteingroup -i qpsms.txt --dbfile psmlookup
  msspsmtable proteingroup -i gpsms -o pgpsms --dbfile psmlookup
  msspsmtable split -i pgpsms --bioset
  mv pgpsms ${td}_psmtable.txt
  """
}

setnames_psm
  .toList()
  .map { it -> [it.sort()]}
  .set { setlist_psm }

setpsmtables
  .map{ it -> [it[0], it[1].sort { a, b -> a.baseName.tokenize('.')[0] <=> b.baseName.tokenize('.')[0] }] } // names are setnames, sort on them then merge with sorted setnames
  .merge(setlist_psm)
  .transpose()
  .set { psm_pep }


process psm2Peptides {
  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  set val(td), file('psms'), val(setname) from psm_pep
  
  output:
  set val(td), val(setname), file('psms'), file('peptides') into prepep
  """
  msspeptable psm2pep -i psms -o peptides --scorecolpattern svm --spectracol 1 ${params.isobaric && td == 'target' ? "--isobquantcolpattern plex" : "" } --ms1quantcolpattern area
  """
}


process shuffleHeaderPeptidesMakeProttables {
 
  container 'ubuntu:latest'

  input:
  set val(td), val(setname), file(psms), file('peptides') from prepep

  output:
  set val(setname), val(td), file(psms), file('peptide_table.txt'), file('proteins') into peptides
  set val(setname), val(td), file(psms), file('proteins'), val('proteins') into proteins
  set val(setname), val(td), file(psms), file('genes'), val('genes') into genes
  set val(setname), val(td), file(psms), file('symbols'), val('assoc') into symbols

  // gene and symbol table are outputted regardless of necessity
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
  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

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

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  set val(setname), val(td), file('psms'), file('peptides'), file('proteins'), file(pratios) from peptable_quant

  output:
  set val(setname), val(td), file("${setname}_linmod"), file(pratios) into pepslinmod
  set val(setname), val('peptides'), val(td), file("${setname}_linmod") into peptides_out

  script:
  if (params.isobaric && td=='target')
  """
  msspsmtable isoratio -i psms -o pepisoquant --targettable peptides --protcol ${accolmap.peptides} --isobquantcolpattern plex --minint 0.1 --denompatterns ${setdenoms[setname].join(' ')} ${normalize ? "--normalize median --norm-ratios $pratios" : ''}
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
  .view()
  .set { prepgs_in }

process prepProteinGeneSymbolTable {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

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
  ${td == 'target' ? 'mssprottable ms1quant -i proteins -o proteintable --psmtable psms --protcol ${accolmap[acctype]}' : 'mv proteins proteintable'}
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
  container 'quay.io/biocontainers/msstitch:2.5--py36_0'
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
  .concat(protfdrout)
  .groupTuple(by: 1)
  .set { ptables_to_merge }

psmlookup
  .filter { it[0] == 'target' }
  .collect()
  .map { it[1] }
  .set { tlookup }

process proteinPeptideSetMerge {

  container 'quay.io/biocontainers/msstitch:2.5--py36_0'

  input:
  set val(setnames), val(acctype), file(tables) from ptables_to_merge
  file(lookup) from tlookup

//  #msslookup proteins --protcol 1 --genecentric genes
//  #msslookup proteins --protcol 1 --genecentric assoc
//  mssprottable build --genecentric genes
//  mssprottable build --genecentric assoc
//  msspeptable build
//  msspeptable build --genecentric
//  msspeptable build --noncentric 
  """
  cp $lookup db.sqlite
  msslookup ${acctype == 'peptides' ? 'peptides --fdrcolpattern \'^q-value\' --peptidecol' : 'proteins --fdrcolpattern \'q-value\' --protcol'} 1 --dbfile db.sqlite -i ${tables.join(' ')} --setnames ${setnames.join(' ')} --ms1quantcolpattern area ${params.isobaric ? '--psmnrcolpattern quanted --isobquantcolpattern plex' : ''} ${acctype in ['genes', 'assoc'] ? "--genecentric ${acctype}" : ''}
              
  ${acctype == 'peptides' ? 'msspeptable build' : 'mssprottable build --mergecutoff 0.01'} --dbfile db.sqlite -o proteintable ${params.isobaric ? '--isobaric' : ''} --precursor --fdr ${acctype in ['genes', 'assoc'] ? "--genecentric ${acctype}" : ''}

  """
}
