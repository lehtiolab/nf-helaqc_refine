include { msconvert; createNewSpectraLookup } from '../modules.nf' 


// AlLibrary can be made like this:
// diann --threads 16 --fasta tdb.fa  --gen-spec-lib --fasta-search --out-lib libfile --dir rawfile_dir/ --var-mod 'UniMod:35,15.994915,M' --var-mod 'UniMod:4,57.021464,C' --var-mods 2


process generateTrackPeptideLibrary {

  container "ghcr.io/lehtiolab/nfhelaqc:${workflow.manifest.version}"
  
  input:
  tuple path(baselib), val(trackedpeptides), path(tdb)
  
  output:
  tuple path('combined.parquet'), path('fulldb.fa')

  script:
  """
  ${trackedpeptides.collect {
  "echo '>${it}' >> tp.fa && echo ${it} >> tp.fa"
  }.join('\n')}
  diann-linux --threads 16 \
     --fasta tp.fa  \
     --gen-spec-lib \
     --fasta-search \
     --out-lib trackpeplibfile \
     --var-mod 'UniMod:35,15.994915,M' \
     --var-mod 'UniMod:4,57.021464,C' \
     --var-mods 2 \
     --predictor
  diann-linux --gen-spec-lib --lib trackpeplibfile.predicted.speclib --out-lib tp.parquet
  diann-linux --gen-spec-lib --lib $baselib --lib tp.parquet --out-lib combined.parquet
  cat tp.fa $tdb > fulldb.fa
  """
}

process DiaNN {

  /*
  q-value cutoff seems 0.01 by default
  */
  
  container "ghcr.io/lehtiolab/nfhelaqc:${workflow.manifest.version}"
  
  input:
  tuple path(raw), path(lib), path(tdb), val(ms1acc), val(ms2acc)
  
  output:
  tuple path('out.txt'), path('out.stats.tsv'), emit: tsv
  path('out.parquet'), emit: pq

  script:
  """
  diann-linux --threads 16 \
    --f $raw \
     --lib $lib \
    --fasta $tdb \
    --window 8 \
    --mass-acc-ms1 $ms1acc \
    --mass-acc $ms2acc \
    --missed-cleavages 2 \
    --var-mods 2 \
    --out out.txt
  """
    //--var-mod 'UniMod:35,15.994915,M' \
    //--var-mod 'UniMod:4,57.021464,C' \
}

process prepareNumbers {
  container 'quay.io/biocontainers/bioawk:1.0--he4a0461_9'

  input:
  tuple path('precursors'), path('stats')

  output:
  tuple path('peptides'), eval('cat nrprots'), eval('cat nrpsms'), eval('wc -l < <(tail -n+2 peptides)'), eval('cat fwhmscans')

  script:
  """
  cut -f ${Utils.get_field_nr("stats", "Proteins.Identified")} stats | tail -n1 > nrprots
  cut -f ${Utils.get_field_nr("stats", "Precursors.Identified")} stats | tail -n1 > nrpsms
  cut -f ${Utils.get_field_nr("stats", "FWHM.Scans")} stats | tail -n1 > fwhmscans
 
  bioawk -v s=${Utils.get_field_nr("precursors", "Modified.Sequence")} -v b=${Utils.get_field_nr("precursors", "Stripped.Sequence")} -v g=${Utils.get_field_nr("precursors", "Genes")} -v m=${Utils.get_field_nr("precursors", "Ms1.Area")} -t \
    '{print \$s,\$b,\$g,\$m }' precursors > pepfields
  head -n1 pepfields > peptides
  # Sort to get highest MS1 (col3, numeric, reverse) for each peptide
  # awk does the uniqueing since BSD? sort in container did not use the column 
  # nr 1 only for some reason to unique. Awk explanation in link
  # (field nr1 !not in _variable, add++)
  # https://stackoverflow.com/questions/1915636/is-there-a-way-to-uniq-by-column
  tail -n+2 pepfields | sort -k1b,1 -k4,4nr | bioawk -t '!_[\$1]++' >> peptides
  """
}

// MS1 area, IM, Precursor quant (ms2) seem quite different between charged peptides, whereas RT is more similar ?


workflow DIAQC {

  take:
  raw
  mzml
  library
  db 
  instrument
  trackedpeptides

  main:
  
  if (raw) {
    channel.fromPath(raw)
    | branch { 
      thermo: it.extension == 'raw' 
      bruker: it.extension == 'd'
      }
    | set { raw_c }

    raw_c.thermo
    | map { [it, instrument, false, false] }
    | msconvert
    | concat(raw_c.bruker)
    | set { diann_in }

    raw_c.bruker
    | set { raw_bruker }

    msconvert.out
    | set { mzml_c }

  } else if (mzml) {
    raw_bruker = channel.empty()

    channel.fromPath(mzml)
    | set { mzml_c }

    mzml_c
    | set { diann_in }

  } else {
    exit 1, 'Must either input a --raw file.raw, a --raw file.d, or an --mzml file.mzML'
  }

  ms1acc = [timstof: 20, orbitrap: 10, velos: 10, qe: 10][instrument]
  ms2acc = 20

  mzml_c
  | map { [it, file('NO__FILE')] }
  | createNewSpectraLookup
  | concat(raw_bruker)
  | set { scandb }
   
  lib_c = channel.fromPath(library)
  db_c = channel.fromPath(db)

  if (trackedpeptides) {
    trackpeps = trackedpeptides.collect { it.tokenize('_')[0] }.toList()
    lib_c
    | map { [it, trackpeps] }
    | combine(db_c)
    | generateTrackPeptideLibrary
    | set { full_lib }
  } else {
    lib_c
    | combine(db_c)
    | set { full_lib }
  }

  diann_in
  | combine(full_lib)
  | map { [it, ms1acc, ms2acc].flatten() }
  | DiaNN
  DiaNN.out.tsv
  | combine(scandb)
  | prepareNumbers

  emit:
  DiaNN.out.pq
  | combine(prepareNumbers.out)
  | combine(scandb)
}
