include { createNewSpectraLookup; getScanNumbers } from '../modules.nf'

/* A library can be made like this:
diann --threads 16 --fasta tdb.fa  --gen-spec-lib --fasta-search --out-lib libfile --var-mod 'UniMod:35,15.994915,M' --var-mod 'UniMod:4,57.021464,C' --var-mods 2 --predictor

Instead of --predictor, one can use raw files:
--dir rawfile_dir/ 
*/

process extractThermoScans {
  /* Used in DIA for Thermo, so we dont need to convert to mzML
  just to get the nr of scans */
  container 'ghcr.io/lehtiolab/kantele_thermoreader:latest'

  input:
  path(raw)

  output:
  env('NRSCANS')

  script:
  outfile = "${raw.baseName}.csv"
  """
  wine /scanheadsman/ScanHeadsman.exe $raw 
NRSCANS=\$(tail -n+2 $outfile | wc -l)
  """
}

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
  tuple path('out.parquet'), path('out.stats.tsv')

  script:
  """
  diann-linux --threads 16 \
    --f \$(realpath $raw) \
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
  container "ghcr.io/lehtiolab/nfhelaqc:${workflow.manifest.version}"

  input:
  tuple path('precursors'), path('stats')

  output:
  tuple path('precursors'), path('peptides'), env('NRPROTS'), env('NRPSMS'), env('NRPEPS'), env('NRUNIPEPS'), env('NRSCANS')

  script:
  """
  dia_prep_numbers.py
  sed -i '1s/Stripped.Sequence_one/Stripped.Sequence/' peptides
  source envvars
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
    | extractThermoScans
    | set { thermonrscans }
    
    raw_c.thermo
    | concat(raw_c.bruker)
    | set { diann_in }

    raw_c.bruker
    | getScanNumbers
    | concat(thermonrscans)
    | set { nrscans }

  } else if (mzml) {
    channel.fromPath(mzml)
    | set { mzml_c }

    mzml_c
    | combine(channel.fromPath('NO__FILE'))
    | createNewSpectraLookup
    | getScanNumbers
    | set { nrscans }

    mzml_c
    | set { diann_in }

  } else {
    exit 1, 'Must either input a --raw file.raw, a --raw file.d, or an --mzml file.mzML'
  }

  ms1acc = [timstof: 20, orbitrap: 10, velos: 10, qe: 10][instrument]
  ms2acc = 20

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
  | prepareNumbers
  | combine(nrscans)
  | set { output }

  emit:
  output
}
