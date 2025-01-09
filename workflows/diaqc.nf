include { msconvert; createNewSpectraLookup } from '../modules.nf' 


// AlLibrary can be made like this:
// diann --threads 16 --fasta tdb.fa  --gen-spec-lib --fasta-search --out-lib libfile --dir rawfile_dir/ --var-mod 'UniMod:35,15.994915,M' --var-mod 'UniMod:4,57.021464,C' --var-mods 2


process DiaNN {

  /*
  q-value cutoff seems 0.01 by default
  */
  
  container 'biocontainers/diann:1.8.1_cv2'
  
  input:
  tuple path(raw), path(lib), path(tdb)
  
  output:
  tuple path('out.txt'), path('out.stats.tsv')
  script:
  """
  diann --threads 16 \
    --f $raw \
     --lib $lib \
    --fasta $tdb \
    --window 8 \
    --mass-acc 15 \
    --mass-acc-ms1 15 \
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
  tuple path('precursors'), path('peptides'), eval('cat nrprots'), eval('cat nrpsms'), eval('wc -l < <(tail -n+2 peptides)')

  script:
  """
  cut -f ${Utils.get_field_nr("stats", "Proteins.Identified")} stats | tail -n1 > nrprots
  cut -f ${Utils.get_field_nr("stats", "Precursors.Identified")} stats | tail -n1 > nrpsms
  
  bioawk -v s=${Utils.get_field_nr("precursors", "Modified.Sequence")} -v g=${Utils.get_field_nr("precursors", "Genes")} -v m=${Utils.get_field_nr("precursors", "Ms1.Area")} -t \
    '{print \$s,\$g,\$m }' precursors > pepfields
  head -n1 pepfields > peptides
  tail -n+2 pepfields | sort -k1b,1 -u >> peptides
  """
}

// MS1 area, IM, Precursor quant (ms2) seem quite different between charged peptides, whereas RT is more similar ?


workflow DIAQC {

  take:
  raw
  library
  db 

  main:
  
  raw_c = channel.fromPath(raw)

  if (file(raw).extension == '.raw') {
    raw_c
    | map { [it, false, false] }
    | msconvert
    | createNewSpectraLookup
    | set { scandb }
  } else {
    raw_c
    | set { scandb }
  }

  raw_c 
  | combine(channel.fromPath(library))
  | combine(channel.fromPath(db))
  | DiaNN
  | combine(scandb)
  | prepareNumbers

  emit:
  prepareNumbers.out
  | combine(raw_c)
}
