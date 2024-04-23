/*
==============================
MZML REFINER PIPELINE
==============================
For when you have discovered you have systematic errors (i.e. global m/z shifts) in your
spectra.

@Authors
Jorrit Boekel @glormph
*/

nextflow.enable.dsl = 2


/* input is a list with filenames to run:
/path/to/file1.mzML
/path/to/file2.mzML
...
*/

process msgfPlus {
  container 'quay.io/biocontainers/msgf_plus:2023.01.1202--hdfd78af_0'

  input:
  tuple path(x), val(sample), path('tdb.fa'), path(mods)

  output:
  tuple path(x), val(sample), path("search.mzid")
  
  script:
  instrument = params.instrument ? params.instrument : false
  msgfprotocol = 0 // automatic protocol
  msgfinstrument = [lowres:0, velos:1, qe:3, qehf: 3, 0:0, qehfx:1, lumos:1, timstof:2][instrument]
  """
  msgf_plus -Xmx16g -d tdb.fa -s "$x" -o search.mzid -thread 12 -mod "${mods}" -tda 0 -t 50.0ppm -ti -1,2 -m 0 -inst ${msgfinstrument} -e 1 -protocol ${msgfprotocol} -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
  rm tdb.c*
  """
}

process mzRefine {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  tuple path(mzml), val(sample), path("${sample}.mzid")

  output:
  path(outfile)

  script:
  outfile = "${sample}_refined.mzML"
  """
  wine msconvert $mzml --outfile "${outfile}" --filter "mzRefiner ${sample}.mzid thresholdValue=-1e-10 thresholdStep=10 maxSteps=2 thresholdScore=SpecEValue"
  """
}


workflow {

  modsfn = [
    itraq8plex: "${baseDir}/data/itraq8mods.txt",
    itraq4plex: "${baseDir}/data/itraq4mods.txt",
    tmt10plex: "${baseDir}/data/tmtmods.txt",
    tmt6plex: "${baseDir}/data/tmtmods.txt",
    tmtpro: "${baseDir}/data/tmt16mods.txt",
    tmt16plex: "${baseDir}/data/tmt16mods.txt",
    tmt18plex: "${baseDir}/data/tmt16mods.txt",
    lf: "${baseDir}/data/labelfreemods.txt",
  ][params.isobaric]

  channel.from(file(params.input).readLines())
  | map { it -> it.tokenize('\t') }
  | map { it -> [file(it[0]), file(it[0]).baseName.replaceFirst(/.*\/(\S+)\.mzML/, "\$1")] } // file, samplename, fn_dbid
  | map { it + [file(params.db), file(modsfn)] }
  | msgfPlus
  | mzRefine
}
