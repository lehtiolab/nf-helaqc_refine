
process msconvert {

  input:
  tuple path(raw), val(instrument), val(usr_filters), val(usr_options)

  output:
  path(outfile)

  script:
  // the string "infile" does not have NF escaping characters like & (e.g. in FAIMS 35&65),
  // which it does to "raw". That would work fine but not if the files are quoted in the 
  // script, then they cant be found when there is \&.
  infile = "${raw.baseName}.${raw.extension}"
  outfile = "${raw.baseName}.mzML"

  std_filters = ['"peakPicking true 2"', '"precursorRefine"']
  def additional_filters = [
    timstof: ['"scanSumming precursorTol=0.02 scanTimeTol=10 ionMobilityTol=0.1"'],
    qe: [],
    velos: [],
    astral: [],
  ]
  filters = usr_filters ? usr_filters.tokenize(';') : std_filters + additional_filters[instrument]
  filters = filters.collect() { x -> "--filter ${x}" }.join(' ')

  std_options = [
    timstof: ['combineIonMobilitySpectra'],
    qe: [],
    velos: [],
    astral: [],
  ]
  options = usr_options ? usr_options.tokenize(';') : std_options[instrument]
  options = options.collect() {x -> "--${x}"}.join(' ')
  """
  # Resolve directory if necessary, pwiz tries to read NF soft links as if they are files, which
  # does not work in case of directory
  ${raw.isDirectory() ?  "mv '${infile}' tmpdir && cp -rL tmpdir '${infile}'" : ''}
  wine msconvert "${infile}" ${filters} ${options}
  """
}


process createNewSpectraLookup {
  container 'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'

  input:
  tuple path(mzml), path(dino)
  
  output:
  path('mslookup_db.sqlite')

  script:
  """
  msstitch storespectra --spectra "${mzml}" --setnames 'QC'
  ${dino.baseName != 'NO__FILE' ? "msstitch storequant --dbfile mslookup_db.sqlite --dinosaur \"${dino}\" --spectra \"${mzml}\" --mztol 20.0 --mztoltype ppm --rttol 5.0" : ''}
  """
}


