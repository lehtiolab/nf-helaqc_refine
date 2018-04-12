from jinja2 import Template
from lxml.html import parse, tostring
import sys
from collections import OrderedDict
import os

main = Template("""<!DOCTYPE html>
<html lang="en">
<head>
    <title>Lehtio proteomics QC report</title>
<link rel="stylesheet" type="text/css" href="https://cdnjs.cloudflare.com/ajax/libs/bulma/0.6.2/css/bulma.min.css">
</head>
<body>
<div class="container">
  <h3 class="title is-3">Protein/peptide level QC</h3>
{% for graphtype in ["featyield", "precursorarea", "isobaric", "nrpsms"] %}
  {% if graphtype in features[features.keys()[0]] %}
  <h4 class="title is-4">{{ titles[graphtype] }}</h4>
  <div class="columns">
    {% for feat in features %}
    <div class="column">
      <h5 class="title is-5">{{ featnames[feat] }}</h3>
      {{ features[feat][graphtype] }}
    </div>
    {% endfor %}
  </div>
<hr>
{% endif %}
{% endfor %}
</div>
<div class="container">
  <h4 class="title is-4">Overall protein coverage</h3>
    {{ features.proteins.coverage }}
</div>
<hr>
<div class="container">
  <h3 class="title is-3">PSM level QC</h3>
  <div class="columns">
    {% for graph in psms %}
    <div class="column">
      <h5 class="title is-5">{{ titles[graph] }}</h3>
      {{ psms[graph] }}
    </div>
    {% endfor %}
  </div>
</div>
{% for graph in ppsms[plates[0]] %}
<div class="container">
  <h4 class="title is-4">{{ titles[graph] }}</h4>
{% for plate in plates %}
<div class="container">
  <h5 class="title is-5">Plate: {{ plate }}</h5>
  {{ ppsms[plate][graph] }}
</div>
{% endfor %}
</div>
{% endfor %}
</body>
</html>
""")


# FIXME
# PSMs
# coverage if protein
# feat yield if multiple sets
# venn diagrams
# isobaric if it is there
titles = {'psm-scans': '# PSMs and scans', 'miscleav': 'Missed cleavages',
          'missing-tmt': 'Isobaric missing values', 'fr-yield': 'Fraction yield', 
          'retentiontime': 'Retention time', 'prec-error': 'Precursor error', 
          'featyield': 'Identifications', 'isobaric': 'Isobaric intensities', 
          'precursorarea': 'Precursor area intensity', 
          'nrpsms': '# PSMs with isobaric quantitation per identification', 
          'coverage': 'Overall protein coverage',
}
featnames = {'assoc': 'Gene symbols', 'peptides': 'Peptides', 'proteins': 'Proteins', 'genes': 'Genes'}
ppsms = {}
for ppsm in sys.argv[1:]:
    plate = os.path.basename(ppsm).replace('_psms.html', '')
    with open(ppsm) as fp:
        ppsms[plate] = {x.attrib['id']: tostring(x) for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'}
        
with open('psms.html') as fp:
    psms = {x.attrib['id']: tostring(x) for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'}

graphs = OrderedDict()
for feat in ['peptides', 'proteins', 'genes', 'assoc']:
    try:
        with open('{}.html'.format(feat)) as fp:
            graphs[feat] = {x.attrib['id']: tostring(x) for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'}
    except IOError as e:
        print(feat, e)

with open('qc.html', 'w') as fp:
    fp.write(main.render(titles=titles, featnames=featnames, psms=psms, plates=ppsms.keys(), ppsms=ppsms, features=graphs))
