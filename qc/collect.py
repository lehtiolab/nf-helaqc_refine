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
  
  <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAA0oAAAB/CAYAAADRhTJuAAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAABDXAAAQ1wBgBiQDQAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAACAASURBVHic7d15fFx1uT/wz3POTCbdW8pWlluhpWxFAQFZReQKiKKoEJvJUipeqqKgFyyZBH6OSpOURQQULVJKM5NJmCL3CrhckUURkJ0KVloKyL6WrrSZzJzv8/sjKdTQdr5nMmdmmnzeLytt83zP90kynZznnO95vgKyMrV9wc8EOLeQsWKc/Zc3z3qm2DkNJ/u0L3gZwB5+xymwfEXT2fsGkBIRERERDWFOuRMgIiIiIiKqNCyUiIiIiIiIBmChRERERERENAALJSIiIiIiogFYKBEREREREQ3AQomIiIiIiGgAFkpEREREREQDsFAiIiIiIiIagIUSERERERHRACyUiIiIiIiIBgiVOwGioWjeM2+PyfZU7Z830Jg3Lz50woslSImIiIiIfGChRBQAb0P1oY6Ye/PFqTjXATg3+IyIiIiIyA8uvSMiIiIiIhqAhRIREREREdEALJSIiIiIiIgGYKFEREREREQ0AAslIiIiIiKiAVgoERERERERDcBCiYiIiIiIaADuo0SBiD/62sgqd9R1+eIEuix2yPi2UuRERERERGSLhRIFIlI1MqIeZuaLM5C/AmChREREREQVhUvviIiIiIiIBqiYO0pzn1xT7ygOyxfnSfjyiw8e+WopciIqh9Yn1pwqwEn54ozj3NDysTFPlyInIiIiouGmYgolUZyiQF3+uMxNAFgo0RAmxyj0/HxRrtF7ALBQIiIiIgoAl94RERERERENwEKJiIiIiIhoABZKREREREREA7BQIiIiIiIiGoCFEhERERER0QAslIiIiIiIiAZgoURERERERDQACyUiIiIiIqIBWCgRERERERENECp3AkRERENNZ2fTx23iPC/0bmPjpS8EnQ/9u/Silt2N6IR8cTlHTX1929JS5ES0PensbNobQN5/Q65xdEZD6+MlSCkQLJSIiIiKTIw8DItVGyHxugHUBp8Rbc5zvUtV5ax8cY4iA6A6yFwSiaapLuRh6wGCn0br235kE5rqiF0PwRm2h87kIrvNmhXvsc6Fhi0xMg/I/9oy0F4AkeAz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%3D" width="500px"><h2 class="title is-2">QC for {{ searchname }}</h2>
  <hr>
  <h3 class="title is-3">Protein/peptide level QC</h3>
{% for graphtype in ["featyield", "precursorarea", "isobaric", "nrpsms", "nrpsmsoverlapping", "percentage_onepsm"] %}
  {% if graphtype in features[features.keys()[0]] %}
  <h4 class="title is-4">{{ titles[graphtype] }}</h4>
  <div class="columns">
    {% for feat in features %}
    <div class="column">
      <h5 class="title is-5">{{ featnames[feat] }}</h5>
      {{ features[feat][graphtype] }}
    </div>
    {% endfor %}
    {% if graphtype == "isobaric" %}
    <div class="column">
      <h5 class="title is-5">Median centering</h5>
        {{ norm }}
    </div>
    {% endif %}
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
  {% if hirief == 'hirief' %}
  <div class="columns">
    {% for graph in psms %}
    <div class="column">
      <h5 class="title is-5">{{ titles[graph] }}</h3>
      {{ psms[graph] }}
    </div>
    {% endfor %}
  </div>
  {% endif %}
</div>
{% for graph in ppsms[firstplate] %}
<div class="container">
  <h4 class="title is-4">{{ titles[graph] }}</h4>
{% for plate, graphs in ppsms|dictsort %}
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
# venn diagrams
ppsms = {}
searchname = sys.argv[1]
frac = sys.argv[2]
if frac == 'hirief':
    fryield = 'Fraction yield'
    for ppsm in sys.argv[3:]:
        plate = os.path.basename(ppsm).replace('_psms.html', '')
        with open(ppsm) as fp:
            ppsms[plate] = {x.attrib['id']: tostring(x) for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'}
    
    with open('psms.html') as fp:
        psms = {x.attrib['id']: tostring(x) for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'}
else:
    fryield = 'Yield'
    psms = {}
    with open('psms.html') as fp:
        ppsms['No plate'] = {x.attrib['id']: tostring(x) for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'}
    
titles = {'psm-scans': '# PSMs and scans', 'miscleav': 'Missed cleavages',
          'missing-tmt': 'Isobaric missing values', 'fr-yield': fryield,
          'retentiontime': 'Retention time', 'prec-error': 'Precursor error',
          'featyield': 'Identifications', 'isobaric': 'Isobaric intensities',
          'precursorarea': 'Precursor area intensity',
          'nrpsms': '# PSMs used for isobaric quantitation per identification',
          'nrpsmsoverlapping': '# PSMs used for isobaric quantitation per identification for only complete overlapping set',
          'percentage_onepsm': 'Percentage of identifications with >1 quantifying PSM in the complete overlapping set',
          'coverage': 'Overall protein coverage',
}
featnames = {'assoc': 'Gene symbols', 'peptides': 'Peptides', 'proteins': 'Proteins', 'genes': 'Genes'}

graphs = OrderedDict()
for feat in ['peptides', 'proteins', 'genes', 'assoc']:
    try:
        with open('{}.html'.format(feat)) as fp:
            graphs[feat] = {x.attrib['id']: tostring(x) for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'}
    except IOError as e:
        print(feat, e)
try:
    with open('norm.html') as fp:
        normgraph = [tostring(x) for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'][0]
except IOError:
    print('No normalization file')
    normgraph = False
except AssertionError:
    print('Normalization file not XML, assuming no normalizing done')
    normgraph = False

with open('qc.html', 'w') as fp:
    fp.write(main.render(hirief=frac, searchname=searchname, titles=titles, featnames=featnames, psms=psms, firstplate=sorted(ppsms.keys())[0], ppsms=ppsms, features=graphs, norm=normgraph))
