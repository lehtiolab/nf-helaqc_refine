#!/usr/bin/env python3

from pyarrow import parquet as pq
from pyarrow import compute as pc
from pyarrow import csv as csv


with open('stats') as fp:
    header = next(fp).strip().split('\t')
    stats = {k: v for k,v in zip(header, next(fp).strip().split('\t'))}



precursors = pq.read_table('precursors')
peptides = precursors.group_by('Modified.Sequence').aggregate([('Ms1.Area', 'max'),
    ('Stripped.Sequence', 'one'), ('Genes', 'one')])
csv.write_csv(peptides, 'peptides', write_options=csv.WriteOptions(delimiter='\t', quoting_style='none', quoting_header='none'))

envstr = f'''
NRPROTS={stats['Proteins.Identified']}
NRPSMS={stats['Precursors.Identified']}
NRSCANS={stats['FWHM.Scans']}
NRPEPS={len(peptides)}
NRUNIPEPS={len(peptides.filter(~pc.match_substring(pc.field("Genes_one"), ';')))}
'''

with open('envvars', 'w') as fp:
    fp.write(envstr)
