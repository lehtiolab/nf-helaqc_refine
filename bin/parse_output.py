#!/usr/bin/env python3

import os
import re
import sys
import json
import argparse
from collections import defaultdict
from sqlite3 import Connection
from pyarrow import compute as pc
from pyarrow import parquet as pq
from pyarrow import csv as pcsv


def calc_boxplot_qs(vals):
    filt = pc.is_finite(vals)
    vals = vals.filter(filt)
    if len(vals):
        q1 = pc.quantile(vals, q=0.25)[0].as_py()
        q2 = pc.quantile(vals, q=0.5)[0].as_py()
        q3 = pc.quantile(vals, q=0.75)[0].as_py()
        return {'q1': q1, 'q2': q2, 'q3': q3, }
    else:
        return False


def parse_wc_output(wc_out):
    return int(wc_out[:wc_out.index(' ')])

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--scandb', dest='scandb')
parser.add_argument('--acquisition', dest='acq')
parser.add_argument('--nrpsms', dest='numpsms', type=int)
parser.add_argument('--nrpeps', dest='numpeps', type=int)
parser.add_argument('--nruni', dest='numuni', type=int)
parser.add_argument('--nrprot', dest='numprot', type=int)
parser.add_argument('--peaks_on_lc', dest='peaks_fwhm', type=float)
parser.add_argument('--trackedpeptides', dest='trackpep', nargs='+', default=[])
args = parser.parse_args(sys.argv[1:])


if args.acq == 'dda':
    numpsms = args.numpsms - 1
    numpeps = args.numpeps - 1
    numuni = args.numuni - 1
    numprot = args.numprot - 1
else:
    numpsms = args.numpsms
    numpeps = args.numpeps
    numuni = args.numuni
    numprot = args.numprot


if os.path.isdir(args.scandb):
    scansql = 'SELECT COUNT(*) FROM Frames'
    scandb = os.path.join(args.scandb, 'analysis.tdf')
else:
    # msstitch sqlite file
    scansql = 'SELECT COUNT(*) FROM mzml'
    scandb = args.scandb

#if os.path.exists(args.nrscans_or_db) and os.path.isfile(args.nrscans_or_db):
con = Connection(scandb)
nrscans = con.execute(scansql).fetchone()[0]

qcout = {
    'nrpsms': numpsms,
    'nrscans': nrscans,
    'nrpeptides': numpeps,
    'nr_unique_peptides': numuni,
    'nrproteins': numprot,
    'missed_cleavages': {},
    }

headers = {
    'dia': {
        'p_error': False,
        'fwhm': False,
		'injtime': False,
        'rt': 'RT',
        'ch': 'Precursor.Charge',
        'score': 'Q.Value',
        'ionmob': 'IM',
        'ms1': 'Ms1.Area',
        'seq': 'Stripped.Sequence',
        'fwhm': 'FWHM',
        },
    
    'dda': {
        'p_error': 'precursor_ppm',
        'fwhm': 'FWHM',
        'seq': 'bareseq',
        'ch': 'charge',
        'rt': 'rt',
        'score': 'sage_discriminant_score',
        'ionmob': 'ion_mobility',
        'ms1': 'MS1 area (highest of all PSMs)',
        'injtime': 'Ion injection time(ms)',
        'matchedpeaks': 'matched_peaks',
        },
    }[args.acq]

# FIXME Q.Value is not a good score thing? in DIA

peps = pcsv.read_csv('peptable.txt', parse_options=pcsv.ParseOptions(delimiter='\t'))

if args.acq == 'dia':
    precursors = pq.read_table('tpsms')
    qcout['peaks_fwhm'] = args.peaks_fwhm
elif args.acq == 'dda':
    precursors = pcsv.read_csv('tpsms', parse_options=pcsv.ParseOptions(delimiter='\t'))
    precursors = precursors.add_column(0, 'bareseq', pc.replace_substring_regex(precursors['peptide'], '[^A-Z]', ''))
    qcout['injtime'] = calc_boxplot_qs(precursors[headers['injtime']])
    qcout['matchedpeaks'] = calc_boxplot_qs(precursors[headers['matchedpeaks']])
    peps = peps.add_column(0, 'bareseq', pc.replace_substring_regex(peps['Peptide sequence'], '[^A-Z]', ''))

miscl = pc.value_counts(pc.count_substring_regex(precursors[headers['seq']], '[KR][^P]'))
qcout['missed_cleavages'] = {x['values'].as_py(): x['counts'].as_py() for x in miscl}
qcout['score'] = calc_boxplot_qs(precursors[headers['score']])
qcout['rt'] = calc_boxplot_qs(precursors[headers['rt']])
if headers['p_error']:
    qcout['p_error'] = calc_boxplot_qs(precursors[headers['p_error']])

qcout['fwhm'] = calc_boxplot_qs(precursors[headers['fwhm']])
ionmob = calc_boxplot_qs(precursors[headers['ionmob']])
if ionmob and ionmob['q2'] != 0.0:
    qcout['ionmob'] = ionmob

# Peptide MS1
qcout['peptide_areas'] = calc_boxplot_qs(peps[headers['ms1']])

qcout['trackedpeptides'] = defaultdict(dict)
precursortrackfields = {headers[x]: x for x in  ['rt', 'score']}
for pep_ch in args.trackpep:
    pep, ch = pep_ch.split('_')
    seq_filter_exp = pc.field(headers['seq']) == pep
    ch_filter_exp = pc.field(headers['ch']) == int(ch)
    for field, val in precursors.filter(seq_filter_exp & ch_filter_exp).to_pydict().items():
        if field in precursortrackfields:
            # TODO median val instead of first, in case of DDA
            qcout['trackedpeptides'][pep_ch][precursortrackfields[field]] = val[0]
    qcout['trackedpeptides'][pep_ch]['ms1'] = peps.filter(seq_filter_exp).to_pydict()[headers['ms1']][0]



with open('qc.json', 'w') as fp:
    json.dump(qcout, fp, indent=2)
