#!/usr/bin/env python3

import sys
import json
from numpy import quantile
from sqlite3 import Connection

def calc_boxplot_qs(vals):
    vals = [float(x) for x in vals if x != 'NA']
    if len(vals):
        q1 = quantile(vals, 0.25)
        q2 = quantile(vals, 0.5)
        q3 = quantile(vals, 0.75)
        iqr = q3 - q1
        return {'q1': q1, 'q2': q2, 'q3': q3, 'upper': q3 + 1.5 * iqr, 'lower': q1 - 1.5 * iqr}
    else:
        return False


def parse_wc_output(wc_out):
    return int(wc_out[:wc_out.index(' ')])


lookup = sys.argv[1]
numpsms = parse_wc_output(sys.argv[2]) - 1
numpeps = parse_wc_output(sys.argv[3]) - 1
numuni = int(sys.argv[4]) - 1
numprot = parse_wc_output(sys.argv[5]) - 1
with Connection(lookup) as con:
    nrscans = con.execute('SELECT COUNT(*) FROM mzml').fetchone()[0]

qcout = {
    'nrpsms': numpsms,
    'nrscans': nrscans,
    'nrpeptides': numpeps,
    'nr_unique_peptides': numuni,
    'nrproteins': numprot,
    'missed_cleavages': {},
    }

with open('tpsms') as fp:
    header = next(fp).strip('\n').split('\t')
    perrorix = header.index('precursor_ppm')
    calc_ms1data = True
    try:
        fwhmix = header.index('FWHM')
    except ValueError:
        print('No FWHM in PSM table, probably --noms1 is specified')
        calc_ms1data = False
    msgfix = header.index('sage_discriminant_score')
    rtix = header.index('rt')
    misclix = header.index('missed_cleavages')
    ionmobix = header.index('ion_mobility')
    use_ionmob = False
    qcpsms = []
    for line in fp:
        line = line.strip('\n').split('\t')
        if line[ionmobix] != 'NA':
            use_ionmob = True
        qcpsms.append(line)
        if int(line[misclix]) < 4:
            try:
                qcout['missed_cleavages'][line[misclix]] += 1
            except KeyError:
                qcout['missed_cleavages'][line[misclix]] = 1
    qcout['precursor_errors'] = calc_boxplot_qs([psm[perrorix] for psm in qcpsms])
    if calc_ms1data:
        qcout['fwhms'] = calc_boxplot_qs([psm[fwhmix] for psm in qcpsms])
    qcout['sagescores'] = calc_boxplot_qs([psm[msgfix] for psm in qcpsms])
    qcout['retention_times'] = calc_boxplot_qs([psm[rtix] for psm in qcpsms])
    if use_ionmob:
        qcout['ionmobilities'] = calc_boxplot_qs([psm[ionmobix] for psm in qcpsms])


if calc_ms1data:
    with open('peptable.txt') as fp:
        header = next(fp).strip('\n').split('\t')
        areaix = header.index('MS1 area (highest of all PSMs)')
        qcout['peptide_areas'] = calc_boxplot_qs([x.strip('\n').split('\t')[areaix] for x in fp])


with open('qc.json', 'w') as fp:
    json.dump(qcout, fp, indent=2)
