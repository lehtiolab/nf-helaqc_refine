#!/usr/bin/env python3

import os
import re
import sys
import json
import argparse
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


def count_missed_cleavage(pepseq, count=0):
    # From msstitch
    '''Regex .*[KR][^P] matches until the end and checks if there is a final
    charachter so this will not match the tryptic residue'''
    #pepseq = re.sub('[\+\-]\d*.\d*', '', full_pepseq)
    match = re.match('.*[KR][^P]', pepseq.upper())
    if match:
        count += 1
        return count_missed_cleavage(match.group()[:-1], count)
    else:
        return count


def parse_wc_output(wc_out):
    return int(wc_out[:wc_out.index(' ')])

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--scandb', dest='scandb')
parser.add_argument('--nrpsms', dest='numpsms', type=int)
parser.add_argument('--nrpeps', dest='numpeps', type=int)
parser.add_argument('--nruni', dest='numuni', type=int)
parser.add_argument('--nrprot', dest='numprot', type=int)
parser.add_argument('--acquisition', dest='acq')
args = parser.parse_args(sys.argv[1:])


# FIXME -1 is only for wc w header! put behind an if statement
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
        'miscl': False,
        'rt': 'RT',
        'score': 'CScore',
        'ionmob': 'IM',
        'ms1': 'Ms1.Area',
        'seq': 'Stripped.Sequence',
        },
    
    'dda': {
        'p_error': 'precursor_ppm',
        'fwhm': 'FWHM',
        'rt': 'rt',
        'score': 'sage_discriminant_score',
        'miscl': 'missed_cleavages',
        'ionmob': 'ion_mobility',
        'ms1': 'MS1 area (highest of all PSMs)',
        },
    }[args.acq]

with open('tpsms') as fp:
    header = next(fp).strip('\n').split('\t')
    calc_ms1data = True
    try:
        # FIXME
        fwhmix = header.index(headers['fwhm'])
    except ValueError:
        fwhmix = False
        if args.acq == 'dda':
            print('No FWHM in PSM table, probably --noms1 is specified, e.g. for TIMS DDA')
            calc_ms1data = False

    use_ionmob = False
    qcpsms = []
    if headers['miscl']:
        misclix = header.index(headers['miscl'])
    else:
        misclix = False
        seqix = header.index(headers['seq'])
    try:
        ionmobix = header.index(headers['ionmob'])
    except ValueError:
        ionmobix = False
    for line in fp:
        line = line.strip('\n').split('\t')
        if ionmobix and line[ionmobix] != 'NA':
            use_ionmob = True
        qcpsms.append(line)
        if not misclix:
            mcnum = count_missed_cleavage(line[seqix])
        else:
            mcnum = int(line[misclix])
        if mcnum < 4:
            try:
                qcout['missed_cleavages'][mcnum] += 1
            except KeyError:
                qcout['missed_cleavages'][mcnum] = 1

    score_ix = header.index(headers['score'])
    qcout['scores'] = calc_boxplot_qs([psm[score_ix] for psm in qcpsms])

    rtix = header.index(headers['rt'])
    qcout['retention_times'] = calc_boxplot_qs([psm[rtix] for psm in qcpsms])

    if fwhmix:
        qcout['fwhms'] = calc_boxplot_qs([psm[fwhmix] for psm in qcpsms])
        
    if headers['p_error']:
        perrorix = header.index(headers['p_error'])
        qcout['precursor_errors'] = calc_boxplot_qs([psm[perrorix] for psm in qcpsms])

    if use_ionmob:
        qcout['ionmobilities'] = calc_boxplot_qs([psm[ionmobix] for psm in qcpsms])


if calc_ms1data:
    # We dont have MS1 for TIMS DDA
    with open('peptable.txt') as fp:
        header = next(fp).strip('\n').split('\t')
        areaix = header.index(headers['ms1'])
        qcout['peptide_areas'] = calc_boxplot_qs([x.strip('\n').split('\t')[areaix] for x in fp])


with open('qc.json', 'w') as fp:
    json.dump(qcout, fp, indent=2)
