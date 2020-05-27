#!/usr/bin/env python3
"""
Script to visualize non-reference reads in a set of files created by VarScan readcounts
"""
import os
import argparse 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import datetime

BASES = ['A', 'C', 'G', 'T']

def build_parser(parser):
    parser.add_argument('readcounts',
        help='Output from VarScan readcounts')
    parser.add_argument('outpath',
        help='Path prefix for output files')
    parser.add_argument('-l', '--label', default='errors',
        help='Label for chart titles')


# older version that created dataframe from output of VarScan readcounts
def VARSCAN_readcounts_to_dataframe(readcounts):
    with open(readcounts, 'r') as f:
        series_list = []
        print(datetime.now(), "reading from readcounts")
        header = next(f)
        for line in f:
            row_fields = line.strip().split('\t')
            s_dict = {}
            s_dict['chrom'] = row_fields[0]
            s_dict['pos'] = int(row_fields[1])
            s_dict['ref'] = row_fields[2]
            s_dict['depth'] = int(row_fields[4])
            
            for entry in row_fields[5:]:
                entry_fields = entry.split(':')
                base = entry_fields[0]
                if base in BASES:
                    s_dict[base] = int(entry_fields[1])

            series_list.append(pd.Series(s_dict))

    print(datetime.now(), "file read")

    df = pd.concat(series_list, axis=1).T
    df.fillna(0, inplace=True)

    for base in BASES:
        df[base + '_freq'] = df[base] / df['depth']

    df.set_index(['chrom', 'pos'], inplace=True)
    print(datetime.now(), "dataframe constructed")
    return df

# current version uses output from mpileup2readcounts binary
def readcounts_to_dataframe(readcounts):
    col_names = ['chrom', 'pos', 'ref', 'depth', 'A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'ins', 'del', 'empty']
    df = pd.read_csv(readcounts, sep='\t', header=0, names=col_names, dtype={'chrom' : str})

    for base in BASES:
        df[base] = df[base] + df[base.lower()]
        df[base + '_freq'] = df[base] / df['depth']

    df.set_index(['chrom', 'pos'], inplace=True)
    print(datetime.now(), "dataframe constructed")
    return df

def create_vaf_histogram(df, title, outpath):
    vaf_list = []
    for ref in BASES:
        ref_df = df[df['ref'] == ref]
        for var in [x for x in BASES if x != ref]:
            var_df = ref_df[ref_df[var + '_freq'] > 0]
            vaf_list.append(var_df[var + '_freq'].ravel())

    vafs = np.concatenate(vaf_list)
    fig, ax = plt.subplots()
    ax.hist(vafs, bins=np.arange(0, 1.001, 0.005), log=True)
    ax.set_title(title)
    ax.set_xlabel('non-reference read franction')
    fig.tight_layout()
    fig.savefig(outpath + '.png')

def create_error_heatmap(df, title, outpath):
    error_rates = pd.DataFrame(columns=BASES, index=pd.Series(BASES, name='ref'))

    for ref in BASES:
        ref_df = df[df['ref'] == ref]
        for var in BASES:
            ref_reads = ref_df['depth'].sum()
            try:
                error_rates.loc[ref, var] = ref_df[var].sum() / ref_reads
            except ZeroDivisionError:
                continue
    
    fig, ax = plt.subplots()
    log_matrix = np.log10(error_rates.to_numpy(dtype=np.float64))
    # np.fill_diagonal(log_matrix, np.nan)

    ref_bases = error_rates.index
    var_bases = error_rates.columns

    ax.imshow(log_matrix, cmap='cool')
    ax.set_title(title)
    ax.set_xlabel('ALT Base')
    ax.set_xticks(range(len(var_bases)))
    ax.set_xticklabels(var_bases)
    ax.set_ylabel('REF Base')
    ax.set_yticks(range(len(ref_bases)))
    ax.set_yticklabels(ref_bases)

    for i in range(len(ref_bases)):
        for j in range(len(var_bases)):
            value = round(log_matrix[i, j], 2)
            text = ax.text(j, i, value, ha="center", va="center", color="black")

    fig.tight_layout()
    fig.savefig(outpath + '.png')
    error_rates.to_csv(outpath + '.tsv', sep='\t')

if __name__ == '__main__':
    # parse the args
    parser = argparse.ArgumentParser()
    build_parser(parser)
    args = parser.parse_args()

    # convert the readcounts file to a dataframe
    df = readcounts_to_dataframe(args.readcounts)
    # df = VARSCAN_readcounts_to_dataframe(args.readcounts)

    # create histogram of non-reference VAFs
    create_vaf_histogram(df, args.label, args.outpath + '.VAFs')
    create_error_heatmap(df, args.label, args.outpath + '.error_rates')

    print(datetime.now(), 'error rates calculated')