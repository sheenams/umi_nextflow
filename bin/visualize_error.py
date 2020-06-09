#!/usr/bin/env python3
"""
Script to visualize non-reference reads in a set of files created by VarScan readcounts
"""
import os
import argparse 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

BASES = ['A', 'C', 'G', 'T']
COLORBAR_LABELS = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0]
COLORBAR_TICKS = np.log10(COLORBAR_LABELS)

def build_parser(parser):
    parser.add_argument('readcounts', nargs='+',
        help='Output from mpileup2readcounts')
    parser.add_argument('-o', '--outpath',
        help='Path prefix for output files')
    parser.add_argument('-l', '--labels', nargs='+',
        help='Labels for readcount files (default: infer from filename)')
    parser.add_argument('-s', '--sample', 
        help='Sample name (default_infer from filename)')
    parser.add_argument('-c', '--min_coverage', type=int, default=1000,
        help='Minimum coverage to consider a site (default: %(default)s)')

# older version that created dataframe from output of VarScan readcounts
def VARSCAN_readcounts_to_dataframe(readcounts):
    with open(readcounts, 'r') as f:
        series_list = []
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


    df = pd.concat(series_list, axis=1).T
    df.fillna(0, inplace=True)

    for base in BASES:
        df[base + '_freq'] = df[base] / df['depth']

    df.set_index(['chrom', 'pos'], inplace=True)
    return df

# current version uses output from mpileup2readcounts binary
def readcounts_to_dataframe(readcounts, min_coverage):
    col_names = ['chrom', 'pos', 'ref', 'depth', 'A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'ins', 'del', 'empty']
    df = pd.read_csv(readcounts, sep='\t', header=0, names=col_names, dtype={'chrom' : str})

    df = df[df['depth'] >= min_coverage]

    for base in BASES:
        df[base] = df[base] + df[base.lower()]
        df[base + '_freq'] = df[base] / df['depth']

    df.set_index(['chrom', 'pos'], inplace=True)
    return df

def dataframe_to_vafs(df):
    # create array of every non-zero VAF
    vaf_list = []
    for ref in BASES:
        ref_df = df[df['ref'] == ref]
        for var in [x for x in BASES if x != ref]:
            var_df = ref_df[ref_df[var + '_freq'] > 0]
            vaf_list.append(var_df[var + '_freq'].ravel())
    # concatenate the arrays
    vafs = np.concatenate(vaf_list)
    return vafs

def create_vaf_histogram(vafs, bins, axis, title=None):   
    # generate histogram
    axis.set_xlabel('non-reference read franction')
    if len(vafs) > 0:
        axis.hist(vafs, bins=bins, log=True)
    if title:
        axis.set_title(title)

def create_error_heatmap(df, axis, title=None):
    # calculate error rates
    error_rates = pd.DataFrame(columns=BASES, index=pd.Series(BASES, name='ref'))
    total_ref_reads = 0
    total_var_reads = 0
    for ref in BASES:
        ref_df = df[df['ref'] == ref]
        ref_reads = ref_df['depth'].sum()
        total_ref_reads += ref_reads
        for var in BASES:
            var_reads = ref_df[var].sum()
            if var != ref:
                total_var_reads += var_reads
            try:
                error_rates.loc[ref, var] =  var_reads / ref_reads
            except ZeroDivisionError:
                continue

    try:
        error_rate = total_var_reads / total_ref_reads
    except ZeroDivisionError:
        error_rate = np.nan

    # plot heatmap of log-adjusted errors for better contrast
    log_matrix = np.log10(error_rates.to_numpy(dtype=np.float64))
    im = axis.imshow(log_matrix, cmap='cool', vmin=COLORBAR_TICKS[0], vmax=COLORBAR_TICKS[-1])

    # add axis ticks and labels
    ref_bases = error_rates.index
    var_bases = error_rates.columns
    axis.set_xlabel('ALT Base')
    axis.set_xticks(range(len(var_bases)))
    axis.set_xticklabels(var_bases)
    axis.set_xlim(-0.5, len(var_bases) - 0.5)
    axis.set_ylabel('REF Base')
    axis.set_yticks(range(len(ref_bases)))
    axis.set_yticklabels(ref_bases)
    axis.set_ylim(len(ref_bases) - 0.5, -0.5)
    if title:
        axis.set_title(title + ' (error: {0:.4f})'.format(error_rate))
    else:
        axis.set_title('average error: {0:.4f}'.format(error_rate))

    # add text to each grid square with the non-log error rate
    for i in range(len(ref_bases)):
        for j in range(len(var_bases)):
            value = round(error_rates.iloc[i, j], 5)
            text = axis.text(j, i, value, ha="center", va="center", color="black")

    return im

if __name__ == '__main__':
    # parse the args
    parser = argparse.ArgumentParser()
    build_parser(parser)
    args = parser.parse_args()

    num_files = len(args.readcounts)

    if args.labels and len(args.labels == num_files):
        labels = args.labels
    else:
        file_names = [os.path.basename(x) for x in args.readcounts]
        labels = [x.split('.')[1] for x in file_names]

    if args.sample:
        sample = args.sample
    else:
        sample = os.path.basename(args.readcounts[0]).split('.')[0]

    data = []
    # parse the readcounts files
    for rc in args.readcounts:
        data.append(readcounts_to_dataframe(rc, args.min_coverage)) 

    # create multi-histogram
    hist_fig, hist_axs = plt.subplots(2, num_files, sharey='row', sharex='row', figsize=(4*num_files, 8), constrained_layout=True)
    
    for i in range(num_files):
        vafs = dataframe_to_vafs(data[i])
        # create complete histogram
        create_vaf_histogram(vafs, np.arange(0, 1.005, 0.005), hist_axs[0, i], title=labels[i])
        # create low-AF histogram
        create_vaf_histogram(vafs, np.arange(0, 0.051, 0.001), hist_axs[1, i])
        
    # save histogram(s) 
    hist_fig.suptitle(sample)   
    hist_fig.savefig(args.outpath + '.VAFs.png')

    # create error rate heatmap
    heat_fig, heat_axs = plt.subplots(1, num_files, figsize=(num_files*4 + 2, 5))

    for df, l, heat_ax in zip(data, labels, heat_axs):
        im = create_error_heatmap(df, heat_ax, title=l)

    heat_fig.suptitle(sample)
    # add color scale
    heat_fig.subplots_adjust(right=0.8)
    cbar_ax = heat_fig.add_axes([0.85, 0.2, 0.05, 0.6])
    heat_fig.colorbar(im, cax=cbar_ax, ticks=COLORBAR_TICKS)
    cbar_ax.set_yticklabels(COLORBAR_LABELS)
    # save heatmap(s)
    heat_fig.savefig(args.outpath + '.error_rates.png')