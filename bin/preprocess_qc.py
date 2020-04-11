#!/usr/bin/env python3
"""
Parse MultiQC output and create UMI-specific CSV report.
Nik Krumm, 2020.
"""
import sys
import glob
import json
import argparse 
import pandas as pd
from collections import defaultdict

picard_format = {
    'TOTAL_READS': ("Total Reads", "{:.0f}"),
    'MEDIAN_TARGET_COVERAGE': ("Median Target Cov.", "{:.0f}"),
    'PCT_ON_BAIT_BASES': ("% on bait", "{:.1%}"),
    'PCT_ON_TARGET_BASES': ("% on target", "{:.1%}"),
    'PCT_TARGET_BASES_100X': ("% bases >100x", "{:.1%}"),
}

def create_header(name):
    return f"""# plot_type: 'table'
# section_name: '{name}'

"""

def write_summary(args):
    j = json.load(args.input)
    df = pd.DataFrame(j["report_saved_raw_data"]["multiqc_picard_HsMetrics"]).T
    df["sample_names"] = df.index.map(lambda i: i.split(".")[0])
    df["bam_type"] = df.index.map(lambda i: i.split(".")[-1])
    df = df.sort_values(["sample_names", "bam_type"], ascending=[True, True])
    
    # custom column generation
    df["PCT_ON_BAIT_BASES"] = df["ON_BAIT_BASES"]/df["PF_UQ_BASES_ALIGNED"]
    df["PCT_ON_TARGET_BASES"] = df["ON_TARGET_BASES"]/df["PF_UQ_BASES_ALIGNED"]

    # limit to only output columns
    cols = ["sample_names", "bam_type"] + list(picard_format.keys())
    df = df[cols]

    # rename and reformat columns
    for key, col_def in picard_format.items():
        col_name, col_fmt = col_def
        df[col_name] = df[key].apply(col_fmt.format)
        del df[key]

    # unstack bam_type rows into columns (so that standard/final appear as columns)
    df = df.set_index(['sample_names', 'bam_type']).unstack(level=-1).reset_index()
    df = df.sort_index(level=[0,1], axis=1, ascending=[False, False])

    # set sample name to index
    df = df.set_index("sample_names")
    
    # Add in fgbio stats
    # these are the UMI histogram percentages
    umi_data = j["report_plot_data"]["fgbio-groupreadsbyumi-plot"]["datasets"][1]
    sample_names = [d["name"].split(".")[0] for d in umi_data]
    umi_df = pd.DataFrame(index = sample_names)
    # calculate non-singleton fraction (i.e., 1.0 - fraction_singletons)
    umi_df[("% UMI ≥ 2", "final")] = [1.0 - d["data"][0][1] for d in umi_data]   
    # format as percentage
    umi_df[("% UMI ≥ 2", "final")] = umi_df[("% UMI ≥ 2", "final")].apply("{:.1%}".format)
    # merge umi_df and df
    df = pd.merge(df, umi_df, left_index=True, right_index=True)

    # rename multicolumn index to columns ("TOTAL_READS", "standard") -> "TOTAL_READS (standard)"
    # also add some formatting to help distinguish. Note MultiQC will respect HTML in these headers!
    df.columns = [f"{i} <br/><span style='color: blue'>({t})</span>" for i, t in df.columns]

    args.output.write(create_header('UMI Summary'))
    df.to_csv(args.output, sep=",", index=True)

def write_counts(args):
    """
    Reads flagstat files created by samtools/sambamba and outputs a CSV.
    Input flagstat files MUST be named by {sample_name}.{bam_type}.flagstat.txt
    """
    readcounts = defaultdict(dict)
    with args.output as out:
        for infile in args.input:
            contents = infile.read()
            readcount = contents.split(" ")[0]
            sample_name = infile.name.split(".")[0]
            bam_type = infile.name.split(".")[1]
            readcounts[bam_type][sample_name] = readcount
            
        df = pd.DataFrame(readcounts)
        df.index.name = 'sample_names'
        out.write(create_header('Readcounts by pipeline stage'))
        df.to_csv(out)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers()

    summary_parser = subparsers.add_parser('summary')
    summary_parser.add_argument('input', type=argparse.FileType(), help="MultiQC JSON file to parse")
    summary_parser.add_argument('output', type=argparse.FileType('w'), nargs='?', 
        default=sys.stdout, help="Output CSV file")
    summary_parser.set_defaults(func=write_summary)

    count_parser = subparsers.add_parser('counts', description=write_counts.__doc__)
    count_parser.add_argument('input', type=argparse.FileType(), nargs='+', help="Flagstat file(s)")
    count_parser.add_argument('--output', type=argparse.FileType('w'), nargs='?', 
        default=sys.stdout, help="Output CSV file")
    count_parser.set_defaults(func=write_counts)


    args = parser.parse_args()


    try:
        args.func(args)
    except AttributeError:
        parser.print_help(sys.stderr)
        sys.exit(1)        

