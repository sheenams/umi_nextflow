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

picard_format = {
    'TOTAL_READS': ("Total Reads", "{:.0f}"),
    'MEDIAN_TARGET_COVERAGE': ("Median Target Cov.", "{:.0f}"),
    'PCT_ON_BAIT_BASES': ("% on bait", "{:.1%}"),
    'PCT_ON_TARGET_BASES': ("% on target", "{:.1%}"),
    'PCT_TARGET_BASES_100X': ("% bases >100x", "{:.1%}"),
}

HEADER = """# plot_type: 'table'
# section_name: 'UMI Summary'
"""

def parse_picard(args):
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
    print(df)
    # rename multicolumn index to columns ("TOTAL_READS", "standard") -> "TOTAL_READS (standard)"
    # also add some formatting to help distinguish. Note MultiQC will respect HTML in these headers!
    df.columns = [f"{i} <br/><span style='color: blue'>({t})</span>" for i, t in df.columns]

    args.output.write(HEADER + "\n")
    df.to_csv(args.output, sep=",", index=False)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers()

    parser_picard = subparsers.add_parser('picard')
    parser_picard.add_argument('input', type=argparse.FileType(), help="MultiQC JSON file to parse")
    parser_picard.add_argument('output', type=argparse.FileType('w'), nargs='?', default=sys.stdout, help="Output CSV file")
    parser_picard.set_defaults(func=parse_picard)

    args = parser.parse_args()


    #try:
    args.func(args)
    #except AttributeError:
        # parser.print_help(sys.stderr)
        # sys.exit(1)        

