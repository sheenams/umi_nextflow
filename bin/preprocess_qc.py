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

picard_fields = ["TOTAL_READS", "MEAN_TARGET_COVERAGE", 
                 "PCT_USABLE_BASES_ON_TARGET", "PCT_OFF_BAIT", "PCT_TARGET_BASES_100X",
                 "PCT_SELECTED_BASES"]

picard_format = {
    'TOTAL_READS': "{:.0f}",
    'MEAN_TARGET_COVERAGE': "{:.0f}",
    'PCT_USABLE_BASES_ON_TARGET': "{:.1%}",
    'PCT_TARGET_BASES_100X': "{:.1%}",
    'PCT_SELECTED_BASES': "{:.1%}",
    'PCT_OFF_BAIT': "{:.1%}",
}

HEADER = """# plot_type: 'table'
# section_name: 'UMI Summary'
"""

def parse_picard(args):
    j = json.load(args.input)
    df = pd.DataFrame(j["report_saved_raw_data"]["multiqc_picard_HsMetrics"]).T
    df["sample_names"] = df.index.map(lambda i: i.split(".")[0])
    df["bam_type"] = df.index.map(lambda i: i.split(".")[-1])
    df = df.sort_values(["sample_names", "bam_type"], ascending=[True, False])
    
    for key, formatter in picard_format.items():
        df[key] = df[key].apply(formatter.format)

    args.output.write(HEADER + "\n")
    cols = ["bam_type"] + picard_fields
    df[cols].to_csv(args.output, sep="\t")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers()

    parser_picard = subparsers.add_parser('picard')
    parser_picard.add_argument('input', type=argparse.FileType(), help="MultiQC JSON file to parse")
    parser_picard.add_argument('output', type=argparse.FileType('w'), help="Output CSV file")
    parser_picard.set_defaults(func=parse_picard)

    args = parser.parse_args()


    try:
        args.func(args)
    except AttributeError:
        parser.print_help(sys.stderr)
        sys.exit(1)        

