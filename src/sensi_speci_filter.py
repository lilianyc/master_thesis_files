
"""
Filter GWAS file with sensitivity and specificity of presence and absence of
genes.

Account for number given as percentages ?

Needs python 3.6+ and pandas
"""

__authors__ = "Lilian Yang-Crosson"
__version__ = "0.1.0"
__maintainer__ = "Lilian Yang-Crosson"


import argparse
from pathlib import Path
import sys

import pandas as pd


DEBUG = False
# Default column names for sensitivity and specificity of presence/absence.
SENSITIVITY_PRESENCE = "sensi_1"
SENSITIVITY_ABSENCE = "sensi_0"

SPECIFICITY_PRESENCE = "speci_1"
SPECIFICITY_ABSENCE = "speci_0"

# Default thresholds, change accordingly.
SE1_THRESHOLD = 0.083
SP1_THRESHOLD = 0.983

SE0_THRESHOLD = 0.083
SP0_THRESHOLD = 0.983
# Add an auto-evaluation for threshold based on number of samples ?

P_VALUE_COLUMN = "lrt-pvalue"

def user_input():
    """Have args ?
    """
    parser = argparse.ArgumentParser(
        description=("Read tab-separated GWAS with sensitivity and "
                     "specificity of mutations and filter out "
                     "according to defined thresholds"
                     ))

    # Casting to path here seems to prevent the bash CLI to search for files ?
    # Make it a positional argument ?
    parser.add_argument("-g", "--gwas_file", required=True,
                        metavar="SENSI_SPECI_GWAS_FILE",
                        type=Path, help="GWAS file to parse (REQUIRED)")

    # Default column names.  As required is False by default, might remove it.
    colnames = parser.add_argument_group(
                    title="Column names for sensitivity and specificity")
    colnames.add_argument("--sensi-presence-col", required=False, type=str,
                          default=SENSITIVITY_PRESENCE,
                          help="column name of sensitivity of presence")
    colnames.add_argument("--speci-presence-col", required=False, type=str,
                          default=SPECIFICITY_PRESENCE,
                          help="column name of specificity of presence")

    colnames.add_argument("--sensi-absence-col", required=False, type=str,
                          default=SENSITIVITY_ABSENCE,
                          help="column name of sensitivity of absence")
    colnames.add_argument("--speci-absence-col", required=False, type=str,
                          default=SPECIFICITY_ABSENCE,
                          help="column name of specificity of absence")

    # Option for percentage or not (bool) ? + for thresh
    # Thresholds of specificity/sensitivity.
    thresholds =  parser.add_argument_group(
                    title="Thresholds for sensitivity and specificity")
    thresholds.add_argument("--sensi-presence-thresh", type=float,
                            default=SE1_THRESHOLD)
    thresholds.add_argument("--speci-presence-thresh", type=float,
                            default=SP1_THRESHOLD)

    thresholds.add_argument("--sensi-absence-thresh", type=float,
                            default=SE0_THRESHOLD)
    thresholds.add_argument("--speci-absence-thresh", type=float,
                            default=SP0_THRESHOLD)

    # Filtering on p-value threshold.
    p_values = parser.add_argument_group(
                    title="Arguments for filtering on p-value")
    p_values.add_argument("-p", "--p-value-thresh", type=float, default=1)
    p_values.add_argument("--p-value-col", type=str, default=P_VALUE_COLUMN)

    # Optional arguments.
    # Allow modification of the 2 outputs ?
    # Or at least give an output dir (see how it works with relative paths)
    # Verify thresholds in correct range.
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {}".format(__version__))

    options = parser.parse_args()

    # Control input + logger ?.
    if not options.gwas_file.is_file():
        print(f"GWAS file {options.gwas_file} not found")
        sys.exit()
    
    # List of thresholds (need to keep order).
    thresholds = [options.sensi_presence_thresh,
                  options.speci_presence_thresh,
                  options.sensi_absence_thresh,
                  options.speci_absence_thresh]
    # Print thresholds
    print(f"Thresholds defined for SE1, SP1, SE0, SP0 : {thresholds}")
    if not all(0 <= threshold <= 1 for threshold in thresholds):
        print("Error: incorrect threshold, must be in [0,1]")
        sys.exit()

    print(f"P-value threshold defined: {options.p_value_thresh}")
    if not 0 <= options.p_value_thresh <= 1:
        print("Error: p-value must be between in [0, 1]")
        sys.exit()

    return options



if __name__ == "__main__":
    options = user_input()

    # print(options)
    df_gwas = pd.read_csv(str(options.gwas_file), sep="\t")

    print(f"Parsed GWAS file with {df_gwas.shape[0]} associated genes")
    # Set of column names
    sensi_speci_columns = {options.sensi_presence_col,
                           options.speci_presence_col,
                           options.sensi_absence_col,
                           options.speci_absence_col}

    # See if columns names are in the GWAS file.
    if sensi_speci_columns - set(df_gwas.columns):
        print(f"Error: Columns {sensi_speci_columns - set(df_gwas.columns)} not found")
        sys.exit()

    # User decided to filter by p-value.
    if options.p_value_thresh != 1:
        if options.p_value_col not in df_gwas.columns:
            print(f"Error: p-value column {df_gwas.columns} not found")
            sys.exit()
        df_gwas = df_gwas[df_gwas[options.p_value_col] < options.p_value_thresh]
        print(f"After filtering out p-values over {options.p_value_thresh}, "
              f"{df_gwas.shape[0]} potentially associated genes remain")


    df_gwas_presence = df_gwas[(df_gwas[options.sensi_presence_col] > options.sensi_presence_thresh) &
                               (df_gwas[options.speci_presence_col] > options.speci_presence_thresh)]
    print(f"After filtering on presence, {df_gwas_presence.shape[0]} "
          "potentially associated genes remain")

    df_gwas_absence = df_gwas[(df_gwas[options.sensi_absence_col] > options.sensi_absence_thresh) &
                              (df_gwas[options.speci_absence_col] > options.speci_absence_thresh)]
    print(f"After filtering on absence, {df_gwas_absence.shape[0]} "
          "potentially associated genes remain")


    presence_file = "presence_filtered_" + options.gwas_file.stem + ".txt"
    absence_file = "absence_filtered_" + options.gwas_file.stem + ".txt"

    print(f"Writing GWAS files : {presence_file}, {absence_file}")

    df_gwas_presence.to_csv(presence_file, sep="\t", index=False)
    df_gwas_absence.to_csv(absence_file, sep="\t", index=False)


    if DEBUG:

        se1_threshold = 0.083
        sp1_threshold = 0.983
        
        se0_threshold = 0.083
        sp0_threshold = 0.983

        gwas_dir = Path("//sas-vp-lsgw1/Homes/l.yang-crosson")
        gwas_file = gwas_dir.joinpath("2020-06-04_sensi_speci.txt")
        gwas_file = gwas_dir.joinpath("human_sensi_speci.txt")
    
    
        gwas_file.is_file()
        # Obsolete
        output_file = "filtered_" + gwas_file.stem + ".txt"


        df_gwas = pd.read_csv(str(gwas_file), sep="\t")

        # Verify if columns exist
    
        # df_gwas[SENSITIVITY_PRESENCE] > se1_threshold
        # df_gwas[SPECIFICITY_PRESENCE] > sp1_threshold
    
        df_gwas_presence = df_gwas[(df_gwas[SENSITIVITY_PRESENCE] > se1_threshold) &
                                   (df_gwas[SPECIFICITY_PRESENCE] > sp1_threshold)]
    
        print(f"After filtering on presence, the shape is {df_gwas_presence.shape}")
        # 1st output
    
        df_gwas_absence = df_gwas[(df_gwas[SENSITIVITY_ABSENCE] > se0_threshold) &
                                   (df_gwas[SPECIFICITY_ABSENCE] > sp0_threshold)]
    
        print(f"After filtering on absence, the shape is {df_gwas_absence.shape}")

        print(f"Writing gwas_file : {output_file}")
        #df_gwas_absence.to_csv(output_file, sep="\t", index=False)
