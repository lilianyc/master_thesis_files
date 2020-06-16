
"""
Filter GWAS file with sensitivity and specificity of presence and absence of
genes.

Account for percentage or not ?
Needs python 3.6+ and pandas
"""

from pathlib import Path

import pandas as pd


# Column names.
SENSITIVITY_PRESENCE = "sensi_1"
SENSITIVITY_ABSENCE = "sensi_0"

SPECIFICITY_PRESENCE = "speci_1"
SPECIFICITY_ABSENCE = "speci_0"

# Percentage
# se1_threshold = 8.3
# sp1_threshold = 98.3

# se0_threshold = 8.3
# sp0_threshold = 98.3

se1_threshold = 0.083
sp1_threshold = 0.983

se0_threshold = 0.083
sp0_threshold = 0.983



if __name__ == "__main__":
    gwas_dir = Path("//sas-vp-lsgw1/Homes/l.yang-crosson")
    gwas_file = gwas_dir.joinpath("2020-06-04_sensi_speci.txt")


    gwas_file.is_file()
    output_file = "filtered_" + gwas_file.stem + ".txt"


    df_gwas = pd.read_csv(str(gwas_file), sep="\t")
    print(f"Parsed GWAS file of shape {df_gwas}")

    df_gwas[SENSITIVITY_PRESENCE] > se1_threshold
    df_gwas[SPECIFICITY_PRESENCE] > sp1_threshold

    df_gwas_filtered = df_gwas[(df_gwas[SENSITIVITY_PRESENCE] > se1_threshold) &
                               (df_gwas[SPECIFICITY_PRESENCE] > sp1_threshold)]

    print(f"After filtering on presence, the shape is {df_gwas_filtered.shape}")
    # 1st output

    # !!!!
    # Output 2 files on presence and absence
    df_gwas_filtered = df_gwas[(df_gwas[SENSITIVITY_ABSENCE] > se0_threshold) &
                               (df_gwas[SPECIFICITY_ABSENCE] > sp0_threshold)]

    print(f"After filtering on absence, the shape is {df_gwas_filtered.shape}")

    print(f"Writing gwas_file : {output_file}")
    df_gwas_filtered.to_csv(output_file, sep="\t", index=False)
