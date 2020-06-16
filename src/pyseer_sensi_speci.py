
"""
Develop a script (PyseerSeSp) to compute sensitivity and specificity of each associated mutations (lrt-pvalue), according to presence (insertion) and absence (deletion) of mutations in the group of genomes of interest (sensitivity) and in the group of compared genomes (specificity), respectively.

In practice, we have to add 4 new columns in the file of associated mutations (gwas.txt) and output a new file (SeSp.txt): 
-	Sensitivity of presence (Se1) = 
Number of mutations present (encoded 1 in R.tab) in the group of genomes of interest (encoded 1 in the phenotype.txt) / Number of genomes of interest (encoded 1 in the phenotype.txt) x 100.
-	Sensitivity of absence (Se0) = Number of mutations absent (encoded 0 in R.tab) in the group of genomes of interest (encoded 1 in the phenotype.txt) / Number of genomes of interest (encoded 1 in the phenotype.txt) x 100.
-	Specificity of presence (Sp1) = Number of mutations absente (encoded 0 in R.tab) in the compared group of genomes (encoded 0 in the phenotype.txt) / Number of genomes in the compared group (encoded 0 in the phenotype.txt) x 100.
-	Specificity of absence (Sp0) = Number of mutations present (encoded 1 in R.tab) in the compared group of genomes (encoded 0 in the phenotype.txt) / Number of genomes in the compared group (encoded 0 in the phenotype.txt) x 100.

Python 3.6 +
"""

__authors__ = "Lilian Yang-Crosson"
__version__ = "0.1.0"
__maintainer__ = "Lilian Yang-Crosson"


import argparse
from pathlib import Path
import sys

import pandas as pd


DEBUG = False


def user_input():
    """Have args ?
    """
    parser = argparse.ArgumentParser(
        description=("Read tab-separated GWAS, phenotype and "
                     "gene_presence_absence files and return sensitivity and "
                     "specificity of mutations"
                     ))
    # Input file is required as a positional argument.
    parser.add_argument("-g", "--gwas_file", required=True,
                        type=Path, help="GWAS file to parse")
    parser.add_argument("-p", "--phenotype_file", required=True,
                        type=Path, help="phenotype file to parse")
    parser.add_argument("-a", "--presence_absence_file", required=True,
                        type=Path, help="gene_presence_absence file to parse")

    parser.add_argument("-c", "--phenotype_col_name", required=True,
                        type=str, help="column name of phenotype")


    # Optional arguments.
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="name of output file")

    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {}".format(__version__))

    options = parser.parse_args()
    
    # TODO: control input + logger ?.
    if not options.gwas_file.is_file():
        print(f"GWAS file {options.gwas_file} not found")
        sys.exit()
    if not options.phenotype_file.is_file():
        print(f"Phenotype file {options.phenotype_file} not found")
        sys.exit()
    if not options.presence_absence_file.is_file():
        print(f"Presence_absence file {options.presence_absence_file} not found")
        sys.exit()
    
    return options


# =============================================================================
# 
# =============================================================================

if DEBUG:
    pyseer_dir = Path("//sas-vp-lsgw1/Homes/l.yang-crosson/data/pyseer")
    # One of the 2 files, as test.
    gwas_file = pyseer_dir.joinpath("pyseer_out/before_2014_gwas.txt")
    phenotype_file = pyseer_dir.joinpath("2020-04-22_dataset_phenotype.txt")
    rtab_file = pyseer_dir.parent.joinpath("roary/result_roary/gene_presence_absence.Rtab")

    print(pyseer_dir.exists())
    phenotype_file.exists()
    rtab_file.exists()

# =============================================================================
# 
# =============================================================================

# try using series filtering ?
# ungrouped_variants = [variant for variant in df_gwas["variant"] if not variant.startswith("group_")]
# len(ungrouped_variants)  # 110 out of 716


# col included
# df_rtab_bfr2014[df_rtab_bfr2014.index.isin(df_gwas["variant"])].sum(axis=1)

# Not it
# df_rtab[df_rtab_bfr2014["Gene"].isin(df_gwas["variant"])].sum(numeric_only=True, axis=1) / 60

# =============================================================================
# 
# =============================================================================

def get_mutation_count(df_rtab, relevant_mutations, mutation=1):
    # return df_rtab[df_rtab.index.isin(relevant_mutations)].sum(numeric_only=True, axis=1)
    df_results = df_rtab[df_rtab.index.isin(relevant_mutations)]
    # return df_rtab[df_rtab.index.isin(relevant_mutations) == mutation].sum(numeric_only=True, axis=1)
    # Do substraction rather ? Have to keep order
    return (df_results == mutation).sum(numeric_only=True, axis=1)

# Unused
# more intelligent way to get phenotype size and specificity ?
def get_binary_classif(df_rtab, relevant_mutations, phenotype_size, specificity=False):
    # Have a genes col arg ? and phenotype ?
    # df_rtab[df_rtab[genes_col].isin(relevant_mutations)].sum(numeric_only=True, axis=1) / phenotype
    if specificity:
        return (phenotype_size - df_rtab[df_rtab.index.isin(relevant_mutations)].sum(numeric_only=True, axis=1)) / phenotype_size
    return df_rtab[df_rtab.index.isin(relevant_mutations)].sum(numeric_only=True, axis=1) / phenotype_size
    


if __name__ == "__main__":
    
    # Make a CLI func ?
    options = user_input()
    # Verif before + cast ?
    df_gwas = pd.read_csv(str(options.gwas_file), sep="\t")
    # see if pertinent to keep index col or not.
    df_pheno = pd.read_csv(str(options.phenotype_file), sep="\t",  index_col=0)
    df_rtab = pd.read_csv(str(options.presence_absence_file), sep="\t", index_col=0)

    # verif on col name. Print the phenotype on which it is done ?
    phenotype_col = options.phenotype_col_name
    output_file = options.output
    
    if not phenotype_col in df_pheno.columns:
        print(f"Column '{phenotype_col}' not found in phenotype file\n"
              f"Columns found: {list(df_pheno.columns)}")
        sys.exit()

    interest_phenotype = df_pheno[df_pheno[phenotype_col] == 1].index
    compared_phenotype = df_pheno[df_pheno[phenotype_col] == 0].index
    
    presence = interest_phenotype
    absence = compared_phenotype


    # Force pval arg, no need.
    if DEBUG:
        df_gwas = pd.read_csv(gwas_file, sep="\t")
        df_pheno = pd.read_csv(phenotype_file, sep="\t",  index_col=0)
        df_rtab = pd.read_csv(rtab_file, sep="\t", index_col=0)
        # Redo pyseer after
        # argument_recu = "before_2014" # nom du phenotype dans le header
        # 
        genomes_bfr_2014 = df_pheno[df_pheno["before_2014"] == 1].index
        genomes_aftr_2014 = df_pheno[df_pheno["before_2014"] == 0].index
        # give better names after. do dataframes on rtab directly ?
        presence = genomes_bfr_2014 
        absence = genomes_aftr_2014
        
        
        # df_rtab.shape[1]
        # substract or conditional ?
        sensi_pres = get_binary_classif(df_rtab[presence], df_gwas["variant"], df_rtab[presence].shape[1])
        sensi_abs = get_binary_classif(df_rtab[absence], df_gwas["variant"], df_rtab[presence].shape[1])
        
        speci_pres =  get_binary_classif(df_rtab[presence], df_gwas["variant"], df_rtab[absence].shape[1], True)
        speci_abs =  get_binary_classif(df_rtab[absence], df_gwas["variant"], df_rtab[absence].shape[1], True)


    # Encapsulate in a function ?
    tmp = df_gwas

    tmp["sensi_1"] = list(get_mutation_count(df_rtab[presence], df_gwas["variant"]) / df_rtab[presence].shape[1])
    tmp["sensi_0"] = list(get_mutation_count(df_rtab[presence], df_gwas["variant"], mutation=0)  / df_rtab[presence].shape[1])
    
    tmp["speci_1"] = list(get_mutation_count(df_rtab[absence], df_gwas["variant"])  / df_rtab[absence].shape[1])
    tmp["speci_0"] = list(get_mutation_count(df_rtab[absence], df_gwas["variant"], mutation=0)  / df_rtab[absence].shape[1])
    
    # Count
    tmp["mutation_in_interest"] = list(get_mutation_count(df_rtab[presence], df_gwas["variant"]))
    tmp["mutation_not_in_interest"] = list(get_mutation_count(df_rtab[presence], df_gwas["variant"], mutation=0))
    
    tmp["mutation_in_compared"] = list(get_mutation_count(df_rtab[absence], df_gwas["variant"]))
    tmp["mutation_not_in_compared"] = list(get_mutation_count(df_rtab[absence], df_gwas["variant"], mutation=0))

    tmp.to_csv(output_file, sep="\t", index=False)
    print(f"Wrote {output_file}")
