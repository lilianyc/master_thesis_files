#### Enterobase metadata filtering ####

# Some steps in this script require some manual fine-tuning and you might want
# to ignore some steps for whatever reason so feel free to adjust the code and/or
# execute only what you need.

# Set working directory.
setwd("~/Lilian_stage_M2_2020/")

# Clean the R environment.
rm(list=ls())

# Used to split strings in a column.
library(tidyr)

# R version 3.6.1 (2019-07-05)
# Enterobase v1.1.2 has by default 41 columns when downloading the metadata table.
# When loading the data into R, spaces and other non-alphabetical characters
# of column names will be converted to dots ".", unless check.names option is F:
#
metadata_columns = c("Uberstrain", "Name",
                     "Data Source(Accession No.;Sequencing Platform;Sequencing Library;Insert Size;Experiment;Bases;Average Length;Status)",
                     "Barcode", "Source Niche", "Source Type", "Source Details",
                     "Collection Year", "Collection Month", "Collection Day", "Collection Time",
                     "Continent", "Country", "Region", "District", "City",
                     "Post Code", "Latitude", "Longitude",
                     "Serovar", "Subspecies", "Disease", "Antigenic Formulas",
                     "Lab Contact", "Phage Type", "Comment",
                     "Bio Project ID", "Project ID", "Sample ID", "Secondary Sample ID",
                     "Date Entered", "Release Date",
                     "Status", "Coverage", "N50", "Length", "Species",
                     "Contig Number(>=200 bp)", "Low Quality Bases", "Version", "Assembly Barcode")
#
# These names will be hardcoded into the following code.


# Data is in a data/ folder. CHANGE accordingly.
METADATA_FILE = "data/2020-01-30_enteritidis_metadata.tsv"

# Get the filename without the path and extension.
FILENAME = tools::file_path_sans_ext(basename(METADATA_FILE))
# Append "_filtered" before extension to the filtered table filename.
OUT_FILE = paste("data/", FILENAME, "_filtered.tsv", sep = "")



#### Data loading ####

# Load dataframe.  Argument fill and quote are used to correctly load the table.
# You can add elements to na.strings if there are other elements considered as NA.
# e.g.: You do not want to consider "ND" entries as a possible factor at all.
df = read.table(METADATA_FILE, header=TRUE, dec=".", sep="\t",
                na.strings=c("NA", " ", "", "NaN") , quote="", fill=TRUE)
# Data exploration
head(df)
colnames(df)
summary(df)

# See if metadata has not changed. If it has, you might need to change the code.
# make.names convert strings to valid R variable names.
make.names(metadata_columns) == colnames(df)

#### General filtering ####

#### _Project ID ####

# No Bio Project ID found. NOTE: this does not test the validity of the ID.
print(paste("Number of entries with no Bio Project ID:",
            sum(is.na(df$Bio.Project.ID))))
# Remove corresponding rows.
df = df[!is.na(df$Bio.Project.ID), ]

# Entries with a Bio Project should have a Sample ID, thus this step should be redundant.
print(paste("Number of entries with no Sample ID:", sum(is.na(df$Sample.ID))))
df = df[!is.na(df$Sample.ID), ]

#### _Collection Date ####

# The entry should at least have a Collection Year.
# NOTE: this only checks for missing dates, not abnormal ones.
print(paste("Number of entries with missing Collection Year:",
            sum(is.na(df$Collection.Year))))
df = df[!is.na(df$Collection.Year), ]

#### _Source ####

# Check existing values
summary(df$Source.Type)

## Source Type

print(paste("Number of entries with missing Source Type:",
            sum(is.na(df$Source.Type))))
df = df[!is.na(df$Source.Type), ]

# !!!: REMOVE the NA BEFORE checking for equality in a column.
# Also check to see if there are other spellings of ND.
print(paste("Number of entries with Not Determined Source Type:",
            sum(df$Source.Type == "ND")))
df = df[df$Source.Type != "ND", ]
# There are ND organisms with known Source Niche or Details but here we ignored them.

## Source Niche

summary(df$Source.Niche)

print(paste("Number of entries with missing Source Niche:",
            sum(is.na(df$Source.Niche))))
df = df[!is.na(df$Source.Niche), ]
# There are entries with missing ecosystem with known Source Type but here we ignored them.

# Uncomment if you need to get rid of rows with ND Niches.
# print(paste("Number of entries with Not Determined Source Type:",
#             sum(df$Source.Type == "ND")))
# df = df[df$Source.Type != "ND", ]

#### _Location ####

## Not missing Continent

# Find if there are entries with Country but missing Continent.
df[is.na(df$Continent) & !is.na(df$Country), ]
# Filtering is based solely on Continent field.
summary(df$Continent)

print(paste("Number of entries with missing Continent:",
            sum(is.na(df$Continent))))
df = df[!is.na(df$Continent), ]

# If you want to manually solve Unresolved continent with a Country:
# df[df$Continent == "Unresolved" & !is.na(df$Country), ]
# The same can be done with Region, District, City ...


#### Custom filtering ####

#### _From Europe ####
print(paste("Number of entries from Europe:", sum(df$Continent == "Europe")))
df = df[df$Continent == "Europe", ]

#### _Only human and food ####

# Source Niche
summary(df$Source.Niche)
# Adjust according to the levels you want to keep.
df = df[df$Source.Niche=="Human" | df$Source.Niche=="Livestock" | 
        df$Source.Niche=="Food" | df$Source.Niche=="Poultry", ]

# Source Type
summary(df$Source.Type)
# Adjust according to the levels you want to keep.
df = df[df$Source.Type!="Air", ]

# Another interesting filter is to split the "Species" column and only take 
# entries where the percentage is over a threshold.  I am not sure what does 
# the percentage represents though, and enterobase documentation is elusive about it. 
# df = tidyr::separate(df, col="Species", sep=";", into=c("Species", "Percentage"), remove=TRUE)
# Remove the "%" and cast to numeric.
# df$Percentage = as.numeric(sub("%", "", df$Percentage))
# df[df$Percentage < 60, ]


#### Other manipulations ####

# No filtering done on whether the source niche and type are coherent.


#### _Relabelling ####
## Change some values label.
# Add "Not Determined" as a possible value for the column.
levels(df$Source.Type) = c(levels(df$Source.Type), "Not Determined")
# Replace the value.
df$Source.Type[df$Source.Type=="ND/Others"] = "Not Determined"
# Check if correct.
# df[df$Source.Type=="Not Determined", ]

# Remove ND/Others after further thought.
# df = df[df$Source.Type!="ND/Others", ]

#### _Splitting columns####

# The problem with the Data Accession column is that while it should have 
# 8 elements, most entries have 6 or 7 elements (often status is missing) and
# some people try to add multiple experiments which make it more difficult for parsing.
# The best way to proceed right now is to only take 7 columns (the 8th one would
# not correspond to status) and fill in for entries with 6 elements.
splitted_cols = c("Accession No.", "Sequencing Platform", "Sequencing Library",
                  "Insert Size", "Experiment", "Bases", "Average Length")




# See frequency of all lengths
sort(table(lengths(strsplit(as.character(df[, 3]), ";"))), decreasing=TRUE)

# Make it so the col can be splitted by separate.
df[, 3] = as.character(df[, 3])

# separate should throw a warning for entries too big BUT you might want
# to remove them beforehand to avoid parsing unwanted data in columns.
# Use remove=FALSE if you want to keep the original column.
# Unfortunately, I could not find a way to not insert the new columns at the 
# position of the splitted column. 
# Otherwise, a Python script could be made for splitting with custom order.
df2 = tidyr::separate(df, col=3, sep=";", into=splitted_cols, remove=TRUE)
# Remove 2nd experiment placed after the Average Length.
df2$`Average Length` = as.numeric(sub(",.*", "", df2$`Average Length`))

# Only get short reads from Illumina.
df2 = df2[df2$`Sequencing Platform` == "ILLUMINA", ]
df2 = df2[df2$`Sequencing Library` == "Paired", ]

#### File writing ####
# Write the table as a tsv.  If using the "default" output file
# make sure a data/ folder is in the current working dir.
# Replace df2 by df if you do not want the splitted columns (might change index).
# Uncomment col.names if you want to have original names I.F.F. the columns are
# the same as the one from Enterobase 1.1.2.
write.table(df2, file=OUT_FILE, sep="\t", row.names=FALSE, quote=FALSE,
            # col.names=metadata_columns
            )

test = read.table(OUT_FILE, sep="\t", header = T)
