#!/usr/bin/env python
# Author: Moe
# Generate a list of read names by using inclusion list file provided by SQ3 filter and read grouping file provide by collapse. Read names are then used to filter mapped bam.

import sys
import pandas as pd

# ../collapse/m84014_pbmc_0.01.collapsed.group.txt 
# ../sqanti3_filter/rules/m84014_pbmc_0.01/m84014_pbmc_0.01_inclusion-list.txt 
# m84014_pbmc_0.01
# ../pbmm2/m84014_pbmc_0.01_read_names.txt

sys.stdout.write(f"{sys.executable}\n")

# Collapse group output
group_file = sys.argv[1]
# List of PB.*** tx names as one column of \n separated values
list_file = sys.argv[2]
# Sample name (added as prefix to the output file)
wld = sys.argv[3]
# Output path
outpath = sys.argv[4]

with open(group_file, "r") as f:
    group_data = pd.read_csv(f, sep="\t", header=None, names=["tx_name", "read_name"])

with open(list_file, "r") as f:
    list_data = pd.read_csv(f, header=None, names=["tx_name"])
# Get the shared subset of reads
merged_data = list_data.merge(group_data, on="tx_name", how="inner")
print(merged_data)

# Test
# Number of rows in the resulting merged dataframe should be equal to `list_data`
assert merged_data.shape[0] == list_data.shape[0], "*** ERROR *** \nListed tx name count does not match collpase tx name count\n*** ERROR ***"

# Extract read names (mix number of values in each row delim=',')
reads_list = merged_data["read_name"].str.split(",").explode().tolist()
with open(outpath, "w") as f:
    for read in reads_list:
        f.write(read + "\n")