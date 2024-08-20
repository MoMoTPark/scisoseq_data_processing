#!/usr/bin/env python
# Author: Moe
# This script takes the abundance output from collapse and calculates the normalised counts values; then generates a SQANTI3 QC input compatible file


import sys
import pandas as pd


sys.stdout.write(f"{sys.executable}\n")
# Collapse abundance file as input
in_file = sys.argv[1]
wld = sys.argv[2]
with open(in_file, 'r') as f:
    data = pd.read_csv(f, sep="\t", header=3)

out_data = data[["pbid", "count_fl"]]
# Add column with normalised count values
out_data["norm_fl"] = out_data["count_fl"] / out_data["count_fl"].sum()
with open(f"../collapse/{wld}.collapsed.abundance.sq3.txt", 'w') as f:
    out_data.to_csv(f, sep="\t", index=False)