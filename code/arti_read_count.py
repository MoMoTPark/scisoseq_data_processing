#!/usr/bin/env python
# author: Moe
# Add an additional column to SQANTI3 Filter reasons output to reflect reads associated with the isoforms, generate a stacked bar plot of counts and artifact tx names

import sys, os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

sys.stdout.write(sys.executable + "\n")
reasons_path = sys.argv[1]
abund_path = sys.argv[2]
wld = sys.argv[3]
outfile = sys.argv[4]

print(reasons_path)
print(abund_path)
print(wld)
print(outfile)

with open(reasons_path, 'r') as f:
    reasons = pd.read_csv(f, sep="\t")

# print(reasons.head())

dedup_reasons = reasons.groupby('isoform').agg({
    'structural_category': 'first',
    'reasons': lambda x: ",".join(x),
}).reset_index()

# print(dedup_reasons.head())

# Sanity check: resulting deduplicated reasons dataframe should have the same number of records as unique isoforms values in the overall reasons dataframe
assert dedup_reasons.shape[0] == reasons['isoform'].unique().shape[0], "*** ERROR: Number of unique rows in reasons do not match deduplicated reasons df"

with open(abund_path, 'r') as f:
    abund = pd.read_csv(f, sep="\t", header=3)
    abund = abund.rename({"pbid":"isoform"}, axis=1)

# print(abund.head())

merged_df = abund.merge(dedup_reasons, on='isoform', how='inner')
# print(merged_df.head())

# Sanity check: resulting merged dataframe should always have the same record count as the deduplicated reasons dataframe
assert merged_df.shape[0] == dedup_reasons.shape[0], "*** ERROR: Rows between merged and deduplicated reasons dataframe do not match ***"

out_df = merged_df[['isoform', 'structural_category', 'reasons', 'count_fl']]

with open(outfile, 'w') as f:
    out_df.to_csv(f, sep='\t', index=False)

# Store transcript names of artifacts
with open(f"../sqanti3_filter/rules/{wld}/{wld}_artifact_tx_names.txt", 'w') as f:
    out_df['isoform'].to_csv(f, index=False, header=None)

plot_df = out_df.groupby(['structural_category', 'reasons'])['count_fl'].sum().unstack(fill_value=0)
ax = plot_df.plot(kind="bar", stacked=True)
plt.legend(title="Filtering reason", loc="right", bbox_to_anchor=[1.8, 0.5])
plt.xticks(rotation=40, ha='right')
def millions(x, pos):
    'The two args are the value and tick position'
    return '%1.1fM' % (x * 1e-6)
formatter = ticker.FuncFormatter(millions)
ax.yaxis.set_major_formatter(formatter)
plt.savefig(f'../sqanti3_filter/rules/{wld}/{wld}_artifact_read_count_stacked_bar.png', bbox_inches='tight')

sys.stdout.write("Finished processing!\n")