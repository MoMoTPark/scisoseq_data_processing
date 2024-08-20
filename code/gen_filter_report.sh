#!/usr/bin/bash


source ~/.bashrc
dir=$1
wld=$2
export R_LIBS_USER="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

./SQANTI3_filter_report_moe.R -d $dir -o $wld -u SQANTI3-5.2.2/utilities/ -f rules