### Filtering rules

#### Filter file description

- `filter_gene_level.json`: Very light filtering steps applied to capture all possible viable transcripts for gene-level quantification. Note that intra-priming is not checked for in most robust categories as these transcripts could still be used in gene-level expression quantitation.
    - `full-splice_match`: apart from models with 1 or 2 exons which are filtered for intrapriming almost all other models with more than 2 exons are included.
    - `incomplete-splice_match`: models with more than 2 exons pass unless percentage of A in downstream sequence is 100%. Models with 1 or 2 exons are mostly filtered due to accumulation of 3' truncated reads due to 10X 3 protocol bias.
    - `novel_in_catalog`: models with up to 85% downstream A sequence enrichment will be included.
    - `novel_not_in_catalog`: models with with 2 or more exons and up to 85% downstream A sequence enrichment will be included. Mono-exonic models are mostly filtered as they do not represent a reliable structure.
    - `intergenic`: mono-exonic models 200-5500b in length are checked for intrapriming. Otherwise, most models are discarded.
    - `rest`: models that pass intra-priming filter, all SJs must be recorded as canonical, and the number of exons must exceed 2 (to remove mono-exon transcripts). 
- `filter_tx_level.json`: Strict filtering applied to capture only informative and reliable transcript models for transcript annotation purposes.
    - `full-splice_match`: tx category would be checked for intra-priming.
    - `incomplete-splice_match`: tx category requires at least 3 read coverage and passing intra-priming filter.
    - `novel_in_catalog`: either tx with at least 2 read coverage OR reads without intra-priming and RT switching that contain all canonical SJs.
    - `novel_not_in_catalog`: tx must not have more than 50(bp) differnece in the immediate upstream or downstream TSS/TTS sites. Either all SJ must be canonical OR at least 2 read coverage recorded for the tx to pass filter.
    - `genic`: preserve models that are spliced but belong to the 
    - `rest`: of tx must pass intra-priming and RT switching filter with all SJ marked as canonical with at least 2 read coverage and 2 exons recorded. Moreover, only coding trancripts pass this filter.
- `AS_event_filter.json`: TO BE COMPLETED

**Note**

`FL` column represents full length read count