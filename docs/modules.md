# Analysis Modules

First, SQUAT will randomly sample 1M entries of input reads and provide the **percentage of poor quality reads** based on overall assessment of the alignment experiments. (You can find how the percentage is calculated [here](seq_map.md#barchart))

It goes on to manipulate various analysis modules to get to the nitty-gritty, listed as follows:
 
- [Basic Statistics](basic_stats.md)
- [Read Labeling](read_labeling.md)
- [Mismatch, Clip, and N Ratio](ratio.md)
- [Alignment Score](aln_score.md)