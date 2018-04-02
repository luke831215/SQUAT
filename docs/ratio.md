# Mismatch, Clip, and N Ratio

<br>

## Summary

In this section, we introduce **mismatch ratio**, **clip ratio**, and **N ratio** for S, C, and N reads respectively. It is worthy of note that most of the time. BWA-backtrack does not report the condition of clips in SAM file. Thus, there is no clip rate analysis for such mapping algorithm.

| Label | Description|
|:-:|:-:|
| P | Perfectly-matched reads |
| S | Reads with substitution errors |
| C | Reads that contain clips |
| O | Reads with other errors|
| M | Multi-mapped reads |
| F | Unmapped reads |
| N | Reads that contain N |

<br>

## Mismatch Ratio

If according to SAM file, a read can be mapped to reference assembly but with substitution errors (type S), we then compute its mismatch ratio as follows,

	Mismatch ratio = No. of mismatch / read length

Reads whose mismatch ratio exceeds a specified threshold (default 20%) will be determined poorly-mapped.

SQUAT will plot the distribution of mismatch ratio for S reads in the report as follow,

![Mismatch Ratio Distribution](imgs/mismatch_ratio.png)
 
<br>

## Clip Ratio

If according to SAM file, a read can be mapped to reference assembly with clips marked on either side of the read (type C), we then compute its clip ratio as follows,

	Clip Ratio = total length of clips / read length

Reads whose clip ratio exceeds a specified threshold (default 30%) will be determined poor quality.

SQUAT will plot the distribution of clip ratio for C reads in the report as follow,

![Clip Ratio Distribution](imgs/clip_ratio.png)

<br>

## N Ratio
If a read contains at least one 'N', we can also compute its N ratio as follows,

	N Ratio = No. of N / read length

Reads whose N ratio exceeds a specified threshold (default 10%) will be determined  poor quality.
