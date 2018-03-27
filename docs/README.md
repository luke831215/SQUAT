# SQUAT manual

## version: 1.0

This manual introduces all the information you need to know about SQUAT.

<br>

SQUAT stands for **Sequencing QUality Assessment Tool**. It allows users to examine their sequencing data and determines if they are truly representative of the original specie before conducting any further assembly experiments. The tool aligns sequencing reads against the reference genome with **[BWA](http://bio-bwa.sourceforge.net/)** and performs detailed analysis according to the SAM file generated from the alignment results.

It will generate an html report at the end of the analysis pipeline. The validity of the sequencing data will also be indicated in the report.

Version 1.0 of SQUAT was released under GPLv3 on XXX.	