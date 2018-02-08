#!/usr/bin/perl
# Merge labels of ids/*.txt to all-reads_class.txt
# By Yu-Jung Chang, last modified: 2017/04/14

# Usage
$numArgs = $#ARGV + 1;
if ($numArgs != 2) {
	print "Usage:\n";
	print "./merge-ids.pl idsPath PrjPrefix\n";
	print "./merge-ids.pl /mnt/biocollab/FEEC/D1-SRR352384 D1\n";
	exit;
}

$idsPath = $ARGV[0];
$PrjPrefix = $ARGV[1];

$cmd = sprintf "awk '{ printf \"%%s\\tN\\n\", \$1 }' %s/ids/%s_ecv_1_contain_N.ids > %s/classify/tmp_class.txt\n",
        $idsPath, $PrjPrefix, $idsPath;
print $cmd;
system $cmd;
$cmd = sprintf "awk '{ printf \"%%s\\tF\\n\", \$1 }' %s/ids/%s_ecv_2_unmappable.ids >> %s/classify/tmp_class.txt\n",
        $idsPath, $PrjPrefix, $idsPath;
print $cmd;
system $cmd;
$cmd = sprintf "awk '{ printf \"%%s\\tM\\n\", \$1 }' %s/ids/%s_ecv_3_mappable_multi.ids >> %s/classify/tmp_class.txt\n",
        $idsPath, $PrjPrefix, $idsPath;
print $cmd;
system $cmd;
$cmd = sprintf "awk '{ printf \"%%s\\tU0\\n\", \$1 }' %s/ids/%s_ecv_4_mappable_unique_noerr.ids >> %s/classify/tmp_class.txt\n",
        $idsPath, $PrjPrefix, $idsPath;
print $cmd;
system $cmd;
$cmd = sprintf "awk '{ printf \"%%s\\tUs\\n\", \$1 }' %s/ids/%s_ecv_5_mappable_unique_suberr.ids >> %s/classify/tmp_class.txt\n",
        $idsPath, $PrjPrefix, $idsPath;
print $cmd;
system $cmd;
$cmd = sprintf "awk '{ printf \"%%s\\tUx\\n\", \$1 }' %s/ids/%s_ecv_6_mappable_unique_othererr.ids >> %s/classify/tmp_class.txt\n",
        $idsPath, $PrjPrefix, $idsPath;
print $cmd;
system $cmd;

$cmd = sprintf "awk '{ printf \"%%s\\tV\\n\", \$1 }' %s/ids/%s_ecv_7_mappable_unique_altsite.ids >> %s/classify/tmp_class.txt\n",
        $idsPath, $PrjPrefix, $idsPath;
print $cmd;
system $cmd;

$cmd = sprintf "awk '{ printf \"%%s\\tR\\n\", \$1 }' %s/ids/%s_ecv_8_mappable_ref_contain_N.ids >> %s/classify/tmp_class.txt\n",
        $idsPath, $PrjPrefix, $idsPath;
print $cmd;
system $cmd;

$cmd = sprintf "sort -n %s/classify/tmp_class.txt -o %s/classify/all-reads_class.txt\n",
        $idsPath, $idsPath;
print $cmd;
system $cmd;
$cmd = sprintf "rm %s/classify/tmp_class.txt\n", $idsPath;
print $cmd;
system $cmd;
