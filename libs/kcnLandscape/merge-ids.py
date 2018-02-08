import sys

def gettype(label):
	return{
		'N': 'N',
		'P': 'u',
		'S': 's',
		'O': 'x',
		'M': 'M',
		'F': 'F',
		'X': 'SP'
	}.get(label, 0)


with open('/home/ph/mnt/biocollab/leon/{0}/ids/{0}_ecv_0_reads.info'.format(sys.argv[1]), 'r') as infile:
	with open('/home/ph/mnt/biocollab/leon/{0}/classify/all-reads_class2.csv'.format(sys.argv[1]), 'w') as w:
		for line in infile:
			[read_id, read_type, rpt_cnt] = line.strip().split('\t')		
			w.write('{0},{1},{2}\n'.format(int(read_id), gettype(read_type), int(rpt_cnt)))
		
