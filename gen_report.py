import sys
import optparse
import os
import re
import numpy as np 
import matplotlib; matplotlib.use('pdf')
import matplotlib.pyplot as plt
import pickle
import shutil

from libs import plotter
from importlib import reload

def build_subset(subset_label, label_list, in_fpath, out_fpath):
	id_list = []
	
	for label in subset_label:
		id_list = id_list + list(np.where(label_list == label)[0])
	id_list.sort()

	if len(id_list) == 0:
		raise ValueError("Subset size is zero.")
		return
	else:
		with open (in_fpath, 'r') as infile:
			with open(out_fpath, 'w') as w:
				current_idx, read_id = 0, 0
				length = len(id_list)
				for line in infile:
					if read_id != id_list[current_idx]:
						next(infile)
						next(infile)
						next(infile)
					else:
						w.write(line)
						for j in range(3):
							w.write(infile.readline())
						if current_idx + 1 == length:
							break
						current_idx += 1
					read_id += 1


def extract_sam_info(data, read_size, sam_file):
	"""Extract most information from sam file."""

	from collections import defaultdict
	#array indexed by read ID
	align_array = np.array([defaultdict(list) for i in range(read_size)])

	with open(sam_file, 'r') as infile:
		crt_id = None
		idx = -1
		for line in infile:
			#skip header
			if re.search('^(?!@)', line):
				line = line.split('\t')
				[_id, flag, scaffold, pos, mapQ, cigar, _, _, _, seq, _] = line[:11]
				if _id != crt_id:
					idx += 1
					crt_id = _id
				#_id = int(_id)
				align_array[idx]['id'].append(_id)
				align_array[idx]['flag'].append(int(flag))
				align_array[idx]['scaffold'].append(scaffold)
				align_array[idx]['pos'].append(pos)
				align_array[idx]['mapQ'].append(int(mapQ))
				align_array[idx]['cigar'].append(cigar)
				align_array[idx]['seq'].append(seq)
				for col in line:
					if re.search('^NM:i:[0-9]+$', col):
						num_NM = int(col.split(':')[-1])
						align_array[idx]['num_mismatch'].append(num_NM)
					if re.search('^AS:i:-*[0-9]+$', col):
						AS = int(col.split(':')[-1])
						align_array[idx]['AS'].append(AS)

	#assert _id + 1 == read_size
	return align_array


def extract_sam_info_primary(data, read_size, sam_file):
	"""Extract most information from sam file."""

	from collections import defaultdict
	#array indexed by read ID
	align_array = np.array([{} for i in range(read_size)])

	with open(sam_file, 'r') as infile:
		crt_id = None
		idx = -1
		for line in infile:
			#skip header
			if re.search('^(?!@)', line):
				line = line.split('\t')
				[_id, flag, scaffold, pos, mapQ, cigar, _, _, _, seq, _] = line[:11]
				if _id != crt_id:
					idx += 1
					crt_id = _id
					#_id = int(_id)
					align_array[idx]['id'] = int(_id)
					align_array[idx]['flag'] = int(flag)
					align_array[idx]['scaffold'] = scaffold
					align_array[idx]['pos'] = pos
					align_array[idx]['mapQ'] = int(mapQ)
					align_array[idx]['cigar'] = cigar
					align_array[idx]['seq'] = seq
					for col in line:
						if re.search('^NM:i:[0-9]+$', col):
							num_NM = int(col.split(':')[-1])
							align_array[idx]['num_mismatch'] = num_NM
						if re.search('^AS:i:-*[0-9]+$', col):
							AS = int(col.split(':')[-1])
							align_array[idx]['AS'] = AS

	#assert _id + 1 == read_size
	return align_array


def get_label_distribution(labels, aln_tool_list, src_dir, data, read_size):
	num_label = len(labels)
	stats = {aln_tool: {key: None for key in labels} for aln_tool in aln_tool_list}

	#stats[aln_tool][label] = #reads / #total: 4*8
	for aln_tool in aln_tool_list:
		with open('{0}/{1}/ids/{2}_0_reads.cnt'.format(src_dir, aln_tool, data), 'r') as infile:
			for i in range(num_label):
				[num_reads, name] = infile.readline().strip().split()
				#if 'backtrack' in aln_tool and labels[i] == 'C':
				#	stats[aln_tool][labels[i]]	= 'N/A'
				stats[aln_tool][labels[i]] = int(num_reads) / read_size

	return stats


def draw_label_dis_bar(label_dict, align_info_dict, src_dir, aln_tool_list, plot_figures):
	cigar_dict, poor_pct_list = {}, np.zeros(len(aln_tool_list))
	fig = plt.figure(figsize=(15, 10))
	ymax, ymin = 0, 0
	ax_list = [None] * len(aln_tool_list)
	for i, aln_tool in enumerate(aln_tool_list):
		ax_list[i] = fig.add_subplot(len(aln_tool_list) / 2, 2, i+1)
		cigar_dict[aln_tool], poor_pct_list[i], ax_ymax, ax_ymin = plotter.do_label_dis_bar(ax_list[i], align_info_dict[aln_tool], aln_tool, label_dict[aln_tool], thre)
		#x-axis
		ax_list[i].axhline(color='black')
		ymax = ax_ymax if ymax < ax_ymax else ymax
		ymin = ax_ymin if ymin > ax_ymin else ymin

	for ax in ax_list:
		ax.set_ylim([ymin - 10, ymax + 10])
	
	#add footnote under the barplot
	footnote = ("1. Bar above the x-axis: portion of properly-mapped reads"
				"\n"
				"2. Bar below the x-axis: portion of poorly-mapped reads"
				)
	plt.figtext(0.1, 0.05, footnote, va="bottom", ha="left")
	fig.savefig('{}/images/label_dis_bar.png'.format(src_dir))
	plot_figures.append(fig)
	plt.close()
	avg_poor_pct = "{:.1%}".format(np.sum(poor_pct_list) / len(poor_pct_list))
	return cigar_dict, avg_poor_pct


def draw_label_dis(label_distribution, aln_tool_list, src_dir, data, read_size, ecv_fpath, plot_figures):
	plotter.do_label_dis_table(label_distribution, src_dir, aln_tool_list, plot_figures)
	plotter.do_label_piechart(label_distribution, src_dir, aln_tool_list, plot_figures)


def draw_genome_eval_table(stats, src_dir, plot_figures):
	new_stats = list(stats)
	new_stats.insert(0, ["Reference sequence assembly", "Value"])
	fig = plt.figure(figsize=(15, 10))
	ax = fig.add_subplot(111)
	#hide axes
	fig.patch.set_visible(False)
	ax.axis('off')
	ax.axis('tight')
	#draw table and savefig
	ax.table(cellText=new_stats, cellLoc='center', loc='center', fontsize=15, bbox=(0.25, 0.25, 0.5, 0.5))
	ax.set_title('Genome evaluation stats', y = 0.8)
	#fig.tight_layout()
	fig.savefig('{}/images/eval_table.png'.format(src_dir))
	plot_figures.insert(0, fig)


def draw_dis_graph(aln_tool, label_array, align_array, src_dir, read_size, cigar_dict, plot_figures):

	#mismatch ratio
	nm_list = cigar_dict['S']
	if len(nm_list):
		plotter.plot_sam_dis(src_dir, nm_list, aln_tool, label_array, read_size, plot_figures, xlabel='Mismatch%', label='S', thre=thre['MR'])

	#clip ratio
	cr_list = cigar_dict['C']
	if len(cr_list):
		plotter.plot_sam_dis(src_dir, cr_list, aln_tool, label_array, read_size, plot_figures, xlabel='Clip%', label='C', thre=thre['CR'])
	
	#aln score
	#bwa-backtrack records no AS field
	if aln_tool == 'bwa-mem':
		for label in (['P', 'S', 'C']):
			as_list = plotter.do_as(align_array, label_array, label)
			if len(as_list):
				plotter.plot_sam_dis(src_dir, as_list, aln_tool, label_array, read_size, plot_figures, xlabel='Alignment Score', label=label)
		


def get_label_dict(data, aln_tool_list, read_size):
	label_dict = {}
	for aln_tool in aln_tool_list:
		label_dict[aln_tool] = np.array([None for i in range(read_size)])
		info_file = src_dir + '/{0}/ids/{1}_0_reads.info'.format(aln_tool, data)
		with open(info_file, 'r') as infile:
			idx = 0
			for line in infile:
				label_dict[aln_tool][idx] = line.split('\t')[1]
				idx += 1
	
	return label_dict


def fastq_examine(fpath, read_size):
	len_min, len_max = float("inf"), 0
	gc_count, total_len = 0, 9
	cnt = 0

	with open(fpath) as infile:
		while True:
			tmp = infile.readline()
			if tmp == '':
				break
			seq = infile.readline().strip()
			seq_len = len(seq)
			len_min = seq_len if len_min > seq_len else len_min
			len_max = seq_len if len_max < seq_len else len_max
			gc_count += seq.count('G') + seq.count('C')
			total_len += seq_len
			next(infile)
			next(infile)
			cnt += 1

	seq_gc = '{:.0%}'.format(gc_count / total_len)
	return len_min, len_max, seq_gc


		
def draw_basic_table(avg_poor_pct, fpath, src_dir, read_size, total_size, plot_figures):
	seq_name = fpath.split('/')[-1]
	len_min, len_max, seq_gc = fastq_examine(fpath, read_size)
	stats = [
			["File name", seq_name], ["No. of sequence", '{:,}'.format(total_size)],
			["Sample size", '{:,}'.format(read_size)],
			["Sequence length", "{0} - {1}".format(len_min, len_max)],
			["Avg. poorly mapped sequence%", avg_poor_pct],
			["GC%", seq_gc]
			]
	fig = plt.figure(figsize=(15, 10))
	ax = fig.add_subplot(111)
	#hide axes
	fig.patch.set_visible(False)
	ax.axis('off')
	ax.axis('tight')

	#draw table and savefig
	new_stats = list(stats)
	new_stats.insert(0, ["Seqeuencing reads information", "Value"])
	the_table = ax.table(cellText=new_stats, 
		colWidths=[0.4, 0.2], cellLoc='center', loc='center', bbox=(0.25, 0.25, 0.5, 0.5))
	the_table.set_fontsize(20)
	ax.set_title('Sequence basic stats', y=0.8)
	fig.set_tight_layout(True)
	fig.savefig('{}/images/basic_stats.png'.format(src_dir))
	#first page of the report
	plot_figures.insert(0, fig)

	return stats


if __name__ == '__main__':
	opt_parser = optparse.OptionParser(usage='usage: %prog [options] \"args\"')
	opt_parser_group = optparse.OptionGroup(opt_parser, "generate post-assembly report based on the output of SQUAT")

	#opt_parser_group.add_option('-o', dest='output_file', help='Output file name.', default = input + '_test/1M-reads.arff')
	opt_parser_group.add_option('-o', dest='out_dir', help = 'output directory')
	opt_parser_group.add_option('-i', dest='ecv_fpath', help = 'the path of read file')
	opt_parser_group.add_option('-d', dest='data', help='name of read data')
	opt_parser_group.add_option('-r', dest='sample_size', help = 'sample size, equals to read size if no random sampling')
	opt_parser_group.add_option('-t', dest='total_size', help = 'the read size of the whole dataset')
	opt_parser_group.add_option('-s', dest='subset', help = 'Return the subset of sequencing reads according to labels (in capitals, e.g. PSCO)', default = '')

	opt_parser.add_option_group(opt_parser_group)
	(options, args) = opt_parser.parse_args()

	if not options.out_dir:
		opt_parser.error('Need to specify the output directory')
	elif not options.ecv_fpath:
		opt_parser.error('Need to specify the path of fastq file')
	elif not options.data:
		opt_parser.error('Need to specify name of the reads')
	elif not options.sample_size:
		opt_parser.error('Need to specify the sample size')
	elif not options.total_size:
		opt_parser.error('Need to specify the read size')

	out_dir, ecv_fpath, data, read_size, total_size = options.out_dir, options.ecv_fpath, options.data, options.sample_size, options.total_size
	src_dir = out_dir + '/' + data

	#aln_tool_list = ['bwa-mem', 'bowtie2-local', 'bwa-backtrack', 'bowtie2-backtrack']
	thre = {}
	aln_tool_list = ['bwa-mem', 'bwa-backtrack']
	labels = ['P', 'S', 'C', 'O', 'M', 'F', 'N']
	read_size, total_size = int(read_size), int(total_size)
	align_info_dict = {}
	plot_figures = []
	#table_figures = []

	with open(out_dir + '/config') as infile:
		thre['PM'] = float(infile.readline().split(':')[1].strip())
		thre['MR'] = float(infile.readline().split(':')[1].strip())
		thre['CR'] = float(infile.readline().split(':')[1].strip())
		thre['OR'] = float(infile.readline().split(':')[1].strip())
		thre['NR'] = float(infile.readline().split(':')[1].strip())

	print("Generate label distribution graph")
	#label distribution table
	label_dict = get_label_dict(data, aln_tool_list, read_size)
	label_distribution = get_label_distribution(labels, aln_tool_list, src_dir, data, read_size)
	draw_label_dis(label_distribution, aln_tool_list, src_dir, data, read_size, ecv_fpath, plot_figures)

	
	#save alignment info for each aligner tool
	print('Extract alignment information from sam files')
	for aln_tool in aln_tool_list:
		sam_file = '{0}/{1}/{2}_ecv_all.sam'.format(src_dir, aln_tool, data)
		path = '/'.join(sam_file.split('/')[:-1])+'/align_info'
		try:
			align_info_dict[aln_tool] = pickle.load(open(path, 'rb'))
		except:
			align_info_dict[aln_tool] = extract_sam_info_primary(data, read_size, sam_file)
			pickle.dump(align_info_dict[aln_tool], open(path, 'wb'))

	#save label distribution bar and return cigar information
	cigar_dict, avg_poor_pct = draw_label_dis_bar(label_dict, align_info_dict, src_dir, aln_tool_list, plot_figures)

	#Plot distribution graph in terms of NM, CR, AS
	print("Plot label distribution graph")
	for aln_tool in aln_tool_list:
		dirname = "{0}/{1}/imgs".format(src_dir, aln_tool)
		if not os.path.isdir(dirname):
			os.makedirs(dirname)
		else:
			shutil.rmtree(dirname)
			os.makedirs(dirname)

		align_array = align_info_dict[aln_tool]
		draw_dis_graph(aln_tool, label_dict[aln_tool], align_array, src_dir, read_size, cigar_dict[aln_tool], plot_figures)

	print('Plot genome evaluation table')
	#draw genome evaluation table
	genome_stats = plotter.get_genome_eval_stat(src_dir)
	draw_genome_eval_table(genome_stats, src_dir, plot_figures)

	print('Plot basic stats table')
	#plot basic stats of sequences, put in the front of plot_figures
	basic_stats = draw_basic_table(avg_poor_pct, ecv_fpath, src_dir, read_size, total_size, plot_figures)

	#make report
	print('Writing report')
	all_pdf_fpath = src_dir+'/report.pdf'
	all_html_fpath = '{0}/{1}.html'.format(out_dir, data)
	#all_html_fpath = '{0}/{1}/post-assembly.html'.format(out_dir, data)
	template_fpath = os.path.dirname(sys.argv[0])+'/template/template.html'
	plotter.save_to_pdf(all_pdf_fpath, plot_figures)
	plotter.save_to_html(all_html_fpath, template_fpath, data, thre, aln_tool_list, label_distribution, basic_stats, genome_stats)

	#output subset reads if specified
	if options.subset != 'NONE':
		print("Building subset")
		dirname='{}/subset'.format(src_dir)
		if not os.path.isdir(dirname):
			os.makedirs(dirname)

		subset_labels = [label for label in sys.argv[5]]
		for aln_tool in aln_tool_list:
			out_fpath = '{0}/{1}_{2}.fastq'.format(dirname, aln_tool, sys.argv[5])
			build_subset(subset_labels, label_dict[aln_tool], ecv_fpath, out_fpath)