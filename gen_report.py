import sys
import os
import re
import numpy as np 
import matplotlib; matplotlib.use('pdf')
import matplotlib.pyplot as plt
import pickle

from libs import plotter
from importlib import reload


def extract_sam_info(data, read_size, sam_file):
	"""Extract most information from sam file."""

	from collections import defaultdict
	#array indexed by read ID
	align_array = np.array([defaultdict(list) for i in range(read_size)])

	with open(sam_file, 'r') as infile:
		for line in infile:
			#skip header
			if re.search('^(?!@)', line):
				line = line.split('\t')
				[_id, flag, scaffold, pos, mapQ, cigar, _, _, _, seq, _] = line[:11]
				_id = int(_id)
				align_array[_id]['flag'].append(int(flag))
				align_array[_id]['scaffold'].append(scaffold)
				align_array[_id]['pos'].append(pos)
				align_array[_id]['mapQ'].append(int(mapQ))
				align_array[_id]['cigar'].append(cigar)
				align_array[_id]['seq'].append(seq)
				for col in line:
					if re.search('^NM:i:[0-9]+$', col):
						num_NM = int(col.split(':')[-1])
						align_array[_id]['num_mismatch'].append(num_NM)
					if re.search('^AS:i:-*[0-9]+$', col):
						AS = int(col.split(':')[-1])
						align_array[_id]['AS'].append(AS)

	#assert _id + 1 == read_size
	return align_array


def get_label_distribution(labels, aln_tool_list, src_dir, data, read_size):
	num_label = len(labels)
	stats = {aln_tool: {key: None for key in labels} for aln_tool in aln_tool_list}

	#stats[aln_tool][label] = #reads / #total: 4*8
	for aln_tool in aln_tool_list:
		with open('{0}/{1}/Pauto/ids/{2}_ecv_0_reads.cnt'.format(src_dir, aln_tool, data), 'r') as infile:
			for i in range(num_label):
				[num_reads, name] = infile.readline().strip().split()
				stats[aln_tool][labels[i]] = int(num_reads) / read_size

	return stats


def get_label_dis_bar(label_distribution, src_dir, labels, aln_tool_list, plot_figures):
	fig = plt.figure(figsize=(15, 10))
	for i in range(len(aln_tool_list)):
		ax = fig.add_subplot(2, 2, i+1)
		plotter.do_label_dis_bar(ax, aln_tool_list[i], label_distribution[aln_tool_list[i]])

	fig.savefig('{}/imgs/dis_bar.png'.format(src_dir))
	plot_figures.append(fig)


def draw_report_tables(label_distribution, aln_tool_list, src_dir, data, read_size, ecv_fpath):
	plotter.do_basic_stats(table_figures, ecv_fpath)
	plotter.do_label_dis_table(label_distribution, src_dir, aln_tool_list, table_figures)


def draw_report_imgs(aln_tool, label_array, src_dir, data, read_size, plot_figures):

	#draw the rest, extract info from sam file first
	sam_file = '{0}/{1}/Pauto/{2}_ecv_all.sam'.format(src_dir, aln_tool, data)
	path = '/'.join(sam_file.split('/')[:-2])+'/align_info'

	align_array = extract_sam_info(data, read_size, sam_file)
	pickle.dump(align_array, open(path, 'wb'))
	#align_array = pickle.load(open(path, 'rb'))

	#no clips
	if aln_tool == 'bowtie2-endtoend':
		nm_list = plotter.do_nm(aln_tool, align_array, label_array, read_size)
		plotter.plot_sam_dis(src_dir, nm_list, aln_tool, label_distribution[aln_tool], plot_figures, xlabel='#Mismatch', label='S')

		for label in (['P', 'S']):
			as_list = plotter.do_as(aln_tool, align_array, label_array, read_size, label)
			plotter.plot_sam_dis(src_dir, as_list, aln_tool, label_distribution[aln_tool], plot_figures, xlabel='Alignment Score', label=label)

	#no AS field, no clips
	elif aln_tool == 'bwa-endtoend':
		nm_list = plotter.do_nm(aln_tool, align_array, label_array, read_size)
		plotter.plot_sam_dis(src_dir, nm_list, aln_tool, label_distribution[aln_tool], plot_figures, xlabel='#Mismatch', label='S')

	else:
		nm_list = plotter.do_nm(aln_tool, align_array, label_array, read_size)
		plotter.plot_sam_dis(src_dir, nm_list, aln_tool, label_distribution[aln_tool], plot_figures, xlabel='#Mismatch', label='S')

		cr_list = plotter.do_cr(aln_tool, align_array, label_array, read_size)
		plotter.plot_sam_dis(src_dir, cr_list, aln_tool, label_distribution[aln_tool], plot_figures, xlabel='Clip%', label='C')
		for label in (['P', 'S', 'C']):
			as_list = plotter.do_as(aln_tool, align_array, label_array, read_size, label)
			plotter.plot_sam_dis(src_dir, as_list, aln_tool, label_distribution[aln_tool], plot_figures, xlabel='Alignment Score', label=label)


def get_label_array(info_file, read_size):
	label_array = np.array([None for i in range(read_size)])
	with open(info_file, 'r') as infile:
		idx = 0
		for line in infile:
			label_array[idx] = line.split('\t')[1]
			idx += 1
	return label_array

		
if __name__ == '__main__':
	src_dir, ecv_fpath, data, read_size = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
	all_pdf_fpath = src_dir+'/report.pdf'
	all_html_fpath = src_dir+'/report.html'
	aln_tool_list = ['bwa-mem', 'bowtie2-local', 'bwa-endtoend', 'bowtie2-endtoend'] #kart
	labels = ['N', 'F', 'M', 'P', 'S', 'C', 'O', 'X']
	read_size = int(read_size)
	plot_figures = []
	table_figures = []

	#label distribution table
	print("Generate label distribution graph")
	label_distribution = get_label_distribution(labels, aln_tool_list, src_dir, data, read_size)
	draw_report_tables(label_distribution, aln_tool_list, src_dir, data, read_size, ecv_fpath)

	#save label distribution bar, one img for each aligner
	get_label_dis_bar(label_distribution, src_dir, labels, aln_tool_list, plot_figures)

	#Draw distribution graph in terms of NM, CR, AS
	for aln_tool in aln_tool_list:
		print("Plot distribution graph from {}".format(aln_tool))
		info_file = src_dir+'/{0}/Pauto/ids/{1}_ecv_0_reads.info'.format(aln_tool, data)
		label_array = get_label_array(info_file, read_size)
		draw_report_imgs(aln_tool, label_array, src_dir, data, read_size, plot_figures)
	
	#regression_plot()
	template_fpath = os.path.dirname(sys.argv[0])+'/template.html'
	plotter.save_to_pdf(all_pdf_fpath, plot_figures, table_figures)
	plotter.save_to_html(src_dir+'/report.html', template_fpath)