import sys
import os
import re
import numpy as np 
import matplotlib; matplotlib.use('pdf')
import matplotlib.pyplot as plt
import pickle
import shutil

from libs import plotter
from importlib import reload

def build_subset(id_list, in_fpath, out_fpath):
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


def get_label_distribution(labels, aln_tool_list, src_dir, data, read_size):
	num_label = len(labels)
	stats = {aln_tool: {key: None for key in labels} for aln_tool in aln_tool_list}

	#stats[aln_tool][label] = #reads / #total: 4*8
	for aln_tool in aln_tool_list:
		with open('{0}/{1}/ids/{2}_0_reads.cnt'.format(src_dir, aln_tool, data), 'r') as infile:
			for i in range(num_label):
				[num_reads, name] = infile.readline().strip().split()
				if 'endtoend' in aln_tool and labels[i] == 'C':
					stats[aln_tool][labels[i]]	= 'N/A'
				else:
					stats[aln_tool][labels[i]] = int(num_reads) / read_size

	return stats


def get_label_dis_bar(label_dict, align_info_dict, src_dir, aln_tool_list, plot_figures):
	cigar_dict, poor_pct_list = {}, np.zeros(len(aln_tool_list))
	fig = plt.figure(figsize=(15, 10))
	for i in range(len(aln_tool_list)):
		ax = fig.add_subplot(len(aln_tool_list) / 2, 2, i+1)
		cigar_dict[aln_tool_list[i]], poor_pct_list[i] = plotter.do_label_dis_bar(ax, align_info_dict[aln_tool_list[i]], aln_tool_list[i], label_dict[aln_tool_list[i]])
		#x-axis
		ax.axhline(color='black')
		
	#add footnote under the barplot
	footnote = ("1. Bar above the x-axis: portion of reads with good quality"
				"\n"
				"2. Bar below the x-axis: portion of reads with bad quality"
				)
	plt.figtext(0.1, 0.05, footnote, va="bottom", ha="left")
	fig.savefig('{}/images/bar.png'.format(src_dir))
	plot_figures.append(fig)
	plt.close()
	avg_poor_pct = "{:.1%}".format(np.sum(poor_pct_list) / len(poor_pct_list))
	return cigar_dict, avg_poor_pct


def draw_report_tables(label_distribution, aln_tool_list, src_dir, data, read_size, ecv_fpath, plot_figures):
	plotter.do_label_dis_table(label_distribution, src_dir, aln_tool_list, plot_figures)


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


def draw_report_imgs(aln_tool, label_array, align_array, src_dir, data, read_size, cigar_dict, plot_figures):

	#no clips
	if aln_tool == 'bowtie2-endtoend':
		nm_list = plotter.do_nm(align_array, label_array)
		if len(nm_list):
			plotter.plot_sam_dis(src_dir, nm_list, aln_tool, label_array, read_size, plot_figures, xlabel='Mismatch%', label='S', cigar_list=cigar_dict['S'])

		for label in (['P', 'S']):
			as_list = plotter.do_as(align_array, label_array, label)
			if len(as_list):
				plotter.plot_sam_dis(src_dir, as_list, aln_tool, label_array, read_size, plot_figures, xlabel='Alignment Score', label=label, cigar_list=cigar_dict[label])

	#no AS field, no clips
	elif aln_tool == 'bwa-endtoend':
		nm_list = plotter.do_nm(align_array, label_array)
		if len(nm_list):
			plotter.plot_sam_dis(src_dir, nm_list, aln_tool, label_array, read_size, plot_figures, xlabel='Mismatch%', label='S', cigar_list=cigar_dict['S'])

	else:
		nm_list = plotter.do_nm(align_array, label_array)
		if len(nm_list):
			plotter.plot_sam_dis(src_dir, nm_list, aln_tool, label_array, read_size, plot_figures, xlabel='Mismatch%', label='S', cigar_list=cigar_dict['S'])

		cr_list = plotter.do_cr(align_array, label_array)
		if len(cr_list):
			plotter.plot_sam_dis(src_dir, cr_list, aln_tool, label_array, read_size, plot_figures, xlabel='Clip%', label='C', cigar_list=cigar_dict['C'])
		for label in (['P', 'S', 'C']):
			as_list = plotter.do_as(align_array, label_array, label)
			if len(as_list):
				plotter.plot_sam_dis(src_dir, as_list, aln_tool, label_array, read_size, plot_figures, xlabel='Alignment Score', label=label, cigar_list=cigar_dict[label])


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


		
def draw_basic_table(avg_poor_pct, fpath, src_dir, read_size, plot_figures):
	seq_name = fpath.split('/')[-1]
	len_min, len_max, seq_gc = fastq_examine(fpath, read_size)
	stats = [
			["File name", seq_name], ["No. of sequence", '{:,}'.format(read_size)],
			["Sequence length", "{0} - {1}".format(len_min, len_max)],
			["Avg. sequence% labeled as poor quality", avg_poor_pct],
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
	out_dir, ecv_fpath, data, read_size = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
	src_dir = out_dir + '/' + data

	aln_tool_list = ['bwa-mem', 'bowtie2-local', 'bwa-endtoend', 'bowtie2-endtoend']
	labels = ['P', 'S', 'C', 'O', 'M', 'F', 'N']
	read_size = int(read_size)
	align_info_dict = {}
	plot_figures = []
	#table_figures = []

	print("Generate label distribution graph")
	#label distribution table
	label_dict = get_label_dict(data, aln_tool_list, read_size)
	label_distribution = get_label_distribution(labels, aln_tool_list, src_dir, data, read_size)
	draw_report_tables(label_distribution, aln_tool_list, src_dir, data, read_size, ecv_fpath, plot_figures)

	
	#save alignment info for each aligner tool
	print('Extract alignment information from sam files')
	for aln_tool in aln_tool_list:
		sam_file = '{0}/{1}/{2}_ecv_all.sam'.format(src_dir, aln_tool, data)
		path = '/'.join(sam_file.split('/')[:-1])+'/align_info'
		try:
			align_info_dict[aln_tool] = pickle.load(open(path, 'rb'))
		except:
			align_info_dict[aln_tool] = extract_sam_info(data, read_size, sam_file)
			#pickle.dump(align_info_dict[aln_tool], open(path, 'wb'))

	#save label distribution bar and return cigar information
	cigar_dict, avg_poor_pct = get_label_dis_bar(label_dict, align_info_dict, src_dir, aln_tool_list, plot_figures)

	#Plot distribution graph in terms of NM, CR, AS
	print("Plot label distribution graph")
	for aln_tool in aln_tool_list:
		dirname = "{0}/{1}/imgs".format(src_dir, aln_tool)
		if not os.path.isdir(dirname):
			#shutil.rmtree(dirname)
			os.makedirs(dirname)

		align_array = align_info_dict[aln_tool]
		draw_report_imgs(aln_tool, label_dict[aln_tool], align_array, src_dir, data, read_size, cigar_dict[aln_tool], plot_figures)

	print('Plot genome evaluation table')
	#draw genome evaluation table
	genome_stats = plotter.get_genome_eval_stat(src_dir)
	draw_genome_eval_table(genome_stats, src_dir, plot_figures)

	print('Plot basic stats table')
	#plot basic stats of sequences, put in the front of plot_figures
	basic_stats = draw_basic_table(avg_poor_pct, ecv_fpath, src_dir, read_size, plot_figures)

	#make report
	print('Writing report')
	all_pdf_fpath = src_dir+'/report.pdf'
	all_html_fpath = '{0}/{1}.html'.format(out_dir, data)
	template_fpath = os.path.dirname(sys.argv[0])+'/template/template.html'
	plotter.save_to_pdf(all_pdf_fpath, plot_figures)
	plotter.save_to_html(all_html_fpath, template_fpath, data, aln_tool_list, label_distribution, basic_stats, genome_stats)

	#output subset reads if specified
	if len(sys.argv) == 6:
		print("Building subset")
		dirname='{}/subset'.format(src_dir)
		if not os.path.isdir(dirname):
			os.makedirs(dirname)

		label_list = [label for label in sys.argv[5]]
		for aln_tool in aln_tool_list:
			id_list = []
			for label in label_list:
				id_list = id_list + list(np.where(label_dict[aln_tool] == label)[0])
			id_list.sort()
			out_fpath = '{0}/{1}_{2}.fastq'.format(dirname, aln_tool, sys.argv[5])
			build_subset(id_list, ecv_fpath, out_fpath)