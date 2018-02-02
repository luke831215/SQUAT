import sys
import os
import re
import numpy as np 
import matplotlib; matplotlib.use('pdf')
import matplotlib.pyplot as plt
import pickle
import mpld3
import glob
import base64

from matplotlib.backends.backend_pdf import PdfPages
from bs4 import BeautifulSoup


plot_figures = []
table_figures = []


def save_to_html(all_html_fpath):
	"""concat all figures into report.html by inserting figures into template.html"""
	template_fpath = os.path.dirname(sys.argv[0])+'/template.html'
	with open(template_fpath) as file:
		soup = BeautifulSoup(file, "lxml")
		main_div = soup.find('div', class_='main')
		cnt = 0
		#for table in table_figures:
		files = glob.glob('{}/imgs/*.png'.format(src_dir))
		files.sort(key=os.path.getmtime)
		for img_name in files:
			data_uri = base64.b64encode(open(img_name, 'rb').read()).decode('utf-8').replace('\n', '')
			new_div = soup.new_tag("div", id="F{}".format(cnt))
			img_src = soup.new_tag("img", src="data:image/png;base64,{0}".format(data_uri))
			new_div.append(img_src)
			main_div.append(new_div)
			cnt += 1
		with open(all_html_fpath, 'w') as w:
			w.write(str(soup))

def save_to_pdf(all_pdf_fpath):
	"""concat all figures into report.pdf"""

	all_pdf_file = PdfPages(all_pdf_fpath)
	for figure in table_figures:
		all_pdf_file.savefig(figure, bbox_inches='tight')
	for figure in plot_figures:
		try:
			all_pdf_file.savefig(figure)
		except:
			print(figure)
			sys.exit(1)
	try:  # for matplotlib < v.1.0
	    d = all_pdf_file.infodict()
	    d['Title'] = 'QUAST full report'
	    d['Author'] = 'QUAST'
	    import datetime
	    d['CreationDate'] = datetime.datetime.now()
	    d['ModDate'] = datetime.datetime.now()
	except AttributeError:
	    pass
	all_pdf_file.close()
	plt.close('all')  # closing all open figures


def plot_sam_dis(data, xlabel, title):
	hist, bins = np.histogram(data, bins=20)
	
	fig, ax = plt.subplots(figsize=(15,10))
	ax.bar(bins[:-1], hist.astype(np.float32) / hist.sum() * 100, width=(bins[1]-bins[0]), alpha=0.5, color='steelblue', linewidth=0)
	plt.xlabel(xlabel)
	plt.ylabel('Percentile')
	plt.title(title)
	plot_figures.append(fig)
	plt.savefig('{0}/imgs/{1}.png'.format(src_dir, title.replace(" ", "")))
	plt.close()


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


def draw_report_tables(aln_tool_list, src_dir, data, read_size):
	read_dict = {}
	labels = ['N (contain N)', 'F (unmapped)', 'M (multi-mapped)', 'P (no error)', 'S (substitution error)', 'C (clips)', 'O (other error)', 'X (special)']
	read_dict = [[] for i in range(10)]
	stats = [[None for j in range(len(aln_tool_list))] for i in range(8)]

	for align_tool in aln_tool_list:
		idx = 0
		with open('{0}/{1}/Pauto/ids/{2}_ecv_0_reads.cnt'.format(src_dir, align_tool, data), 'r') as infile:
			for line in infile:
				[num_reads, name] = line.strip().split()
				read_dict[idx].append(int(num_reads))
				idx += 1
			total = int(num_reads)

	read_dict = np.array(read_dict[:len(labels)])
	for i in range(7):
		for j in range(read_dict.shape[1]):
			stats[i][j] = '{:.2%}'.format(read_dict[i][j] / total)

	fig, ax = plt.subplots(figsize=(15,10))
	ax.axis('off')
	the_table = ax.table(cellText=stats,
						 rowLabels=labels,
						 colLabels=aln_tool_list,
						 loc='center'
						  )
	the_table.scale(1, 4)
	the_table.set_fontsize(14)
	plt.subplots_adjust(0.25)

	title = 'Label Distribution Table'
	plt.title(title, fontdict={'fontsize': 30})
	table_figures.append(fig)
	plt.savefig('{0}/imgs/{1}.png'.format(src_dir, title.replace(" ", "")))
	plt.close()


def regression_plot():
	pass

def do_nm(aln_tool, align_array, label_array, read_size):
	"""plot mismatch in unique reads with sub. error"""

	cnt = 0
	nm_list = np.zeros(read_size)

	for _id in np.where(label_array == 'S')[0]:
		nm_list[cnt] = align_array[_id]['num_mismatch'][0]
		cnt += 1

	plot_sam_dis(nm_list[:cnt], '#Mismatch', '#Mismatch Distribution in S - {}'.format(aln_tool))


def do_cr(aln_tool, align_array, label_array, read_size):
	"""clip rate, excluding bowtie2-endtoend (no clip labels)"""

	cnt = 0
	cr_list = np.zeros(read_size)

	for _id in np.where(label_array == 'C')[0]:
		length, clip_size = 0, 0
		cigar = align_array[_id]['cigar'][0]
		m = re.findall('\d+[MIDNSHP=X]', cigar)
		for value in m:
			length += int(value[:-1])
			if value[-1] == 'S':
				clip_size += int(value[:-1])

		cr_list[cnt] = clip_size / length
		cnt += 1
	
	plot_sam_dis(cr_list[:cnt], 'Clip Rate', 'Clip Rate Distribution in C - {}'.format(aln_tool))


def do_as(aln_tool, align_array, label_array, read_size, label):
	""""alignment score, excluding bwa-endtoend (no AS field) and bowtie2-endtoend for clipped reads"""

	as_list = np.zeros(read_size)
	align_subset = align_array[np.where(label_array == label)]
	cnt = 0

	for alignment in align_subset:
		as_list[cnt] = alignment['AS'][0]
		cnt += 1

	plot_sam_dis(as_list[:cnt], '#Alignment Score - {}'.format(label), 'Alignment Score Distribution in {0} - {1}'.format(label, aln_tool))


def draw_report_imgs(aln_tool, label_array, src_dir, data, read_size):
	sam_file = '{0}/{1}/Pauto/{2}_ecv_all.sam'.format(src_dir, aln_tool, data)
	path = '/'.join(sam_file.split('/')[:-2])+'/align_info'

	align_array = extract_sam_info(data, read_size, sam_file)
	pickle.dump(align_array, open(path, 'wb'))
	#align_array = pickle.load(open(path, 'rb'))

	#no clips
	if aln_tool == 'bowtie2-endtoend':
		do_nm(aln_tool, align_array, label_array, read_size)
		for label in (['P', 'S']):
			do_as(aln_tool, align_array, label_array, read_size, label)

	#no AS field, no clips
	elif aln_tool == 'bwa-endtoend':
		do_nm(aln_tool, align_array, label_array, read_size)

	else:
		do_nm(aln_tool, align_array, label_array, read_size)
		do_cr(aln_tool, align_array, label_array, read_size)
		for label in (['P', 'S', 'C']):
			do_as(aln_tool, align_array, label_array, read_size, label)


def get_label_array(info_file, read_size):
	label_array = np.array([None for i in range(read_size)])
	with open(info_file, 'r') as infile:
		idx = 0
		for line in infile:
			label_array[idx] = line.split('\t')[1]
			idx += 1
	return label_array

		
if __name__ == '__main__':
	src_dir, data, read_size = sys.argv[1], sys.argv[2], sys.argv[3]
	#src_dir = os.path.dirname(all_pdf_fpath)
	all_pdf_fpath = src_dir+'/report.pdf'
	all_html_fpath = src_dir+'/report.html'
	aln_tool_list = ['kart', 'bwa-mem', 'bowtie2-local', 'bwa-endtoend', 'bowtie2-endtoend']
	read_size = int(read_size)

	#label distribution table
	print("Generate label distribution graph")
	draw_report_tables(aln_tool_list, src_dir, data, read_size)

	#plots
	for aln_tool in aln_tool_list:
		print("Plot distribution graph from {}".format(aln_tool))
		info_file = src_dir+'/{0}/Pauto/ids/{1}_ecv_0_reads.info'.format(aln_tool, data)
		label_array = get_label_array(info_file, read_size)
		draw_report_imgs(aln_tool, label_array, src_dir, data, read_size)
	
	#regression_plot()
	save_to_pdf(all_pdf_fpath)
	save_to_html(src_dir+'/report.html')