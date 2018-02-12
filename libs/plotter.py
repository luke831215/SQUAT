import os
import re
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import gridspec
import glob
import base64
from matplotlib.backends.backend_pdf import PdfPages
from bs4 import BeautifulSoup
import xlsxwriter
import csv

labels = ['P', 'S', 'C', 'O', 'M', 'F', 'N']

def save_to_csv(stats, src_dir, aln_tool_list):
	with open('{}/label_dis/label_dis.csv'.format(src_dir), 'w') as w:
		cnt = 0
		writer = csv.writer(w, delimiter=',')
		writer.writerow(['']+aln_tool_list)
		for label in labels:
			writer.writerow([label] + list(stats[cnt]))
			cnt += 1


def fill_in_table(stats, main_div):
	for i in range(len(stats)):
		for j in range(len(stats[0])):
			cell = main_div.find('td', id='{0}{1}'.format(labels[i], j+1))
			try:
				cell.string = stats[i][j]
			except:
				print(i, j)
				return


def save_to_html(all_html_fpath, template_fpath, aln_tool_list, label_distribution):
	"""concat all figures into report.html by inserting figures into template.html"""

	src_dir = os.path.dirname(all_html_fpath)
	with open(template_fpath) as file:
		soup = BeautifulSoup(file, "lxml")
		main_div = soup.find('div', class_='main')
		#set up label dis. table
		fill_in_table(flatten(label_distribution, aln_tool_list), main_div)

		cnt = 0
		for aln_tool in ['label_dis'] + aln_tool_list:
			files = glob.glob('{0}/{1}/imgs/*.png'.format(src_dir, aln_tool))
			files.sort(key=os.path.getmtime)
			for img_name in files:
				data_uri = base64.b64encode(open(img_name, 'rb').read()).decode('utf-8').replace('\n', '')
				new_div = soup.new_tag("div", id="inner", class_="{0}-{1}".format(aln_tool, cnt))
				img_src = soup.new_tag("img", src="data:image/png;base64,{0}".format(data_uri))
				new_div.append(img_src)
				main_div.append(new_div)
				cnt += 1

		with open(all_html_fpath, 'w') as w:
			w.write(str(soup))


def save_to_pdf(all_pdf_fpath, plot_figures, table_figures):
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


def plot_sam_dis(src_dir, data, aln_tool, label_dis, plot_figures, xlabel='', label=''):
	hist, bins = np.histogram(data, bins=20)
	
	fig = plt.figure(figsize=(15, 10))
	gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 

	#first plot
	ax = fig.add_subplot(gs[0])
	#second plot: dis bar
	ax2 = fig.add_subplot(gs[1])

	ax.bar(bins[:-1], hist.astype(np.float32) / hist.sum() * 100, width=(bins[1]-bins[0]), alpha=0.5, color='steelblue', linewidth=0, align='edge')
	ax.set_xlabel(xlabel)
	ax.set_ylabel('Percentile')
	
	if label:
		title = '{0} distribution of {1} reads\n{2}'.format(xlabel, label, aln_tool)
		ax.set_title(title)
		do_label_dis_bar(ax2, aln_tool, label_dis, label)

	#im = plt.imread(insert_fig)
	#ax.imshow(im)
	#newax.axis('off')

	plot_figures.append(fig)

	title = '{0} distribution of {1} reads {2}'.format(xlabel, label, aln_tool)
	plt.savefig('{0}/{1}/imgs/{2}.png'.format(src_dir, aln_tool, title.replace(" ", "_")))
	plt.close()


def do_label_dis_bar(ax, aln_tool, label_dis, crt_label=''):
	colors ='rgbymc'
	patch_handles = {}
	x_labels = ['Mappable', 'Repeat', 'Unmapped']
	bottom = np.zeros(len(x_labels))
	value_dict = {}
	for key in label_dis:
		value_dict[key] = label_dis[key]*100

	#P, S
	patch_handles['P'] = ax.bar(0, value_dict['P'], color=colors[0],
	align='center', bottom=bottom[0], width = 0.5, label='P')
	bottom[0] += value_dict['P']
	patch_handles['S'] = ax.bar(0, value_dict['S'], color=colors[1],
	align='center', bottom=bottom[0], width = 0.5, label='S')
	bottom[0] += value_dict['S']

	#M
	patch_handles['M'] = ax.bar(1, value_dict['M'], color=colors[2],
	align='center', bottom=bottom[1], width = 0.5, label='M')
	bottom[1] += value_dict['M']

	#F, C, O
	patch_handles['F'] = ax.bar(2, value_dict['F'], color=colors[3],
	align='center', bottom=bottom[2], width = 0.5, label='F')
	bottom[2] += value_dict['F']

	patch_handles['C'] = ax.bar(2, value_dict['C'], color=colors[4],
	align='center', bottom=bottom[2], width = 0.5, label='C')
	bottom[2] += value_dict['C']

	patch_handles['O'] = ax.bar(2, value_dict['O'], color=colors[5],
	align='center', bottom=bottom[2], width = 0.5, label='O')
	bottom[2] += value_dict['O']

	#Add label in the middle of bar
	if crt_label:
		patch = patch_handles[crt_label].get_children()[0]
		coo_ax = patch.get_xy()
		x = 0.5*patch.get_width() + coo_ax[0]
		y = 0.5*patch.get_height() + coo_ax[1]
		ax.text(x, y, '*', ha='center', fontsize=30)

	ax.set_xticks(np.arange(3))
	ax.set_xticklabels(x_labels)
	ax.set_ylabel('Percentage')
	ax.legend(loc=1)
	ax.set_title('{}'.format(aln_tool))
	#return ax.get_figure()


def add_text(ax, text, fontsize='x-large', weight='medium'):
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.xaxis.set_ticks_position('none')
	ax.yaxis.set_ticks_position('none')
	ax.text(0.5, 0.5, text, fontsize=fontsize, va="center", ha='center', weight=weight)


def flatten(label_distribution, aln_tool_list):
	stats = np.array([None] * len(aln_tool_list) * len(labels))
	any_key = list(label_distribution.keys())[0]
	idx = 0

	for label in labels:
		for aln_tool in aln_tool_list:
			stats[idx] = '{:.1%}'.format(label_distribution[aln_tool][label])
			idx += 1

	#7*4 array
	return np.ravel(stats).reshape(len(labels), len(aln_tool_list))


def do_label_dis_table(label_dis, src_dir, aln_tool_list, table_figures):
	descriptions = ['P (no error)', 'S (substitution error)', 'C (contain clips)', 'O (other error)', 'M (multi-mapped)', 'F (unmapped)', 'N (contain N)']
	axes = []
	num_col = len(aln_tool_list)+3
	num_row = 8
	
	fig = plt.figure(figsize=(15, 10))
	gs = gridspec.GridSpec(num_row, num_col, width_ratios=[0.8, 1.1, 1.1] + [1]*(num_col-3))
	# set zero spacing between axes.
	gs.update(wspace=0, hspace=0)

	#Header 
	add_text(plt.subplot(gs[0, :3]), 'Label Distribution Table', weight='bold')
	
	#Unique reads
	add_text(plt.subplot(gs[1:5, 0]), 'Uniquely\nmapped\nReads')
	#P S C O
	add_text(plt.subplot(gs[1, 1:3]), '(P) Perfectly-mapped Reads')
	add_text(plt.subplot(gs[2, 1:3]), '(S) Reads with Substitution Error')
	add_text(plt.subplot(gs[3, 1:3]), '(C) Reads that Contain Clips')
	add_text(plt.subplot(gs[4, 1:3]), '(O) Reads that contain other error')
	#M U F
	add_text(plt.subplot(gs[5, 0:3]), '(M) Multi-mapped Reads')
	add_text(plt.subplot(gs[6, 0:3]), '(F) Unmapped Reads')
	add_text(plt.subplot(gs[7, 0:3]), '(N) Reads that contain N')

	#columns of stats from each aligner
	for j in range(3, num_col):
		aln_tool = aln_tool_list[j-3].replace('-', '\n')
		add_text(plt.subplot(gs[0, j]), aln_tool, weight='bold')
		for i in range(1, num_row):
			add_text(plt.subplot(gs[i, j]), '{:.1%}'.format(label_dis[aln_tool_list[j-3]][labels[i-1]]))
			

	
	#title = 'Label Distribution Table'
	#plt.title(title, fontdict={'fontsize': 30})
	table_figures.append(fig)
	plt.savefig('{0}/label_dis/imgs/table.png'.format(src_dir))
	plt.close()

	save_to_csv(flatten(label_dis, aln_tool_list), src_dir, aln_tool_list)


def do_basic_stats(table_figures, ecv_fpath):
	pass


def do_nm(aln_tool, align_array, label_array, read_size):
	"""plot mismatch in unique reads with sub. error"""

	cnt = 0	
	nm_list = np.zeros(read_size)

	for _id in np.where(label_array == 'S')[0]:
		nm_list[cnt] = align_array[_id]['num_mismatch'][0]
		cnt += 1

	#plot_sam_dis(nm_list[:cnt], '#Mismatch', '#Mismatch Distribution in S - {}'.format(aln_tool), insert_fig)
	return nm_list[:cnt]


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
	
	#plot_sam_dis(cr_list[:cnt], 'Clip Rate', 'Clip Rate Distribution in C - {}'.format(aln_tool), insert_fig)
	return cr_list[:cnt]


def do_as(aln_tool, align_array, label_array, read_size, label):
	""""alignment score, excluding bwa-endtoend (no AS field) and bowtie2-endtoend for clipped reads"""

	as_list = np.zeros(read_size)
	align_subset = align_array[np.where(label_array == label)]
	cnt = 0

	for alignment in align_subset:
		as_list[cnt] = alignment['AS'][0]
		cnt += 1

	#plot_sam_dis(as_list[:cnt], '#Alignment Score - {}'.format(label), 'Alignment Score Distribution in {0} - {1}'.format(label, aln_tool), insert_fig)
	return as_list[:cnt]


def regression_plot():
	pass

