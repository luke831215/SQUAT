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

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
labels = ['P', 'S', 'C', 'O', 'M', 'F', 'N']
crt_thre = 0.3
N_thre = 0.1


def save_to_csv(stats, src_dir, aln_tool_list):
	with open('{}/label_dis/label_dis.csv'.format(src_dir), 'w') as w:
		cnt = 0
		writer = csv.writer(w, delimiter=',')
		writer.writerow(['']+aln_tool_list)
		for label in labels:
			writer.writerow([label] + list(stats[cnt]))
			cnt += 1


def fill_in_table(stats, main_div, soup, template_fpath):
	label_list = ['P', 'S', 'C', 'O', 'M', 'F', 'N']
	icon_dirpath = "{}/icons".format(os.path.dirname(template_fpath))
	label_div = main_div.findAll('div', {"class": 'label'})

	for i in range(len(label_div)):
		#remove original img tag
		tmp = label_div[i].find('img')
		tmp.replace_with('')

		#generate new icon
		data_uri = base64.b64encode(open("{0}/{1}.png".format(icon_dirpath, label_list[i]), 'rb').read()).decode('utf-8').replace('\n', '')
		img_src = soup.new_tag("img", src="data:image/png;base64,{0}".format(data_uri), **{'class': 'icon'})
		label_div[i].append(img_src)

	for i in range(len(stats)):
		for j in range(len(stats[0])):
			cell = main_div.find('td', id='{0}{1}'.format(labels[i], j+1))
			cell.string = stats[i][j]


def save_to_html(all_html_fpath, template_fpath, aln_tool_list, label_distribution):
	"""concat all figures into report.html by inserting figures into template.html"""

	src_dir = os.path.dirname(all_html_fpath)
	with open(template_fpath) as file:
		soup = BeautifulSoup(file, "lxml")
		main_div = soup.find('div', {"class": 'main'})
		
		#set up label dis. table
		fill_in_table(flatten(label_distribution, aln_tool_list), main_div, soup, template_fpath)

		#plot label dis bar first
		img_fpath = '{}/label_dis/bar.png'.format(src_dir)
		data_uri = base64.b64encode(open(img_fpath, 'rb').read()).decode('utf-8').replace('\n', '')
		new_div = soup.new_tag("div", id="label-dis-barchart", **{'class': 'inner'})
		img_src = soup.new_tag("img", src="data:image/png;base64,{0}".format(data_uri))
		new_div.append(img_src)
		main_div.append(new_div)

		#the rest of distribution graph for each aligner
		
		for aln_tool in aln_tool_list:
			cnt = 1
			files = glob.glob('{0}/{1}/imgs/*.png'.format(src_dir, aln_tool))
			files.sort(key=os.path.getmtime)
			for img_fpath in files:
				data_uri = base64.b64encode(open(img_fpath, 'rb').read()).decode('utf-8').replace('\n', '')
				new_div = soup.new_tag("div", id="{0}-{1}".format(aln_tool, cnt), **{'class': 'inner'})
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


def plot_sam_dis(src_dir, data, aln_tool, label_array, read_size, plot_figures, xlabel='', label='', cigar_list={}):
	hist, bins = np.histogram(data, bins=20)
	
	fig = plt.figure(figsize=(15, 10))
	#gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1]) 

	ax = fig.add_subplot(111)
	ax.bar(bins[:-1], hist.astype(np.float32) / hist.sum() * 100, width=(bins[1]-bins[0]), alpha=0.5, color='steelblue', linewidth=0, align='edge')
	ax.set_xlabel(xlabel)
	ax.set_ylabel('Percentile')

	#add legend
	num_read = np.sum(label_array == label)
	if 'Alignment Score' not in xlabel:
		ax.axvline(crt_thre, color='red', zorder=20)
		ratio_good = np.sum(data < crt_thre)/num_read
		ratio_bad = 1 - ratio_good
		ax.plot(1, 1, label='Below threshold (eligible): {:.1%}'.format(ratio_good), marker='', ls='')
		ax.plot(1, 1, label='Above threshold (poor): {:.1%}'.format(ratio_bad), marker='', ls='')

		#zoom in the graph if necessary
		if np.max(data) < 0.5:
			ax.set_xlim(0, 0.5)
			ax.text(crt_thre*2, 0.5, 'threshold', transform=ax.transAxes)
		else:
			ax.text(crt_thre, 0.5, 'threshold', transform=ax.transAxes)
	
	label_dis = num_read / read_size
	ax.plot(1, 1, label='No. of {0} reads: {1} ({2:.1%})'.format(label, num_read, label_dis), marker='', ls='')
	ax.legend(loc=1, prop={'size': 15})

	title = '{0} distribution of {1} reads\n{2}'.format(xlabel, label, aln_tool)
	ax.set_title(title)

	plot_figures.append(fig)

	#title = '{0} distribution of {1} reads {2}'.format(xlabel, label, aln_tool)
	plt.savefig('{0}/{1}/imgs/{2}_dis_{3}.png'.format(src_dir, aln_tool, xlabel, label))
	plt.close()


def draw_bar_bi(ax, field_list, read_size, idx, label='', thre=crt_thre):
	upper_cnt, lower_cnt = 0, 0
	for field in field_list:
		if field < thre:
			upper_cnt += 1
		else:
			lower_cnt += 1
	
	#print(label, upper_cnt, lower_cnt)
	value = round(upper_cnt / read_size * 100)
	ax.bar(idx, value, color=colors[idx], align='center', width=0.5)#, label=label)#label=r'${}_+$'.format(label))
	if value:
		ax.text(idx, value+10, '{}'.format(value), ha='center')
	
	value = round(-1 * lower_cnt / read_size * 100)	
	ax.bar(idx, value, color=colors[idx], align='center', width=0.5)#label=r'${}_-$'.format(label))
	if value:
		ax.text(idx, value-10, '{}'.format(value), ha='center')

	return lower_cnt


def draw_bar(ax, label_array, read_size, idx, label=''):
	constant = -1 if label == 'F' else 1

	value = round(int(constant * np.sum(label_array == label) / read_size * 100))
	ax.bar(idx, value, color=colors[idx], align='center', width=0.5)#, label=label)
	if value:
		ax.text(idx, value+10*constant, '{}'.format(value), ha='center')


def do_label_dis_bar(ax, align_array, aln_tool, label_array, crt_label=''):
	read_size = len(align_array)
	patch_handles = {}
	x_labels = ['P', 'S', 'C', 'O', 'M', 'F', 'N', ]
	below_cnt = np.sum(label_array == 'F')

	#P
	draw_bar(ax, label_array, read_size, 0, label='P')

	#S
	nm_list = do_nm(align_array, label_array)
	lower_cnt = draw_bar_bi(ax, nm_list, read_size, 1, label='S')
	below_cnt += lower_cnt
	
	#C
	cr_list = do_cr(align_array, label_array)
	lower_cnt = draw_bar_bi(ax, cr_list, read_size, 2, label='C')
	below_cnt += lower_cnt

	#O
	o_list = do_others(align_array, label_array)
	lower_cnt = draw_bar_bi(ax, o_list, read_size, 3, label='O')
	below_cnt += lower_cnt

	#M
	draw_bar(ax, label_array, read_size, 4, label='M')

	#F
	draw_bar(ax, label_array, read_size, 5, label='F')

	#N
	n_list = do_N(align_array, label_array)
	lower_cnt = draw_bar_bi(ax, n_list, read_size, 6, label='N', thre=N_thre)
	below_cnt += lower_cnt

	below_pct = below_cnt / read_size
	above_pct = 1 - below_pct

	ax.plot(1, 1, label='Above: {:.1%}'.format(above_pct), marker='', ls='')
	ax.plot(1, 1, label='Below: {:.1%}'.format(below_pct), marker='', ls='')
	ax.set_ylim(-100, 100)
	ax.set_xticks(np.arange(len(x_labels)))
	ax.set_xticklabels(x_labels)
	ax.set_ylabel('Percentage')
	ax.legend(loc=1, prop={'size': 10}, frameon=False)
	ax.set_title('{}'.format(aln_tool))

	#record barplot input for future use
	return {'P': None, 'S': nm_list, 'C': cr_list, 'O': o_list, 'M': None, 'F': None, 'N': n_list}


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
	#plt.savefig('{0}/label_dis/imgs/table.png'.format(src_dir))
	plt.close()

	save_to_csv(flatten(label_dis, aln_tool_list), src_dir, aln_tool_list)


def do_basic_stats(table_figures, ecv_fpath):
	pass


def do_nm(align_array, label_array):
	"""plot mismatch in unique reads with sub. error"""

	cnt = 0
	label_id_list = np.where(label_array == 'S')[0]
	nm_list = np.zeros(len(label_id_list))

	for _id in label_id_list:
		nm_list[cnt] = align_array[_id]['num_mismatch'][0]/len(align_array[_id]['seq'][0])
		cnt += 1

	return nm_list


def do_cr(align_array, label_array):
	"""clip rate, excluding bowtie2-endtoend (no clip labels)"""

	cnt = 0
	label_id_list = np.where(label_array == 'C')[0]
	cr_list = np.zeros(len(label_id_list))

	for _id in label_id_list:
		length, clip_size = 0, 0
		cigar = align_array[_id]['cigar'][0]
		m = re.findall('\d+[MIDNSHP=X]', cigar)
		for value in m:
			length += int(value[:-1])
			if value[-1] == 'S':
				clip_size += int(value[:-1])

		cr_list[cnt] = clip_size / length
		cnt += 1
	
	return cr_list


def do_as(align_array, label_array, label):
	""""alignment score, excluding bwa-endtoend (no AS field) and bowtie2-endtoend for clipped reads"""

	label_id_list = np.where(label_array == label)
	as_list = np.zeros(len(label_id_list[0]))
	align_subset = align_array[label_id_list]
	cnt = 0

	for alignment in align_subset:
		as_list[cnt] = alignment['AS'][0]
		cnt += 1

	return as_list


def do_others(align_array, label_array):
	cnt = 0
	label_id_list = np.where(label_array == 'O')[0]
	o_list = np.zeros(len(label_id_list))

	for _id in label_id_list:
		length, o_size = 0, 0
		cigar = align_array[_id]['cigar'][0]
		m = re.findall('\d+[MIDNSHP=X]', cigar)
		for value in m:
			length += int(value[:-1])
			if value[-1] in ['I', 'D', 'P', '=', 'X']:
				o_size += int(value[:-1])

		o_list[cnt] = o_size / length
		cnt += 1
	
	return o_list


def do_N(align_array, label_array):
	cnt = 0
	label_id_list = np.where(label_array == 'N')[0]
	n_list = np.zeros(len(label_id_list))
	
	for _id in label_id_list:
		seq = align_array[_id]['seq'][0]
		if 'N' in seq:
			n_size = seq.count('N')
			n_list[cnt] = n_size / len(seq)
			cnt += 1

	return n_list


def regression_plot():
	pass

