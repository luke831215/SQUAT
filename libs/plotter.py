import os
import re
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import gridspec
import glob
import base64
from matplotlib.backends.backend_pdf import PdfPages
from bs4 import BeautifulSoup

def save_to_html(all_html_fpath, template_fpath):
	"""concat all figures into report.html by inserting figures into template.html"""

	src_dir = os.path.dirname(all_html_fpath)
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


def plot_sam_dis(src_dir, data, align_tool, label_dis, plot_figures, xlabel='', label=''):
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
		title = '{0} distribution of {1} reads\n{2}'.format(xlabel, label, align_tool)
		ax.set_title(title)
		do_label_dis_bar(ax2, align_tool, label_dis, label)

	#im = plt.imread(insert_fig)
	#ax.imshow(im)
	#newax.axis('off')

	plot_figures.append(fig)

	title = '{0} distribution of {1} reads {2}'.format(xlabel, label, align_tool)
	plt.savefig('{0}/imgs/{1}.png'.format(src_dir, title.replace(" ", "_")))
	plt.close()


def do_label_dis_bar(ax, align_tool, label_dis, crt_label=''):
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
	ax.set_title('{}'.format(align_tool))
	return ax.get_figure()


def do_label_dis_table(label_distribution, src_dir, aln_tool_list, table_figures):
	descriptions = ['N (contain N)', 'F (unmapped)', 'M (multi-mapped)', 'P (no error)', 'S (substitution error)', 'C (contain clips)', 'O (other error)', 'X (special)']
	stats = np.array([None] * len(aln_tool_list) * len(descriptions))
	idx = 0
	
	#transform label_distribution from dict to list of shape (4, 8)
	any_key = list(label_distribution.keys())[0]
	for label in label_distribution[any_key]:
		for align_tool in label_distribution:
			stats[idx] = '{:.1%}'.format(label_distribution[align_tool][label])
			idx += 1

	stats = np.ravel(stats).reshape(len(descriptions), len(aln_tool_list))

	fig, ax = plt.subplots(figsize=(15,10))
	ax.axis('off')
	the_table = ax.table(cellText=stats,
						 rowLabels=descriptions,
						 colLabels=aln_tool_list,
						 colWidths=[0.15]*4,
						 loc='center',
						 cellLoc='center'
						)
	the_table.scale(1, 4)
	the_table.set_fontsize(15)
	plt.subplots_adjust(left=0.15)

	title = 'Label Distribution Table'
	#plt.title(title, fontdict={'fontsize': 30})
	table_figures.append(fig)
	plt.savefig('{0}/imgs/{1}.png'.format(src_dir, title.replace(" ", "_")))
	plt.close()


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

