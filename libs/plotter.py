import os
import re
import datetime
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import gridspec
import glob
import base64
from matplotlib.backends.backend_pdf import PdfPages
from bs4 import BeautifulSoup
import csv

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
labels = ['P', 'S', 'C', 'O', 'M', 'F', 'N']


def save_to_csv(stats, src_dir, aln_tool_list):
	dirname='{}/label_dis'.format(src_dir)
	if not os.path.isdir(dirname):
		os.makedirs(dirname)
	with open('{}/label_dis/label_dis.csv'.format(src_dir), 'w') as w:
		cnt = 0
		writer = csv.writer(w, delimiter=',')
		writer.writerow(['']+aln_tool_list)
		for label in labels:
			writer.writerow([label] + list(stats[cnt]))
			cnt += 1


def fill_in_label(stats, soup, template_fpath):
	label_list = ['P', 'S', 'C', 'O', 'M', 'F', 'N']
	icon_dirpath = "{}/icons".format(os.path.dirname(template_fpath))
	label_div = soup.findAll('div', {"class": 'label'})

	for i in range(len(label_div)):
		#remove original img tag
		tmp = label_div[i].find('img')
		tmp.replace_with('')

		#generate new icon
		data_uri = base64.b64encode(open("{0}/{1}.png".format(icon_dirpath, label_list[i]), 'rb').read()).decode('utf-8').replace('\n', '')
		img_src = soup.new_tag("img", src="data:image/png;base64,{0}".format(data_uri), **{'class': 'icon'})
		label_div[i].append(img_src)

	label_dis_table = soup.find('table', class_="label-dis-table")
	for i in range(len(stats)):
		for j in range(len(stats[0])):
			cell = label_dis_table.find('td', id='{0}{1}'.format(labels[i], j+1))
			cell.string = stats[i][j]


def save_to_html(out_dir, template_fpath, data, thre, aln_tool_list, label_distribution, basic_stats, gen_stats):
	"""concat all figures into report.html by inserting figures into template.html"""
	toc_fpath = '{}/toc.html'.format(os.path.dirname(template_fpath))
	all_html_fpath = '{0}/{1}/post-assembly_report.html'.format(out_dir, data)
	index_fpath = '{0}/{1}.html'.format(out_dir, data)
	src_dir = os.path.dirname(all_html_fpath)
	icon_dirpath = "{}/icons".format(os.path.dirname(template_fpath))

	with open(template_fpath) as file:
		soup = BeautifulSoup(file, "lxml")
		time_sec = soup.find('p', id='time')
		time_sec.string = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')

		#fill in thre table
		names = ['Overall PM%', 'S ratio', 'C ratio', 'O ratio', 'N ratio']
		thre_table = soup.find('table', **{"class": 'threshold-table'})
		keys = list(thre.keys())
		for i in range(len(names)):
			new_tr= soup.new_tag("tr")
			td1, td2 = soup.new_tag('td'), soup.new_tag('td')
			td1.string = names[i]
			td2.string = '{:.0%}'.format(thre[keys[i]])
			new_tr.append(td1)
			new_tr.append(td2)
			thre_table.append(new_tr)


		main = soup.find('div', **{"class": 'main'})
		
		data_sec = main.find('h2', **{"class": 'data'})
		data_sec.string = data

		#Overall review depends on avg. poor ratio
		poor_ratio = None
		for row in basic_stats:
			if row[0] == "Avg. poorly mapped sequence%":
				poor_ratio = float(row[1].strip('%'))/100
				break

		#poor_ratio = float(basic_stats[3][1].strip('%'))/100
		if poor_ratio < thre['PM']:
			span_class = ["result"]
			review = "Percentage of poorly-mapped reads: {:.1%}".format(poor_ratio)
			icon_fpath = '{}/check.png'.format(icon_dirpath)

		else:
			span_class = ["result", "far", "fa-times-circle"]
			review = "Percentage of poorly-mapped reads: {:.1%}".format(poor_ratio)
			icon_fpath = '{}/cross.png'.format(icon_dirpath)

		result = main.find('span', class_="result")
		data_uri = base64.b64encode(open(icon_fpath, 'rb').read()).decode('utf-8').replace('\n', '')
		img_src = soup.new_tag("img", src="data:image/png;base64,{0}".format(data_uri), **{'class': "title-icon"})
		result.insert(0, img_src)
		result.h3.string = review
		
		#section I
		#set up basic sequence stats table
		#seq information
		seq_info = soup.find('table', class_="seq-info")
		for i in range(len(basic_stats)):
			new_tr = soup.new_tag('tr')
			stat, value = soup.new_tag('td'), soup.new_tag('td')
			stat.string = basic_stats[i][0]
			value.string = basic_stats[i][1]
			new_tr.append(stat)
			new_tr.append(value)
			seq_info.append(new_tr)

		#genome evaluation table
		gen_eval_table = soup.find('table', class_="gen-eval-table")
		for i in range(len(gen_stats)):
			new_tr = soup.new_tag('tr')
			stat, value = soup.new_tag('td'), soup.new_tag('td')
			stat.string = gen_stats[i][0]
			value.string = gen_stats[i][1]
			new_tr.append(stat)
			new_tr.append(value)
			gen_eval_table.append(new_tr)

		#set up label dis. table
		fill_in_label(flatten(label_distribution, aln_tool_list), soup, template_fpath)

		#plot label dis piechart
		img_fpath = '{}/images/piechart.png'.format(src_dir)
		data_uri = base64.b64encode(open(img_fpath, 'rb').read()).decode('utf-8').replace('\n', '')
		label_dis_bar = soup.find("div", id="label-dis-piechart")
		img_src = soup.new_tag("img", src="data:image/png;base64,{0}".format(data_uri))
		label_dis_bar.append(img_src)

		#plot label dis bar
		img_fpath = '{}/images/label_dis_bar.png'.format(src_dir)
		data_uri = base64.b64encode(open(img_fpath, 'rb').read()).decode('utf-8').replace('\n', '')
		label_dis_bar = soup.find("div", id="label-dis-barchart")
		img_src = soup.new_tag("img", src="data:image/png;base64,{0}".format(data_uri))
		label_dis_bar.append(img_src)

		#section II
		#the rest of distribution graph for each aligner
		i = 0
		#title_list = ['BWA - MEM', 'Bowtie2 - local', 'BWA - backtrack', 'Bowtie2 - backtrack']
		title_list = ['BWA - MEM', 'BWA - backtrack']
		sec_two = soup.find('div', {"class": 'sec-2'})
		for aln_tool in aln_tool_list:
			title = soup.new_tag("h2", **{'class': 'plot-name'})
			title.string = title_list[i]
			sec_two.append(title)
			i += 1

			cnt = 1
			files = glob.glob('{0}/{1}/plot/*.png'.format(src_dir, aln_tool))
			files.sort(key=os.path.getmtime)
			for img_fpath in files:
				data_uri = base64.b64encode(open(img_fpath, 'rb').read()).decode('utf-8').replace('\n', '')
				new_div = soup.new_tag("div", id="{0}-{1}".format(aln_tool, cnt), **{'class': 'inner'})
				img_src = soup.new_tag("img", src="data:image/png;base64,{0}".format(data_uri))
				new_div.append(img_src)
				sec_two.append(new_div)
				cnt += 1

		#write report.html
		with open(all_html_fpath, 'w') as w:
			w.write(str(soup))

	#fill in the href in toc.html
	with open(toc_fpath) as file:
		soup = BeautifulSoup(file, "lxml")
		link = soup.find('a', id="post-assembly")
		link['href'] = '{}/post-assembly_report.html'.format(data)
		link = soup.find('a', id="pre-assembly")
		link['href'] = '{}/pre-assembly_report.htm'.format(data)

		with open(index_fpath, 'w') as w:
			w.write(str(soup))


def save_to_pdf(all_pdf_fpath, plot_figures):
	"""concat all figures into report.pdf"""

	all_pdf_file = PdfPages(all_pdf_fpath)
	for figure in plot_figures:
		all_pdf_file.savefig(figure)
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


def get_genome_eval_stat(src_dir, ref_name):
	if os.path.isfile('{}/quast/gage_report.txt'.format(src_dir)):
		GAGE = True
		report_list = ['gage_report.txt', 'report.txt']
	else:
		GAGE = False
		report_list = ['report.txt']

	stats = {}
	if GAGE:
		stats_name = [
				"# Contigs", "Assembly Size (Mbp)", "N50 (Kbp)",
				"Max Contigs (Kbp)", "NG25 (Kbp)", "NG50 (Kbp)", "NG75 (Kbp)", 
				"LG80", "LG90", "LG99", "# N's per 100 kbp",
				"GC (%)", "Genome Fraction (%)", "Indels >= 5", "Inversions", "Relocation", "Translocation"
		]
		cor_stats_name = [
				"Contigs #", "Assembly size", "N50",
				"Largest contig", "NG25", "NG50", "NG75", "LG80", "LG90", "LG99", "# N's per 100 kbp",
				"GC (%)", "Genome fraction (%)", "Indels >= 5", "Inversions", "Relocation", "Translocation"
		]
	else:
		stats_name = [
			"# Contigs", "Assembly Size (Mbp)", "Max Contigs (Kbp)", "N25 (Kbp)", "N50 (Kbp)", "N75 (Kbp)",
			"L80", "L90", "L99", "# N's per 100 kbp", "GC (%)"
		]
		cor_stats_name = [
			"# contigs", "Total length", "Largest contig", "N25", "N50", "N75", 
			"L80", "L90", "L99", "# N's per 100 kbp", "GC (%)"
		]

	#retrieve all stats in quast report
	for file in report_list:
		with open("{0}/quast/{1}".format(src_dir, file)) as infile:
			next(infile)
			next(infile)
			next(infile)
			for line in infile:
				try:
					[name, stat] = re.split(r'\s{2,}', line.strip())[:2]
					stat = re.search(r"^\d+[.\d+]*", stat.split()[0]).group(0)
				except:
					#print(line.strip())
					continue
					#sys.exit(1)
				stats[name] = stat

	#extract the wanted stats
	rows = [['File name', ref_name]]
	for i in range(len(stats_name)):
		if cor_stats_name[i] not in stats:
			print(cor_stats_name[i]+' not found')
			if cor_stats_name[i] == 'LG99':
				rows.append(['L99', stats['L99']])
		else:
			value = stats[cor_stats_name[i]]
			if '.' not in value:
				value = int(value)
				if 'Kbp' in stats_name[i]:
					value = value // 1000

				if 'Mbp' in stats_name[i]:
					value = value // 1000000
				value = '{:,}'.format(round(value, 1))
			rows.append([stats_name[i], value])

	return rows


def plot_sam_dis(src_dir, data, aln_tool, label_array, read_size, plot_figures, xlabel='', label='', thre=''):
	hist, bins = np.histogram(data, bins=20)

	fig = plt.figure(figsize=(15, 10))
	#gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1]) 

	ax = fig.add_subplot(111)
	ax.bar(bins[:-1], hist.astype(np.float32) / hist.sum() * 100, width=(bins[1]-bins[0]), alpha=0.5, color='steelblue', linewidth=0, align='edge')
	ax.set_xlabel(xlabel)
	ax.set_ylabel('Percentile')

	#add legend
	num_read = np.sum(label_array == label)
	#if 'Alignment Score' not in xlabel:
	if thre:
		ax.axvline(thre, color='red', zorder=20)
		ratio_good = np.sum(data < thre) / num_read if num_read else 0
		ratio_bad = 1 - ratio_good
		ax.plot(1, 1, label='Below threshold (eligible): {:.1%}'.format(ratio_good), marker='', ls='')
		ax.plot(1, 1, label='Above threshold (poor): {:.1%}'.format(ratio_bad), marker='', ls='')

		#zoom in the graph if necessary
		if np.max(data) < 0.5:
			ax.set_xlim(0, 0.5)
			ax.text(thre*2, 0.5, '{} ratio'.format(label), transform=ax.transAxes, zorder=21)
		else:
			ax.text(thre, 0.5, '{} ratio'.format(label), transform=ax.transAxes, zorder=21)
	
	label_dis = num_read / read_size
	ax.plot(1, 1, label='No. of {0} reads: {1} ({2:.1%})'.format(label, num_read, label_dis), marker='', ls='')
	ax.legend(loc=1, prop={'size': 15})

	title = '{0} distribution of {1} reads\n{2}'.format(xlabel, label, aln_tool)
	ax.set_title(title)

	plot_figures.append(fig)

	if xlabel == 'Mismatch%':
		fig_name = 'mismatch_ratio'
	elif xlabel == 'Clip%':
		fig_name = 'clip_ratio'
	else:
		fig_name = 'aln_score_{}'.format(label)
	#title = '{0} distribution of {1} reads {2}'.format(xlabel, label, aln_tool)
	plt.savefig('{0}/{1}/plot/{2}.png'.format(src_dir, aln_tool, fig_name))
	plt.close()


def draw_bar_bi(ax, field_list, read_size, idx, thre, label=''):
	upper_cnt, lower_cnt = 0, 0
	for field in field_list:
		if field < thre:
			upper_cnt += 1
		else:
			lower_cnt += 1
	
	#print(label, upper_cnt, lower_cnt)
	value1 = round(upper_cnt / read_size * 100)
	ax.bar(idx, value1, color=colors[idx], align='center', width=0.5)#, label=label)#label=r'${}_+$'.format(label))
	if upper_cnt:
		num = '~0' if value1 == 0 else value1
	else:
		num = 0
	ax.text(idx, value1+5, '{}'.format(num), ha='center')
	
	value2 = round(-1 * lower_cnt / read_size * 100)	
	ax.bar(idx, value2, color=colors[idx], align='center', width=0.5)#label=r'${}_-$'.format(label))
	if lower_cnt:
		num = '~0' if value2 == 0 else value2
	else:
		num = 0
	ax.text(idx, value2-5, '{}'.format(num), ha='center')

	return value1, value2, lower_cnt


def draw_bar(ax, label_array, read_size, idx, label=''):
	constant = -1 if label == 'F' else 1
	label_cnt = np.sum(label_array == label)
	value = round(int(constant * label_cnt / read_size * 100))
	ax.bar(idx, value, color=colors[idx], align='center', width=0.5)#, label=label)
	if label_cnt:
		num = '~0' if value == 0 else value
	else:
		num = 0
	ax.text(idx, value+5*constant, '{}'.format(num), ha='center')

	return value


def do_label_dis_bar(ax, align_array, aln_tool, label_array, thre):
	read_size = len(align_array)
	patch_handles = {}
	#x_labels = ['P', 'S', 'C', 'O', 'M', 'F', 'N'] if 'backtrack' not in aln_tool else ['P', 'S', 'O', 'M', 'F', 'N']
	x_labels = ['P', 'S', 'C', 'O', 'M', 'F', 'N']
	idx = 0
	below_cnt = 0
	ymax, ymin = 0, 0

	#P
	value = draw_bar(ax, label_array, read_size, idx, label='P')
	idx += 1
	ymax = value if ymax < value else ymax

	#S
	nm_list = do_nm(align_array, label_array)
	pos_value, neg_value, lower_cnt = draw_bar_bi(ax, nm_list, read_size, idx, thre['MR'], label='S')
	below_cnt += lower_cnt
	idx += 1
	ymax = pos_value if ymax < pos_value else ymax
	ymin = neg_value if ymin > pos_value else ymin
	
	#C
	#if 'backtrack' not in aln_tool:
	cr_list = do_cr(align_array, label_array)
	pos_value, neg_value, lower_cnt = draw_bar_bi(ax, cr_list, read_size, idx, thre['CR'], label='C')
	below_cnt += lower_cnt
	idx += 1
	ymax = pos_value if ymax < pos_value else ymax
	ymin = neg_value if ymin > pos_value else ymin

	#O
	o_list = do_others(align_array, label_array)
	pos_value, neg_value, lower_cnt = draw_bar_bi(ax, o_list, read_size, idx, thre['OR'], label='O')
	below_cnt += lower_cnt
	idx += 1
	ymax = pos_value if ymax < pos_value else ymax
	ymin = neg_value if ymin > pos_value else ymin

	#M
	value = draw_bar(ax, label_array, read_size, idx, label='M')
	idx += 1
	ymax = value if ymax < value else ymax

	#F
	value = draw_bar(ax, label_array, read_size, idx, label='F')
	below_cnt += np.sum(label_array == 'F')
	idx += 1
	ymin = value if ymin > value else ymin

	#N
	n_list = do_N(align_array, label_array)
	pos_value, neg_value, lower_cnt = draw_bar_bi(ax, n_list, read_size, idx, thre['NR'], label='N')
	below_cnt += lower_cnt
	idx += 1
	ymax = pos_value if ymax < pos_value else ymax
	ymin = neg_value if ymin > pos_value else ymin


	below_pct = below_cnt / read_size
	above_pct = 1 - below_pct

	ax.plot(1, 1, label='Above: {:.1%}'.format(above_pct), marker='', ls='')
	ax.plot(1, 1, label='Below (PM%): {:.1%}'.format(below_pct), marker='', ls='')
	ax.set_ylim(-100, 100)
	ax.set_xticks(np.arange(len(x_labels)))
	ax.set_xticklabels(x_labels)
	ax.set_ylabel('Percentage')
	ax.legend(loc=1, prop={'size': 10}, frameon=False)
	title = aln_tool + ' (local)' if aln_tool == 'bwa-mem' else aln_tool + ' (end2end)'
	ax.set_title('{}'.format(title))

	#record barplot input for future use
	return {'P': None, 'S': nm_list, 'C': cr_list, 'O': o_list, 'M': None, 'F': None, 'N': n_list}, below_pct, ymax, ymin


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
			try:
				stats[idx] = '{:.1%}'.format(label_distribution[aln_tool][label])
			except ValueError:
				stats[idx] = '{}'.format(label_distribution[aln_tool][label])
			idx += 1

	#7*4 array
	return np.ravel(stats).reshape(len(labels), len(aln_tool_list))


def do_label_dis_table(label_dis, src_dir, aln_tool_list, plot_figures):
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
			try:
				add_text(plt.subplot(gs[i, j]), '{:.1%}'.format(label_dis[aln_tool_list[j-3]][labels[i-1]]))
			except ValueError:
				add_text(plt.subplot(gs[i, j]), '{}'.format(label_dis[aln_tool_list[j-3]][labels[i-1]]))

	#title = 'Label Distribution Table'
	#plt.title(title, fontdict={'fontsize': 30})
	plot_figures.append(fig)
	#plt.savefig('{0}/label_dis/imgs/table.png'.format(src_dir))
	plt.close()

	#save_to_csv(flatten(label_dis, aln_tool_list), src_dir, aln_tool_list)


def do_label_piechart(label_dis, src_dir, aln_tool_list, plot_figures):
	labels = ['P', 'C', 'S', 'O', 'F', 'M', 'N']
	fracs = [None] * len(labels)

	fig = plt.figure(figsize=(15, 8))
	ax_list = [None] * len(aln_tool_list)
	for i, aln_tool in enumerate(aln_tool_list):
		for j, label in enumerate(labels):
			fracs[j] = label_dis[aln_tool][label] * 100

		ax = fig.add_subplot(len(aln_tool_list) / 2, 2, i+1)
		patches, texts, autotexts = ax.pie(fracs, labels=labels, autopct="%.1f%%", radius=0.8, pctdistance=1.25, labeldistance=1.05)
		patches[0].set_edgecolor('white')
		for text in texts:
			text.set_fontsize(15)
		title = aln_tool + ' (local)' if aln_tool == 'bwa-mem' else aln_tool + ' (end2end)'
		ax.set_title(title)

	fig.savefig('{}/images/piechart.png'.format(src_dir))
	plot_figures.append(fig)
	plt.close()


def do_nm(align_array, label_array):
	"""plot mismatch in unique reads with sub. error"""

	cnt = 0
	label_id_list = np.where(label_array == 'S')[0]
	nm_list = np.zeros(len(label_id_list))

	for _id in label_id_list:
		try:
			nm_list[cnt] = align_array[_id]['num_mismatch']/len(align_array[_id]['seq'])
		except:
			print(_id)
		cnt += 1


	return nm_list


def do_cr(align_array, label_array):
	"""clip rate, excluding bowtie2-backtrack (no clip labels)"""

	cnt = 0
	label_id_list = np.where(label_array == 'C')[0]
	cr_list = np.zeros(len(label_id_list))

	for _id in label_id_list:
		length, clip_size = 0, 0
		cigar = align_array[_id]['cigar']
		m = re.findall('\d+[MIDNSHP=X]', cigar)
		for value in m:
			length += int(value[:-1])
			if value[-1] == 'S':
				clip_size += int(value[:-1])

		cr_list[cnt] = clip_size / length
		cnt += 1
	
	return cr_list


def do_as(align_array, label_array, label):
	""""alignment score, excluding bwa-backtrack (no AS field) and bowtie2-backtrack for clipped reads"""

	label_id_list = np.where(label_array == label)[0]
	as_list = np.zeros(len(label_id_list))
	align_subset = align_array[label_id_list]
	cnt = 0

	for alignment in align_subset:
		as_list[cnt] = alignment['AS']
		cnt += 1

	return as_list


def do_others(align_array, label_array):
	cnt = 0
	label_id_list = np.where(label_array == 'O')[0]
	o_list = np.zeros(len(label_id_list))

	for _id in label_id_list:
		length, o_size = 0, 0
		cigar = align_array[_id]['cigar']
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
		seq = align_array[_id]['seq']
		if 'N' in seq:
			n_size = seq.count('N')
			n_list[cnt] = n_size / len(seq)
			cnt += 1

	return n_list
