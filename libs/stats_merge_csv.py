import sys, csv, re

#subset output_file, sizedir

if len(sys.argv) < 5:
	print("Usage: stats_merge_csv.py subset output_path size_record mode")
	sys.exit(1)

QUASTDIR = '/home/luke831215/biocollab/leon/pipeline/quast/out/'
stats = {}
stats_name = [
		"Subset Size",
		"# Scaffolds", "Scaffold NG50 (Kbp)", "# c. Scaffolds", "c. Scaffold Assembly Size (Mbp)",
		"Max c. Scaffold (Kbp)", "Scaffold c. NG25 (Kbp)", "Scaffold c. NG50 (Kbp)", "Scaffold c. NG75 (Kbp)", 
		"Scaffold LG80", "Scaffold LG90", "Scaffold LG99", "# N's per 100 kbp",
		"# Contigs", "Contig NG50 (Kbp)", "# c. Contigs", "Assembly Size (Mbp)", 
		"Max c. Contig (Kbp)", "Contig c. NG25 (Kbp)", "Contig c. NG50 (Kbp)", "Contig c. NG75 (Kbp)", "Contig LG80", "Contig LG90", "Contig LG99",
		"GC (%)", "Genome Fraction (%)", "Indels >= 5", "Inversions", "Relocation", "Translocation"
]
cor_stats_name = [
		"Subset Size",
		"Contigs #", "Not corrected N50", "Corrected contig #", "Corrected assembly size",
		"Max correct contig", "Corrected N25", "Corrected N50", "Corrected N75", "LG80", "LG90", "LG99", "# N's per 100 kbp",
		"Contigs #", "Not corrected N50", "Corrected contig #", "Corrected assembly size",
		"Max correct contig", "Corrected N25", "Corrected N50", "Corrected N75", "LG80", "LG90", "LG99",
		"GC (%)", "Genome fraction (%)", "Indels >= 5", "Inversions", "Relocation", "Translocation"
]

for file in ['gage_report.txt', 'report.txt']:
	with open(QUASTDIR+sys.argv[1]+'_scf/'+file, 'r') as r1:
		next(r1)
		next(r1)
		next(r1)
		for line in r1:
			try:
				[name, stat] = re.split(r'\s{2,}', line.strip())[:2]
				stat = re.search(r"^\d+[.\d+]*", stat.split()[0]).group(0)
				if '.' not in stat:
					stat = '{:,}'.format(int(stat))
			except:
				print(line.strip())
				continue
				#sys.exit(1)
			stats[name] = stat

if '%' not in sys.argv[3]:
	with open(sys.argv[3], 'r') as r2:
		stats['Subset Size'] = re.search(r"\(\d+.\d%", r2.readline().strip()).group(0)[1:]
else:
	stats['Subset Size'] = sys.argv[3]

rows = []
if sys.argv[4] == 'w':
	rows.append(["", sys.argv[1]])
	for i in range(len(stats_name)):
		if cor_stats_name[i] not in stats:
			print(cor_stats_name[i])
			if cor_stats_name[i] == 'LG99':
				rows.append(['L99', stats['L99']])
		else:
			rows.append([stats_name[i], stats[cor_stats_name[i]]])

else:
	with open(sys.argv[2], 'r') as r3:
		reader = csv.reader(r3)
		row = next(reader)
		row.append(sys.argv[1])
		rows.append(row)
		i = 0
		while i < len(stats_name):
			if cor_stats_name[i] not in stats:
				print(cor_stats_name[i])
				if cor_stats_name[i] == 'LG99':
					row = next(reader)		
					row.append(stats['L99'])
					rows.append(row)
			else:
				row = next(reader)
				while row[0] != stats_name[i]:
					if stats_name.index(row[0]) > i:
						col = len(row) - 1
						rows.append([stats_name[i]]+[""]*col+[stats[cor_stats_name[i]]])
						i+=1
					else:
						rows.append(row+[""])
						row = next(reader)

				row.append(stats[cor_stats_name[i]])
				rows.append(row)
			i+=1

with open(sys.argv[2], 'w') as w:
	writer = csv.writer(w)
	writer.writerows(rows)




