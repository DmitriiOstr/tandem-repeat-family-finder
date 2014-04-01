import os
import sys
import csv
import re
from my_classes import Write_CSV, Read_CSV
import argparse


def trfdat_to_csv(input_file, output_file):
	re_contscaf = re.compile('(contig[0-9]{1,10}|scaffold\w\
		|ctg|Con\w|Contig\w|ctg\w|chromosome\w[0-9][0-9]|chromosome|LGE)')
	re_str = re.compile('\d')
	output = open(output_file, 'w')
	with open(input_file, 'r') as csvfile:
		trf_dat = csv.reader(csvfile, delimiter=' ')
		for row in trf_dat:
			strrow = str(row[0:1])
			if strrow[2:3] == '':
				pass
			elif strrow[2:3] == 'S':
				row_trf_dat = list(row)
				for i in range(len(row_trf_dat)):
					if row_trf_dat[i].find('gi|') >= 0:
						genbank_id = row_trf_dat[i]
						break
				for i in range(len(row_trf_dat)):
					if re_contscaf.search(row_trf_dat[i]) is not None:
						contig_name = re_contscaf.search(row_trf_dat[i])
						break
			elif re_str.search(strrow[2:3]) is not None:
				raw = '\t'.join(row)
				output.write('%s\t%s\t%s\n' % (raw, genbank_id, contig_name.group()))
 

def sorting_csv(infile, outfile):

	'''
	sort by TR monomer lenght
	'''

	Rd = Read_CSV(infile, '\t')
	table = Rd.read_csv()
	table = sorted(table, key=lambda a: int(a[2]), reverse=True)
	Wr = Write_CSV(table, outfile, '\t')
	Wr.write_csv()


def cycl_permutation(infile, outfile):

	'''
	Find and replace cyclic permutation for TR with monomer longer than 4 bp
	'''	

	def cycl_repl(seq):
		seq_repl = []
		for i in range(len(seq)):
			seq_new = seq + seq[0]
			seq_new = seq_new[1:len(seq_new)]
			seq_new_repl = seq_new.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]
			seq_repl.append(seq_new)
			seq_repl.append(seq_new_repl)
			i += 1
			seq = seq_new
		return seq_repl

	Rd = Read_CSV(infile, '\t')
	table = Rd.read_csv()
	for r in range(len(table)):
		try:
			if int(table[r][2]) > 4:
				for t in range(r+1, len(table)):
					if int(table[t][2]) < int(table[r][2]):
						break
					else:
						if len(str(table[r][13])) == len(str(table[t][13])):
							if 	(table[r][13].count('G') == (table[t][13].count('G')
								and table[r][13].count('C') == table[t][13].count('C')
								and table[r][13].count('A') == table[t][13].count('A')
								and table[r][13].count('T') == table[t][13].count('T'))
								or (table[r][13].count('G') == table[t][13].count('C')
								and table[r][13].count('C') == table[t][13].count('G')
								and table[r][13].count('T') == table[t][13].count('A')
								and table[r][13].count('A') == table[t][13].count('T'))):
								cycl_repl_seq = cycl_repl(table[r][13])
								for i in range(len(cycl_repl_seq)):
									if cycl_repl_seq[i] == table[t][13]:
										table[t][14] = table[t][14].replace(table[t][13], cycl_repl_seq[-2])
										table[t][13] = str(cycl_repl_seq[-2])												
		except(IndexError):
			pass
	Wr = Write_CSV(table, outfile, '\t')
	Wr.write_csv()
	

def trf_filter(infile, outfileSAT, outfileSMALL, outfileSSR):

	'''
	TR filter. Big TR: GC 20-80%, pole enthropy > 1.76, pole length > 1000
	create short SAT (pole lenght < 1000)
	and file with SSR repete
	'''

	Rd = Read_CSV(infile ,'\t')
	table = Rd.read_csv()
	table_SAT = []
	table_SSR = []
	table_SMALL = []
	for r in range(len(table)):
		try:
			gc_cont = float(table[r][9]) + float(table[r][11])
			pole_lenght = float(table[r][2]) * float(table[r][3])
			if (gc_cont in range (20,81)) and float(table[r][12]) > 1.76 and float(table[r][2]) > 4 and float(table[r][3]) > 4:
				if pole_lenght >= 1000:
					table_SAT.append(table[r])
				else:
					table_SMALL.append(table[r])
			elif float(table[r][12]) < 1.76 and gc_cont in range (20, 81):
				table_SSR.append(table[r])
		except:
			pass
	Wr = Write_CSV(table_SAT,outfileSAT, '\t')
	Wr.write_csv()
	Wr = Write_CSV(table_SMALL,outfileSMALL, '\t')
	Wr.write_csv()
	Wr = Write_CSV(table_SSR,outfileSSR, '\t')
	Wr.write_csv()

parent_parser = argparse.ArgumentParser(add_help = False)
child_parser1 = argparse.ArgumentParser(parents = [parent_parser])
child_parser1.add_argument('-i', action = 'store')
child_parser1.add_argument('-o', action = 'store')
args = child_parser1.parse_args()
inp_file = args.i
out_file = args.o
trfdat_to_csv(inp_file, out_file)
out_file_sort = out_file + '.sort'
sorting_csv(out_file, out_file_sort)
out_file_perm = out_file_sort + '.perm'
cycl_permutation(out_file_sort, out_file_perm)
outfileSAT = out_file + '.SAT'
outfileSSR = out_file + '.SSR'
outfileSMALL = out_file + 'SAT.short'
trf_filter(out_file_perm, outfileSAT, outfileSMALL, outfileSSR)

'''
set duplication items flag 1
'''

infile = outfileSAT
outfile = outfileSAT + '.clear' + '.csv'
Rd = Read_CSV(infile, '\t')
table = Rd.read_csv()
for i in range(len(table)):
	table[i].append(0)
for i in range(len(table)):
	if table[i][17] == 0:
		for j in range(len(table)):
			if i != j:
				if table[j][17] == 0 and table[i][16] == table[j][16]:
					if int(table[i][0]) <= int(table[j][0]) and int(table[i][1]) >= int(table[j][1]):
						if int(table[i][2]) > int(table[j][2]) and len(table[i][14]) <= len(table[j][14]):
							table[i][17] = 1
						elif int(table[i][2]) > int(table[j][2]) and len(table[j][14]) < len(table[i][14]):
							table[j][17]=1
						elif int(table[i][2]) < int(table[j][2]) and len(table[j][14]) < len(table[i][14]):
							table[j][17] = 1
						else:
							table[i][17] = 1
Wr = Write_CSV(table, outfile, '\t')
Wr.write_csv()

'''
Create fasta file with all TR pole
'''

Rd = Read_CSV(outfile, '\t')
table = Rd.read_csv()
outfile_fasta = args.o + '_all_TR.fasta'
outfile_fasta_all_TR = open(outfile_fasta, 'w')

for i in range(len(table)):
	try:
		if table[i][17] == '0':
			outfile_fasta_all_TR.write('>%s,%s,%s,%s\n%s\n' % (table[i][2], table[i][16], table[i][0], table[i][1], table[i][14]))
	except(IndexError):
		pass