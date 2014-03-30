import subprocess
import os
import argparse
import csv
from my_classes import Read_CSV, Write_CSV, Graph
import time

if __name__ == "__main__":
	parser = argparse.ArgumentParser(add_help = False)
	parser.add_argument('-o',action = 'store')
	parser.add_argument('-dir',action ='store')
	args = parser.parse_args()
	infile = args.o + '.SAT' + '.clear' + '.csv'
	Rd = Read_CSV(infile,'\t')
	table = Rd.read_csv()
	for i in range(len(table)):
		if table[i][17]=='0':
			file_name = args.dir + '/' + table[i][2] + str(i) +'.fasta'
			outfile = open(file_name,'w')
			outfile.write('>'+table[i][2] + ',' + table[i][16] + ',' + table[i][0] + ',' + table[i][1] + '\n' + table[i][14])
			outfile.close()
	filename = os.listdir(args.dir)
	for j in range(len(filename)):
		inp_f=args.dir + '/' + filename[j]
		a = open(inp_f)
		b = a.readline().split(',')[0][1:]
		db_name = args.o+'_all_TR.fasta'
		name = args.dir+'/'+str(b)+'_'+str(j)
		pipe = subprocess.Popen(["blastn","-task","blastn","-query",inp_f,'-db',db_name,'-out',name,'-outfmt','10', '-evalue','10e-16', '-dust','no'],bufsize=0)
		pipe.wait()
	os.chdir(args.dir)
	os.system("rm *.fasta")

	filename=os.listdir()
	def blast_fam_clear(infile):
		Rd = Read_CSV(infile,',')
		table = Rd.read_csv()
		a=True
		table2=[]
		for r in range(len(table)):
			for t in range(len(table2)):
				if table2[t][4]==table[r][4] and table2[t][5]==table[r][5] and table2[t][6]==table[r][6]: #сравнение на длину мономера, номер скаффолда/контига/хромосомы, начало поля
					a=False
					break
			if a==True:
				table2.append(table[r])
			a=True
		Wr = Write_CSV(table2,infile,',')
		Wr.write_csv()
	for i in range(len(filename)):
		infile=str(filename[i])
		blast_fam_clear(infile)
	table2=[]
	for files in filename:
			Rd = Read_CSV(files,',')
			table = Rd.read_csv()
			for i in range(len(table)):
				table2.append(table[i])
			subprocess.call(['rm',files])
	Wr = Write_CSV(table2,'TR_blast_all.txt',',')
	Wr.write_csv()

	inp_f = 'TR_blast_all.txt'
	table = Read_CSV(inp_f,',')
	table = table.read_csv()
	out_f = open('temp.txt','w')
	dct1 = {}
	dct1[(table[0][0],table[0][1],table[0][2],table[0][3])]=1
	a = 1
	for i in range(len(table)):
		if (table[i][0],table[i][1],table[i][2],table[i][3]) in dct1.keys():
			out_f.write('%s '%dct1[(table[i][0],table[i][1],table[i][2],table[i][3])])
		if (table[i][0],table[i][1],table[i][2],table[i][3]) not in dct1.keys():
			a += 1
			dct1[(table[i][0],table[i][1],table[i][2],table[i][3])] = a
			out_f.write('%s '%dct1[(table[i][0],table[i][1],table[i][2],table[i][3])])
		if (table[i][4],table[i][5],table[i][6],table[i][7]) in dct1.keys():
			out_f.write('%s\n'%dct1[(table[i][4],table[i][5],table[i][6],table[i][7])])
		if (table[i][4],table[i][5],table[i][6],table[i][7]) not in dct1.keys():
			a += 1
			dct1[(table[i][4],table[i][5],table[i][6],table[i][7])] = a
			out_f.write('%s\n'%dct1[(table[i][4],table[i][5],table[i][6],table[i][7])]) 
	out_f.close()
	out_f = 'temp.txt'
	G = Graph()
	with open(out_f) as inp_file:
		for lines in inp_file:
			line = lines.strip().split()
			G.add_edge(int(line[0]),int(line[1]))
	inp_file.close()
	subprocess.call(['rm','temp.txt'])
	n = a
	for i in range(1,n+1):
		G.add_node(i)
	lst_visited = []		#list of visited node
	lst_connections = []		#list of connected node
	dct_connections_comp = {}		# Num of conn.comp:[list of node]
	def sv(p,lst_visited,lst_connections,dct_connections_comp,j):		# Depth-first search
		lst_visited.append(p)
		lst_connections.append(p)
		neib_node = G.neighbor[p]
		for node in neib_node:
			if node not in lst_visited:
				sv(node,lst_visited,lst_connections,dct_connections_comp,j)
		dct_connections_comp[j] = lst_connections
		return lst_visited,lst_connections,dct_connections_comp
	j = 1
	for p in range(1,n+1):		# Depth-first search for all node
		if p not in lst_visited and p in list(G.neighbor):
			lst_connections = []
			sv(p,lst_visited,lst_connections,dct_connections_comp,j)
			j += 1
		elif p not in lst_visited and p not in list(G.neighbor):
			dct_connections_comp[j] = [p]
			j += 1
	node_dict = {}		# Num of conn.comp:[list of node] to node: Num of conn.comp
	for key in dct_connections_comp:
		node = dct_connections_comp[key]
		if type(node) == list:
			for i in range(len(node)):
				node_dict[node[i]] = key
		else:	
			node_dict[node]=key
	for item in dct1:
		dct1[item]=node_dict[dct1[item]]
	os.chdir('../')
	infile = args.o + '.SAT' + '.clear' + '.csv'
	table = Read_CSV(infile,'\t')
	table = table.read_csv()
	for i in range(len(table)):
		table[i].append(' ')
		if table[i][17] == '0':
			table[i][18]=dct1[(table[i][2]),table[i][16],table[i][0],table[i][1]]
	table = sorted(table, key=lambda a: str(a[18]), reverse=False)
	Wr = Write_CSV(table,infile,'\t')
	Wr.write_csv()

	table_out = []
	for i in range(len(table)-1):
		if table[i][18] == table[i+1][18] and table[i][18]!=' ':
			table_out.append(table[i])
		elif table[i][18] != table[i+1][18] and table[i][18] != ' ':
			table_out.append(table[i])
			outfile=args.dir + '/' + str(table[i][18]) + '.csv'
			Wr = Write_CSV(table_out,outfile,'\t')
			Wr.write_csv()
			table_out = []