#!/usr/bin/python

import subprocess
import sys
import os
import argparse
from my_classes import Read_CSV, Write_CSV

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-r',action='store')
parser.add_argument('-o',action='store')
parser.add_argument('-dir',action='store')
args = parser.parse_args()

out=open('logname.txt','w')
os.chdir (args.dir)
csv_name = os.listdir()
for item in csv_name:
	infile = item
	Rd = Read_CSV(infile,'\t')
	table = Rd.read_csv()
	outfile = item + '.fasta'
	outfile = open(outfile,'w')
	for i in range(len(table)):
		outfile.write('>%s%s%s\n%s\n'%(table[i][2],table[i][0],table[i][1],table[i][14]))
	outfile.close()

subprocess.call(['mkdir','bowtie'])
subprocess.call(['/bin/sh','-c','mv *.fasta bowtie/'])

os.chdir ('bowtie/')
name = os.listdir()
reads = args.r
for i in range(len(name)):
	out.write(name[i]+'\n')
	pipe = subprocess.Popen(["bowtie2-build","-f",name[i],name[i]],stdout=subprocess.DEVNULL)
	pipe.wait()
	pipe = subprocess.Popen(["bowtie2","-x",name[i],"-U",reads,"-p","3","--local",">","/dev/null"])
	pipe.wait()
out.close()

os.chdir('../')
subprocess.call(['/bin/sh','-c','rm -r bowtie/'])

os.chdir('../')
inp_f = open('log.txt','r')
lst1 = [line.strip().split() for line in inp_f]
inp_f2 = open('logname.txt','r')
lst2 = [line.strip() for line in inp_f2]
dct={}
for i in range(len(lst2)):
	a = float(lst1[i*6+1][0])-float(lst1[i*6+2][0])
	b = float(lst1[i*6+1][0])
	dct[lst2[i]] = '%1.5f'%(100*a/b)
lst3 = [(k,v) for k,v in dct.items()]
lst3 = sorted(lst3, key=lambda a: a[0], reverse=False)
out_f=open('TR_fam_genome_content.csv','w')
writer = csv.writer(out_f,delimiter='\t')
for row in lst3:
	writer.writerow(row)