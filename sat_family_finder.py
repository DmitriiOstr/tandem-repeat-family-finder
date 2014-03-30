#!/usr/bin/python
import subprocess
import sys
import argparse


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-i',action='store')
parser.add_argument('-o',action='store')
parser.add_argument('-dir',action='store')
args = parser.parse_args()

pipe = subprocess.Popen(["python","SAT_fam_ident.py","-i",args.i,"-o",args.o])
pipe.wait()

fasta_all_TR = args.o + '_all_TR.fasta'
pipe = subprocess.Popen(["formatdb","-i",fasta_all_TR,"-t",args.o,"-p","F"])
pipe.wait()

pipe = subprocess.Popen(["mkdir",args.dir])
pipe.wait()

pipe = subprocess.Popen(["python","blast_TR.py","-o",args.o,"-dir",args.dir])
pipe.wait()