#!/usr/bin/env python
# iterate over all files in pos/ and neg/ folder and generate read-sized fragements

import glob, os, random, sys
from termcolor import colored
import sys, os
from Bio import SeqIO

try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
except ImportError:
	print "This script requires BioPython to be installed!"

class Printer():
	"""Print things to stdout on one line dynamically"""
	def __init__(self,data):
		sys.stdout.write("\r\x1b[K"+data.__str__())
		sys.stdout.flush()

def chunks(l, n):
	for i in xrange(0, len(l), n):
		yield l[i:i+n]

def renameheader(prefix, path):
	newpath = path + ".corrected.fasta"
	with open(path) as original, open(newpath, 'w') as corrected:
		records = SeqIO.parse(path, 'fasta')
		i = 0
		for record in records:
			corrected.write(">" + prefix + "-" + str(i) + "\n" + str(record.seq) + "\n")
		#	record.description = prefix + ";" + str(i) 
		#	SeqIO.write(record, corrected, 'fasta')
			i += 1

def balancesize(a,b, num):
	inFiles =  [a,b]
	outFiles = []
	for i in range(len(inFiles)):
		for k in range(3):
			fNames = []
			fSeqs = []
			outFiles.append(file(str(inFiles[i])+'_'+'rand_'+str(num)+'-'+str(k+1)+'.fasta', 'wt'))
			# random subsampling without replacement
			for line in open(inFiles[i]):
				if (line[0] == '>'):
					fNames.append(line)
				else:
					fSeqs.append(line)
				curr = (len(outFiles)-1)
				for j in range(num):
					a = random.randint(0, (len(fNames)-1))
					outFiles[curr].write(fNames.pop(a))
					outFiles[curr].write(fSeqs.pops(a))

def split(length, path, export, chunkspath):
	Printer(colored('create chunks of ' + path, 'green'))
	for filename in glob.iglob(path):
		(prefix, sep, suffix) = filename.rpartition('.')
		new_filename = prefix + '.frag.fasta'	
		handle = open(filename, 'r')
		records = list(SeqIO.parse(handle, "fasta"))
		record = records[0]
		records = []
		for (pos, chunk) in enumerate(chunks(str(record.seq), length)):
				records.append(SeqRecord(Seq(chunk, record.seq.alphabet), id=str(pos), name=record.name, description=record.description))    	
		SeqIO.write(records, open(new_filename, 'w'), "fasta")
		# create multiple sequence fasta file
	Printer(colored('exporting: ' + export, 'green'))
	#sprint("exporting " + export)
	with open(export, 'w') as outfile:
		for fname in glob.iglob(chunkspath):
			with open(fname) as infile:
				for line in infile:
					outfile.write(line)
	Printer(colored('cleanup', 'green'))
	(prefix, sep, suffix) = path.rpartition('.')
	filelist = glob.glob(prefix + '.frag.fasta') 
	for f in filelist:
		os.remove(f)

def readsim(chunksize, simnum, input_file, output_file, method):
	Printer(colored('simulate reads of ' + input_file, 'green'))
	output_file_fq = output_file + '.fq'
	output_file2_q = output_file + 'fq.l2'
	s = " "
	if method == "wgsim":
		cmd = ("src/wgsim -d 100 -e 0 -r 0 -R 0 -X 0 -A 0 -N", str(simnum), "-1", str(chunksize) , str(input_file), str(output_file_fq), str(output_file2_q), '>/dev/null')
		os.system(s.join( cmd ))
		confertfastqtofasta(output_file_fq, output_file)
	if method == "art":
		cmd = ("src/art_illumina -ss HS25 -c", str(simnum), "-1", str(chunksize) ,"-i", str(input_file), "-o",str(output_file_fq), '>/dev/null')
		os.system(s.join( cmd ))
		confertfastqtofasta(output_file_fq, output_file)


def confertfastqtofasta(fq_path, fa_path):
	# conert fastq to fasta
	input_handle = open(fq_path, "rU")
	output_handle = open(fa_path, "w")
	sequences = SeqIO.parse(input_handle, "fastq")
	count = SeqIO.write(sequences, output_handle, "fasta")
	output_handle.close()
	input_handle.close()

if __name__ == "__main__":
	sys.exit(main(sys.argv))
	