#!/usr/bin/env python
# iterate over all files in pos/ and neg/ folder and generate read-sized fragements

import glob, os

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print "This script requires BioPython to be installed!"

def chunks(l, n):
	for i in xrange(0, len(l), n):
		yield l[i:i+n]

def renameheader(prefix, path):
	newpath = path + ".corrected.fasta"
	with open(path) as original, open(newpath, 'w') as corrected:
		records = SeqIO.parse(path, 'fasta')
		i = 0
		for record in records:
			record.description = prefix + ";" + str(i) 
			SeqIO.write(record, corrected, 'fasta')
			i += 1

def split(length, path, export, chunkspath):
	print("create chunks of " + path)
	for filename in glob.iglob(path):
		(prefix, sep, suffix) = filename.rpartition('.')
		new_filename = prefix + '.frag.fasta'	
		#if __name__ == '__main__':
		# extract chunks
		handle = open(filename, 'r')
		records = list(SeqIO.parse(handle, "fasta"))
		record = records[0]
		records = []
		for (pos, chunk) in enumerate(chunks(str(record.seq), length)):
	    		records.append(SeqRecord(Seq(chunk, record.seq.alphabet), id=str(pos), name=record.name, description=record.description))    	
	 	SeqIO.write(records, open(new_filename, 'w'), "fasta")
	 	# create multiple sequence fasta file
	print("exporting " + export)
    	with open(export, 'w') as outfile:
        	for fname in glob.iglob(chunkspath):
        		with open(fname) as infile:
        			for line in infile:
        				outfile.write(line)
if __name__ == "__main__":
    sys.exit(main(sys.argv))
    