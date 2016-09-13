# iterate over all files in pos/ and neg/ folder and generate read-sized fragements

import glob, os

try:
    from progress.bar import Bar
except ImportError:
    print "This script requires progess to be installed!"
try:
    from Bio import SeqIO
except ImportError:
    print "This script requires BioPython to be installed!"
try:
    from Bio.Seq import Seq
except ImportError:
    print "This script requires BioPython to be installed!"
try:
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print "This script requires BioPython to be installed!"

def chunks(l, n):
	for i in xrange(0, len(l), n):
		yield l[i:i+n]

def split(length, path, export):
	print("create chunks")
	for filename in glob.iglob(path):
		(prefix, sep, suffix) = filename.rpartition('.')
		new_filename = prefix + '.frag.fasta'	
		if __name__ == '__main__':
			# extract chunks
			handle = open(filename, 'r')
			records = list(SeqIO.parse(handle, "fasta"))
			record = records[0]
			records = []
			for (pos, chunk) in enumerate(chunks(str(record.seq), splitlength)):
	    			records.append(SeqRecord(Seq(chunk, record.seq.alphabet), id=str(pos), name=record.name, description=record.description))    	
	 		SeqIO.write(records, open(new_filename, 'w'), "fasta")
	 		# create multiple sequence fasta file
	print("exporting " + export)
    	with open(export, 'w') as outfile:
        	for fname in glob.iglob(path):
        		with open(fname) as infile:
        			for line in infile:
        				outfile.write(line)

pla_path = 'pla/*.fasta'
pla_path_export = "plasmid_200.fasta"
chr_path = 'chr/*.fasta'
chr_path_export = "chromosomes_200.fasta"

splitlength = 200 # set length of outpuf fragments
split(splitlength, pla_path, pla_path_export)
split(splitlength, chr_path, chr_path_export)