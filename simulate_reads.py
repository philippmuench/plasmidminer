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

def split(length, path):
	files =  glob.iglob(path)
	#bar = Bar('split plasmid sequences', max=len(files))
	for filename in files:
		print(filename)
		(prefix, sep, suffix) = filename.rpartition('.')
		new_filename = prefix + '.frag.fasta'	
	if __name__ == '__main__':
		handle = open(filename, 'r')
		records = list(SeqIO.parse(handle, "fasta"))
		record = records[0]
		records = []
		for (pos, chunk) in enumerate(chunks(str(record.seq), splitlength)):
	    		records.append(SeqRecord(Seq(chunk, record.seq.alphabet), id=str(pos), name=record.name, description=record.description))    	
	 	SeqIO.write(records, open(new_filename, 'w'), "fasta")

pla_path = 'pla/*.fasta'
chr_path = 'chr/*.fasta'

splitlength = 200 # set length of outpuf fragments
split(splitlength, pla_path)
split(splitlength, chr_path)

#print("merge files")
#os.system('cat pla/*.frag.fasta > plasmids_200.fasta')
#os.system('cat chr/*.frag.fasta > chromosomes_200.fasta')

