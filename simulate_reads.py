# iterate over all files in pos/ and neg/ folder and generate read-sized fragements

import glob os

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

splitlength = 200 # set length of outpuf fragments

# process pos dataset
for filename in glob.iglob('pla/*.fasta'):
	print filename
	(prefix, sep, suffix) = filename.rpartition('.')
 	new_filename = prefix + '.frag.fasta'
	try:
		from Bio import SeqIO
	except ImportError:
		print "This script requires BioPython to be installed!"

	def chunks(l, n):
	 	"""Yield successive n-sized chunks from l."""
		for i in xrange(0, len(l), n):
			yield l[i:i+n]

	if __name__ == '__main__':
		handle = open(filename, 'r')
		records = list(SeqIO.parse(handle, "fasta"))
		record = records[0]
		records = []
		for (pos, chunk) in enumerate(chunks(record.seq.tostring(), splitlength)):
	    		records.append(SeqRecord(Seq(chunk, record.seq.alphabet), id=str(pos), name=record.name, description=record.description))
	 	SeqIO.write(records, open(new_filename, 'w'), "fasta")

# process neg dataset
for filename in glob.iglob('chr/*.fasta'):
	print filename
	(prefix, sep, suffix) = filename.rpartition('.')
 	new_filename = prefix + '.frag.fasta'
	try:
		from Bio import SeqIO
	except ImportError:
		print "This script requires BioPython to be installed!"

	def chunks(l, n):
	 	"""Yield successive n-sized chunks from l."""
		for i in xrange(0, len(l), n):
			yield l[i:i+n]

	if __name__ == '__main__':
		handle = open(filename, 'r')
		records = list(SeqIO.parse(handle, "fasta"))
		record = records[0]
		records = []
		for (pos, chunk) in enumerate(chunks(record.seq.tostring(), splitlength)):
	    		records.append(SeqRecord(Seq(chunk, record.seq.alphabet), id=str(pos), name=record.name, description=record.description))
	 	SeqIO.write(records, open(new_filename, 'w'), "fasta")

# merge to single file
os.system('cat pla/*.frag.fasta > plasmids_200.fasta')
os.system('cat chr/*.frag.fasta > chromosomes_200.fasta')

