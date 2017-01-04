#!/usr/bin/python
import os, csv, sys
import download, simulate, features
import argparse

def loaddata(args):
	#email = "philipp.muench@helmholtz-hzi.de"
	if not os.path.exists('dat/chr'):
		download.downloadChr(args.email)
	else:
		print("WARNING: dat/chr already exists. Skipping this step.")
	if not os.path.exists('dat/pla'):
		download.downloadPla(args.email)
	else:
		print("WARNING: dat/pla already exists. Skipping this step.")

def getchunks():
	if not os.path.exists('dat/plasmid_chunks.fasta'):
		simulate.split(int(results.chunksize), 'dat/pla/*.fasta', 'dat/plasmid_chunks.fasta', 'dat/pla/*.frag.fasta')
		print('rewrite fasta headers')
		simulate.renameheader('positive','dat/plasmid_chunks.fasta')
	else:
		print("WARNING: dat/plasmid_chunks.fasta already exists. Skipping this step.")
	if not os.path.exists('dat/chromosome_chunks.fasta'):
		simulate.split(int(results.chunksize), 'dat/chr/*.fasta', 'dat/chromosome_chunks.fasta', 'dat/chr/*.frag.fasta')
		print('rewrite fasta headers')
		simulate.renameheader('negative','dat/chromosome_chunks.fasta')
	else:
		print("WARNING: dat/chromosome_chunks.fasta already exists. Skipping this step.")
	print('merge pos/negative set to train.txt')
	if not os.path.exists('dat/train.fasta'):
		filenames = ['dat/plasmid_chunks.fasta.corrected.fasta', 'dat/chromosome_chunks.fasta.corrected.fasta']
		with open('dat/train.fasta', 'w') as outfile:
			for fname in filenames:
				with open(fname) as infile:
					for line in infile:
						outfile.write(line)
	else:
		print("WARNING: dat/train.fasta already exists. Skipping this step.")

def getstatfeatures():
	if not os.path.exists('dat/train.features'):
	    print('export features')
	    os.system('python plasmidminer/features.py -I dat/train.fasta -s length -s na -s cpg > dat/train.features')
	    print ('export csv')
	    with open("dat/train.features", "r") as inp, open("dat/train.features.csv", "w") as out:
		w = csv.writer(out, delimiter=",")
		w.writerows(x for x in csv.reader(inp, delimiter="\t"))
		features.clearit("dat/train.features.csv", "dat/train.features.clear.csv")
		os.system("tail -n +2 dat/train.features.clear.csv > dat/train.features.clear2.csv")

def compress():
	if not os.path.exists('dat/train.fasta.gz'):
		print('compressing fasta file')
		os.system("gzip --keep dat/train.fasta")
	if not os.path.exists('dat/train.features.kmer'):
	    print('get kmer profile')
	    os.system('src/fasta2kmers2 -i dat/train.fasta -f dat/train.features.kmer -j 4 -k 5 -s 0')
	    with open("dat/train.features.kmer", "r") as inp, open("dat/train.features.kmer.csv", "w") as out:
		    w = csv.writer(out, delimiter=",")
		    w.writerows(x for x in csv.reader(inp, delimiter="\t"))

def getkmerprofile():
	# get feature matrix
	print('export features')
	if not os.path.exists('dat/train.features'):
	    os.system('python plasmidminer/features.py -I dat/train.fasta -s length -s na -s cpg > dat/train.features')
	    print ('export csv')
	    with open("dat/train.features", "r") as inp, open("dat/train.features.csv", "w") as out:
		w = csv.writer(out, delimiter=",")
		w.writerows(x for x in csv.reader(inp, delimiter="\t"))
		features.clearit("dat/train.features.csv", "dat/train.features.clear.csv")
		os.system("tail -n +2 dat/train.features.clear.csv > dat/train.features.clear2.csv")
	else:
		print("WARNING: dat/train.features already exists. Skipping.")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-c', action='store', dest='chunksize', help='Chunk size in nt')
 	parser.add_argument("--email", help='email adress needed for ncbi file download', dest="email")
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	args = parser.parse_args()
	print 'chunk size     =', args.chunksize
	print 'email     =', args.email

	loaddata(args)
	getchunks()
	getstatfeatures()
	compress()
	getkmerprofile()

