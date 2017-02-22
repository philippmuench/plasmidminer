#!/usr/bin/python
import os, csv, sys
import download, simulate, features
import os.path
import argparse
from termcolor import colored
import pandas as pd
import glob
import ntpath
from Bio import SearchIO

try:
	from Bio import Entrez
	from Bio import SeqIO
except ImportError:
	print "This script requires BioPython to be installed!"

class Printer():
	"""Print things to stdout on one line dynamically"""
	def __init__(self,data):
		sys.stdout.write("\r\x1b[K"+data.__str__())
		sys.stdout.flush()

def loaddata(args):
	"""Download FASTA files from NCBI"""
	download.downloadChr(args)
	download.downloadPla(args)

def creatematrix(features, kmer, args):
	"""Generates feature matrix from raw data"""
	stat = pd.read_csv(features, sep=",")
	kmer = pd.read_csv(kmer, sep="\t", header=None)
	kmer = kmer.iloc[:, :-1]
	id2 = stat.id.str.split("-", expand=True)  # split the string to get label
	id2 = id2.iloc[:, :-1]
	stat2 = stat.iloc[:, 1:]
	df = pd.concat([stat2.reset_index(drop=True), kmer],
				   axis=1)  # concat kmerand stat matrix
	df = pd.concat([id2, df], axis=1)
	df.columns.values[0] = "label"
	# encoding class labels as integers
	df.loc[df.label == 'positive', 'label'] = 1
	df.loc[df.label == 'negative', 'label'] = 0
	# get number of instances per group
	y = df['label'].tolist()  # extract label
	X = df.drop(df.columns[[0]], 1)  # remove label
	return X, y

def getchunks(args):
	"""Split FASTA files into even-sized fragments and exort them as multiple FASTA file"""
	Printer(colored('(processing) ', 'green') + 'rewrite fasta headers (plasmids)')
	if args.readsim:
		if os.path.isfile('dat/plasmid.fasta'):
			simulate.readsim(int(args.chunksize), int(args.simnum), 'dat/plasmid.fasta', 'dat/plasmid_chunks.fasta')
		else:
			print("Error: Cannot find dat/plasmid.fasta. \n")
			sys.exit(1)
	else:
		simulate.split(int(args.chunksize), 'dat/pla/*.fasta', 'dat/plasmid_chunks.fasta', 'dat/pla/*.frag.fasta')

	if os.path.isfile('dat/plasmid_chunks.fasta'):
		simulate.renameheader('positive','dat/plasmid_chunks.fasta')
	else:
		print("Error: Cannot find dat/plasmid_chunks.fasta. \n")
		sys.exit(1)

	if args.readsim:
		simulate.readsim(int(args.chunksize), int(args.simnum), 'dat/chromosome.fasta', 'dat/chromosome_chunks.fasta')
	else:
		simulate.split(int(args.chunksize), 'dat/chr/*.fasta', 'dat/chromosome_chunks.fasta', 'dat/chr/*.frag.fasta')

	Printer(colored('(processing) ', 'green') + 'rewrite fasta headers (chromosomes)')
	simulate.renameheader('negative','dat/chromosome_chunks.fasta')

def createmerged():
	"""writes chr/plasmid sequences to one file"""
	filenames = ['dat/plasmid_chunks.fasta.corrected.fasta.clear', 'dat/chromosome_chunks.fasta.corrected.fasta.clear']
	with open('dat/train.fasta', 'w') as outfile:
		for fname in filenames:
			with open(fname) as infile:
				for line in infile:
					outfile.write(line)

def getstatfeatures():
	"""export statistical properties of fasta files as features"""
	Printer(colored('(processing) ', 'green') + 'generte feature matrix')
	os.system('python plasmidminer/features.py -I dat/train.fasta -s length -s na -s cpg > dat/train.features')
	Printer(colored('(processing) ', 'green') + 'export csv')
	with open("dat/train.features", "r") as inp, open("dat/train.features.csv", "w") as out:
		w = csv.writer(out, delimiter=",")
		w.writerows(x for x in csv.reader(inp, delimiter="\t"))
	features.clearit("dat/train.features.csv", "dat/train.features.clear.csv")
	os.system("tail -n +2 dat/train.features.clear.csv > dat/train.features.clear2.csv")

def compress():
	"""compress train.fasta files"""
	Printer(colored('(processing) ', 'green') + 'compress files')
	os.system("gzip --keep dat/train.fasta")

def extractkmers():
	Printer(colored('(processing) ', 'green') + 'extract kmer profile')
	os.system('src/fasta2kmers2 -i dat/train.fasta -f dat/train.features.kmer -j 4 -k 6 -s 0')
	with open("dat/train.features.kmer", "r") as inp, open("dat/train.features.kmer.csv", "w") as out:
		w = csv.writer(out, delimiter=",")
		w.writerows(x for x in csv.reader(inp, delimiter="\t"))

def sequence_cleaner(fasta_file, args):
	"""removes fragments smaller than chunksize from fasta file"""
	sequences={}
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		sequence = str(seq_record.seq).upper()
		if (len(sequence) >= int(args.chunksize)):
			#if sequence not in sequences:
			sequences[sequence] = seq_record.id
			#else:
			#	sequences[sequence] += "_" + seq_record.id
	output_file = open(fasta_file + ".clear", "w+")
	for sequence in sequences:
		output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")
	output_file.close()

def savepickl(args):
	"""saves dataset as pickl object"""
	Printer(colored('(processing) ', 'green') + 'save dataset as pickl object')
	X, y = creatematrix('dat/train.features.clear2.csv','dat/train.features.kmer', args)
	f = open(args.save, "w")
	pickle.dump(X, f)
	pickle.dump(y, f)
	f.close()

def savemsg(args):
	"""saves dataset as msgpack object"""
	Printer(colored('(processing) ', 'green') + 'save dataset as msg object')
	X, y = creatematrix('dat/train.features.clear2.csv','dat/train.features.kmer', args)
	y = pd.DataFrame(y)
	X.to_msgpack(args.save)
	y.to_msgpack(args.save, append=True)

def runHmmer(args, list_path, file_path, f):
	"""run prodigal and hmmsearch on chr files"""
	if not os.path.exists('dat/tmp'):
		os.makedirs('dat/tmp')
	# get the sample group
	head, group = os.path.split(os.path.split(file_path)[0])
	basename = os.path.splitext(str(ntpath.basename(str(file_path))))[0]
	exportpath = 'dat/tmp/' + ntpath.basename(str(file_path))
	hmmpath = 'dat/tmp/' + ntpath.basename(str(file_path)) + '.out'
	print('Processing %s of group %s' % (basename, group))
	s = " "
	cmd = ("prodigal -p meta -i",  str(file_path), "-a", exportpath, '-d /dev/null > /dev/null 2> /dev/null')
	os.system(s.join( cmd ))
	# run hmmsearch on faa ORF files
	s = " "
	cmd = ("hmmsearch -E 0.001 --domtblout", hmmpath, 'resources/remove.hmm', exportpath, '> /dev/null 2> /dev/null')
	os.system(s.join( cmd ))
	# write it to output file if there is a hit
	with open(hmmpath, 'rU') as input:
		try:
			for qresult in SearchIO.parse(input, 'hmmscan3-domtab'):
				query_id = qresult.id
				hits = qresult.hits
				num_hits = len(hits)
				acc = qresult.accession
				if num_hits > 0:
					f.write(''.join((basename, '\t', str(file_path),'\n')))
		except ValueError:
			print('parsing error on %s' % basename)

def screenhmm(args):
	"""walkes through the files we want to screen for hmm model"""
	remove_list = 'dat/chr_to_remove.txt'
	Printer(colored('(processing) ', 'green') + 'screen downloaded chr for plasmid domains')
	fastalist = filter(os.path.isfile, glob.glob('dat/chr/*.fasta'))
	if not fastalist:
		sys.stderr.write("No fasta files found.\n")
		exit(1)
	else:
		f = open(remove_list,'w')
		for filename in fastalist:
			runHmmer(args, remove_list, filename, f)
		f.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--save', action='store', dest='save', help='Save dataset as msg pack object', default='dat/dataset.msg')
	parser.add_argument('-t', '--taxa', action='store', dest='taxa', help='Taxonomic name for downloaded samples', default='Escherichia coli')
	parser.add_argument('-a', '--planum', action='store', dest='planum', help='Number of plasmids to be downloaded', default=100)
	parser.add_argument('-b', '--chrnum', action='store', dest='chrnum', help='Number of chromosomes to be downloaded', default=20)
	parser.add_argument('-c', '--chunksize', action='store', dest='chunksize', help='Chunk size in nt', default=200)
	parser.add_argument('-s', '--readsim', dest='readsim', action='store_true', help='use read simluation based on wgsim instead of sliding window')
	parser.add_argument('-N', '--simnum', action='store', dest='simnum', help='number of reads to simulate with wgsim', default=1000)
	parser.add_argument('-e', "--email", help='Email adress needed for ncbi file download', dest="email", default="pmu15@helmholtz-hzi.de")
	parser.add_argument('--multi', dest='gpf', action='store_true', help='Prepare data in multi feature format (taxonomic column in train/test samples)')
	parser.add_argument('--no_download', dest='no_download', action='store_true', help='use dataset stored in dat/')
	parser.add_argument('--screen', dest='screen', action='store_true', help='Screen downloaded chr for plasmid specific domains and remove them')
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	args = parser.parse_args()

	# do not download data
	if not args.no_download:
		Entrez.email = args.email # setting Entrez email
		loaddata(args)
	# screen data
	if args.screen:
		screenhmm(args)
	# split data into small read sized subset using sliding window approach
	getchunks(args)
	# remove sequences that are too short
	sequence_cleaner('dat/chromosome_chunks.fasta.corrected.fasta', args)
	sequence_cleaner('dat/plasmid_chunks.fasta.corrected.fasta', args)
	# merge genomes to one file
	createmerged()
	# export features such as GC content
	getstatfeatures()
	# compress fasta file
	compress()
	# extract kmer content
	extractkmers()
	# save feature data as msgpack object
	savemsg(args)
