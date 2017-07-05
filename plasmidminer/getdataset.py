#!/usr/bin/python
# loads pos/neg dataset from NCBI and simulates fragments and outputs the feature matrix
# philipp.muench@helmholtz-hzi.de
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
		if os.path.isfile(str(args.data) + '/plasmid.fasta'):
			simulate.readsim(int(args.chunksize), int(args.simnum), str(args.data) + '/plasmid.fasta', str(args.data) + '/plasmid_chunks.fasta')
		else:
			print("Error: Cannot find plasmid.fasta. \n")
			sys.exit(1)
	else:
		simulate.split(int(args.chunksize), str(args.data) + '/pla/*.fasta', str(args.data) + '/plasmid_chunks.fasta', str(args.data) + '/pla/*.frag.fasta')

	if os.path.isfile(str(args.data) + '/plasmid_chunks.fasta'):
		simulate.renameheader('positive',str(args.data) + '/plasmid_chunks.fasta')
	else:
		print("Error: Cannot find plasmid_chunks.fasta. \n")
		sys.exit(1)

	if args.readsim:
		simulate.readsim(int(args.chunksize), int(args.simnum), str(args.data) + '/chromosome.fasta', str(args.data) + '/chromosome_chunks.fasta')
	else:
		simulate.split(int(args.chunksize), str(args.data) + '/chr/*.fasta', str(args.data) + '/chromosome_chunks.fasta', str(args.data) + '/chr/*.frag.fasta')

	Printer(colored('(processing) ', 'green') + 'rewrite fasta headers (chromosomes)')
	simulate.renameheader('negative', str(args.data) + '/chromosome_chunks.fasta')

def createmerged(args):
	"""writes chr/plasmid sequences to one file"""
	filenames = [str(args.data) + '/plasmid_chunks.fasta.corrected.fasta.clear', str(args.data) + '/chromosome_chunks.fasta.corrected.fasta.clear']
	fname = str(args.data) + '/train.fasta'
	with open(fname, 'w') as outfile:
		for fname in filenames:
			with open(fname) as infile:
				for line in infile:
					outfile.write(line)

def getstatfeatures(args):
	"""export statistical properties of fasta files as features"""
	Printer(colored('(processing) ', 'green') + 'generate feature matrix')
	s = ""
	cmd = ('python plasmidminer/features.py -I ', str(args.data), '/train.fasta -s length -s na -s cpg > ', str(args.data), '/train.features')
	print s.join( cmd )
	os.system(s.join( cmd ))

	Printer(colored('(processing) ', 'green') + 'export csv')
	fname = str(args.data) + "/train.features"
	fname2 = str(args.data) + "/train.features.csv"
	with open(fname, "r") as inp, open(fname2, "w") as out:
		w = csv.writer(out, delimiter=",")
		w.writerows(x for x in csv.reader(inp, delimiter="\t"))
	features.clearit(str(args.data) + "/train.features.csv", str(args.data) + "/train.features.clear.csv")
	s = ""
	cmd = ('tail -n +2 ',str(args.data), '/train.features.clear.csv > ', str(args.data), '/train.features.clear2.csv')
	os.system(s.join( cmd ))


def compress(args):
	"""compress train.fasta files"""
	Printer(colored('(processing) ', 'green') + 'compress files')
	s = ""
	cmd = ('gzip --keep ', str(args.data), '/train.fasta')
	os.system(s.join( cmd ))

def extractkmers(args):
	Printer(colored('(processing) ', 'green') + 'extract kmer profile')
	s = ""
	cmd = ('src/fasta2kmers2 -i ', str(args.data), '/train.fasta -f ', str(args.data), '/train.features.kmer -j 4 -k 6 -s 0')
	os.system(s.join( cmd ))
	fname = str(args.data) + "/train.features.kmer"
	fname2 = str(args.data) + "/train.features.kmer.csv"
	with open(fname, "r") as inp, open(fname2 , "w") as out:
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
	X, y = creatematrix( str(args.data) + '/train.features.clear2.csv', str(args.data) + '/train.features.kmer', args)
	f = open(args.save, "w")
	pickle.dump(X, f)
	pickle.dump(y, f)
	f.close()

def savemsg(args):
	"""saves dataset as msgpack object"""
	Printer(colored('(processing) ', 'green') + 'save dataset as msg object')
	featurepath = str(args.data) + '/train.features.clear2.csv'
	kmerpath = str(args.data) + '/train.features.kmer'
	X, y = creatematrix(featurepath, kmerpath , args)
	y = pd.DataFrame(y)
	X.to_msgpack(args.save)
	y.to_msgpack(args.save, append=True)

def runHmmer(args, list_path, file_path, f):
	"""run prodigal and hmmsearch on chr files"""
	if not os.path.exists( str(args.data) + '/tmp'):
		os.makedirs(str(args.data) + '/tmp')
	# get the sample group
	head, group = os.path.split(os.path.split(file_path)[0])
	basename = os.path.splitext(str(ntpath.basename(str(file_path))))[0]
	exportpath =  str(args.data) + '/tmp/' + ntpath.basename(str(file_path))
	hmmpath =  str(args.data) + '/tmp/' + ntpath.basename(str(file_path)) + '.out'
	print('Processing %s of group %s' % (basename, group))
	s = ""
	cmd = ("prodigal -p meta -i ", str(file_path), " -a ", exportpath, ' -d /dev/null > /dev/null 2> /dev/null')
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
	remove_list =  str(args.data) + '/chr_to_remove.txt'
	Printer(colored('(processing) ', 'green') + 'screen downloaded chr for plasmid domains')
	fastalist = filter(os.path.isfile, glob.glob( str(args.data) + '/chr/*.fasta'))
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
	parser.add_argument('--save', action='store', dest='save', help='Save dataset as msg pack object', default='dataset.msg')
	parser.add_argument('--output', action='store', dest='data', help='path to output folder for raw data without tailing backslash', default='dat')
	parser.add_argument('-t', '--taxa', action='store', dest='taxa', help='Taxonomic name for downloaded samples', default='Escherichia coli')
	parser.add_argument('-a', '--planum', action='store', dest='planum', help='Number of plasmids to be downloaded', default=100)
	parser.add_argument('-b', '--chrnum', action='store', dest='chrnum', help='Number of chromosomes to be downloaded', default=20)
	parser.add_argument('-c', '--chunksize', action='store', dest='chunksize', help='Chunk size in nt', default=200)
	parser.add_argument('-s', '--readsim', dest='readsim', action='store_true', help='use read simluation based on wgsim instead of sliding window')
	parser.add_argument('-N', '--simnum', action='store', dest='simnum', help='number of reads to simulate with wgsim', default=1000)
	parser.add_argument('-e', "--email", help='Email adress needed for ncbi file download', dest="email", default="pmu15@helmholtz-hzi.de")
	parser.add_argument('--multi', dest='gpf', action='store_true', help='Prepare data in multi feature format (taxonomic column in train/test samples)')
	parser.add_argument('--no_download', dest='no_download', action='store_true', help='use dataset stored in data folder')
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
	sequence_cleaner(str(args.data) + '/chromosome_chunks.fasta.corrected.fasta', args)
	sequence_cleaner(str(args.data) + '/plasmid_chunks.fasta.corrected.fasta', args)
	# merge genomes to one file
	createmerged(args)
	# export features such as GC content
	getstatfeatures(args)
	# compress fasta file
	compress(args)
	# extract kmer content
	extractkmers(args)
	# save feature data as msgpack object
	savemsg(args)
