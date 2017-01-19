#!/usr/bin/python
import os, csv, sys
import download, simulate, features
import argparse
from termcolor import colored

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
	"""
	Download FASTA files from NCBI
	args.taxa: str, part of the NCBI search command that contain the taxonomic information.
	args.chrnum: int, max. number of chromosom genomes that will be downloaded
	args.planum: int, max. number of plasmid genomes that will be downloaded
	args.gpf: boolean, If True, gpf files will be downloaded for each NCBI genome. False on deault.
	"""
	download.downloadChr(args.taxa, args.chrnum, args.gpf)
	download.downloadPla(args.taxa, args.planum, args.gpf)

def getchunks(args):
	"""Split FASTA files into same-sized fragments and exort them as multiple FASTA file"""
	Printer(colored('(processing) ', 'green') + 'rewrite fasta headers (plasmids)')
	simulate.split(int(args.chunksize), 'dat/pla/*.fasta', 'dat/plasmid_chunks.fasta', 'dat/pla/*.frag.fasta')
	simulate.renameheader('positive','dat/plasmid_chunks.fasta')
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
	os.system('src/fasta2kmers2 -i dat/train.fasta -f dat/train.features.kmer -j 4 -k 5 -s 0')
	with open("dat/train.features.kmer", "r") as inp, open("dat/train.features.kmer.csv", "w") as out:
		w = csv.writer(out, delimiter=",")
		w.writerows(x for x in csv.reader(inp, delimiter="\t"))

def sequence_cleaner(fasta_file, args):
	"""removes fragments smaller than chunksize from fasta file"""
	sequences={}
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		sequence = str(seq_record.seq).upper()
		if (len(sequence) >= args.chunksize):
			if sequence not in sequences:
				sequences[sequence] = seq_record.id
			else:
				sequences[sequence] += "_" + seq_record.id
	output_file = open(fasta_file + ".clear", "w+")
	for sequence in sequences:
		output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")
	output_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--taxa', action='store', dest='taxa', help='Taxonomic name for downloaded samples', default='Escherichia coli')
    parser.add_argument('-a', '--planum', action='store', dest='planum', help='Number of plasmids to be downloaded', default=100)
    parser.add_argument('-b', '--chrnum', action='store', dest='chrnum', help='Number of chromosomes to be downloaded', default=20)
    parser.add_argument('-c', '--chunksize', action='store', dest='chunksize', help='Chunk size in nt', default=200)
    parser.add_argument('-e', "--email", help='Email adress needed for ncbi file download', dest="email", default="pmu15@helmholtz-hzi.de")
    parser.add_argument('--multi', dest='gpf', action='store_true', help='Prepare data in multi feature format (taxonomic column in train/test samples)')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()

    Entrez.email = args.email # setting Entrez email
    loaddata(args)
    getchunks(args)
    sequence_cleaner('dat/chromosome_chunks.fasta.corrected.fasta', args)
    sequence_cleaner('dat/plasmid_chunks.fasta.corrected.fasta', args)
    createmerged()
    getstatfeatures()
    compress()	
    extractkmers()
