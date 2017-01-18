#!/usr/bin/env python
# philipp.muench@helmholtz-hzi.de
# http://www.ncbi.nlm.nih.gov/genome/browse/
import time, os, glob, fileinput, sys
from termcolor import colored

try:
    from progress.bar import Bar
except ImportError:
    print "This script requires progess to be installed!"

try:
    from Bio import SeqIO
except ImportError:
    print "This script requires BioPython to be installed!"
try:
    from Bio import Entrez
except ImportError:
    print "This script requires BioPython to be installed!"

class Printer():
    """Print things to stdout on one line dynamically"""
    def __init__(self,data):
        sys.stdout.write("\r\x1b[K"+data.__str__())
        sys.stdout.flush()

def downloadChr(taxa, num):
	if not os.path.exists('dat/chr'):
		os.makedirs('dat/chr')
	search_term = '"' + str(taxa) + '"[Organism] AND complete genome[title]'
	Printer(colored('(download chromosomes) ', 'green') + 'searching for chromosomes (search term: ' + search_term + ')')
	time.sleep(5)
	handle = Entrez.esearch(db='nucleotide', term=search_term, retmax =num)
	genome_ids = Entrez.read(handle)['IdList']
	records = len(genome_ids)
	i = 1
	#bar = Bar('download chromosom sequences', max=records)
	for genome_id in genome_ids:
		Printer(colored('(download chromosomes) ', 'green')+ colored('['+str(i) +'/' + str(len(genome_ids)) + '] ', 'blue')+ 'fetching database (ID:' + str(genome_id) + ')') 
		record = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
		#record_gpf = Entrez.efetch(db="protein", id=genome_id, format="gpf") ## comment out if you want to download gfp files
		#filename_gpf = 'dat/chr/geneBankRecord_{}.gpf'.format(genome_id) ## comment out if you want to download gfp files
		filename = 'dat/chr/genBankRecord_{}.fasta'.format(genome_id)
		if not os.path.exists(filename):
			Printer(colored('(download chromosomes) ', 'green') + colored('['+str(i) +'/' + str(len(genome_ids))+ '] ', 'blue') + 'add genome (ID: ' + str(genome_id) + ') to collection')
			with open(filename, 'w') as f:
				f.write(record.read())
			#with open(filename_gpf, 'w') as f: ## comment out if you want to download gfp files
				#f.write(record_gpf.read()) ## comment out if you want to download gfp files
			time.sleep(1) # to make sure not many requests go per second to ncbi
		else:
			Printer(colored('(download plasmids) ', 'green') + colored('[' + str(i) + '/' + str(len(genome_ids))+ '] ', 'blue') + ' (ID: ' + str(genome_id) + ') is in collection. Skipping.')
		i = i + 1

    #    bar.next()
    #bar.finish()
    # create multiple sequence fasta file
	filenames = glob.glob("dat/chr/*.fasta")
	Printer(colored("merging files...", 'green'))
	with open('dat/chromosome.fasta', 'w') as outfile:
		for fname in filenames:
			with open(fname) as infile:
				for line in infile:
					outfile.write(line)

def downloadPla(taxa, num):
	# imports all E.Coli plasmids as positive samples
	# create download folders
	if not os.path.exists('dat/pla'):
		os.makedirs('dat/pla')
	# get the genome ids
	search_term = "'" + str(taxa) + '"[Organism] AND (bacteria[filter] AND plasmid[filter])'
	Printer(colored('(download plasmids) ', 'green') + 'searching for chromosomes (search term: ' + search_term + ')')
	handle = Entrez.esearch(db='nucleotide', term=search_term, retmax =num)
	genome_ids = Entrez.read(handle)['IdList']
  	records = len(genome_ids)
	#bar = Bar('download plasmid sequences', max=records)
	# iterate over genome ids and download informations
	i = 1
	for genome_id in genome_ids:
		Printer(colored('(download plasmids) ', 'green')+ colored('['+str(i) +'/' + str(len(genome_ids)) + '] ', 'blue')+ 'fetching database (ID:' + str(genome_id) + ')') 
		record = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
		record_gpf = Entrez.efetch(db="protein", id=genome_id, format="gpf")
		# write gb record
		filename = 'dat/pla/genBankRecord_{}.fasta'.format(genome_id)
		#filename_gpf = 'dat/pla/geneBankRecord_{}.gpf'.format(genome_id) ## comment out if you want to download gfp files
		if not os.path.exists(filename):
			Printer(colored('(download plasmids) ', 'green') + colored('[' + str(i) + '/' + str(len(genome_ids))+ '] ', 'blue') + 'add genome (ID: ' + str(genome_id) + ') to collection')
			with open(filename, 'w') as f:
				f.write(record.read())
        #with open(filename_gpf, 'w') as f: ## comment out if you want to download gfp files
        	#f.write(record_gpf.read()) ## comment out if you want to download gfp files
		else:
			Printer(colored('(download plasmids) ', 'green') + colored('[' + str(i) + '/' + str(len(genome_ids))+ '] ', 'blue') + ' (ID: ' + str(genome_id) + ') is in collection. Skipping.')
		i = i + 1
		time.sleep(1) # to make sure not many requests go per second to ncbi
	#	bar.next()
	#bar.finish()
    # create multiple sequence fasta file
	filenames = glob.glob("dat/pla/*.fasta")
	Printer(colored("merging files...", 'green'))
	with open('dat/plasmid.fasta', 'w') as outfile:
		for fname in filenames:
			with open(fname) as infile:
				for line in infile:
					outfile.write(line)
if __name__ == "__main__":
    sys.exit(main(sys.argv))
