#!/usr/bin/python
import os
import download, simulate, features

# download training/learning dataset
email = "philipp.muench@helmholtz-hzi.de"
if not os.path.exists('chr'):
	download.downloadChr(email)
if not os.path.exists('pla'):
	download.downloadPla(email)

# get read sized chunks
splitlength = 200 # set length of outpuf fragments
if not os.path.exists('plasmid_200.fasta'):
	simulate.split(splitlength, 'pla/*.fasta', 'plasmid_200.fasta')
if not os.path.exists('chromosomes_200.fasta'):
	simulate.split(splitlength, 'chr/*.fasta', 'chromosomes_200.fasta')

# get feature matrix
#if not os.path.exists('plasmid_200.fasta.features'):
print('export features for plasmid_200.fasta')
os.system('plasmidminer/features.py -I plasmid_200.fasta -s length -s na -s cpg > plasmid_200.fasta.features')
print('export features for chromosomes_200.fasta')
os.system('plasmidminer/features.py -I chromosomes_200.fasta -s length -s na -s cpg > chromosomes_200.fasta.features')