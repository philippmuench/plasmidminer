#!/usr/bin/python
import os
import download, simulate, features

# download training/learning dataset
email = "philipp.muench@helmholtz-hzi.de"
if not os.path.exists('pla'):
	download.downloadChr(email)
if not os.path.exists('chr'):
	download.downloadPla(email)

# get read sized chunks
splitlength = 200 # set length of outpuf fragments
if not os.path.exists('plasmid_200.fasta'):
	simulate.split(splitlength, 'pla/*.fasta', 'plasmid_200.fasta')
if not os.path.exists('chromosomes_200.fasta'):
	simulate.split(splitlength, 'chr/*.fasta', 'chromosomes_200.fasta')

# get feature matrix
#if not os.path.exists('plasmid_200.fasta.features'):
#	features.generate('plasmid_200.fasta')