#!/usr/bin/python
import os
import glob
import ntpath

def runHmmer(folder, path):
	for filename in glob.iglob(folder):
		print('search ORF in: %s' % filename)
		# find ORF in fasta files
		exportpath = str(path) + ntpath.basename(str(filename))
		hmmpath = str(path) + ntpath.basename(str(filename)) + '.out'

		s = " "
		cmd = ("prodigal -i",  str(filename), "-a", exportpath, '-d /dev/null > /dev/null 2> /dev/null')
		os.system(s.join( cmd ))
		# run hmmsearch on faa ORF files
		s = " "
		cmd = ("hmmsearch --pfamtblout", hmmpath, "resources/plasmid.hmm", exportpath, '> /dev/null 2> /dev/null')
		os.system(s.join( cmd ))

if not os.path.exists('dat/chr_faa'):
	os.makedirs('dat/chr_faa')
if not os.path.exists('dat/pla_faa'):
	os.makedirs('dat/pla_faa')

runHmmer('dat/chr/*.fasta', 'dat/chr_faa/')

