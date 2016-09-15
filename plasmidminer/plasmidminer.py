#!/usr/bin/python
import os, csv
import download, simulate, features

# download training/learning dataset
email = "philipp.muench@helmholtz-hzi.de"
if not os.path.exists('chr'):
	download.downloadChr(email)
if not os.path.exists('pla'):
	download.downloadPla(email)

# get read sized chunks
if not os.path.exists('plasmid_200.fasta'):
	simulate.split(200, 'pla/*.fasta', 'plasmid_200.fasta', 'pla/*.frag.fasta')

if not os.path.exists('chromosome_200.fasta'):
	simulate.split(200, 'chr/*.fasta', 'chromosome_200.fasta', 'chr/*.frag.fasta')

print('rewrite fasta headers')
simulate.renameheader('positive','plasmid_200.fasta')
simulate.renameheader('negative','chromosome_200.fasta')

print('merge pos/negative set to train.txt')
filenames = ['plasmid_200.fasta.corrected.fasta', 'chromosome_200.fasta.corrected.fasta']
with open('train.fasta', 'w') as outfile:
	for fname in filenames:
		with open(fname) as infile:
			for line in infile:
				outfile.write(line)

# get feature matrix
if not os.path.exists('train.features'):
	print('export features')
	os.system('python plasmidminer/features.py -I train.fasta -s length -s na -s cpg > train.features')
	print ('export csv')
	with open("train.features", "r") as inp, open("train.features.csv", "w") as out:
		w = csv.writer(out, delimiter=",")
		w.writerows(x for x in csv.reader(inp, delimiter="\t"))

# compress
if not os.path.exists('train.fasta.gz'):
	print('compressing fasta file')
	os.system("gzip --keep train.fasta")

if not os.path.exists('train.features.kmer'):
	print('get kmer profile')
	os.system('python plasmidminer/kmer.py -I train.fasta.gz --kmer-size 4 > train.features.kmer')