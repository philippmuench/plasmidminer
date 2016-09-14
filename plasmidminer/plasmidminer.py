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

# get feature matrix
if not os.path.exists('plasmid_200.fasta.features'):
	print('export features for plasmid_200.fasta')
	os.system('python plasmidminer/features.py -I plasmid_200.fasta -s length -s na -s cpg > plasmid_200.fasta.features')

if not os.path.exists('chromosome_200.fasta.features'):
	print('export features for chromosome_200.fasta')
	os.system('python plasmidminer/features.py -I chromosome_200.fasta -s length -s na -s cpg > chromosome_200.fasta.features')

with open("plasmid_200.fasta.features", "r") as inp, open("plasmid_200.fasta.features.clean.csv", "w") as out:
	w = csv.writer(out, delimiter=",")
	w.writerows(x for x in csv.reader(inp, delimiter="\t"))

with open("chromosome_200.fasta.features", "r") as inp, open("chromosome_200.fasta.features.clean.csv", "w") as out:
	w = csv.writer(out, delimiter=",")
	w.writerows(x for x in csv.reader(inp, delimiter="\t"))

#if not os.path.exists('pos.train.csv'):
#	print('create pos.csv')
#	

#if not os.path.exists('neg.train.csv'):
#	with open('chromosome_200.fasta.features.clean.csv','r') as csvinput:
#		print('create neg.csv')
#		with open('neg.train.csv', 'w') as csvoutput:
#			writer = csv.writer(csvoutput, lineterminator='\n')
#			reader = csv.reader(csvinput)
#			all = []
#			row = next(reader)
#			row.append('Type')
#			all.append(row)
#			for row in reader:
#				row.append(row[0])
#				all.append(row)
#			writer.writerows(all)