#!/usr/bin/python
import os, csv
import download, simulate, features

# download training/testing dataset
email = "philipp.muench@helmholtz-hzi.de"
if not os.path.exists('dat/chr'):
	download.downloadChr(email)
if not os.path.exists('dat/pla'):
	download.downloadPla(email)

# get read sized chunks
if not os.path.exists('dat/plasmid_200.fasta'):
	simulate.split(200, 'dat/pla/*.fasta', 'dat/plasmid_200.fasta', 'dat/pla/*.frag.fasta')

if not os.path.exists('dat/chromosome_200.fasta'):
	simulate.split(200, 'dat/chr/*.fasta', 'dat/chromosome_200.fasta', 'dat/chr/*.frag.fasta')

print('rewrite fasta headers')
simulate.renameheader('positive','dat/plasmid_200.fasta')
simulate.renameheader('negative','dat/chromosome_200.fasta')

print('merge pos/negative set to train.txt')
filenames = ['dat/plasmid_200.fasta.corrected.fasta', 'dat/chromosome_200.fasta.corrected.fasta']
with open('dat/train.fasta', 'w') as outfile:
	for fname in filenames:
		with open(fname) as infile:
			for line in infile:
				outfile.write(line)

# get feature matrix
if not os.path.exists('dat/train.features'):
    print('export features')
    os.system('python plasmidminer/features.py -I dat/train.fasta -s length -s na -s cpg > dat/train.features')
    print ('export csv')
    with open("dat/train.features", "r") as inp, open("dat/train.features.csv", "w") as out:
        w = csv.writer(out, delimiter=",")
        w.writerows(x for x in csv.reader(inp, delimiter="\t"))        
    features.clearit("dat/train.features.csv", "dat/train.features.clear.csv")
    os.system("tail -n +2 dat/train.features.clear.csv > dat/train.features.clear2.csv")

# compress
if not os.path.exists('dat/train.fasta.gz'):
	print('compressing fasta file')
	os.system("gzip --keep dat/train.fasta")

if not os.path.exists('dat/train.features.kmer'):
    print('get kmer profile')
    os.system('src/fasta2kmers2 -i dat/train.fasta -f dat/train.features.kmer -j 4 -k 5 -s 0')
    with open("dat/train.features.kmer", "r") as inp, open("dat/train.features.kmer.csv", "w") as out:
            w = csv.writer(out, delimiter=",")
            w.writerows(x for x in csv.reader(inp, delimiter="\t"))

