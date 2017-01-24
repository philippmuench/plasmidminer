#!/usr/bin/python

import argparse, os, csv, sys, re
import download, simulate, features
import pandas as pd
import cPickle

def createpredictmatrix(features, kmer):
	stat = pd.read_csv(features, sep=",")
	kmer = pd.read_csv(kmer, sep="\t", header = None) 
	df = pd.concat([stat.reset_index(drop=True), kmer], axis=1) # concat kmer and stat matrix
	df_complete = df.dropna(axis=1) # remove columns with NAN
	return df_complete

def predict(args):
	print ('predicting %s ') % args.fasta_file

	# generate compostition statistic features	
	s = " "
	cmd = ("python plasmidminer/features.py -I", args.fasta_file, "-s length -s na -s cpg > dat/test.features")
	os.system(s.join( cmd ))
	
	with open("dat/test.features", "r") as inp, open("dat/test.features.csv", "w") as out:
		w = csv.writer(out, delimiter=",")
		w.writerows(x for x in csv.reader(inp, delimiter="\t"))
	features.clearit("dat/test.features.csv", "dat/test.features.clear.csv")
	os.system("tail -n +2 dat/test.features.clear.csv > dat/test.features.clear2.csv")

	# generate kmer features
	s = " "
	cmd = ("src/fasta2kmers2 -i", args.fasta_file, "-f dat/test.features.kmer -j 4 -k 5 -s 0")
	os.system(s.join( cmd ))
	with open("dat/test.features.kmer", "r") as inp, open("dat/test.features.kmer.csv", "w") as out:
		w = csv.writer(out, delimiter=",")
		w.writerows(x for x in csv.reader(inp, delimiter="\t"))
	
	# predict data use pickle model
	dat = createpredictmatrix('dat/test.features.clear2.csv', 'dat/test.features.kmer')
	X = dat.drop(dat.columns[[0]], 1) # remove label

	print(X)
	# load classifier
	with open('dat/rf.pkl', 'rb') as fid:
		rf_mod = cPickle.load(fid)
	print('classifier loaded')
	label = {0:'chromosomal', 1:'plasmid'}
	print('Prediction: %s\nProbability: %.2f%%' %\
		(label[rf_mod.predict(X)[0]], rf_mod.predict_proba(X)[0].max()*100))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta_file', type=str, help='path to input fasta file')
	parser.add_argument('--save_dir', type=str, default='save', help='directory to store output files')
	
	parser.add_argument("--model", type=str, default="dat/rf.pkl", help='path to model .pkl file')
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	args = parser.parse_args()

	print 'input folder     =', args.fasta_file
	print 'output folder     =', args.save_dir
	print 'model file     =', args.model

	predict(args)
