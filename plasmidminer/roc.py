#!/usr/bin/python
# script to compare model performances
# philipp.muench@helmholtz-hzi.de

import os, glob, sys
import argparse
import findparameters
from termcolor import colored
import gzip
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.externals import joblib
from sklearn.model_selection import cross_val_score, train_test_split

class Printer():
	"""Print things to stdout on one line dynamically"""
	def __init__(self,data):
		sys.stdout.write("\r\x1b[K"+data.__str__())
		sys.stdout.flush()
def showinfo(args):
	"""outputs basic info about input files"""
	Printer(colored('(processing) ', 'green') + 'listing input file')
	modellist = filter(os.path.isfile, glob.glob(str(args.model) + '/*.pkl'))
	if not modellist:
		sys.stderr.write("Error: No models found in models_folder. Please point '--model_folder' parameter to the input folder that contains one or multiple pkl models. This is usually the output folder of findparameters.py\n")
		exit(1)
	print('\n-------------------------------------------------------------')
	print "models found: ", len(modellist), '\t', str(modellist)
	print('-------------------------------------------------------------\n')
	return modellist

def drawroc(clf, clf_labels, X_train, y_train, X_test, y_test):
	"""draw a roc curve for each model in clf, save as roc.png"""
	Printer(colored('(running) ', 'green') + 'draw ROC')
	colors = ['black', 'orange', 'blue', 'green']
	linestyles = [':', '--', '-.', '-']
	for clf, label, clr, ls \
			in zip(clf, clf_labels, colors, linestyles):
		scores = cross_val_score(estimator=clf, X=X_train, y=y_train, cv=int(
			args.cv), scoring='roc_auc', n_jobs=-1)
		#print("ROC AUC: %0.2f (+/- %0.2f) [%s]" %
		#	  (scores.mean(), scores.std(), label))
		y_pred = clf.fit(X_train, y_train).predict_proba(X_test)[:, 1]
		fpr, tpr, thresholds = metrics.roc_curve(y_true=y_test, y_score=y_pred)
		roc_auc = metrics.auc(x=fpr, y=tpr)
		plt.plot(fpr, tpr, color=clr, linestyle=ls,
				 label='%s (auc = %0.2f)' % (label, roc_auc))
	plt.legend(loc='lower right')
	plt.plot([0, 1], [0, 1], linestyle='--', color='gray', linewidth=2)
	plt.xlim([-0.1, 1.1])
	plt.ylim([-0.1, 1.1])
	plt.grid()
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	filename = 'roc_all.png'
	# plt.tight_layout()
	plt.savefig(filename, dpi=300)
	plt.close()

def loadmodels(modellist, args):
	"""iterates over pkl files in foulder and loads pkl objects"""
	clf = []
	clf_labels = []
	for filename in modellist:
		Printer(colored('(processing) ', 'green') + 'load model %s' % filename)
		with open(filename, 'rb') as fid:
			# extract model name from file name
			head, tail = os.path.split(os.path.split(filename)[1])
			basename = os.path.splitext(tail)[0]
			modelname = basename.split('_', 1)[0]
			#modelname
			clf_labels.append(modelname)
			clf.append(joblib.load(fid))
	return clf, clf_labels


def showinfoclf(clf, clf_labels):
	"""outputs basic info about input files"""
	Printer(colored('(processing) ', 'green') + 'listing cls')
	print('\n-------------------------------------------------------------')
	print "models loaded:\t", len(clf)
	print "annotations loaded: ", len(clf_labels), ':\t', str(clf_labels)
	print('-------------------------------------------------------------\n')

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--model_folder', action='store', dest='model',
						help='path to folder where .pkl models are located', default='cv')
	parser.add_argument('-d', '--data_folder', action='store', dest='data',
						help='path to file foulder', default='dat')
	parser.add_argument('-t', '--test_size', action='store', dest='test_size',
						help='size of test set from whole dataset in percent', default=0.33)
	parser.add_argument('-c', '--cv', action='store', dest='cv',
						help='cross validation size (e.g. 10 for 10-fold cross validation)', default=3)
	parser.add_argument('--version', action='version', version='%(prog)s 0.1')
	args = parser.parse_args()

	# import model
	modellist = showinfo(args)
	clf, clf_labels = loadmodels(modellist, args)

	# print model information to screen
	showinfoclf(clf, clf_labels)

	# import samples
	X, y = findparameters.creatematrix(str(args.data) + '/train.features.clear2.csv', str(args.data) + '/train.features.kmer', args)

	# generate subsample with equal amount of pos/neg instances
	X_sub, y_sub = findparameters.balanced_subsample(X, y)

	# split train/testset
	#Printer(colored('(preprocessing) ', 'green') + 'generate train/test set')
	X_train, X_test, y_train, y_test = train_test_split(X_sub, y_sub, test_size=args.test_size)

	# show input data statistics
 	findparameters.showinput(X, y, 'raw dataset')
 	findparameters.showinput(X_train, y_train, 'train subset')
 	findparameters.showinput(X_test, y_test, 'test subset')
	
	# get figure
	drawroc(clf, clf_labels, X_train, y_train, X_test, y_test)
