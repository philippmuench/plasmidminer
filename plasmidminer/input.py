#!/usr/bin/python
import os, csv, sys
import download, simulate, features
import argparse
from termcolor import colored
import cPickle
import pandas as pd
import numpy as np
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.base import BaseEstimator
from sklearn.base import ClassifierMixin
from sklearn.preprocessing import LabelEncoder
from sklearn.externals import six
from sklearn.base import clone
from sklearn.pipeline import _name_estimators
try:
	from Bio import Entrez
	from Bio import SeqIO
except ImportError:
	print "This script requires BioPython to be installed!"

class MajorityVoteClassifier(BaseEstimator,
                             ClassifierMixin):
    """ A majority vote ensemble classifier
    Parameters
    ----------
    classifiers : array-like, shape = [n_classifiers]
      Different classifiers for the ensemble
    vote : str, {'classlabel', 'probability'} (default='label')
      If 'classlabel' the prediction is based on the argmax of
        class labels. Else if 'probability', the argmax of
        the sum of probabilities is used to predict the class label
        (recommended for calibrated classifiers).
    weights : array-like, shape = [n_classifiers], optional (default=None)
      If a list of `int` or `float` values are provided, the classifiers
      are weighted by importance; Uses uniform weights if `weights=None`.
    """
    def __init__(self, classifiers, vote='classlabel', weights=None):

        self.classifiers = classifiers
        self.named_classifiers = {key: value for key, value
                                  in _name_estimators(classifiers)}
        self.vote = vote
        self.weights = weights

    def fit(self, X, y):
        """ Fit classifiers.
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Matrix of training samples.
        y : array-like, shape = [n_samples]
            Vector of target class labels.
        Returns
        -------
        self : object
        """
        if self.vote not in ('probability', 'classlabel'):
            raise ValueError("vote must be 'probability' or 'classlabel'"
                             "; got (vote=%r)"
                             % self.vote)

        if self.weights and len(self.weights) != len(self.classifiers):
            raise ValueError('Number of classifiers and weights must be equal'
                             '; got %d weights, %d classifiers'
                             % (len(self.weights), len(self.classifiers)))

        # Use LabelEncoder to ensure class labels start with 0, which
        # is important for np.argmax call in self.predict
        self.lablenc_ = LabelEncoder()
        self.lablenc_.fit(y)
        self.classes_ = self.lablenc_.classes_
        self.classifiers_ = []
        for clf in self.classifiers:
            fitted_clf = clone(clf).fit(X, self.lablenc_.transform(y))
            self.classifiers_.append(fitted_clf)
        return self

    def predict(self, X):
        """ Predict class labels for X.
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Matrix of training samples.
        Returns
        ----------
        maj_vote : array-like, shape = [n_samples]
            Predicted class labels.
        """
        if self.vote == 'probability':
            maj_vote = np.argmax(self.predict_proba(X), axis=1)
        else:  # 'classlabel' vote

            #  Collect results from clf.predict calls
            predictions = np.asarray([clf.predict(X)
                                      for clf in self.classifiers_]).T

            maj_vote = np.apply_along_axis(
                lambda x:
                np.argmax(np.bincount(x,
                                      weights=self.weights)),
                axis=1,
                arr=predictions)
        maj_vote = self.lablenc_.inverse_transform(maj_vote)
        return maj_vote

    def predict_proba(self, X):
        """ Predict class probabilities for X.
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.
        Returns
        ----------
        avg_proba : array-like, shape = [n_samples, n_classes]
            Weighted average probability for each class per sample.
        """
        probas = np.asarray([clf.predict_proba(X)
                             for clf in self.classifiers_])
        avg_proba = np.average(probas, axis=0, weights=self.weights)
        return avg_proba

    def get_params(self, deep=True):
        """ Get classifier parameter names for GridSearch"""
        if not deep:
            return super(MajorityVoteClassifier, self).get_params(deep=False)
        else:
            out = self.named_classifiers.copy()
            for name, step in six.iteritems(self.named_classifiers):
                for key, value in six.iteritems(step.get_params(deep=True)):
                    out['%s__%s' % (name, key)] = value
        return out
    
class Printer():
	"""Print things to stdout on one line dynamically"""
	def __init__(self,data):
		sys.stdout.write("\r\x1b[K"+data.__str__())
		sys.stdout.flush()

def getstatfeatures():
	"""export statistical properties of fasta files as features"""
	Printer(colored('(processing) ', 'green') + 'generte feature matrix')
	os.system('python plasmidminer/features.py -I dat/input.fasta -s length -s na -s cpg > dat/input.features')

	Printer(colored('(processing) ', 'green') + 'export csv')
	with open("dat/input.features", "r") as inp, open("dat/input.features.csv", "w") as out:
		w = csv.writer(out, delimiter=",")
		w.writerows(x for x in csv.reader(inp, delimiter="\t"))
	features.clearit("dat/input.features.csv", "dat/input.features.clear.csv")
	os.system("tail -n +2 dat/input.features.clear.csv > dat/input.features.clear2.csv")

def extractkmers():
	Printer(colored('(processing) ', 'green') + 'extract kmer profile')
	os.system('src/fasta2kmers2 -i dat/input.fasta -f dat/input.features.kmer -j 4 -k 5 -s 0 -n 1')
	with open("dat/input.features.kmer", "r") as inp, open("dat/input.features.kmer.csv", "w") as out:
		w = csv.writer(out, delimiter=",")
		w.writerows(x for x in csv.reader(inp, delimiter="\t"))

def createpredictmatrix(features, kmer):
	stat = pd.read_csv(features, sep=",")
	kmer = pd.read_csv(kmer, sep="\t", header = None) 
	df = pd.concat([stat.reset_index(drop=True), kmer], axis=1) # concat kmer and stat matrix
	df_complete = df.dropna(axis=1) # remove columns with NAN
	return df_complete

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	args = parser.parse_args()

	getstatfeatures()	
	extractkmers()
	print('create matrix')
	dat = createpredictmatrix('dat/input.features.clear2.csv', 'dat/input.features.kmer')
	X = dat.drop(dat.columns[[0]], 1) # remove label
	with open('model.pkl', 'rb') as fid:
		pipe = cPickle.load(fid)
	print('classifier loaded')
	label = {0:'chromosomal', 1:'plasmid'}

	for index, row in X.iterrows():
		print('Prediction: %s\nProbability: %.2f%%' %\
			(label[pipe.predict(X)[index]], pipe.predict_proba(X)[index].max()*100))
