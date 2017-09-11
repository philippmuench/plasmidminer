#!/usr/bin/python
import os, csv, sys
import download, simulate, features
import argparse
from termcolor import colored
import cPickle
import pandas as pd
import numpy as np
import tempfile
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.base import BaseEstimator
from sklearn.base import ClassifierMixin
from sklearn.preprocessing import LabelEncoder
from sklearn.externals import six
from sklearn.base import clone
from sklearn.pipeline import _name_estimators
from sklearn.externals import joblib
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

def getstatfeatures(args):
    """export statistical properties of fasta files as features"""
    Printer(colored('(processing) ', 'green') + 'generate feature matrix')
    s = " "
    cmd = ("python plasmidminer/features.py -I ", str(args.input), "-s length -s na -s cpg > dat_tmp/input.features")
    os.system(s.join( cmd ))
    Printer(colored('(processing) ', 'green') + 'export csv')
    with open("dat_tmp/input.features", "r") as inp, open("dat_tmp/input.features.csv", "w") as out:
        w = csv.writer(out, delimiter=",")
        w.writerows(x for x in csv.reader(inp, delimiter="\t"))
    features.clearit("dat_tmp/input.features.csv", "dat_tmp/input.features.clear.csv")
    os.system("tail -n +2 dat_tmp/input.features.clear.csv > dat_tmp/input.features.clear2.csv")

def extractkmers(args):
    Printer(colored('(processing) ', 'green') + 'extract kmer profile')
    s = " "
    cmd = ("src/fasta2kmers2 -i ", str(args.input), "-f dat_tmp/input.features.kmer -j 4 -k 6 -s 0 -n 1")
    os.system(s.join( cmd ))
    with open("dat_tmp/input.features.kmer", "r") as inp, open("dat_tmp/input.features.kmer.csv", "w") as out:
        w = csv.writer(out, delimiter=",")
        w.writerows(x for x in csv.reader(inp, delimiter="\t"))

def createpredictmatrix(features, kmer):
    stat = pd.read_csv(features, sep=",")
    kmer = pd.read_csv(kmer, sep="\t", header = None) 
    df = pd.concat([stat.reset_index(drop=True), kmer], axis=1) # concat kmer and stat matrix
    df_complete = df.dropna(axis=1) # remove columns with NAN
    return df_complete

def sepseq(fasta_file, predictions, probabilities, args):
    """seperates plasmid and chromosomal fragments based on prediction"""
    plasmid_sequences={}
    chromsom_sequences={}
    i = 0
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(seq_record.seq).upper()
        if (predictions[i] == 1):
            if (args.probability):
                plasmid_sequences[sequence] = seq_record.description + "-" + str(probabilities[i])
            else:
                plasmid_sequences[sequence] = seq_record.description 
        else:
            if (args.probability):
                chromsom_sequences[sequence] = seq_record.description + "-" +  str(probabilities[i])
            else:
                chromsom_sequences[sequence] = seq_record.description
        i += 1

    output_file = open(fasta_file + ".plasmids", "w+")
    for sequence in plasmid_sequences:
        output_file.write(">" + plasmid_sequences[sequence] + "\n" + sequence + "\n")
    output_file.close()

    output_file = open(fasta_file + ".chromosomes", "w+")
    for sequence in chromsom_sequences:
        output_file.write(">" + chromsom_sequences[sequence] + "\n" + sequence + "\n")
    output_file.close()
    return plasmid_sequences, chromsom_sequences

def slidewindowfragments(args):
    """Generate sliding window fragments from fasta file"""
    with open('dat_tmp/window_fragments.fasta',"w") as f:
        for seq_record in SeqIO.parse(args.input, "fasta"):
            for i in range(len(seq_record.seq) - int(args.window) - 1) :
               f.write(str(">" + seq_record.id) + "\n")
               f.write(str(seq_record.seq[i:i + int(args.window)]) + "\n")  #first 5 base positions

def showresults(args, plasmid_sequences, chromsom_sequences):
	"""outputs basic statistics of input files"""
	print('\n-------------------------------------------------------------')
	print len(plasmid_sequences), "\t plasmid reads found"
	print len(chromsom_sequences), "\t chromosome reads found"
	print "chromosome fasta file:\t", args.input+".chromosomes"
	print "plasmid fasta file:\t", args.input+".plasmids"
	print "probability file:\t", "probabilities.txt"
	print "prediction file:\t", "predictions.txt"
	print('-------------------------------------------------------------\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input', help='Path to input FASTA file', default='dat/input.fasta')
    parser.add_argument('-m', '--model', action='store', dest='model', help='Path to model (.pkl)', default='model.pkl')
    parser.add_argument('-s', '--split', dest='split', action='store_true', help='split data based on prediction')
    parser.add_argument('-p', '--probability', dest='probability', action='store_true', help='add probability information to fasta header')
    parser.add_argument('--sliding', dest='sliding', action='store_true', help='Predict on sliding window')
    parser.add_argument('--window', action='store', dest='window', help='size of sliding window', default='100')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()
    if (args.sliding):
        Printer(colored('(preprocessing) ', 'green') + 'run sliding window')
        slidewindowfragments(args)
        args.input = 'dat_tmp/window_fragments.fasta'

    getstatfeatures(args)
    extractkmers(args)
    Printer(colored('(preprocessing) ', 'green') + 'create matrix')
    dat = createpredictmatrix('dat_tmp/input.features.clear2.csv', 'dat_tmp/input.features.kmer')
    X = dat.drop(dat.columns[[0]], 1) # remove label
    Printer(colored('(preprocessing) ', 'green') + 'import model')
    with open(args.model, 'rb') as pkl_source:
        pipe = joblib.load(pkl_source)
    label = {0:'chromosomal', 1:'plasmid'}
    Printer(colored('(running) ', 'green') + 'save predictions to file')
    predictions = pipe.predict(X)
    np.savetxt('predictions.txt', predictions)
    if (args.sliding):
    	Printer(colored('(running) ', 'green') + 'make prediction within sliding window')
        probabilities = pipe.predict_proba(X)[:,1]# probabilitiesility that this is a plasmid
        np.savetxt('window_probabilities.txt', probabilities)
		#for index, row in X.iterrows():
		#	print('Prediction: %s\nProbability: %.2f%%' %\
		#		(label[pipe.predict(X)[index]], pipe.predict_proba(X)[index].max()*100))

    else:
    	Printer(colored('(running) ', 'green') + 'predict input fasta file')
        probabilities = np.amax(pipe.predict_proba(X), axis=1) # max probability over two classes
        np.savetxt('probabilities.txt', probabilities)
    if (args.split):
    	Printer(colored('(running) ', 'green') + 'generate fasta files')
        plasmid_sequences, chromsom_sequences = sepseq(args.input, predictions, probabilities, args)
        Printer(colored('(running) ', 'green') + 'finished')
    	showresults(args, plasmid_sequences, chromsom_sequences)

