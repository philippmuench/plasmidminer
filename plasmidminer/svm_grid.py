#!/usr/bin/python
from __future__ import print_function
try:
	import pandas as pd
except ImportError:
    print "This script requires pandas to be installed!"

try:
	import numpy as np
except ImportError:
    print "This script requires numpy to be installed!"

try:
	import operator
except ImportError:
    print "This script requires numpy to be installed!"
try:
	import cPickle
except ImportError:
    print "This script requires numpy to be installed!"
try:
	import argparse
except ImportError:
    print "This script requires numpy to be installed!"
try:
	from sklearn import datasets
	from sklearn.model_selection import train_test_split
	from sklearn.model_selection import GridSearchCV
	from sklearn.metrics import classification_report
	from sklearn.svm import SVC
	from sklearn.cross_validation import train_test_split
	from sklearn.preprocessing import StandardScaler
	from sklearn.base import BaseEstimator
	from sklearn.base import ClassifierMixin
	from sklearn.preprocessing import LabelEncoder
	from sklearn.externals import six
	from sklearn.base import clone
	from sklearn.pipeline import _name_estimators
	from sklearn.pipeline import Pipeline
	from sklearn import svm
	from sklearn.svm import SVC
	from sklearn.model_selection import cross_val_score
except ImportError:
    print "This script requires sklearn to be installed!"
try:
	from skrvm import RVC
except ImportError:
    print "This script requires skrvm to be installed!"


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
    
def creatematrix(features, kmer):
    stat = pd.read_csv(features, sep=",")
    kmer = pd.read_csv(kmer, sep="\t", header = None) 
    df = pd.concat([stat.reset_index(drop=True), kmer], axis=1) # concat kmerand stat matrix
    id2 = df.id.str.split("-",expand=True) # split the string to get label
    df2 = pd.concat([id2, df], axis=1)
    all = df2.drop(df2.columns[[1,2]], 1)
    all_complete = all[(all.length == 200)] # keep only complete fragments
    all_complete.columns.values[0] = "label"
    # encoding class labels as integers
    all_complete.loc[all_complete.label == 'positive', 'label'] = 1
    all_complete.loc[all_complete.label == 'negative', 'label'] = 0
    all_complete = all_complete.dropna(axis=1) # remove columns with NAN
    return all_complete

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--model', action='store', dest='model', help='Patho to model file.')
	parser.add_argument('--features', action='store', dest='features', help='Path to features matrix', default='dat/train.features.clear2.csv')
	parser.add_argument('--kmers', action='store', dest='kmers', help='Path to features matrix', default='dat/train.features.kmer')
	parser.add_argument('--balanced', action='store_true')
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	args = parser.parse_args()
	print 'model =', args.model
	print 'features =', args.features
	print 'kmers =', args.kmers
	print 'balanced = ', args.balanced 

	# read data
	dat = creatematrix(args.features, args.kmers)

	# check if the dataset is inbalanced
	dat.pos = dat.loc[dat['label'] == 1]
	dat.neg = dat.loc[dat['label'] == 0]

	print('number of negative instances: %d ' % (dat.neg.shape[0]))
	print('number of positive instances: %d ' % (dat.pos.shape[0]))

	if args.balanced:
		# generate a random subset of both with the size of $num
		num = min(dat.neg.shape[0],dat.pos.shape[0])
		print('limit dataset size to %d' % (num))
		posrand = dat.pos.sample(n=num)
		negrand = dat.neg.sample(n=num)
		df = posrand.copy()
		dat = df.append(negrand)

	# split to test and training set
	y = dat['label'].tolist() # extract label
	X = dat.drop(dat.columns[[0]], 1) # remove label
	print("data procesed")

	# generate cross validation datasets, split x,y arrays into 30percent test data and 70 percent training data
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=0)
	print("train/test set generated")

	scores = ['precision', 'recall']

	# Set the parameters by cross-validation
	tuned_parameters = [{'kernel': ['rbf'], 'gamma': [100, 10, 1, 0.1, 1e-2 ,1e-3, 1e-4, 1e-5],
		             'C': [1, 10, 100, 1000, 10000]},
		            {'kernel': ['linear'], 'C': [0.01, 0.1, 1, 10, 100, 1000, 10000]},
		            {'kernel': ['poly'], 'degree': [1, 3, 5, 10, 50]]

	for score in scores:
	    print("# Tuning hyper-parameters for %s" % score)
	    print()

	    clf = GridSearchCV(SVC(C=1), tuned_parameters, cv=5,
		               scoring='%s_macro' % score, n_jobs=40)
	    clf.fit(X_train, y_train)

	    print("Best parameters set found on development set:")
	    print()
	    print(clf.best_params_)
	    print()
	    print("Grid scores on development set:")
	    print()
	    means = clf.cv_results_['mean_test_score']
	    stds = clf.cv_results_['std_test_score']
	    for mean, std, params in zip(means, stds, clf.cv_results_['params']):
		print("%0.3f (+/-%0.03f) for %r"
		      % (mean, std * 2, params))
	    print()

	    print("Detailed classification report:")
	    print()
	    print("The model is trained on the full development set.")
	    print("The scores are computed on the full evaluation set.")
	    print()
	    y_true, y_pred = y_test, clf.predict(X_test)
	    print(classification_report(y_true, y_pred))
	    print()
