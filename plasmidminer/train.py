#!/usr/bin/python

import operator

from time import time
from scipy.stats import randint as sp_randint
import sys
import os
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cPickle
import argparse
from termcolor import colored
from skbayes.rvm_ard_models import RVR,RVC
try:
	from sklearn.model_selection import GridSearchCV
	from sklearn.model_selection import RandomizedSearchCV
	from sklearn.datasets import load_digits
	from sklearn.ensemble import RandomForestClassifier
	from sklearn.linear_model import LogisticRegression
	from sklearn.tree import DecisionTreeClassifier
	from sklearn.neighbors import KNeighborsClassifier
	from sklearn.pipeline import Pipeline
	from sklearn.ensemble import VotingClassifier, RandomForestClassifier
	from sklearn import svm
	from sklearn.svm import SVC
	from sklearn.model_selection import cross_val_score, train_test_split
	from sklearn.preprocessing import StandardScaler
	from sklearn.base import BaseEstimator
	from sklearn.base import ClassifierMixin
	from sklearn.preprocessing import LabelEncoder
	from sklearn.externals import six
	from sklearn.base import clone
	from sklearn.pipeline import _name_estimators
	from sklearn.metrics import roc_curve
	from sklearn.metrics import auc
	from sklearn.linear_model import LogisticRegression
	from sklearn.pipeline import Pipeline
except ImportError:
	print "This script requires sklearn to be installed!"

# pd.set_option('display.mpl_style', 'default') # jupyter

class Printer():
	"""Print things to stdout on one line dynamically"""
	def __init__(self,data):
		sys.stdout.write("\r\x1b[K"+data.__str__())
		sys.stdout.flush()

def creatematrix(features, kmer, args):
	stat = pd.read_csv(features, sep=",")
	kmer = pd.read_csv(kmer, sep="\t", header = None)
	kmer = kmer.iloc[:, :-1]
	id2 = stat.id.str.split("-",expand=True) # split the string to get label
	id2 = id2.iloc[:, :-1]
	stat2 = stat.iloc[:, 1:]
	df = pd.concat([stat2.reset_index(drop=True), kmer], axis=1) # concat kmerand stat matrix
	df = pd.concat([id2, df], axis=1)
	df.columns.values[0] = "label"
	# encoding class labels as integers
	df.loc[df.label == 'positive', 'label'] = 1
	df.loc[df.label == 'negative', 'label'] = 0
	# get number of instances per group
#	dat_pos = df.loc[df['label'] == 1]
#	dat_neg = df.loc[df['label'] == 0]
#	print('number of negative instances: %d ' % (dat_neg.shape[0]))
#	print('number of positive instances: %d ' % (dat_pos.shape[0]))
#
#	if(args.balance):
#		num = min(dat_neg.shape[0],dat_pos.shape[0])
#		print('limit dataset size to %d' % (num))
#		# generate a random subset of both with the size of $num
#		posrand = dat_pos.sample(n=num)
#		negrand = dat_neg.sample(n=num)
#		df = posrand.copy()
#		dat = df.append(negrand)

	y = df['label'].tolist() # extract label
	X = df.drop(df.columns[[0]], 1) # remove label
	return X, y

def drawroc(clf, clf_labels, X_train, y_train, X_test, y_test):
	""" draw a roc curve for each model in clf, save as roc.png"""
	colors = ['black', 'orange', 'blue', 'green']
	linestyles = [':', '--', '-.', '-']
	for clf, label, clr, ls \
									in zip(clf, clf_labels, colors, linestyles):
					scores = cross_val_score(estimator=clf,X=X_train,y=y_train,cv=int(args.cv), scoring='roc_auc', n_jobs=-1, verbose=3)
					print("ROC AUC: %0.2f (+/- %0.2f) [%s]" % (scores.mean(), scores.std(), label))
					y_pred = clf.fit(X_train, y_train).predict_proba(X_test)[:, 1]
					fpr, tpr, thresholds = roc_curve(y_true= y_test, y_score=y_pred)
					roc_auc = auc(x=fpr, y=tpr)
					plt.plot(fpr, tpr, color=clr, linestyle=ls, label='%s (auc = %0.2f)' % (label, roc_auc))
	plt.legend(loc='lower right')
	plt.plot([0, 1], [0, 1], linestyle='--', color='gray', linewidth=2)
	plt.xlim([-0.1, 1.1])
	plt.ylim([-0.1, 1.1])
	plt.grid()
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	# plt.tight_layout()
	plt.savefig('roc.png', dpi=300)
	plt.show()

class MajorityVoteClassifier(BaseEstimator, ClassifierMixin):
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


def balanced_subsample(x,y,subsample_size=1.0):

	class_xs = []
	min_elems = None

	for yi in np.unique(y):
		elems = x[(y == yi)]
		class_xs.append((yi, elems))
		if min_elems == None or elems.shape[0] < min_elems:
			min_elems = elems.shape[0]

	use_elems = min_elems
	if subsample_size < 1:
		use_elems = int(min_elems*subsample_size)

	xs = []
	ys = []

	for ci,this_xs in class_xs:
		if len(this_xs) > use_elems:
			this_xs = this_xs.reindex(np.random.permutation(this_xs.index))
			
		x_ = this_xs[:use_elems]
		y_ = np.empty(use_elems)
		y_.fill(ci)

		xs.append(x_)
		ys.append(y_)

	xs = pd.concat(xs)
	ys = pd.Series(data=np.concatenate(ys),name='target')

	return xs,ys

def report(results, n_top=3):
	for i in range(1, n_top + 1):
		candidates = np.flatnonzero(results['rank_test_score'] == i)
		for candidate in candidates:
			print("Model with rank: {0}".format(i))
			print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
				  results['mean_test_score'][candidate],
				  results['std_test_score'][candidate]))
			print("Parameters: {0}".format(results['params'][candidate]))
			print("")


def build_randomForest(X, y, args):
	Printer(colored('(training) ', 'green') + 'random forest')
	# specify parameters and distributions to sample from
	param_dist = {"max_depth": [5, 4, 3, None],
			  "n_estimators": [500, 2000],
			  "max_features": sp_randint(1, 50),
			  "min_samples_split": sp_randint(2, 50),
			  "min_samples_leaf": sp_randint(1, 50),
			  "bootstrap": [True, False],
			  "criterion": ["gini", "entropy"]}
	clf = RandomForestClassifier()
	random_search = RandomizedSearchCV(clf, param_distributions=param_dist, scoring = 'accuracy',
								   n_iter=args.iter, n_jobs=-1, refit=True)
	random_search.fit(X, y)
	acc = random_search.cv_results_['mean_test_score']
	filename= 'cv/randomforest_' + str(np.amax(acc))
	# save model
	with open(filename, 'wb') as fid:
		cPickle.dump(random_search, fid) 
	return random_search

def build_logisticregression(X, y, args):
	Printer(colored('(training) ', 'green') + 'logistic regression')
	# specify parameters and distributions to sample from
	param_dist = {"C": np.logspace(-9, 3, 13),
			  "solver": ['newton-cg', 'lbfgs', 'liblinear', 'sag'],
			  "dual":[False],
			  "tol": np.logspace(-9, 3, 13)
			}
	clf = LogisticRegression(penalty='l2')
	random_search = RandomizedSearchCV(clf, param_distributions=param_dist, scoring = 'accuracy',
								   n_iter=args.iter, n_jobs=-1, refit=True)
	random_search.fit(X, y)
#	report(random_search.cv_results_)

	acc = random_search.cv_results_['mean_test_score']
	filename= 'cv/logisticregression_' + str(np.amax(acc))
	# save model
	with open(filename, 'wb') as fid:
		cPickle.dump(random_search, fid) 
	return random_search

def build_svc(X, y, args):
	Printer(colored('(training) ', 'green') + 'SVC')
	# specify parameters and distributions to sample from
	param_dist = {'C': pow(2.0, np.arange(-10, 11, 0.1)),
	'gamma': pow(2.0, np.arange(-10, 11, 0.1)),'kernel': ['linear', 'rbf']}
	clf = SVC(probability=True)
	random_search = RandomizedSearchCV(clf, param_distributions=param_dist, scoring = 'accuracy',
								   n_iter=args.iter, n_jobs=-1, refit=True)
	random_search.fit(X, y)
	
	acc = random_search.cv_results_['mean_test_score']
	filename= 'cv/svc_' + str(np.amax(acc))
	# save model
	with open(filename, 'wb') as fid:
		cPickle.dump(random_search, fid) 
	return random_search
#	report(random_search.cv_results_)

def build_rvc(X, y, args):
	Printer(colored('(training) ', 'green') + 'RVC')
	# specify parameters and distributions to sample from
	param_dist = {'gamma': pow(2.0, np.arange(-10, 11, 0.1)),
'kernel': ['linear', 'rbf']}
	clf = RVC(kernel = 'rbf', gamma = 1)
	random_search = RandomizedSearchCV(clf, param_distributions=param_dist, scoring = 'accuracy',
								   n_iter=args.iter, n_jobs=-1, refit=True)
	random_search.fit(X, y)

#	report(random_search.cv_results_)

	acc = random_search.cv_results_['mean_test_score']
	filename= 'cv/rvc_' + str(np.amax(acc))
	# save model
	with open(filename, 'wb') as fid:
		cPickle.dump(random_search, fid) 
	return random_search
	
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-t', '--test_size', action='store', dest='test_size', help='size of test set from whole dataset in percent', default=30)
	parser.add_argument('-r', '--random_size', action='store', dest='random_size', help='size of balanced random subset of data in percent', default=10)
	parser.add_argument('-i', '--iterations', action='store', dest='iter', help='number of random iterations or hyperparameter optimization', default=10)
	parser.add_argument('-c', '--cv', action='store', dest='cv', help='cross validation size (e.g. 10 for 10-fold cross validation)', default=3)
	parser.add_argument('--lhs', dest='lhs', action='store_true', help='optimize parameters')
	parser.add_argument('--roc', dest='roc', action='store_true', help='plot ROC curve')
	parser.add_argument('--balance', dest='balance', action='store_true', help='balance dataset')
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	args = parser.parse_args()

	### load/preprocess data
	Printer(colored('(preprocessing) ', 'green') + 'import data')
	X, y = creatematrix('dat/train.features.clear2.csv', 'dat/train.features.kmer', args)

	# generate a random subset
	Printer(colored('(preprocessing) ', 'green') + 'generate a random subset')
	X_sub, y_sub = balanced_subsample(X, y, subsample_size=0.1)

	# split train/testset
	Printer(colored('(preprocessing) ', 'green') + 'generate train/test set')
	X_train, X_test, y_train, y_test = train_test_split(X_sub, y_sub, test_size=0.3)

	if not os.path.exists('cv'):
		os.makedirs('cv')

	### build model
	rf_model = build_randomForest(X_train, y_train, args)
	lg_model = build_logisticregression(X_train, y_train, args)
	svc_model = build_svc(X_train, y_train, args)
	rvc_model = build_rvc(X_train, y_train, args)

	### save model
#'	Printer(colored('(training) ', 'green') + 'save model')
#'	with open('model.pkl', 'wb') as fid:
#'		cPickle.dump(pipe, fid) 

	### draw ROC
#'	if(args.roc):
#'		Printer(colored('(training) ', 'green') + 'draw ROC curve')
#'		drawroc(all_clf, clf_labels, X_train, y_train, X_test, y_test)