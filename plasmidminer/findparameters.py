#!/usr/bin/python
# script to get individual model
# philipp.muench@helmholtz-hzi.de
from scipy.stats import randint as sp_randint
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cPickle
import pickle
import scipy
import collections
from scipy import stats
import argparse
from scipy.stats import expon
from termcolor import colored
import sobol_seq
from skbayes.rvm_ard_models import RVR,RVC
try:
	from sklearn.externals import joblib
	from sklearn.pipeline import Pipeline
	from sklearn.model_selection import cross_validate
	from sklearn import preprocessing
	from sklearn.preprocessing import StandardScaler
	from sklearn.ensemble import VotingClassifier
	from sklearn.preprocessing import MaxAbsScaler
	from sklearn.model_selection import RandomizedSearchCV
	from sklearn.ensemble import RandomForestClassifier
	from sklearn.ensemble import GradientBoostingClassifier
	from sklearn.linear_model import LogisticRegression
	from sklearn.svm import SVC
	from sklearn.model_selection import cross_val_score, train_test_split
	from sklearn.metrics import roc_curve
	from sklearn.metrics import auc
except ImportError:
	print 'This script requires sklearn to be installed!'

def loaddataset(filename):
	"""Loads feature matrix from a msgpack obj"""
	Printer(colored('(preprocessing) ', 'green') + 'import data (from binary file)')
	X, y= pd.read_msgpack('dat/dataset.msg')
	return(X, y)

def creatematrix(features, kmer, args):
	"""Generates feature matrix from raw data"""
	Printer(colored('(preprocessing) ', 'green') + 'import data (from raw filess)')
	stat = pd.read_csv(features, sep=",")
	kmer = pd.read_csv(kmer, sep="\t", header=None)
	kmer = kmer.iloc[:, :-1]
	id2 = stat.id.str.split("-", expand=True)  # split the string to get label
	id2 = id2.iloc[:, :-1]
	stat2 = stat.iloc[:, 1:]
	df = pd.concat([stat2.reset_index(drop=True), kmer],axis=1)  # concat kmerand stat matrix
	df = pd.concat([id2, df], axis=1)
	df.columns.values[0] = "label"
	# encoding class labels as integers
	df.loc[df.label == 'positive', 'label'] = 1
	df.loc[df.label == 'negative', 'label'] = 0
	# get number of instances per group
	y = df['label'].tolist()  # extract label
	X = df.drop(df.columns[[0]], 1)  # remove label
	return X, y

def savemodel(model, filename):
	"""Saves the model to the cv/ folder"""
	print('saved model to ' + filename)
	with open(filename, 'wb') as fid:
		joblib.dump(model, fid, compress=9)

def saveparams(model_best_params, filename):
	with open(filename, 'wb') as fid:
		joblib.dump(model_best_params, fid, compress=9)

class Printer():
	"""Print things to stdout on one line dynamically"""
	def __init__(self, data):
		sys.stdout.write("\r\x1b[K" + data.__str__())
		sys.stdout.flush()

def drawroc(clf, clf_labels, X_train, y_train, X_test, y_test):
	"""draw a roc curve for each model in clf, save as roc.png"""
	colors = ['black', 'orange', 'blue', 'green']
	linestyles = [':', '--', '-.', '-']
	for clf, label, clr, ls \
			in zip(clf, clf_labels, colors, linestyles):
		scores = cross_val_score(estimator=clf, X=X_train, y=y_train, cv=int(
			args.cv), scoring='roc_auc', n_jobs=-1)
	#	print("ROC AUC: %0.2f (+/- %0.2f) [%s]" %
	#		  (scores.mean(), scores.std(), label))
		y_pred = clf.fit(X_train, y_train).predict_proba(X_test)[:, 1]
		fpr, tpr, thresholds = roc_curve(y_true=y_test, y_score=y_pred)
		roc_auc = auc(x=fpr, y=tpr)
		plt.plot(fpr, tpr, color=clr, linestyle=ls,
				 label='%s (auc = %0.2f)' % (label, roc_auc))
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

def balanced_subsample(x, y, subsample_size=1):
	"""Balances the dataset in  a 1:1 fashion"""
	class_xs = []
	min_elems = None
	for yi in np.unique(y):
		elems = x[(y == yi)]
		class_xs.append((yi, elems))
		if min_elems == None or elems.shape[0] < min_elems:
			min_elems = elems.shape[0]
	use_elems = min_elems
	if subsample_size < 1:
		use_elems = int(min_elems * subsample_size)
	xs = []
	ys = []
	for ci, this_xs in class_xs:
		if len(this_xs) > use_elems:
			this_xs = this_xs.reindex(np.random.permutation(this_xs.index))
		x_ = this_xs[:use_elems]
		y_ = np.empty(use_elems)
		y_.fill(ci)
		xs.append(x_)
		ys.append(y_)
	xs = pd.concat(xs)
	ys = pd.Series(data=np.concatenate(ys), name='target')
	return xs, ys

def report(results, n_top=3):
	"""Prints the model validation scores and best parameter to screen"""
	for i in range(1, n_top + 1):
		candidates = np.flatnonzero(results['rank_test_score'] == i)
		for candidate in candidates:
			print("Model with rank: {0}".format(i))
			print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
				  results['mean_test_score'][candidate],
				  results['std_test_score'][candidate]))
			print("Parameters: {0}".format(results['params'][candidate]))
			print("")

def build_randomForest(X_loc, y_loc, args):
	"""finds best parameters for random forest"""
	Printer(colored('(training) ', 'green') +
			'searching for best parameters for random forest')
	# specify parameters and distributions to sample from
	param_dist = {"clf__max_depth": [100, 50, 20 , 10, 5, 4, 3, 2, None],
				  "clf__max_features": sp_randint(1, 50),
				  "clf__min_samples_split": sp_randint(2, 100),
				  "clf__min_samples_leaf": sp_randint(1, 100),
				  "clf__bootstrap": [True, False],
				  "clf__criterion": ["gini", "entropy"]}
	clf = RandomForestClassifier(n_estimators = 2000)
	pipe = Pipeline([['sc', MaxAbsScaler()],['clf', clf]])
	random_search = RandomizedSearchCV(pipe, param_distributions=param_dist, scoring='accuracy', n_iter=int(args.iter), n_jobs=-1, refit=True, cv=3)
	random_search.fit(X_loc, y_loc)
	acc = random_search.cv_results_['mean_test_score']
	filename = 'cv/randomforest_' + str(np.mean(acc)) + '.pkl'
	# save model
	savemodel(random_search, filename)
	# save best params
	filename_param = 'cv/randomforest_param_' + str(np.mean(acc)) + '.pkl'
	saveparams(random_search.best_params_, filename_param)
	return random_search, acc

def build_logisticregression(X_loc, y_loc, args):
	"""finds best parameters for logistic regression"""
	Printer(colored('(training) ', 'green') +
			'searching for best parameters for logistic regression')
	# specify parameters and distributions to sample from
	param_dist = {"C": np.logspace(-9, 3, 13),
				  "solver": ['newton-cg', 'lbfgs', 'liblinear', 'sag'],
				  "dual": [False],
				  "tol": np.logspace(-9, 3, 13)
				  }
	clf = LogisticRegression(penalty='l2')
	random_search = RandomizedSearchCV(clf, param_distributions=param_dist, scoring='accuracy', n_iter=int(args.iter), n_jobs=-1, refit=True, cv=3)
	random_search.fit(X_loc, y_loc)
	acc = random_search.cv_results_['mean_test_score']
	filename = 'cv/logisticregression_' + str(np.mean(acc)) + '.pkl'
	# save model
	savemodel(random_search, filename)
	# save best params
	filename_param = 'cv/logisticregression_param_' + str(np.mean(acc)) + '.json'
	saveparams(random_search.best_params_, filename_param)
	return random_search

def build_svc(X_loc, y_loc, args):
	"""finds best parameters for support vector machine"""
	Printer(colored('(training) ', 'green') +
			'searching for best parameters for SVC')
	# specify parameters and distributions to sample from
	if (args.sobol):
		param_dist = {'clf__C': sobol_seq.i4_sobol_generate(1, int(args.iter)) * 2** (15),
				  'clf__gamma': sobol_seq.i4_sobol_generate(1, int(args.iter)), 'clf__kernel': ['linear', 'rbf', 'poly'], 'clf__degree': [1,2,4,6,8]}
	else:
		param_dist = {'clf__C': scipy.stats.expon(scale=100),
		'clf__gamma': pow(2.0, np.arange(-10, 11, 0.1)), 'clf__kernel': ['linear', 'rbf', 'poly'], 'clf__degree': [1,2,4,6,8]}
	clf = SVC(probability=True)
	pipe = Pipeline([['sc', MaxAbsScaler()],['clf', clf]])
	random_search = RandomizedSearchCV(pipe, param_distributions=param_dist, scoring='accuracy', n_iter=int(args.iter), n_jobs=-1, refit=True, cv=3)
	random_search.fit(X_loc, y_loc)
	acc = random_search.cv_results_['mean_test_score']
	# save model
	filename = 'cv/svc_' + str(np.mean(acc)) + '.pkl'
	savemodel(random_search, filename)
	# save best params
	filename_param = 'cv/svc_param_' + str(np.mean(acc)) + '.pkl'
	saveparams(random_search.best_params_, filename_param)
	return random_search, acc
#   report(random_search.cv_results_)

def build_rvc(X_loc, y_loc, args):
	"""finds best parameters for relevance vector machine"""
	Printer(colored('(training) ', 'green') +
			'searching for best parameters for RVC')
	# specify parameters and distributions to sample from
	if (args.sobol):
		param_dist = {'clf__gamma': sobol_seq.i4_sobol_generate(1, int(args.iter)),
		'clf__kernel': ['linear', 'rbf', 'poly'],
		'clf__degree': [1,2,4,6,8]}
	else:
		param_dist = {'clf__gamma': pow(2.0, np.arange(-10, 11, 0.1)),
		'clf__kernel': ['linear', 'rbf', 'poly'],
		'clf__degree': [1,2,4,6,8]}
	clf = RVC()
	pipe = Pipeline([['sc', MaxAbsScaler()],['clf', clf]])
	random_search = RandomizedSearchCV(pipe, param_distributions=param_dist, scoring='accuracy', n_iter=int(args.iter), n_jobs=-1, refit=True, cv=3)
	random_search.fit(X_loc, y_loc)
	acc = random_search.cv_results_['mean_test_score']
	filename = 'cv/rvc_' + str(np.mean(acc)) + '.pkl'
	# save model
	savemodel(random_search, filename)
	# save best params
	filename_param = 'cv/rvc_param_' + str(np.mean(acc)) + '.pkl'
	saveparams(random_search.best_params_, filename_param)
	return random_search, acc

def build_gbc(X_loc, y_loc, args):
	"""finds best parameters for gradient boosting classifier"""
	Printer(colored('(training) ', 'green') +
			'searching for best parameters for GBC forest')
	# specify parameters and distributions to sample from
	param_dist = {"clf__n_estimators": [100, 500, 1000],
	'clf__learning_rate': [0.01, 0.1, 0.1, 0.5, 1.0],
	'clf__max_depth':[1, 3, 5, 7, 9],
	'clf__max_features':['auto', 'sqrt', 'log2'],
	'clf__loss':['deviance','exponential'],
	'clf__criterion':['friedman_mse', 'mae', 'mse']
	}
	clf = GradientBoostingClassifier()
	pipe = Pipeline([['sc', MaxAbsScaler()],['clf', clf]])
	random_search = RandomizedSearchCV(pipe, param_distributions=param_dist, scoring='accuracy', n_iter=int(args.iter), n_jobs=-1, refit=True, cv=3)
	random_search.fit(X_loc, y_loc)
	acc = random_search.cv_results_['mean_test_score']
	filename = 'cv/gbc_' + str(np.mean(acc)) + '.pkl'
	# save model
	savemodel(random_search, filename)
	# save best params
	filename_param = 'cv/gcb_param_' + str(np.mean(acc)) + '.pkl'
	saveparams(random_search.best_params_, filename_param)
	return random_search, acc

def build_voting(rf_model, svc_model, rvc_model, gbc_model, X_loc, y_loc, args):
	Printer(colored('(training) ', 'green') +
		'searching for voting sheme')
	clf = VotingClassifier(estimators=[('rf', rf_model), ('svc', svc_model),
	 ('rvc', rvc_model), ('gbc', gbc_model)], voting='hard')
	scoring = {'accuracy': 'accuracy'}

	scores = cross_validate(clf, X_loc, y_loc, scoring=scoring, cv=3, return_train_score=False)
	acc = scores["test_accuracy"].mean()
	filename = 'cv/voting_' + str(acc) + '.pkl'
	# save model
	savemodel(clf, filename)
	return clf, acc


def showinput(X_loc,y_loc,string):
	"""outputs basic statistics of input files"""
	print('\n-------------------------------------------------------------')
	print(str(string))
	print "number of instances:\t", X_loc.shape[0]
	print "number of features:\t", X_loc.shape[1]
	print "response:\t", collections.Counter(y_loc)
	print('-------------------------------------------------------------\n')

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--dataset', action='store', dest='dataset',
						help='path to data folder', default='dat')
	parser.add_argument('-i', '--iterations', action='store', dest='iter',
						help='number of random iterationsfor hyperparameter optimization', default=50)
	parser.add_argument('-c', '--cv', action='store', dest='cv',
						help='cross validation size (e.g. 10 for 10-fold cross validation)', default=3)
	# binary
	parser.add_argument('--sobol', dest='sobol',
						action='store_true', help='use sobol sequence for random search')
	parser.add_argument('--version', action='version', version='%(prog)s 0.2')
	args = parser.parse_args()

	# load data from msgpack object
#	X, y = loaddataset(args.dataset)

	X, y = creatematrix(args.dataset + '/train.features.clear2.csv', args.dataset +'/train.features.kmer', args)

	# print input data
	showinput(X, y, 'imported data')

	# generate a random subset
	Printer(colored('(preprocessing) ', 'green') + 'generate a random balanced subset')

	X_sub, y_sub = balanced_subsample(X, y)

	# split train/testset
	#Printer(colored('(preprocessing) ', 'green') + 'generate train/test set')
	#X_train, X_test, y_train, y_test = train_test_split(X_sub, y_sub, test_size=0.33)

	# print input data
	showinput(X_sub, y_sub, 'training data')

	# print input data
	#showinput(X_test, y_test, 'testing data')

	# create output folder
	if not os.path.exists('cv'):
		os.makedirs('cv')

	# optimize hyperparameters model
	rf_model, rf_acc = build_randomForest(X_sub, y_sub, args)
	svc_model, svc_acc = build_svc(X_sub, y_sub, args)
	rvc_model, rvc_acc = build_rvc(X_sub, y_sub, args)
	gbc_model, gbc_acc = build_gbc(X_sub, y_sub, args)

	# create voting classifier
	voting_model, voting_acc = build_voting(rf_model, svc_model, rvc_model, 
		gbc_model, X_sub, y_sub, args)

	print('\n')
	print('--------------------------------------------')
	if 'svc_model' in locals():
		print('SVC accuracy: %0.2f' % np.mean(svc_acc))
	if 'rvc_model' in locals():
		print('RVC accuracy: %0.2f' % np.mean(rvc_acc))
	if 'rf_model' in locals():
		print('RF accuracy: %0.2f' % np.mean(rf_acc))
	if 'gbc_model' in locals():
		print('GBC accuracy: %0.2f' % np.mean(gbc_acc))
	if 'voting_model' in locals():
		print('voting accuracy: %0.2f' % np.mean(voting_acc))
	print('--------------------------------------------')
