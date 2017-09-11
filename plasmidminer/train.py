import os, glob, sys
import argparse
#import cPickle
from termcolor import colored
import matplotlib.pyplot as plt
import numpy as np
import gzip
from sklearn import metrics
from sklearn.externals import joblib
from sklearn.model_selection import cross_val_score, train_test_split
from skbayes.rvm_ard_models import RVR,RVC
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import LinearSVC
from sklearn.model_selection import learning_curve
from sklearn.model_selection import ShuffleSplit
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate
try:
	from sklearn.externals import joblib
	from sklearn.pipeline import Pipeline
	from sklearn import preprocessing
	from sklearn.preprocessing import StandardScaler
	from sklearn.preprocessing import MaxAbsScaler
	from sklearn.model_selection import RandomizedSearchCV
	from sklearn.ensemble import RandomForestClassifier
	from sklearn.ensemble import GradientBoostingClassifier
	from sklearn.linear_model import LogisticRegression
	from sklearn.svm import SVC
	from sklearn.model_selection import cross_val_score, train_test_split
	from sklearn.metrics import roc_curve
	from sklearn.metrics import auc
	from sklearn.metrics import (brier_score_loss, precision_score, recall_score,
                             f1_score)
	from sklearn.calibration import CalibratedClassifierCV, calibration_curve
	from sklearn.model_selection import cross_validate
except ImportError:
	print 'This script requires sklearn to be installed!'
import findparameters
import roc

class Printer():
	"""Print things to stdout on one line dynamically"""
	def __init__(self,data):
		sys.stdout.write("\r\x1b[K"+data.__str__())
		sys.stdout.flush()

def plot_learning_curve(estimator, title, X, y, ylim=None, cv=None,
                        n_jobs=-1, train_sizes=np.linspace(.1, 1.0, 5)):
    """
    Generate a simple plot of the test and training learning curve.

    Parameters
    ----------
    estimator : object type that implements the "fit" and "predict" methods
        An object of that type which is cloned for each validation.

    title : string
        Title for the chart.

    X : array-like, shape (n_samples, n_features)
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.

    y : array-like, shape (n_samples) or (n_samples, n_features), optional
        Target relative to X for classification or regression;
        None for unsupervised learning.

    ylim : tuple, shape (ymin, ymax), optional
        Defines minimum and maximum yvalues plotted.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
          - None, to use the default 3-fold cross-validation,
          - integer, to specify the number of folds.
          - An object to be used as a cross-validation generator.
          - An iterable yielding train/test splits.

        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`StratifiedKFold` used. If the estimator is not a classifier
        or if ``y`` is neither binary nor multiclass, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validators that can be used here.

    n_jobs : integer, optional
        Number of jobs to run in parallel (default 1).
    """
    plt.figure()
    plt.title(title)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=cv, n_jobs=-1, train_sizes=train_sizes)
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    plt.grid()

    plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1,
                     color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r",
             label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g",
             label="Cross-validation score")

    plt.legend(loc="best")
    plt.savefig('evaluation/learning_curve.png')   # save the figure to file
    plt.close()

def drawsingleroc(clf, clf_label, X_train, y_train, X_test, y_test):
	"""draw a roc curve for each model in clf, save as roc.png"""
	scores = cross_val_score(estimator=clf, X=X_train, y=y_train, cv=int(
		args.cv), scoring='roc_auc', n_jobs=-1)
	#print("\nROC AUC: %0.2f (+/- %0.2f) [%s]" %
#		(scores.mean(), scores.std(), clf_label))
	y_pred = clf.fit(X_train, y_train).predict_proba(X_test)[:, 1]
	fpr, tpr, thresholds = metrics.roc_curve(y_true=y_test, y_score=y_pred)
	roc_auc = metrics.auc(x=fpr, y=tpr)
	plt.plot(fpr, tpr, color='black', linestyle='-',
			label='%s (auc = %0.2f)' % (clf_label, roc_auc))
	plt.legend(loc='lower right')
	plt.plot([0, 1], [0, 1], linestyle='--', color='gray', linewidth=2)
	plt.xlim([-0.1, 1.1])
	plt.ylim([-0.1, 1.1])
	plt.grid()
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	filename = 'evaluation/roc_auc_' + str(np.amax(roc_auc)) + '.png'
	# plt.tight_layout()
	plt.savefig(filename, dpi=300)
	plt.close()

def evaluate_randomForest(X, y, X_val, y_val, args, param_dist):
	"""finds best parameters for random forest"""
	clf = RandomForestClassifier(n_estimators = 2000)
	pipe = Pipeline([['sc', MaxAbsScaler()],['clf', clf]])
	scoring = {'accuracy': 'accuracy', 'precision': 'precision', 'recall': 'recall',
	 'roc_auc' : 'roc_auc', 'average_precision' : 'average_precision', 'f1' : 'f1',
	 'f1_micro' : 'f1_micro','f1_macro' : 'f1_macro','f1_weighted' : 'f1_weighted'}
	scores = cross_validate(pipe, X, y, scoring=scoring, cv=3, return_train_score=False)
	label = 'final model'
	print("\ntest_accuracy: %0.2f (+/- %0.2f) [%s]" %
		(scores["test_accuracy"].mean(), scores["test_accuracy"].std(), label))
	print("test_recall: %0.2f (+/- %0.2f) [%s]" %
		(scores["test_recall"].mean(), scores["test_recall"].std(), label))
	print("test_precision: %0.2f (+/- %0.2f) [%s]" %
		(scores["test_precision"].mean(), scores["test_precision"].std(), label))
	print("test_roc_auc: %0.2f (+/- %0.2f) [%s]" %
		(scores["test_roc_auc"].mean(), scores["test_roc_auc"].std(), label))
	print("test_average_precision: %0.2f (+/- %0.2f) [%s]" %
		(scores["test_average_precision"].mean(), scores["test_average_precision"].std(), label))
	print("test_f1: %0.2f (+/- %0.2f) [%s]" %
		(scores["test_f1"].mean(), scores["test_f1"].std(), label))
	print("test_f1_micro: %0.2f (+/- %0.2f) [%s]" %
		(scores["test_f1_micro"].mean(), scores["test_f1_micro"].std(), label))
	print("test_f1_macro: %0.2f (+/- %0.2f) [%s]" %
		(scores["test_f1_macro"].mean(), scores["test_f1_macro"].std(), label))
	print("test_f1_weighted: %0.2f (+/- %0.2f) [%s]" %
		(scores["test_f1_weighted"].mean(), scores["test_f1_weighted"].std(), label))

	filename = 'evaluation/randomforest_final_' + str(scores["test_accuracy"].mean()) + '.pkl'
	findparameters.savemodel(pipe, filename)
	return pipe

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--model_param_file', action='store', dest='model_param',
						help='path to file where .pkl models is located')
	parser.add_argument('-d', '--data_folder', action='store', dest='data',
						help='path to file foulder', default='dat')
	parser.add_argument('-t', '--test_size', action='store', dest='test_size',
						help='size of test set from whole dataset in percent', default=33)
	parser.add_argument('-c', '--cv', action='store', dest='cv',
						help='cross validation size (e.g. 10 for 10-fold cross validation)', default=3)
	parser.add_argument('--version', action='version', version='%(prog)s 0.1')
	args = parser.parse_args()

	# setting up log file
	old_stdout = sys.stdout
	log_file = open("evaluation/log.txt", "w")

	# print input data
	X, y = findparameters.creatematrix(args.data + '/train.features.clear2.csv', args.data +'/train.features.kmer', args)
	sys.stdout = log_file

	findparameters.showinput(X, y, 'imported data')
	
	# generate a random subset
	sys.stdout = old_stdout
	Printer(colored('(preprocessing) ', 'green') + 'generate a random balanced subset')
	sys.stdout = log_file

	# generate subsample with equal amount of pos/neg instances
	X_sub, y_sub = findparameters.balanced_subsample(X, y)

	# split train/testset
	sys.stdout = old_stdout
	Printer(colored('(preprocessing) ', 'green') + 'generate train/test set')
	sys.stdout = log_file

	X_train, X_val, y_train, y_val = train_test_split(X_sub, y_sub, test_size=0.33)

	# print input data
	findparameters.showinput(X_train, y_train, 'training data')

	# print input data
	findparameters.showinput(X_val, y_val, 'validation data')

	# create output folder
	if not os.path.exists('evaluation'):
		os.makedirs('evaluation')

	# load parameters
	sys.stdout = old_stdout
	Printer(colored('(processing) ', 'green') + 'load model' )
	sys.stdout = log_file

	param_dist = joblib.load('cv/randomforest_param_0.774145616642.pkl')

	# tain model
	sys.stdout = old_stdout
	Printer(colored('(processing) ', 'green') + 'train model' )
	sys.stdout = log_file

	estimator = evaluate_randomForest(X_train, y_train, X_val, y_val, args, param_dist)

	# plot learning curve
	sys.stdout = old_stdout
	Printer(colored('(processing) ', 'green') + 'plot learning curve' )
	sys.stdout = log_file
	cv_shuffle = ShuffleSplit(n_splits=5, test_size=0.2, random_state=0)
	plot_learning_curve(estimator, 'Learning Curves', X, y, ylim=(0.7, 1.01), cv=cv_shuffle, n_jobs=4)

	# draw ROC curve
	sys.stdout = old_stdout
	Printer(colored('(processing) ', 'green') + 'draw ROC curve' )
	sys.stdout = log_file

	drawsingleroc(estimator,'final model', X_train, y_train, X_val, y_val)
	log_file.close()