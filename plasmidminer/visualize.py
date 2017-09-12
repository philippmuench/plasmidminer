
import os, glob, sys
import argparse
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
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from time import time
from sklearn.decomposition import PCA
import numpy as np
from numpy  import array
import matplotlib.pyplot as plt
from matplotlib import offsetbox
from sklearn import (manifold, datasets, decomposition, ensemble,
                     discriminant_analysis, random_projection)
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
	from sklearn.decomposition import PCA as sklearnPCA
	from sklearn.manifold import TSNE
except ImportError:
	print 'This script requires sklearn to be installed!'
import findparameters
import roc



if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--data_folder', action='store', dest='data',
						help='path to file foulder', default='dat')
	parser.add_argument('--version', action='version', version='%(prog)s 0.1')
	args = parser.parse_args()

	# print input data
	X, y = findparameters.creatematrix(args.data + '/train.features.clear2.csv', args.data +'/train.features.kmer', args)
	findparameters.showinput(X, y, 'imported data')

	X = X.as_matrix(columns=None)
	y = array(y)
	print(X)
	print(type(X))
	print(y)
	print(type(y))

	sc = StandardScaler()
	X = sc.fit_transform(X)
	pca = PCA(n_components=2)
	X_r = pca.fit(X).transform(X)

	X_tsne = TSNE(learning_rate=100).fit_transform(X)


	# Percentage of variance explained for each components
	print('explained variance ratio (first two components): %s'
		  % str(pca.explained_variance_ratio_))

	plt.figure()
	#plt.scatter(X_r[y == 1, 0], X_r[y == 1, 1], color='red', alpha=.5, lw=1,label='plasmid')
	#plt.scatter(X_r[y == 0, 0], X_r[y == 0, 1], color='blue', alpha=.5, lw=1,label='chromosome')
	#plt.title('PCA')
	#plt.legend(loc='best', shadow=False, scatterpoints=1)
	
	plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=y)
	plt.show()
