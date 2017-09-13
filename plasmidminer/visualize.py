
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
from sklearn.decomposition import TruncatedSVD
from sklearn.decomposition import PCA
from sklearn import manifold, datasets
from sklearn.random_projection import sparse_random_matrix
import numpy as np
from numpy  import array
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
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


class Printer():
	"""Print things to stdout on one line dynamically"""
	def __init__(self,data):
		sys.stdout.write("\r\x1b[K"+data.__str__())
		sys.stdout.flush()

def column(matrix, i):
    return [row[i] for row in matrix]

def density_estimation(m1, m2):
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--data_folder', action='store', dest='data',
						help='path to file foulder', default='dat')
	parser.add_argument('--version', action='version', version='%(prog)s 0.1')
	args = parser.parse_args()

	Printer(colored('(processing) ', 'green') + 'import data')
	# print input data
	X, y = findparameters.creatematrix(args.data + '/train.features.clear2.csv', args.data +'/train.features.kmer', args)
	findparameters.showinput(X, y, 'imported data')

	
	X = X.as_matrix(columns=None)
	y = array(y)

	Printer(colored('(processing) ', 'green') + 'scale data')
	# resaling
	sc = StandardScaler()
	X_sc = sc.fit_transform(X)
	np.savetxt('x_scaled.csv', X_sc, delimiter='\t')
	np.savetxt('y.csv', y, delimiter='\t')
	# data for http://projector.tensorflow.org/
	# plot PCA
	Printer(colored('(processing) ', 'green') + 'plot PCA')
	pca = PCA(n_components=2)
	X_r = pca.fit(X_sc).transform(X_sc)
	X_r_no_scale = pca.fit(X).transform(X)
	# Percentage of variance explained for each components
	print('explained variance ratio (first two components): %s'
		  % str(pca.explained_variance_ratio_))

	plt.figure()
	plt.scatter(X_r[y == 1, 0], X_r[y == 1, 1], color='red', edgecolor='', alpha=.5, label='plasmid', marker='o')
	plt.scatter(X_r[y == 0, 0], X_r[y == 0, 1], color='blue', edgecolor='', alpha=.5, label='chromosome', marker='o')
	plt.title('PCA')
	plt.legend(loc='best', shadow=False, scatterpoints=1)
	plt.savefig('pca.pdf')
	plt.close()

	plt.figure()
	plt.scatter(X_r_no_scale[y == 1, 0], X_r_no_scale[y == 1, 1], color='red', edgecolor='', alpha=.5, label='plasmid', marker='o')
	plt.scatter(X_r_no_scale[y == 0, 0], X_r_no_scale[y == 0, 1], color='blue', edgecolor='', alpha=.5, label='chromosome', marker='o')
	plt.title('PCA - no scale')
	plt.legend(loc='best', shadow=False, scatterpoints=1)
	plt.savefig('pca_unscaled.pdf')
	plt.close()

	perplexity = 50
	Printer(colored('(processing) ', 'green') + 'plot tSNE')
	# plot tSNE
	marker = y
	marker = np.where(marker > 0, '+', '.')

	cols = array(column(X, 5))	
	scaled_z = (cols - cols.min()) / cols.ptp()
	colors = plt.cm.coolwarm(scaled_z)
	col_2 = np.where(y > 0, 'black', 'red')
	print(col_2)


#	lda_model = lda.LDA(n_topics=20, n_iter=500)
#	X_topics = lda_model.fit_transform(X_sc)

	X_tsne = TSNE(n_components = 2, perplexity=perplexity, init='pca', verbose=1).fit_transform(X_sc)

	# color by label
	plt.figure()
	plt.scatter( X_tsne[:, 0], X_tsne[:, 1], marker='o', c=col_2, edgecolor='', alpha = 0.2)
	plt.title('tSNE')
	plt.savefig('tsne_binary_colored.pdf')

	plt.figure()
	for _col, _m, _x, _y in zip(colors, marker, X_tsne[:, 0], X_tsne[:, 1]):
		plt.scatter(_x, _y, marker=_m, c=_col, edgecolor='')
	plt.title('tSNE')
	plt.savefig('tsne_4_colored.pdf')
	plt.close()



	# color
#	plt.figure()
#	plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=array(column(X, 3)), alpha=0.1, marker='o', edgecolor='')
#	filename = 'tsne_by_3.pdf'
#	plt.savefig(filename, dpi=300)
#	plt.close()

	# color
#	plt.figure()
#	plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=array(column(X, 4)), alpha=0.1, marker='o', edgecolor='')
#	filename = 'tsne_by_4.pdf'
#	plt.savefig(filename, dpi=300)
#	plt.close()

	# color 
#	plt.figure()
#	plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=array(column(X, 5)), alpha=0.1, marker='o', edgecolor='')
#	plt.savefig('tsne_by_5.png', dpi=300)
#	plt.close()

