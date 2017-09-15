
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
	parser.add_argument('--pca', dest='pca', action='store_true', help='plot PCA')
	parser.add_argument('--tsne', dest='tsne', action='store_true', help='plot tSNE')

	parser.add_argument('--version', action='version', version='%(prog)s 0.1')
	args = parser.parse_args()

	# create output folder
	if not os.path.exists('visualization'):
		os.makedirs('visualization')

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
	np.savetxt('visualization/x_scaled.csv', X_sc, delimiter='\t')
	np.savetxt('visualization/y.csv', y, delimiter='\t')

	marker = y
	marker = np.where(marker > 0, '+', '.')

	cols_gc = array(column(X, 16))	
	cols_met = array(column(X, 18))	

	scaled_z_gc = (cols_gc - cols_gc.min()) / cols_gc.ptp()
	colors_gc = plt.cm.coolwarm(scaled_z_gc)

	scaled_z_met = (cols_met - cols_met.min()) / cols_met.ptp()
	colors_met = plt.cm.coolwarm(scaled_z_met)
	col_2 = np.where(y > 0, 'black', 'red')

	if (args.pca):
		Printer(colored('(processing) ', 'green') + 'plot PCA')
		pca = PCA(n_components=2)
		X_r = pca.fit(X_sc).transform(X_sc)
		print('scaled - explained variance ratio (first two components): %s'
			  % str(pca.explained_variance_ratio_))
		X_r_no_scale = pca.fit(X).transform(X)
		# Percentage of variance explained for each components
		print('unscaled - explained variance ratio (first two components): %s'
			  % str(pca.explained_variance_ratio_))

		plt.figure()
		sc = plt.scatter(X_r[y == 1, 0], X_r[y == 1, 1], marker='.', c=colors_gc, edgecolor='')
		plt.savefig('visualization/pca_gc.pdf', transparent=True)
		plt.close()

#		plt.figure()
#		for _lab, _col, _x, _y in zip(y, colors_gc, X_r[y == 1, 0], X_r[y == 1, 1]):
#			plt.scatter(_x, _y, marker='.', c=_col, edgecolor='')
#		plt.colorbar(colors_gc)
#		plt.title('PCA colored by GC content')
#		plt.savefig('visualization/pca_gc.pdf', transparent=True)
#		plt.close()

	#	plt.figure()
#		for _col, _x, _y in zip(colors_met, X_r[y == 1, 0], X_r[y == 1, 1]):
#			plt.scatter(_x, _y, marker='.', c=_col, edgecolor='')
#		plt.title('PCA colored by GpG content')
#		plt.colorbar(colors_met)
#		plt.savefig('visualization/pca_met.pdf', transparent=True)
#		plt.close()

		plt.figure()
		plt.scatter(X_r[y == 1, 0], X_r[y == 1, 1], color='red', edgecolor='', alpha=.5, label='plasmid', marker='o')
		plt.scatter(X_r[y == 0, 0], X_r[y == 0, 1], color='blue', edgecolor='', alpha=.5, label='chromosome', marker='o')
		plt.title('PCA')
		plt.legend(loc='best', shadow=False, scatterpoints=1)
		plt.savefig('visualization/pca.pdf')
		plt.close()

		plt.figure()
		plt.scatter(X_r[y == 1, 0], X_r[y == 1, 1], color='red', edgecolor='', alpha=.2, label='plasmid', marker='o')
		plt.scatter(X_r[y == 0, 0], X_r[y == 0, 1], color='blue', edgecolor='', alpha=.2, label='chromosome', marker='o')
		plt.title('PCA')
		plt.legend(loc='best', shadow=False, scatterpoints=1)
		plt.savefig('visualization/pca_2.pdf')
		plt.close()

		plt.figure()
		plt.scatter(X_r[y == 1, 0], X_r[y == 1, 1], color='red', edgecolor='', alpha=.2, label='plasmid', marker='.')
		plt.scatter(X_r[y == 0, 0], X_r[y == 0, 1], color='blue', edgecolor='', alpha=.2, label='chromosome', marker='.')
		plt.title('PCA')
		plt.legend(loc='best', shadow=False, scatterpoints=1)
		plt.savefig('visualization/pca_3.pdf')
		plt.close()

		plt.figure()
		plt.scatter(X_r_no_scale[y == 1, 0], X_r_no_scale[y == 1, 1], color='red', edgecolor='', alpha=.5, label='plasmid', marker='o')
		plt.scatter(X_r_no_scale[y == 0, 0], X_r_no_scale[y == 0, 1], color='blue', edgecolor='', alpha=.5, label='chromosome', marker='o')
		plt.title('PCA - no scale')
		plt.legend(loc='best', shadow=False, scatterpoints=1)
		plt.savefig('visualization/pca_unscaled.pdf')
		plt.close()

	if (args.tsne):
		Printer(colored('(processing) ', 'green') + 'plot tSNE (with perplexity set to 10)')
		X_tsne_10 = TSNE(n_components = 2, perplexity=10, init='pca', verbose=1).fit_transform(X_sc)
		Printer(colored('(processing) ', 'green') + 'plot tSNE (with perplexity set to 25)')
		X_tsne_25 = TSNE(n_components = 2, perplexity=25, init='pca', verbose=1).fit_transform(X_sc)
		Printer(colored('(processing) ', 'green') + 'plot tSNE (with perplexity set to 50)')
		X_tsne_50 = TSNE(n_components = 2, perplexity=50, init='pca', verbose=1).fit_transform(X_sc)

		# 10
		plt.figure()
		plt.scatter( X_tsne_10[:, 0], X_tsne_10[:, 1], marker='o', c=col_2, edgecolor='', alpha = 0.2)
		plt.title('tSNE (perplexity 10)')
		plt.savefig('visualization/tsne_10_binary_colored.pdf')

		plt.figure()
		for _col, _m, _x, _y in zip(colors_gc, marker, X_tsne_10[:, 0], X_tsne_10[:, 1]):
			plt.scatter(_x, _y, marker=_m, c=_col, edgecolor='')
		plt.title('tSNE (perplexity 10)')
		plt.savefig('visualization/tsne_10_gc_colored.pdf')
		plt.close()

		plt.figure()
		for _col, _m, _x, _y in zip(colors_met, marker, X_tsne_10[:, 0], X_tsne_10[:, 1]):
			plt.scatter(_x, _y, marker=_m, c=_col, edgecolor='')
		plt.title('tSNE (perplexity 10)')
		plt.savefig('visualization/tsne_10_met_colored.pdf')
		plt.close()

		# 25
		plt.figure()
		plt.scatter( X_tsne_25[:, 0], X_tsne_25[:, 1], marker='o', c=col_2, edgecolor='', alpha = 0.2)
		plt.title('tSNE (perplexity 25)')
		plt.savefig('visualization/tsne_25_binary_colored.pdf')

		plt.figure()
		for _col, _m, _x, _y in zip(colors_gc, marker, X_tsne_25[:, 0], X_tsne_25[:, 1]):
			plt.scatter(_x, _y, marker=_m, c=_col, edgecolor='')
		plt.title('tSNE (perplexity 25)')
		plt.savefig('visualization/tsne_25_gc_colored.pdf')
		plt.close()

		plt.figure()
		for _col, _m, _x, _y in zip(colors_met, marker, X_tsne_25[:, 0], X_tsne_25[:, 1]):
			plt.scatter(_x, _y, marker=_m, c=_col, edgecolor='')
		plt.title('tSNE (perplexity 25)')
		plt.savefig('visualization/tsne_25_met_colored.pdf')
		plt.close()

		# 50
		plt.figure()
		plt.scatter( X_tsne_50[:, 0], X_tsne_50[:, 1], marker='o', c=col_2, edgecolor='', alpha = 0.2)
		plt.title('tSNE (perplexity 50)')
		plt.savefig('visualization/tsne_50_binary_colored.pdf')

		plt.figure()
		for _col, _m, _x, _y in zip(colors_gc, marker, X_tsne_50[:, 0], X_tsne_50[:, 1]):
			plt.scatter(_x, _y, marker=_m, c=_col, edgecolor='')
		plt.title('tSNE (perplexity 50)')
		plt.savefig('visualization/tsne_50_gc_colored.pdf')
		plt.close()

		plt.figure()
		for _col, _m, _x, _y in zip(colors_met, marker, X_tsne_50[:, 0], X_tsne_50[:, 1]):
			plt.scatter(_x, _y, marker=_m, c=_col, edgecolor='')
		plt.title('tSNE (perplexity 50)')
		plt.savefig('visualization/tsne_50_met_colored.pdf')
		plt.close()