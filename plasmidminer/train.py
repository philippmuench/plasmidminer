#!/usr/bin/python

try:
import pandas as pd
from sklearn.cross_validation import train_test_split
except ImportError:
    print "This script requires pandasand sklearn to be installed!"

# pd.set_option('display.mpl_style', 'default') # jupyter

# features are located at ../dat/train.features.clear2.csv
# kmers are located at ../dat/train.features.kmer
def creatematrix(features, kmer):
	stat = pd.read_csv(features, sep=",")
	kmer = pd.read_csv(kmer, sep="\t", header = None)

	# concat
	df = pd.concat([stat.reset_index(drop=True), kmer], axis=1)

	# prepare matrix
	id2 = df.id.str.split("-",expand=True) # split the string to get label
	df2 = pd.concat([id2, df], axis=1)
	all = df2.drop(df2.columns[[1,2]], 1)
	all_complete = all[(all.length == 200)] # keep only complete fragments
	all_complete.columns.values[0] = "label"
	# generate test/train set
	y = all_complete['label'].tolist() # extract label
	X = all_complete.drop(all_complete.columns[[0]], 1) # remove label
	return X,y


# generate balanced subsets 
#rand= all_complete.sample(n=3)

# generate cross validation datasets, split x,y arrays into 30percent test data and 70 percent training data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0)


