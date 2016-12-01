#!/usr/bin/python

# function to prepare the dataset for training process

try:
	import pandas as pd
	from sklearn.cross_validation import train_test_split
    from sklearn.preprocessing import StandardScaler
except ImportError:
    print "This script requires pandasand sklearn to be installed!"

# pd.set_option('display.mpl_style', 'default') # jupyter

def creatematrix(features, kmer):
	stat = pd.read_csv(features, sep=",")
	kmer = pd.read_csv(kmer, sep="\t", header = None) 
	df = pd.concat([stat.reset_index(drop=True), kmer], axis=1) # concat kmerand stat matrix
	id2 = df.id.str.split("-",expand=True) # split the string to get label
	df2 = pd.concat([id2, df], axis=1)
	all = df2.drop(df2.columns[[1,2]], 1)
	all_complete = all[(all.length == 200)] # keep only complete fragments
	all_complete.columns.values[0] = "label"
	return all_complete

print("load data")
dat = creatematrix('dat/train.features.clear2.csv', 'dat/train.features.kmer')

print("split data")
# split to test and training set
y = dat['label'].tolist() # extract label
X = dat.drop(dat.columns[[0]], 1) # remove label

# TODO(pmuench): generate balanced subsets
# get max occurence or labels and use it to draw random subsample from it 
#rand= all_complete.sample(n=3)

# generate cross validation datasets, split x,y arrays into 30percent test data and 70 percent training data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0)

# scale data
sc = StandardScaler()
sc.fit(X_train) # estimate sample mean and standard deviation for each feature dimension from the traing data
X_train_std = sc.transform(X_train) # use same scaling parameters to standardize the test set so that both the values in the trainign and test dataset are comparable to each other
X_test_std = sc.transform(X_test)


print("finished!")
