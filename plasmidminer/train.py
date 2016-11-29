#!/usr/bin/python

from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt, savetxt, recfromcsv
import re

# remove a column from a structured numpy array based on column name
def remove_field_name(a, name):
    names = list(a.dtype.names)
    if name in names:
        names.remove(name)
    b = a[names]
    return b

d = recfromcsv('/home/pmuench/github.com/philippmuench/plasmidminer/features.small')
d.shape

# iterate over first column and extract the label and id
labels=[]
for x in np.nditer(d["lab"]):
    str=np.array_str(x)
    x = str.split("-")[0]
    labels.append(x)
    
train = remove_field_name(d, "lab")
train
#print d.dtype.fields
