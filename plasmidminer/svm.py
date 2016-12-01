import numpy as np
from sklearn.cross_validation import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression

pipe_lr = Pipeline([('scl', StandardScaler()), ('svm', SVC(kernel='rbf', random_state=0))])


pipe_lr.fit(X_train, y_train)

print('Test Accuracy: %.3f' % pipe_lr.score(X_test, y_test))


kfold = StratifiedKFold(y=y_train, n_folds=2, random_state=1)
scores = []


for k, (train, test) in enumerate(kfold):
            pipe_lr.fit(X_train[train], y_train[train])
                score = pipe_lr.score(X_train[test], y_train[test])
                    scores.append(scroe)
                        #print('Fold: %s, Class dist.: %s, Acc: %.3f' % (k+1, np.bincount(y_train[train]), score))
