# plasmidminer

detection of plasmid fragments in metagenomic samples

## usuage

download the dataset from ncbi and train/evaluate model

```
python plasmidminer/plasmidminer.py -h
usage: plasmidminer.py [-h] [-c CHUNKSIZE] [--email EMAIL] [--version]

optional arguments:
  -h, --help     show this help message and exit
  -c CHUNKSIZE   Chunk size in nt
  --email EMAIL  email adress needed for ncbi file download
  --version      show program's version number and exit
```

output:

```
ROC AUC: 0.86 (+/- 0.01) [Logistic Regression, saved to logreg.pcl] 
ROC AUC: 0.93 (+/- 0.01) [Random Forest, saved to rf.pcl]
```
plus ROC curves under roc.png



## install
you may want use virtualenv:

`virtualenv env`
`source env/bin/activate`

to install plasmidminer just type: `python setup.py install`

tested with python 2.7.12


`sudo apt-get install libmysqlclient-dev libpq-dev libboost-all-dev`

### train
`python plasmidminer/plasmidminer`

![alt text](index.png "oob")
