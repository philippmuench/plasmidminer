#!/usr/bin/python

import download, simulate

# download training/learning dataset
# TODO(pmuench): test if files exist
email = "philipp.muench@helmholtz-hzi.de"
download.downloadChr(email)
download.downloadPla(email)

# get read sized chunks
# TODO(pmuench): test if files exist
splitlength = 200 # set length of outpuf fragments
simulate.split(splitlength, 'pla/*.fasta', 'plasmid_200.fasta')
simulate.split(splitlength, 'chr/*.fasta', 'chromosomes_200.fasta')

#TODO(pmuench): get feature matrix