#!/usr/bin/python
# resample fragments from raw input data
# philipp.muench@helmholtz-hzi.de

import os
import getdataset
import argparse
import shutil
import errno


def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:  # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--data_folder', action='store', dest='data_folder', help='Input data folder which should be processed without tailing backslash', default='dat')
	parser.add_argument('--output', action='store', dest='data', help='output data folder without tailing backslash', default='dat_resampled')
	parser.add_argument('--save', action='store', dest='save', help='Save dataset as msg pack object', default='dataset_resampled.msg')
	parser.add_argument('-c', '--chunksize', action='store', dest='chunksize', help='Chunk size in nt', default=200)
	parser.add_argument('-s', '--readsim', dest='readsim', action='store_true', help='use read simluation based on wgsim instead of sliding window')
	parser.add_argument('-N', '--simnum', action='store', dest='simnum', help='number of reads to simulate with wgsim', default=1000)
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	args = parser.parse_args()

	# copy raw data to new dir
	#if not os.path.exists(str(args.data)):
	#	os.makedirs(str(args.data))
	copyanything(str(args.data_folder), str(args.data))

	# split data into small read sized subset using sliding window approach
	getdataset.getchunks(args)

	# remove sequences that are too short
	getdataset.sequence_cleaner(str(args.data) + '/chromosome_chunks.fasta.corrected.fasta', args)
	getdataset.sequence_cleaner(str(args.data) + '/plasmid_chunks.fasta.corrected.fasta', args)

	# merge genomes to one file
	getdataset.createmerged(args)

	# export features such as GC content
	getdataset.getstatfeatures(args)

	# compress fasta file
	getdataset.compress(args)

	# extract kmer content
	getdataset.extractkmers(args)

	# save feature data as msgpack object
	getdataset.savemsg(args)
