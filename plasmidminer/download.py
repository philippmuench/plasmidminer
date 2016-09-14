#!/usr/bin/env python
# philipp.muench@helmholtz-hzi.de
# this script downloades plasmid/chromosome sequences from entrez

import time, os, glob, fileinput

try:
    from progress.bar import Bar
except ImportError:
    print "This script requires progess to be installed!"

try:
    from Bio import SeqIO
except ImportError:
    print "This script requires BioPython to be installed!"
try:
    from Bio import Entrez
except ImportError:
    print "This script requires BioPython to be installed!"

def downloadChr(email):
    # imports all complete E.Coli genomes as negative samples
    Entrez.email ="philipp.muench@helmholtz-hzi.de"
    if not os.path.exists('chr'):
        os.makedirs('chr')
    # get the genome ids 
    search_term = '"Escherichia coli"[Organism] AND complete genome[title]'
    print("searching in entrez for " + search_term)
    handle = Entrez.esearch(db='nucleotide', term=search_term, retmax =10)
    genome_ids = Entrez.read(handle)['IdList']
    records = len(genome_ids)
    bar = Bar('download chromosom sequences', max=records)
    for genome_id in genome_ids:
        record = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
        record_gpff = Entrez.efetch(db="protein", id=genome_id, format="gp")
        # write record
        filename_gpff = 'chr/geneBankRecord_{}.gpff'.format(genome_id)
        filename = 'chr/genBankRecord_{}.fasta'.format(genome_id)
        with open(filename, 'w') as f:
            f.write(record.read())
        with open(filename_gpff, 'w') as f:
            f.write(record_gpff.read())
        time.sleep(1) # to make sure not many requests go per second to ncbi
        bar.next()
    bar.finish()
    # create multiple sequence fasta file
    filenames = glob.glob("chr/*.fasta")
    print("creating chromosomes.fasta")
    with open('chromosomes.fasta', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

def downloadPla(email):
    # imports all E.Coli plasmids as positive samples
    # create download folders
    if not os.path.exists('pla'):
        os.makedirs('pla')
    # get the genome ids
    search_term = '"Escherichia coli"[Organism] AND (bacteria[filter] AND plasmid[filter])'
    print("searching in entrez for " + search_term)
    handle = Entrez.esearch(db='nucleotide', term=search_term, retmax =10)
    genome_ids = Entrez.read(handle)['IdList']
    records = len(genome_ids)
    bar = Bar('download plasmid sequences', max=records)
    # iterate over genome ids and download informations
    for genome_id in genome_ids:
        record = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
        record_gpff = Entrez.efetch(db="protein", id=genome_id, format="gp")
        # write gb record
        filename = 'pla/genBankRecord_{}.fasta'.format(genome_id)
        filename_gpff = 'pla/geneBankRecord_{}.gpff'.format(genome_id)
        with open(filename, 'w') as f:
            f.write(record.read())
        with open(filename_gpff, 'w') as f:
            f.write(record_gpff.read())
        time.sleep(1) # to make sure not many requests go per second to ncbi
        bar.next()
    bar.finish()
    # create multiple sequence fasta file
    filenames = glob.glob("pla/*.fasta")
    print("creating plasmid .fasta file")
    with open('plasmids.fasta', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
