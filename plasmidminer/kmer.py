#!/usr/bin/python

import sys, re, itertools

try:
    from progress.bar import Bar
except ImportError:
    print "This script requires progess to be installed!"

try:
    import CGAT.FastaIterator as FastaIterator
    import CGAT.Experiment as E
except ImportError:
    print "This script requires GCAT to be installed!"

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """
    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-k", "--kmer-size", dest="kmer", type="int",
                      help="supply kmer length")

    parser.add_option(
        "-p", "--output-proportion", dest="proportion", action="store_true",
        help="output proportions - overides the default output")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # do not allow greater than octonucleotide
    assert options.kmer <= 8, "cannot handle kmer of length %i" % options.kmer

    # how we deal with the nucleotides depends on the kmer length
    nucleotides = []
    for nucleotide in ["A", "C", "T", "G"]:
        nucleotides = nucleotides + \
            [x for x in itertools.repeat(nucleotide, options.kmer)]

    E.info("retrieving %imer sequences" % options.kmer)
    # get all kmer sequences to query
    kmers = set()
    for kmer in itertools.permutations(nucleotides, options.kmer):
        kmers.add(kmer)

    E.info("matching %imers in file" % options.kmer)
    # count the number of kmers in each sequence

    result = {}

    # NB assume that non fasta files are caught by FastaIterator
    total_entries = 0
   # bar = Bar('extract kmer information', max=len(FastaIterator.iterate(options.stdin)))
    for fasta in FastaIterator.iterate(options.stdin):
        total_entries += 1
        result[fasta.title] = {}
        for kmer in kmers:
            counts = [m.start()
                      for m in re.finditer("".join(kmer), fasta.sequence)]
            result[fasta.title][kmer] = len(counts)
    #    bar.next()
    #bar.finish()

    E.info("writing results")
    # write out the results
    headers = result.keys()
    rows = set()
    for kmer_counts in result.values():
        for kmer, count in kmer_counts.iteritems():
            rows.add("".join(kmer))
    # write header row
    options.stdout.write("kmer\t" + "\t".join(headers) + "\n")

    # output proportions if required - normalises by
    # sequence length
    E.info("computing total counts")
    totals = {}
    for header in headers:
        totals[header] = sum([result[header][tuple(row)] for row in rows])

    for row in rows:
        if options.proportion:
            options.stdout.write("\t".join(
                [row] + [str(float(result[header][tuple(row)]) / totals[header]) for header in headers]) + "\n")
        else:
            options.stdout.write(
                "\t".join([row] + [str(result[header][tuple(row)]) for header in headers]) + "\n")

    E.info("written kmer counts for %i contigs" % total_entries)
    # write footer and output benchmark information.
    E.Stop()
if __name__ == "__main__":
  sys.exit(main(sys.argv))
