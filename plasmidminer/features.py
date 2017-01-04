#!/usr/bin/python
import sys, re, math

try:
    import CGAT.Experiment as E
    import CGAT.Genomics as Genomics
    import CGAT.IOTools as IOTools
    import CGAT.SequenceProperties as SequenceProperties
    import CGAT.FastaIterator as FastaIterator
except ImportError:
    print "This script requires GCAT to be installed!"

def clearit(inputFileName, outputFileName):
    input = open(inputFileName, "r")
    output = open(outputFileName, "w")
    output.write(input.readline())
    for line in input:
        if not line.lstrip().startswith("#"):
            if not line.lstrip().startswith('"#'):
                    output.write(line)
    input.close()
    output.close()

def main(argv=None):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-w", "--weights-tsv-file", dest="filename_weights",
        type="string",
        help="filename with codon frequencies. Multiple filenames "
        "can be separated by comma.")

    parser.add_option(
        "-s", "--section", dest="sections", type="choice", action="append",
        choices=("length", "sequence", "hid", "na", "aa", "cpg", "dn",
                 "degeneracy", "gaps",
                 "codons", "codon-usage", "codon-translator", "codon-bias"),
        help="which sections to output [%default]")

    parser.add_option(
        "-t", "--sequence-type", dest="seqtype", type="choice",
        choices=("na", "aa"),
        help="type of sequence: na=nucleotides, aa=amino acids [%default].")

    parser.add_option(
        "-e", "--regex-identifier", dest="regex_identifier", type="string",
        help="regular expression to extract identifier from fasta "
        "description line.")

    parser.add_option(
        "--split-fasta-identifier", dest="split_id",
        action="store_true",
        help="split fasta description line (starting >) and use "
        "only text before first space")

    parser.add_option(
        "--add-total", dest="add_total", action="store_true",
        help="add a row with column totals at the end of the table"
        "[%default]")

    parser.set_defaults(
        filename_weights=None,
        pseudocounts=1,
        sections=[],
        regex_identifier="(.+)",
        seqtype="na",
        gap_chars='xXnN',
        split_id=False,
        add_total=False,
    )

    (options, args) = E.Start(parser, argv=argv)
    rx = re.compile(options.regex_identifier)

    reference_codons = []
    if options.filename_weights:
        options.filename_weights = options.filename_weights.split(",")
        for filename in options.filename_weights:
            if filename == "uniform":
                reference_codons.append(Genomics.GetUniformCodonUsage())
            else:
                reference_codons.append(
                    IOTools.ReadMap(IOTools.openFile(filename, "r"),
                                    has_header=True,
                                    map_functions=(str, float)))

        # print codon table differences
        options.stdlog.write(
            "# Difference between supplied codon usage preferences.\n")
        for x in range(0, len(reference_codons)):
            for y in range(0, len(reference_codons)):
                if x == y:
                    continue
                # calculate KL distance
                a = reference_codons[x]
                b = reference_codons[y]
                d = 0
                for codon, p in a.items():
                    if Genomics.IsStopCodon(codon):
                        continue
                    d += b[codon] * math.log(b[codon] / p)
                options.stdlog.write("# tablediff\t%s\t%s\t%f\n" %
                                     (options.filename_weights[x],
                                      options.filename_weights[y],
                                      d))
    iterator = FastaIterator.FastaIterator(options.stdin)

    def getCounter(section):
        if options.seqtype == "na":
            if section == "length":
                s = SequenceProperties.SequencePropertiesLength()
            elif section == "sequence":
                s = SequenceProperties.SequencePropertiesSequence()
            elif section == "hid":
                s = SequenceProperties.SequencePropertiesHid()
            elif section == "na":
                s = SequenceProperties.SequencePropertiesNA()
            elif section == "gaps":
                s = SequenceProperties.SequencePropertiesGaps(
                    options.gap_chars)
            elif section == "cpg":
                s = SequenceProperties.SequencePropertiesCpg()
            elif section == "dn":
                s = SequenceProperties.SequencePropertiesDN()
            # these sections requires sequence length to be a multiple of 3
            elif section == "aa":
                s = SequenceProperties.SequencePropertiesAA()
            elif section == "degeneracy":
                s = SequenceProperties.SequencePropertiesDegeneracy()
            elif section == "codon-bias":
                s = SequenceProperties.SequencePropertiesBias(reference_codons)
            elif section == "codons":
                s = SequenceProperties.SequencePropertiesCodons()
            elif section == "codon-usage":
                s = SequenceProperties.SequencePropertiesCodonUsage()
            elif section == "codon-translator":
                s = SequenceProperties.SequencePropertiesCodonTranslator()
            else:
                raise ValueError("unknown section %s" % section)
        elif options.seqtype == "aa":
            if section == "length":
                s = SequenceProperties.SequencePropertiesLength()
            elif section == "sequence":
                s = SequenceProperties.SequencePropertiesSequence()
            elif section == "hid":
                s = SequenceProperties.SequencePropertiesHid()
            elif section == "aa":
                s = SequenceProperties.SequencePropertiesAminoAcids()
            else:
                raise ValueError("unknown section %s" % section)
        return s

    # setup totals
    totals = {}
    for section in options.sections:
        totals[section] = getCounter(section)
    options.stdout.write("id")
    for section in options.sections:
        options.stdout.write("\t" + "\t".join(totals[section].getHeaders()))

    options.stdout.write("\n")
    options.stdout.flush()
    s = getCounter("hid")
    s.loadSequence("AAAAAAAAA", "na")
    for cur_record in iterator:
        sequence = re.sub(" ", "", cur_record.sequence).upper()
        if len(sequence) == 0:
            raise ValueError("empty sequence %s" % cur_record.title)
        id = rx.search(cur_record.title).groups()[0]
        if options.split_id is True:
            options.stdout.write("%s" % id.split()[0])
        else:
            options.stdout.write("%s" % id)
        options.stdout.flush()
        for section in options.sections:
            s = getCounter(section)
            s.loadSequence(sequence, options.seqtype)
            totals[section].addProperties(s)
            options.stdout.write("\t" + "\t".join(s.getFields()))
        options.stdout.write("\n")
    if options.add_total:
        options.stdout.write("total")
        for section in options.sections:
            options.stdout.write("\t" + "\t".join(totals[section].getFields()))
        options.stdout.write("\n")
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
