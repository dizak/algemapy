#! /usr/bin/env python


import sys
import warnings
import re
import argparse
import glob
import pathos.multiprocessing as ptmp
from Bio import SeqIO
from Bio import Phylo as ph
from tqdm import tqdm


__author__ = "Dariusz Izak IBB PAS"
__version = "testing"


def sanitize_names(input_file_name,
                   output_file_name,
                   leading_char="@",
                   unwanted_char=":",
                   wanted_char="_"):
    """
    Remove unwanted characters from read names. Reads are recognized as
    lines starting with <>> character. Uses simplests possible built-in methods.

    Parameters
    -------
    input_file_name: str
        Path to input file.
    output_file_name: str
        Path to output file.
    leading_char: str or int
        Character to identify lines of interest. Using empty str causes\
        removal of unwanted_char everywhere.
    unwanted_char: str
        Charater to remove.
    wanted_char: str
        Charater to replace unwanted_char with.

    Returns
    -------
    list of str
        Sanitized input file content as list of lines.

    Examples
    -------
    >>> sanitized = sanitize_names("/path/to/your/file",
                                   unwanted_char=":",
                                   wanted_char="_")
    >>> sanitized[0]
    '>M00967_43_000000000-A3JHG_1_1101_10551_7682 1_N_0_188\n'
    """
    corrected_file = []
    with open(input_file_name) as fin:
        for i in fin.readlines():
            if i.startswith(leading_char):
                i = i.replace(unwanted_char, wanted_char)
            else:
                pass
            corrected_file.append(i)
    with open(output_file_name, "w") as fout:
        fout.writelines(corrected_file)


def dots4names(input_file_name,
               output_file_name,
               file_format="newick",
               wanted_char=".",
               rec_lim=10000):
    """
    Replace read names with dot or other desired character. Uses Bio.Phylo
    module. Write results to file.

    Parameters
    -------
    input_file_name: str
        Path to input file.
    output_file_name: str
        Path to output file.
    file_format: str
        Input and output file format. Default <newick>.
    wanted_char: str
        Charater to replace read names with.
    """
    tree = ph.read(input_file_name, file_format)
    warnings.warn("Recursion limit set to {}".format(rec_lim),
                  UserWarning)
    sys.setrecursionlimit(rec_lim)
    for i in tqdm(tree.find_clades()):
        i.name = wanted_char
    ph.write(tree, output_file_name, file_format)


def find_stop_codons(records,
                     threshold,
                     stop_sign="*",
                     below_threshold=False):
    """
    Filter out (high- or lowpass) Bio.SeqRecord.SeqRecord depending on
    specified number of stop codons found in the translation of all potential\
    ORFs.

    Parameters
    -------
    records: list of Bio.SeqRecord.SeqRecord
        Records to processes.
    threshold: into
        Maximum number of stop codons allowed.
    stop_sign: str
        Representation of stop codon. Default: <*>.
    below_threshold: bool
        Save records below or above threshold. Default: <False>.

    Returns
    -------
    below_thr: list of Bio.SeqRecord.SeqRecord
        Records with smaller number of stop codons than threshold.
    above_thr: list of Bio.SeqRecord.SeqRecord
        Records with greater number of stop codons than threshold.
    """
    below_thr = []
    above_thr = []
    for i in tqdm(records):
        fr_ORFs = [i.seq[x:].translate(table=11) for x in range(3)]
        rv_ORFs = [i.seq.reverse_complement()[x:].translate(table=11) for x in range(3)]
        all_ORFs = fr_ORFs + rv_ORFs
        if all(all_ORFs[x].count(stop_sign) > threshold for x in range(6)):
            below_thr.append(i)
        else:
            above_thr.append(i)
    if below_threshold is True:
        print "{0} reads left".format(len(below_thr))
        print "{0} reads removed ({1}%)".format(len(above_thr),
                                                (len(above_thr) * 100 /
                                                 len(below_thr)))
        return below_thr
    else:
        print "{0} reads left".format(len(above_thr))
        print "{0} reads removed ({1}%)".format(len(below_thr),
                                                (len(below_thr) * 100 /
                                                 len(above_thr)))
        return above_thr


def conv_n_filter(input_file_name,
                  output_file_name,
                  input_format="fastq",
                  output_format="fasta",
                  max_stop_codons=3,
                  multiprocessing=False):
    """
    Convert sequences from file to another format, replace <:> from read names\
    with <_> and filter out those with number of stop codons above maximum\
    value. Write results to file.

    Parameters
    -------
    input_file_name: str
        Path to input file.
    output_file_name: str
        Path to output file.
    input_format: str
        Input file format. Default <fastq>.
    output_format: str
        Output file format. Default <fasta>.
    max_stop_codons: int
        Maximum value of stop codons allowed. Default: <3>.
    multiprocessing: bool
        Experimental. Use all the cores on a given machine. Default <False>.
    """
    def f(i):
        print "Processing {}...".format(i)
        records_list = list(SeqIO.parse(i, format=input_format))
        filtered = find_stop_codons(threshold=max_stop_codons,
                                    records=records_list)
        with open(output_file_name, "w") as fout:
            SeqIO.write(filtered, fout, format=output_format)
        print "DONE!"
    if multiprocessing is True:
        ptmp.ProcessPool().map(sanitize_names, input_file_name)
        ptmp.ProcessingPool().map(f, input_file_name)
    else:
        sanitize_names(input_file_name, input_file_name)
        f(input_file_name)


def main():
    parser = argparse.ArgumentParser(prog="agmfi",
                                     usage="agmfi.py [FILE] [OPTION]",
                                     description="Part of \
                                     ALternativeGEnomicMAppingPYpeline.\
                                     Holds algemapy built-in filtering.",
                                     version="testing")
    parser.add_argument(action="store",
                        dest="input_file",
                        metavar="",
                        help="Input file.")
    parser.add_argument("-o",
                        "--output",
                        action="store",
                        dest="output_file_name",
                        metavar="",
                        required=True,
                        help="Output file name.")
    parser.add_argument("-s",
                        "--sanitize-only",
                        action="store",
                        dest="sanitize_only",
                        metavar="",
                        default=None,
                        help="Use if you just want to remove specified\
                        character from the file.")
    parser.add_argument("-d",
                        "--dotize-only",
                        action="store_true",
                        dest="dotize_only",
                        default=None,
                        help="Use if you just want to replace read names\
                        with dots in the tree.")
    args = parser.parse_args()

    if args.sanitize_only is True:
        sanitize_names(input_file_name=args.input_file,
                       output_file_name=args.output_file_name,
                       leading_char="",
                       unwanted_char=args.sanitize_only,
                       wanted_char="")
        exit()
    else:
        pass
    if args.dotize_only is not None:
        dots4names(input_file_name=args.input_file,
                   output_file_name=args.output_file_name,
                   file_format="newick",
                   wanted_char=".",
                   rec_lim=10000)
        exit()
    else:
        pass
    conv_n_filter(input_file_name=args.input_file,
                  output_file_name=args.output_file_name,
                  input_format="fastq",
                  output_format="fasta",
                  max_stop_codons=3,
                  multiprocessing=False)


if __name__ == '__main__':
    main()
