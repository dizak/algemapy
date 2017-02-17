#! /usr/bin/env python


import argparse
from Bio import SeqIO


__author__ = "Dariusz Izak IBB PAS"
__version = "testing"


def id_reform(input_file_name,
              output_file_name):
    output_ids = []
    output_recs = []
    fasta = list(SeqIO.parse(input_file_name, "fasta"))
    for i in fasta:
        num = i.id.split(".")[0]
        tax = i.description.split(";")[-1]
        i.id = "{0}.{1}".format(tax, num)
    with open(output_file_name, "w") as fout:
        SeqIO.write(fasta, fout, "fasta")


def main():
    parser = argparse.ArgumentParser(prog="agmdbf",
                                     usage="agmdbf.py [FILE] [OPTION]",
                                     description="Part of \
                                     ALternativeGEnomicMAppingPYpeline.\
                                     Holds algemapy built-in database\
                                     formatting.",
                                     version="testing")
    parser.add_argument(action="store",
                        dest="files_directory",
                        metavar="",
                        help="Input file path.")
    parser.add_argument("-o",
                        "--output",
                        action="store",
                        dest="output_path",
                        metavar="",
                        default=".",
                        help="Output file path. Default: working\
                        directory")
    args = parser.parse_args()
    id_reform(args.files_directory, args.output_path)


if __name__ == '__main__':
    main()
