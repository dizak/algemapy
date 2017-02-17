#! /usr/bin/env python


import argparse
from Bio import SeqIO


__author__ = "Dariusz Izak IBB PAS"
__version = "testing"


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


if __name__ == '__main__':
    main()
