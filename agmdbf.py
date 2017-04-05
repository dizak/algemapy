#! /usr/bin/env python


import argparse
from Bio import Entrez, SeqIO, SeqRecord, Phylo
import pandas as pd
from tqdm import tqdm


__author__ = "Dariusz Izak IBB PAS"
__version = "testing"


def db_dwn(db,
           term,
           output_file_name,
           email="dariusz.izak@ibb.waw.pl"):
    Entrez.email = email
    db_handle = Entrez.esearch(db=db, term=term)
    db_record = Entrez.read(db_handle)
    db_full_size = db_record["Count"]
    db_handle = Entrez.esearch(db=db, term=term, retmax=db_full_size)
    db_record = Entrez.read(db_handle)
    with open(output_file_name, "w") as fout:
        for i in tqdm(db_record["IdList"]):
            hand = Entrez.efetch(db=db, id=i, rettype="fasta", retmode="text")
            rec = SeqIO.read(hand, "fasta")
            SeqIO.write(rec, fout, "fasta")


def id_reform(input_file_name,
              output_file_name):
    output_ids = []
    output_recs = []
    with open(output_file_name, "w") as fout:
        for i in tqdm(SeqIO.parse(input_file_name, "fasta")):
            record = SeqRecord(id="{0}.{1}".format(i.description.split(";")[-1],
                                                   i.id.split(".")[0]),
                               seq=i.seq,
                               description="")
            SeqIO.write(record, fout, "fasta")


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
    parser.add_argument("--download",
                        action="store_true",
                        dest="download",
                        default=False,
                        help="Use if you want to download database.")
    args = parser.parse_args()
    if args.download is True:
        db_dwn(db="nucleotide",
               term="pheS[Gene] AND Lactobacillaceae [Orgn]",
               output_file_name=args.output_path)
    id_reform(args.files_directory, args.output_path)


if __name__ == '__main__':
    main()
