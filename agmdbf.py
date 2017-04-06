#! /usr/bin/env python


import re
import argparse
from Bio import Entrez, SeqIO, SeqRecord, Phylo
import pandas as pd
from tqdm import tqdm


__author__ = "Dariusz Izak IBB PAS"
__version = "testing"


def db_dwn(output_file_name,
           db,
           term,
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


def gene_tab_sum_reform(input_file_name,
                        sep="\t",
                        cols2rename={"genomic_nucleotide_accession.version":
                                     "ACC_NO",
                                     "start_position_on_the_genomic_accession":
                                     "START",
                                     "end_position_on_the_genomic_accession":
                                     "END"},
                        essenatial_cols=["ACC_NO", "START", "END"]):
    tab_sum_df = pd.read_csv(input_file_name, sep=sep)
    tab_sum_ren = tab_sum_df.rename(columns=cols2rename)
    tab_sum_nona = tab_sum_ren[essenatial_cols].dropna()
    tab_sum_nona.ACC_NO = tab_sum_nona.ACC_NO.map(lambda x: re.sub("\.\d+",
                                                                   "",
                                                                   x))
    tab_sum_nona.START = tab_sum_nona.START.map(lambda x: re.sub("\.\d+",
                                                                 "",
                                                                 str(x)))
    tab_sum_nona.END = tab_sum_nona.END.map(lambda x: re.sub("\.\d+",
                                                             "",
                                                             str(x)))
    return tab_sum_nona


def gene_seq_dwn(output_file_name,
                 sum_df,
                 rettype="fasta"):
    with open(output_file_name, "w") as fout:
        for i in tqdm(sum_df.itertuples()):
            id = getattr(i, "ACC_NO")
            start = getattr(i, "START")
            end = getattr(i, "END")
            hand = Entrez.efetch(db="nucleotide",
                                 id=id,
                                 seq_start=start,
                                 seq_stop=end,
                                 rettype=rettype,
                                 retmode="text")
            rec = SeqIO.read(hand, rettype)
            SeqIO.write(rec, fout, rettype)


def gene_sanit_desc(input_file_name,
                    output_file_name,
                    file_format):
    def f(x):
        x_split = x.split(" ")[1:3]
        x_split_str = " ".join(x_split)
        return x_split_str
    with open(input_file_name) as fin:
        with open(output_file_name, "w") as fout:
            for i in tqdm(SeqIO.parse(fin, file_format)):
                i.id = ""
                i.description = f(i.description)
                SeqIO.write(i, fout, file_format)


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
