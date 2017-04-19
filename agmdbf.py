#! /usr/bin/env python


import os
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


def sanitize_ref_clade(line,
                       bad_words=["bacterium",
                                  "miscellaneous",
                                  "sp."]):
    it = iter(line.split(" "))
    genus = None
    species = None
    bad_words = bad_words
    for i in it:
        if i.istitle() is True and any([x.isdigit() for x in i]) is False:
            genus = i
            species = next(it, None)
    if species and genus is not None:
        if any(i in species for i in bad_words) is True:
            return genus
        else:
            return "{} {}\n".format(genus, species)
    elif genus is not None and species is None:
        return genus
    else:
        return "unknown"


def sanitize_ref_tree(input_file_name,
                      output_file_name,
                      file_format="newick"):
    tree = Phylo.read(input_file_name, file_format)
    for i in tree.find_clades():
        i.name = sanitize_ref_clade(i.name)
    Phylo.write(tree, output_file_name, file_format)


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


def gene_seq_dwn(sum_df,
                 output_file_name,
                 rettype="fasta",
                 email="headnode.notify@gmail.com"):
    Entrez.email = email
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
                    file_format,
                    sep="_",
                    increment=True):
    count = 0
    with open(input_file_name) as fin:
        with open(output_file_name, "w") as fout:
            for i in tqdm(SeqIO.parse(fin, file_format)):
                if increment is True:
                    count += 1
                i.id = ""
                desc_split = i.description.split(" ")[1:3]
                desc_split_str = sep.join(desc_split)
                if increment is True:
                    i.description = "{}_{}".format(desc_split_str, count)
                else:
                    i.description = "{}".format(desc_split_str)
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
                        dest="input_file_name",
                        metavar="",
                        help="Input file path.")
    parser.add_argument("-o",
                        "--output",
                        action="store",
                        dest="output_file_name",
                        metavar="",
                        required=True,
                        help="Output file name.")
    parser.add_argument("--download-from-tab-summary",
                        action="store_true",
                        dest="dwn_from_tab_sum",
                        default=False,
                        help="Download genes from nucleotide database by\
                        entries from Gene tabular summary.")
    parser.add_argument("--sanitize-ref-tree",
                        action="store_true",
                        dest="sanitize_ref_tree",
                        default=False,
                        help="Sanitize reference tree by leaving just\
                        taxonomical names.")
    parser.add_argument("--leave-raw",
                        action="store_true",
                        dest="leave_raw",
                        default=False,
                        help="Do not remove raw file downloaded from\
                        nucleotide db.")
    args = parser.parse_args()

    if args.dwn_from_tab_sum is True:
        raw_seqs_file_name = "raw.{}".format(args.output_file_name)
        print "Reformatting Gene database tabular summary..."
        gene_tab_sum = gene_tab_sum_reform(input_file_name=args.input_file_name,
                                           sep="\t",
                                           cols2rename={"genomic_nucleotide_accession.version":
                                                        "ACC_NO",
                                                        "start_position_on_the_genomic_accession":
                                                        "START",
                                                        "end_position_on_the_genomic_accession":
                                                        "END"},
                                           essenatial_cols=["ACC_NO", "START", "END"])
        print "Downloading sequences by ID and coordinates from Gene database tabular summary..."
        gene_seq_dwn(sum_df=gene_tab_sum,
                     output_file_name=raw_seqs_file_name,
                     rettype="fasta")
        print "Removing unwanted part of sequences names..."
        gene_sanit_desc(input_file_name=raw_seqs_file_name,
                        output_file_name=args.output_file_name,
                        file_format="fasta",
                        sep="_",
                        increment=True)
        if args.leave_raw is False:
            os.remove(raw_seqs_file_name)
        print "Done!"
        exit()
    if args.sanitize_ref_tree is True:
        sanitize_ref_tree(input_file_name=args.input_file_name,
                          output_file_name=args.output_file_name,
                          file_format="newick")
        print "Done!"
        exit()


if __name__ == '__main__':
    main()
