#! /usr/bin/env python


import argparse
import glob
import pathos.multiprocessing as ptmp
from Bio import SeqIO
from tqdm import tqdm


__author__ = "Dariusz Izak IBB PAS"
__version = "testing"


def find_stop_codons(threshold,
                     records,
                     below_threshold=False):
    below_thr = []
    above_thr = []
    for i in tqdm(records):
        fr_ORFs = [i.seq[x:].translate(table=11) for x in range(3)]
        rv_ORFs = [i.seq.reverse_complement()[x:].translate(table=11) for x in range(3)]
        all_ORFs = fr_ORFs + rv_ORFs
        if all(all_ORFs[x].count("*") > threshold for x in range(6)):
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


def conv_n_filter(files_directory,
                  glob_path="*extendedFrags.fastq",
                  max_stop_codons=3):
    for i in glob.glob("{0}/{1}".format(files_directory, glob_path)):
        print "Processing {}...".format(i)
        records_list = list(SeqIO.parse(i, format="fastq"))
        filtered = find_stop_codons(threshold=max_stop_codons,
                                    records=records_list)
        file_name_fasta = "{0}.fasta".format(".".join(i.split("/")[-1].split(".")[:-1]))
        with open(file_name_fasta, "w") as fout:
            SeqIO.write(filtered,
                        fout,
                        format="fasta")
        print "DONE!"
