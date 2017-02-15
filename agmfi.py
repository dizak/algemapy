#! /usr/bin/env python


import argparse
import glob
import pathos.multiprocessing as ptmp
from Bio import SeqIO
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
                  output_path=".",
                  max_stop_codons=3,
                  multiprocessing=False):
    def f(i):
        print "Processing {}...".format(i)
        records_list = list(SeqIO.parse(i, format="fastq"))
        filtered = find_stop_codons(threshold=max_stop_codons,
                                    records=records_list)
        file_name_fasta = "{0}.fasta".format(".".join(i.split("/")[-1].split(".")[:-1]))
        with open("{0}/{1}".format(output_path, file_name_fasta), "w") as fout:
            SeqIO.write(filtered, fout, format="fasta")
        print "DONE!"
    input_files = glob.glob("{0}/{1}".format(files_directory, glob_path))
    if multiprocessing is True:
        ptmp.ProcessPool().map(sanitize_names, input_files)
        ptmp.ProcessingPool().map(f, input_files)
    for i in input_files:
        sanitize_names(i, i)
        f(i)


def main():
    parser = argparse.ArgumentParser(prog="agmfi",
                                     usage="agmfi.py [OPTION]",
                                     description="Part of \
                                     ALternativeGEnomicMAppingPYpeline.\
                                     Holds algemapy built-in filtering.",
                                     version="testing")
    parser.add_argument(action="store",
                        dest="files_directory",
                        metavar="",
                        help="Input directory path.")
    parser.add_argument("-o",
                        "--output",
                        action="store",
                        dest="output_path",
                        metavar="",
                        default=".",
                        help="Output directory path. Default: working\
                        directory")
    args = parser.parse_args()

    conv_n_filter(files_directory=args.files_directory,
                  output_path=args.output_path)


if __name__ == '__main__':
    main()
