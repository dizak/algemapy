#! /usr/bin/env python


import jinja2 as jj2
import argparse
import os
import sys
import glob
from Bio import SeqIO
from tqdm import tqdm


__author__ = "Dariusz Izak IBB PAS"
__version = "testing"


def load_template_file(template_file):
    """
    Load jinja2.environment.Template from HTML_Generator.template_file.
    Search path absolute.

    Parameters
    -------
    template_file: str
        Name of the template file to load.

    Returns
    -------
    jinja2.environment.Template
        Loaded template, ready to be rendered.
    """
    template_Loader = jj2.FileSystemLoader(searchpath="/")
    template_Env = jj2.Environment(loader=template_Loader)
    template = template_Env.get_template(template_file)
    return template


def render_template(template_loaded,
                    files_directory=".",
                    notify_email=None,
                    job_name="algemapy.job",
                    run=None,
                    node_type=None,
                    processors=12,
                    name=None,
                    left=None,
                    right=None,
                    reads=None,
                    ml_software="iqtree-omp"):
    """
    Render previosuly loaded jinja2.environment.Template into str with passed
    vars expanded.

    Parameters
    -------
    template_loaded: loaded jinja2.environment.Template
        jinja2 template, loaded by correctly into jinja2 environment.
    job_name: str, default <algemapy.job>
        Name that identifies your job in slurm. Also used for the project
        identification by algemapy.
    reads: list of two-element tuple/list, default <None>
        Reads file names that will be assembled by flash.

    Returns
    -------
    str
        Ready to be saved with regular file handle.
    """
    template_vars = {"files_directory": files_directory,
                     "notify_email": notify_email,
                     "job_name": job_name,
                     "run": run,
                     "node_type": node_type,
                     "processors": processors,
                     "name": name,
                     "left": left,
                     "right": right,
                     "reads": reads,
                     "ml_software": ml_software}
    template_rendered = template_loaded.render(template_vars)
    return template_rendered


def save_template(out_file_name,
                  template_rendered):
    with open(out_file_name, "w") as fout:
        fout.write(template_rendered)
    """
    Save jinja2 template rendered into str with regular file handle.

    Parameters
    -------
    out_file_name: str
        Output file name.
    template_rendered: str
        jinja2 template rendered into str.
    """


def get_dir_path(file_name=""):
    """
    Find out what is the script system path and return its location. Optionally
    put desired file name at the end of the path. Facilitates access to files
    stored in the same directory as executed script. Requires the executed
    script being added to the system path

    Parameters
    --------
    file_name: str, default <"">
        File name to put at the end of the path. Use empty string if want just
        the directory.

    Returns
    --------
    str
        System path of the executable.

    Examples
    -------
    >>> get_dir_path()
    '/home/user/program/bin/'

    >>> get_dir_path("foo")
    '/home/user/program/bin/foo'
    """
    prog_path = sys.argv[0].replace(sys.argv[0].split("/")[-1],
                                    file_name)
    return prog_path


def left_n_right_generator(files_directory=".",
                           split_sign="_",
                           files_extension="fastq",
                           left_reads_sign="R1",
                           right_reads_sign="R2",
                           return_only=""):
    """
    Align corresponding file names containing the extension of interest in a
    given directory. Extract parts of the file names to use as IDs.

    Parameters
    -------
    files_directory: str, default <.>
        Directory to read the file names from.
    split_sign: str, default <_>
        Characters before this are recognized as IDs.
    files_extension: str, default <fastq>
        Only files with this extension are recognized.
    left_reads_sign: str, default <R1>
        File names with this are recognized as left.
    right_reads_sign: str, default <R2>
        File names with this are recognized as right.
    return_only: str, default <"">
        Return file names recognized as left only if <left>. Return file names
        recognized as right only if <right>. Return IDs only if <name>

    Returns
    -------
    dict of lists of dicts of str
        If containing IDs.
    list of str
        If return_only set to <left>, <right> or <name>.

    Examples
    -------
    >>> left_n_right_generator("/home/user/data/")
    {'left': [{'left_reads': 'F3D1_S189_L001_R1_001.fastq', 'name': 'F3D1'},
              {'left_reads': 'F3D0_S188_L001_R1_001.fastq', 'name': 'F3D0'},
              {'left_reads': 'F3D3_S191_L001_R1_001.fastq', 'name': 'F3D3'},
     'right': [{'name': 'F3D1', 'right_reads': 'F3D1_S189_L001_R2_001.fastq'},
               {'name': 'F3D0', 'right_reads': 'F3D0_S188_L001_R2_001.fastq'},
               {'name': 'F3D3', 'right_reads': 'F3D3_S191_L001_R2_001.fastq'}]}

    >>> left_n_right_generator("/home/user/data/", return_only="left")
    ['F3D1_S189_L001_R1_001.fastq',
     'F3D0_S188_L001_R1_001.fastq',
     'F3D3_S191_L001_R1_001.fastq']
    """
    left_name_reads_list = []
    right_name_reads_list = []
    files_list = os.listdir(files_directory)
    files_list = [i for i in files_list if files_extension == i.split(".")[-1]]
    sample_names_list = [i.split(split_sign)[0] for i in files_list]
    sample_names_list = list(set(sample_names_list))
    for i in sample_names_list:
        for ii in files_list:
            if i == ii.split(split_sign)[0] and left_reads_sign in ii:
                left_name_reads_list.append({"name": i, "left_reads": ii})
            elif i == ii.split(split_sign)[0] and right_reads_sign in ii:
                right_name_reads_list.append({"name": i, "right_reads": ii})
            else:
                pass
    name_reads = {"left": left_name_reads_list,
                  "right": right_name_reads_list}
    if return_only == "left":
        return [i["left_reads"] for i in name_reads["left"]]
    elif return_only == "right":
        return [i["right_reads"] for i in name_reads["right"]]
    elif return_only == "name":
        return [i["name"] for i in name_reads["left"]]
    else:
        return name_reads


def main():
    parser = argparse.ArgumentParser(prog="algemapy",
                                     usage="algemapy.py [OPTION]",
                                     description="ALternativeGEnomicMAppingPYpeline.\
                                                  Facilitates non-16S markers\
                                                  analysis.",
                                     version="testing")
    headnode = parser.add_argument_group("headnode options")
    parser.add_argument(action="store",
                        dest="files_directory",
                        metavar="",
                        help="Input directory path.")
    parser.add_argument("--dry-run",
                        action="store_true",
                        dest="dry_run",
                        default=False,
                        help="Prevents output script execution.")
    parser.add_argument("-n",
                        "--job-name",
                        action="store",
                        dest="job_name",
                        metavar="",
                        default="algemapy.job",
                        help="job name. Used for naming scripts, queued job\
                        and output. Default <algemapy.job>.")
    parser.add_argument("-o",
                        "--output",
                        action="store",
                        dest="output_file_name",
                        metavar="",
                        default="preproc.sh",
                        help="Output file name. Default <preproc.sh>.")
    parser.add_argument("--notify-email",
                        action="store",
                        dest="notify_email",
                        metavar="",
                        default=None,
                        help="Email address you want to notify when job is\
                              done.")
    parser.add_argument("-r",
                        "--run",
                        action="store",
                        dest="run",
                        metavar="",
                        default=None,
                        help="Shell call. Use if you want to run the mothur\
                                script immediately, in current directory.\
                                eg -r sh for regular bash or -r sbatch for\
                                slurm.")
    parser.add_argument("--processors",
                        action="store",
                        dest="processors",
                        metavar="",
                        default=1,
                        help="number of logical processors. Default: <1>")
    parser.add_argument("--ML-software",
                        action="store",
                        dest="ml_software",
                        metavar="",
                        default="iqtree-omp",
                        help="Maximum Likelihood computation software to use.\
                        Use same invocation as when calling the program on its\
                        own. At the moment, only RAxML and iqtree are\
                        recognized and properly set. Anything and everything\
                        can go wrong if using something else. Default\
                        <iqtree-omp>.")
    headnode.add_argument("--node-type",
                          action="store",
                          dest="node_type",
                          metavar="",
                          default=None,
                          help="Node type to use on headnode. N and PHI are\
                          available. Default <None>")
    args = parser.parse_args()

    files_directory_abs = "{0}/".format(os.path.abspath(args.files_directory))
    reads = zip(left_n_right_generator(files_directory_abs,
                                       files_extension="fastq",
                                       return_only="name"),
                left_n_right_generator(files_directory_abs,
                                       return_only="left"),
                left_n_right_generator(files_directory_abs,
                                       return_only="right"))
    if len(reads) == 0:
        print "No fastq files found in {0}. Quitting...".format(files_directory_abs)
        quit()
    if args.run == "sbatch":
        if args.node_type is not None:
            node_type = args.node_type.upper()
        else:
            node_type = "N"
        for name, left, right in reads:
            loaded_templ = load_template_file(get_dir_path("subscript.sh.jj2"))
            rendered_templ = render_template(loaded_templ,
                                             files_directory=files_directory_abs,
                                             job_name=args.job_name,
                                             name=name,
                                             left=left,
                                             right=right,
                                             notify_email=args.notify_email,
                                             node_type=node_type,
                                             processors=args.processors,
                                             ml_software=args.ml_software)
            save_template("{0}.sh".format(name),
                          rendered_templ)
            if args.dry_run is False:
                os.system("sbatch {0}".format("{0}.sh".format(name)))
    elif args.run == "sh":
        if args.node_type is not None:
            node_type = args.node_type.upper()
        else:
            node_type = "N"
        loaded_templ = load_template_file(get_dir_path("sequential.sh.jj2"))
        rendered_templ = render_template(loaded_templ,
                                         files_directory=files_directory_abs,
                                         job_name=args.job_name,
                                         reads=reads,
                                         notify_email=args.notify_email,
                                         node_type=node_type,
                                         processors=args.processors,
                                         ml_software=args.ml_software)
        save_template("{0}.sh".format(args.job_name),
                      rendered_templ)
        if args.dry_run is False:
            os.system("sh {0}".format(args.output_file_name))
    else:
        print "Unknow command for submitting the job. Quitting..."


if __name__ == '__main__':
    main()
