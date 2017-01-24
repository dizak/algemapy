#! /usr/bin/env python


import jinja2 as jj2
import argparse
import os
import sys


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
                    job_name="algemapy.job",
                    partition="long",
                    nodes=1,
                    ntasks_per_node=6,
                    mem_per_cpu=24,
                    node_list=None):
    """
    Render previosuly loaded jinja2.environment.Template into str with passed
    vars expanded.

    Parameters
    -------
    template_loaded: loaded jinja2.environment.Template
        jinja2 template, loaded by correctly into jinja2 environment.
    job_name: str
        Name that identifies your job in slurm. Also used for the project
        identification by algemapy.
    partition: str
        Headnode's partition. Values: test, short, big, long, accel. Accel
        necessary for phi/gpu nodes Default <long>.
    nodes: int
        Number of nodes. Default: <1>.
    ntasks_per_node: int
        Number of tasks to invoke on each node. Default <6>.
    mem_per_cpu: int:
        Maximum amount of real memory per node in gigabytes. Default <24>.
    node_list: str
        Request a specific list of nodes.

    Returns
    -------
    str
        Ready to be saved with regular file handle.
    """
    template_vars = {"job_name": job_name,
                     "partition": "long",
                     "nodes": 1,
                     "ntasks_per_node": 6,
                     "mem_per_cpu": 24,
                     "node_list": None}
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
    file_name: str
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


def main():
    parser = argparse.ArgumentParser(prog="algemapy",
                                     usage="algemapy.py [OPTION]",
                                     description="ALternativeGEnomicMAppingPYpeline.\
                                                  Facilitates non-16S markers\
                                                  analysis.",
                                     version="testing")
    headnode = parser.add_argument_group("headnode options")
    parser.add_argument("-o",
                        "--output",
                        action="store",
                        dest="output_file_name",
                        metavar="",
                        default="preproc.sh",
                        help="output file name. Default <preproc.sh>.")
    headnode.add_argument("--partition",
                          action="store",
                          dest="partition",
                          metavar="",
                          default="long",
                          help="headnode's partition. Values: test, short, big,\
                                  long, accel. Accel necessary for phi/gpu nodes\
                                  Default <long>.")
    headnode.add_argument("--nodes",
                          action="store",
                          dest="nodes",
                          metavar="",
                          default=1,
                          help="number of nodes. Default: <1>.")
    headnode.add_argument("--ntasks-per-node",
                          action="store",
                          dest="ntasks_per_node",
                          metavar="",
                          default=6,
                          help="number of tasks to invoke on each node.\
                                Default <6>")
    headnode.add_argument("--mem-per-cpu",
                          action="store",
                          dest="mem_per_cpu",
                          metavar="",
                          default=24,
                          help="maximum amount of real memory per node in\
                                gigabytes. Default <24>.")
    headnode.add_argument("--node-list",
                          action="store",
                          dest="node_list",
                          metavar="",
                          default=None,
                          help="request a specific list of nodes")
    args = parser.parse_args()

    loaded_templ = load_template_file("/home/darek/Pulpit/algemapy/preproc_template.sh.jj2")
    rendered_templ = render_template(loaded_templ,
                                     partition=args.partition,
                                     nodes=args.nodes,
                                     ntasks_per_node=args.ntasks_per_node,
                                     mem_per_cpu=args.mem_per_cpu,
                                     node_list=args.node_list)
    save_template(args.output_file_name,
                  rendered_templ)


if __name__ == '__main__':
    main()
