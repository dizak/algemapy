#! /usr/bin/env python


import jinja2 as jj2
import argparse
import os


__author__ = "Dariusz Izak IBB PAS"
__version = "testing"


def load_template_file(template_file):
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


def main():
    parser = argparse.ArgumentParser(prog="algemapy",
                                     usage="algemapy.py [OPTION]",
                                     description="ALternativeGEnomicMAppingPYpeline.\
                                                  Facilitates non-16S markers\
                                                  analysis.",
                                     version="testing")
    headnode = parser.add_argument_group("headnode options")
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
    headnode.add_argument("--processors",
                          action="store",
                          dest="processors",
                          metavar="",
                          default=24,
                          help="number of logical processors. Default: <24>")
    args = parser.parse_args()


if __name__ == '__main__':
    main()
