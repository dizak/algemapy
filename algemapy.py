#! /usr/bin/env python


import jinja2 as jj2
import argparse
import os


def load_template_file(template_file):
    template_Loader = jj2.FileSystemLoader(searchpath="/")
    template_Env = jj2.Environment(loader=template_Loader)
    template = template_Env.get_template(template_file)
    return template


def render_template(template_loaded):
    template_vars = {}
    template_rendered = template_loaded.render(template_vars)
    return template_rendered


def save_template(out_file_name,
                  template_rendered):
    with open(out_file_name, "w") as fout:
        fout.write(template_rendered)


__author__ = "Dariusz Izak IBB PAS"
__version = "testing"


def main():
    pass


if __name__ == '__main__':
    main()
