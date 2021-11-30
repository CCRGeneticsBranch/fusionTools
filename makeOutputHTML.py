#!/usr/bin/env python
import os,sys,subprocess
import argparse
import json
from classes import *
        
parser = argparse.ArgumentParser(description='Convert fusion output file to html.')
parser.add_argument("--in_file", "-i", metavar="input text file", required=True, help="fusion text file")
parser.add_argument("--out_file", "-o", metavar="output html file", required=True, help="fusion html file")
parser.add_argument("--template_file", "-t", metavar="template html file", required=True, help="template html file")
parser.add_argument("--cyto_file", "-c", metavar="cytoband file", required=True, help="cytoband file")
try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

in_file = args.in_file.strip()
out_file = args.out_file.strip()
temp_file = args.template_file.strip()
cyto_file = args.cyto_file.strip()

makeOutputHTML(in_file, out_file, temp_file, cyto_file)