"""Individual virus-based forward simulator

By: Shiyi Wang
"""

_epilog="""

python -m gdml --help

NOTES ON __main__.py FILE

The __main__.py file is special and allows one to run the 
files contents from the command line if the program is 
correctly installed.  When running the program remember to
use the '-m' flag:

python -m PROGRAM_NAME [options] arguments

The benefit is that the code is installed on the users system
and managed by pip.  Consequently, the user doesn't need to 
keep track of script files or software versions using a
version control system.

"""

import os
import argparse
import numpy as np
import pandas as pd
import pysam
import math
from pathlib import Path
import time
from joblib import Parallel, delayed
import multiprocessing
import sys
from math import e
from pyfaidx import Fasta, Faidx

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, help='Give the random seed (optional).')
    parser.add_argument('-run_number', type=int, help='define the number of simulations to run.')
    parser.add_argument('--start_number', type=int, help='define the repeat number of the first simulation to run.')
    parser.add_argument('-mode', type=str, help='select mode between initiation (init) and continue (cont)')
    parser.add_argument('--disc_simu', type=str, help='TXT file for target simulations to run.')   
    parser.add_argument('-i', type=str, help='the directory of Fasta files with starting materials')
    parser.add_argument('-ref', type=str, help='import the reference sequence.')
    parser.add_argument('--tag', type=str, help='define the shared tag for individual genome.')
    parser.add_argument('-o', type=str, help='define the folder to store the output files.')
    parser.add_argument('-snv', type=float, help='define the mutation rate.')
    parser.add_argument('-rec', type=float, help='define the recombination rate.')
    parser.add_argument('-g', type=float, help='define the number of generations.')
    parser.add_argument('-R', type=float, help='define the basic reproductive number.')
    parser.add_argument('--sample_time', type=float, help='store viral population information every x generations.')
    parser.add_argument('-treatment', type=str, help='type the treatment name')
    parser.add_argument('--redo_number', type=int, help='set the maximum redo times for each treatment.')
    parser.add_argument('-settings', type=str, help='TXT file with the settings information.')
    parser.add_argument('-score_info', type=str, help='The CSV file with scoring information.')
    parser.add_argument('-rebound_size', type=int, help='Define the VL to define a rebound.')
    parser.add_argument('--cores', type=int, help='Core number to use.')
    args = parser.parse_args()
    
    if args.run_number:
        args.run_number = int(args.run_number)
    else:
        args.run_number = 1
 
    if args.start_number:
        args.start_number = int(args.start_number)
    else:
        args.start_number = 1

    if args.mode:
        args.mode = str(args.mode)
    else:
        args.mode = 'init' # default mode is mode initiation 
        
    if args.snv:
        args.snv = float(args.snv)
    else:
        args.snv = '0.000036' # HIV mutations rete is set to 3.6e-5 per base per generation
        
    if args.rec:
        args.rec = float(args.rec)
    elif args.rec == 0:
        args.rec = float(args.rec)
    else:
        args.rec = '0.002' # HIV average recombination rete is set to 2e-3 per base per generation, but actually not used in the simu
    
    if args.g:
        args.g = int(args.g)
    else:
        args.g = 5 # by default run the program for 5 generation
    
    if args.R:
        args.R = float(args.R)
    else:
        args.R = 2  # by default set the reproductive number of the best strain to 2
                    # also the mean for poisson distribution
    
    if args.redo_number:
        args.redo_number = int(args.redo_number)
    else:
        args.redo_number = 5

    if args.sample_time:
        args.sample_time = int(args.sample_time)
    else:
        args.sample_time = 50 # by default store the viral sequence information every 50 generations
    
    if args.score_info:
        args.score_info = args.score_info
    else:
        pass
    
    if args.treatment:
        args.treatment = str(args.treatment)
    else:
        args.treatment = None
        
    if args.rebound_size:
        args.rebound_size = int(args.rebound_size)
    else:
        args.rebound_size = 500 # by default defining the rebound size to be 5000