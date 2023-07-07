import argparse


def _parse_args(args, **kwargs):

    parser = argparse.ArgumentParser(
                        prog='Simulator',
                        description='An individual virus-based forward simulator for HIV evolution',
                        epilog='')
    
    parser.add_argument('--seed',
                        type=int,
                        help='Give the random seed (optional).')    
    parser.add_argument('--mode',
                        type=str,
                        choices=['init', 'cont'],
                        help='select mode between initiation (init) and continue (cont)')
    parser.add_argument('-i',
                        dest = "input_dir",
                        type=str,
                        help='the directory of Fasta files with starting materials')
    parser.add_argument('--run_number',
                        default = 1,
                        type=int, help='define the number of simulations to run.')
    parser.add_argument('--start_number',
                        default = 0,
                        type=int,
                        help='define the repeat number of the first simulation to run.')
    parser.add_argument('--disc_simu',
                        type=str,
                        help='TXT file for target simulations to run.')   
    parser.add_argument('--ref',
                        type=str,
                        help='import the reference sequence.')
    parser.add_argument('-g',
                        default = 5,
                        type=int,
                        help='define the number of generations.')
    parser.add_argument('--sample_time',
                        default = 5,
                        type=int,
                        help='store viral population information every x generations.')
    parser.add_argument('--score_info',
                        type=str,
                        help='The CSV file with scoring information.')
    parser.add_argument('--treatment',
                        type=str,
                        help='type the treatment name')
    parser.add_argument('--redo_number',
                        default = 1,
                        type=int,
                        help='set the maximum redo times for each treatment.')
    parser.add_argument('--settings',
                        type=str,
                        help='TXT file with the settings information.')
    parser.add_argument('--snv',
                        default = 0.000036,
                        type=float,
                        help='define the mutation rate.')
    parser.add_argument('--rec',
                        default = 0.002,
                        type=float,
                        help='define the recombination rate.')
    parser.add_argument('-R',
                        default = 1.25,
                        type=float,
                        help='define the basic reproductive number.')
    parser.add_argument('--rebound_size',
                        default = 500,
                        type=int,
                        help='Define the VL to define a rebound.')
    parser.add_argument('-o',
                        type=str,
                        help='define the folder to store the output files.')
    parser.add_argument('--tag',
                        type=str,
                        help='define the shared tag for individual genome.')
    parser.add_argument('--cores',
                        default = 1,
                        type=int,
                        help='Core number to use.')
    
    return parser.parse_args(args)
