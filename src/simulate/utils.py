import argparse

def _parse_args(args, **kwargs):
    # argparse.ArgumentParser arguments
    # prog: use this if the name is different then that used in running
    #    application.
    # description: text above packages command line argument list
    # epilog: text below command line argument list
    # formatter_class: maintain text layout in string literals.

    parser = argparse.ArgumentParser(**kwargs,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

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
                        
    return parser.parse_args(args)