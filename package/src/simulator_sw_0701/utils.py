import argparse
def _parse_args(args, **kwargs):

    parser = argparse.ArgumentParser(
                        prog='Simulator',
                        description='An individual virus-based forward simulator for HIV evolution',
                        epilog='')
    
    parser.add_argument('--seed', type=int, help='provide a random seed.')    
    parser.add_argument('-m', '--mode', type=str, choices=['init', 'cont'], help='select mode between initiation (init) and continue (cont)')
    parser.add_argument('-i', type=str, help='provides the directory where Fasta files with starting materials are stored')
    parser.add_argument('--run_number', default = 1, type=int, help='define the number of simulations to run')
    parser.add_argument('--start_number', default = 1, type=int, help='provide the index of the first simulation to run')
    parser.add_argument('--disc_simu', type=str, help='provide the index of seperate simulation to run')   
    parser.add_argument('-r', '--ref', type=str, help='provide the reference sequence')
    parser.add_argument('-g', default = 5, type=int, help='the maximum number of generations in each simulation')
    parser.add_argument('--sample_time', default = 5, type=int, help='save the intermediate simulated population every x generations')
    parser.add_argument('--score_info', type=str, help='provide synthetic mutation pairs act against each synthetic drug pressure ')
    parser.add_argument('--treatment', type=str, help='select one synthetic treatment from A to H')
    parser.add_argument('--redo_number', default = 1, type=int, help='set the maximum redo times for each simulation.')
    parser.add_argument('--settings', type=str, help='assign values to variables used')
    parser.add_argument('--snv', default = 0.000036, type=float, help='define the mutation rate')
    parser.add_argument('--rec', default = 0.002, type=float, help='define the recombination rate')
    parser.add_argument('-R', default = 1.25, type=float, help='define the basic reproductive number')
    parser.add_argument('--rebound_size', default = 500, type=int, help='defines the population size threshold to define a simulated viral rebound')
    parser.add_argument('-o', type=str, help='define the folder to store the output files')
    parser.add_argument('--tag', type=str, help='define the tag to name output files')
    parser.add_argument('--cores', default = 1, type=int, help='numbers of cores to use.')
    
    return parser.parse_args(args)
