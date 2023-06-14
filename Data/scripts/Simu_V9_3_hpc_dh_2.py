#!/usr/bin/python

# 2022-05-23, Simulation version 5, from V4_12D3
# Designed for continuable running
# start the run treatment after treatment manually
# must finish all repeats before moving to the next treatment
# A simplified system

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, help='Give the random seed (optional).')
    parser.add_argument('-run_number', type=int, help='define the number of simulations to run.')
    parser.add_argument('--start_number', type=int, help='define the repeat number of the first simulation to run.')
    parser.add_argument('-mode', type=str, help='select mode between initiation (init) and continue (cont)')
    parser.add_argument('--disc_simu', type=str, help='TXT file for target simulations to run.')   
    parser.add_argument('-seed_pop', type=str, help='import genomes in the population in the FASTA file format.')
    parser.add_argument('-ref', type=str, help='import the reference sequence.')
    parser.add_argument('--n', type=int, help='define the number of Ns inserted between repeats in the concatemer.')
    parser.add_argument('--tag', type=str, help='define the shared tag for individual genome.')
    parser.add_argument('-o', type=str, help='define the folder to store the output files.')
    parser.add_argument('-snv', type=float, help='define the mutation rate.')
    parser.add_argument('-rec', type=float, help='define the recombination rate.')
    parser.add_argument('-g', type=float, help='define the number of generations.')
    parser.add_argument('-R', type=float, help='define the basic reproductive number.')
    parser.add_argument('-kmb', type=str, help='import the csv file with the probability for mutational burden.')
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

    if args.n:
        args.n = int(args.n)
    else:
        args.n = 5
        
    if args.snv:
        args.snv = float(args.snv)
    else:
        args.snv = '0.000036' # 10 fold, HIV mutations rete is set to 3.6e-5 per base per generation
        
    if args.rec:
        args.rec = float(args.rec)
    elif args.rec == 0:
        args.rec = float(args.rec)
    else:
        args.rec = '0.002' # 10 fold, HIV average recombination rete is set to 2e-3 per base per generation, but actually not used in the simu
    
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


def main(simulation_time):
    """The main function."""
    # for parallel repeats
    
    # global settings
    global R0
    R0 = args.R
    
    global codon_lib
    codon_lib = codon_lib_generator()
    
    global df_scores
    df_scores = pd.read_csv(args.score_info)
    
    global df_MB
    df_MB = pd.read_csv(args.kmb)
    
    global treatment
    treatment = args.treatment
    
    if args.mode == 'init':
    
        output_folder = args.o + 'Simulation_time_' + str(simulation_time + 1) + '/' # name the result folder for each repeat.
        Path(output_folder).mkdir(parents=True, exist_ok=True) # create the result folder if not existed.

        # create metadata file:
        Metadata_file = output_folder + 'Metadata' + '.csv'
        # write all args to metadata file.
        write_args(Metadata_file, args.R, args.run_number, args.seed_pop, args.ref, args.snv, args.rec, args.rebound_size)
    
        # check initial viral population
        if args.seed_pop:
            #print(f'Seed population found. Grab 5 reads: {simulation_time*5} - {simulation_time*5 + 5}')
            #print(f'Generating starting viral population...')
            #initial_pop_size = args.rebound_size
            initial_pop_size = 300000
            concatemer, line_numbers = starting_pop_sepNs(simulation_time, args.seed_pop, initial_pop_size, args.n)
            note = 'Initial population size: ' + str(initial_pop_size)
            metadata(Metadata_file, note)
            note = 'Reads selected:' + ', '.join([str(i) for i in line_numbers])
            metadata(Metadata_file, note)
            concatemer_output = output_folder + args.tag + '_simu' + '.fa'
            write_file(seq = concatemer, output_file_path = concatemer_output, tag = args.tag)
        else:
            print('Error: missing initial viral population.')
            sys.exit()
            
    elif args.mode == 'cont':
        output_folder = args.o + 'Simulation_time_' + str(simulation_time + 1) + '/' # direct the result folder to continue.
        Metadata_file = output_folder + 'Metadata' + '.csv'

        note = '---*---*---*---*---*---*---'
        metadata(Metadata_file, note)  
        note = 'Continue the simulation on' + time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
        metadata(Metadata_file, note) # take notes for tracking

        concatemer_output = output_folder + args.tag + '_simu' + '.fa'
        concatemer, aside_pop_size = concatemer_sepNs(concatemer_output, args.n)
        initial_pop_size = len(split_concatamer(output_folder + args.tag + '_simu', args.n))
        note = 'Initial population size: ' + str(initial_pop_size)
        metadata(Metadata_file, note)
         
    if args.settings:
        with open(args.settings,'r') as f:
            p, r, c, MB_DRM = [float(i) for i in f.readline().rstrip().rsplit(', ')]
        #print(f'Settings: {p, r, c, MB_DRM}')
        note = 'Settings: p, r, c, MB_DRM'
        metadata(Metadata_file, note)       
        note = 'Settings: ' + str(p) + ', ' +  str(r) + ', ' + str(c) + ', ' + str(MB_DRM)
        metadata(Metadata_file, note)
    else:
        print('Error: missing setting information.')
        sys.exit()
        
    if treatment:
        note = 'Treatment:' + str(args.treatment)
        metadata(Metadata_file, note)  
        note = '---*---*---*---*---*---*---'
        metadata(Metadata_file, note)     
    else:
        print('Error: missing treatment information.')
        sys.exit()
        
    
    # create VL tracking file:
    tracking_file = output_folder + 'VL_tracking_file' + '.csv'
    
    tracking_treatment_generation_VL = {}
    
    # In each repeat:

    switch = False
        
    concatemer_side = concatemer
    redo_count = 0
        
    while switch == False:
        rng = stochastic_function(child_states[simulation_time_list.index(simulation_time) * args.redo_number + redo_count])
        progeny_pool_size_list = []            

        # start viral evolution until reboud or die out
        for generation in range(1, args.g+1, 1):
            # mutation
            # file name e.g., HXB2_1000_ms.fa
            mutation(input_file = concatemer_output, ori_seed = ori_seed, run_number = (args.run_number * args.redo_number), simulation_time = (simulation_time_list.index(simulation_time) * args.redo_number + redo_count), mutation_rate = args.snv)

            # split the mutated concatemer into individual genomes
            input_mutated_concatemer =  concatemer_output.split('.fa')[0] + '_ms' # ./HXB2_simu_ms
            input_mutated_split = input_mutated_concatemer + '_split' # ./HXB2_simu_ms_split
            write_ms_sp_strain_file(seq_list = split_concatamer(input_file = input_mutated_concatemer, n = args.n), output_file = input_mutated_split)

            #recombination
            if args.rec != 0:
                input_ms_sp_file = input_mutated_split + '.fa' # ./HXB2_simu_ms_split.fa
                recombination(input_file = input_ms_sp_file, ori_seed = ori_seed, run_number = (args.run_number * args.redo_number), simulation_time = (simulation_time_list.index(simulation_time) * args.redo_number + redo_count), recombination_rate = args.rec)
                input_ms_sp_ms_it_file = input_mutated_split + '_ms_it.fa'# ./HXB2_simu_ms_split_ms_it.fa
            elif args.rec == 0:
                input_ms_sp_file = input_mutated_split + '.fa' # ./HXB2_simu_ms_split.fa
                input_ms_sp_ms_it_file = input_ms_sp_file

            # fitness calculation and replication
            progeny_list, progeny_pool_size = replication(input_ms_sp_ms_it_file, args.ref, treatment, df_scores, p, r, c, MB_DRM)
            progeny_pool_size_list.append(progeny_pool_size)
            note = str(generation) + ':' + str(progeny_pool_size)
            metadata(Metadata_file, note)

            if progeny_pool_size < 100 and len(progeny_pool_size_list) >= 3 and progeny_pool_size_list[-2] < 100 and progeny_pool_size_list[-3] < 100:
                break
            
            elif progeny_pool_size > args.rebound_size and len(progeny_pool_size_list) > 1 and progeny_pool_size_list[-1] > progeny_pool_size_list[-2] and progeny_pool_size_list[-1]/progeny_pool_size_list[-2] < 1.17:
                break
                    
            elif progeny_pool_size > args.rebound_size and len(progeny_pool_size_list) > 1 and progeny_pool_size_list[-1] > progeny_pool_size_list[-2] and progeny_pool_size_list[-1]/progeny_pool_size_list[-2] >= 1.17:
                #print(f'In treatment {treatment}: Virus rebound at generation {generation}.')
                note = 'In treatment ' + str(treatment) + ': Virus rebound at generation ' + str(generation) + ' with ' + str(progeny_pool_size_list[-1]/progeny_pool_size_list[-2])
                metadata(Metadata_file, note)

                # store the rebound viral population
                concatemer = progeny_concatemer_sepNs(input_file = input_ms_sp_ms_it_file, progeny_list = progeny_list, n = args.n)
                rebound_pop = output_folder + args.tag + '_simu_' + str(treatment) + '_g' + str(generation) + '_rebound.fa' # ./HXB2_simu_g100_rebound.fa
                write_file(seq = concatemer, output_file_path = rebound_pop, tag = args.tag)
                tracking_treatment_generation_VL[treatment] = progeny_pool_size_list

                # generate the progeny population for the next treatment
                concatemer = progeny_concatemer_sepNs(input_file = input_ms_sp_ms_it_file, progeny_list = progeny_list, n = args.n)
                write_file(seq = concatemer, output_file_path = concatemer_output, tag = args.tag)
                switch = True
                break

            else:
                #if generation == 14 and progeny_pool_size < 5000:
                if generation == 14:
                    concatemer, new_progeny_pool_size = progeny_concatemer_sepNs_dh(input_file = input_ms_sp_ms_it_file, progeny_list = progeny_list, n = args.n, progeny_pool_size = progeny_pool_size)
                    write_file(seq = concatemer, output_file_path = concatemer_output, tag = args.tag)  
                    progeny_pool_size_list[-1] = new_progeny_pool_size  
                else:
                    # generate the progeny population
                    # ATTENTION: output file will cover the previous concatemer_output file
                    concatemer = progeny_concatemer_sepNs(input_file = input_ms_sp_ms_it_file, progeny_list = progeny_list, n = args.n)
                    write_file(seq = concatemer, output_file_path = concatemer_output, tag = args.tag)

                # save intermediate timepoint data
                if generation % args.sample_time == 0:
                    put_aside = output_folder + args.tag + '_simu' + str(treatment) + '_g' + str(generation) + '.fa' # ./HXB2_simu_g100.fa
                    write_file(seq = concatemer, output_file_path = put_aside, tag = args.tag)
                else: pass

        if switch == False:
            note = 'In Round ' + str(redo_count) + ': Virus is completely suppressed in treatment ' + str(treatment) + ' with ' + str(progeny_pool_size_list[-1]/progeny_pool_size_list[-2])
            metadata(Metadata_file, note)
            write_file(seq = concatemer_side, output_file_path = concatemer_output, tag = args.tag)
            redo_count += 1
            if redo_count >= args.redo_number:
                note = 'Virus fails to rebound in all attempts under treatment ' + str(treatment)
                metadata(Metadata_file, note)
                return
            else: pass

        else:
            pass
        
    tracking(tracking_file, treatment, initial_pop_size, tracking_treatment_generation_VL)
    return
     
## -------------------------------------------------------------------------
## STEP0 preparation steps
## -------------------------------------------------------------------------

def stochastic_function(seed): # generate RNG
    rng = np.random.default_rng(seed)
    return rng

def concatemer_sepNs(input_file_path, n):
    # concatenate all seq in a multi-FASTA file
    con_seq = ''
    count_seq = 0
    
    with open(input_file_path, 'r') as f:
        for line in f:
            if '>' not in line:
                con_seq += line.rstrip()
                con_seq += 'N'*n
                count_seq += 1
    return (con_seq[:-n]), count_seq

def starting_pop_sepNs(simulation_time, input_file_path, pop_size, n):
    # generate starting viral population
    read_list = []
    line_numbers = list(range(simulation_time * 10 + 1, simulation_time * 10 + 10, 2))
    with open(input_file_path, 'r') as f:
        for i, line in enumerate(f):
            if i in line_numbers:
                read_list.append(line.strip())
            elif i > (simulation_time * 10 + 10):
                break

    con_seq = ''
    repeats = int(pop_size/len(read_list))
    for read in read_list:
        for i in range(repeats):
            con_seq += read
            con_seq += 'N'*n
    return (con_seq[:-n]), line_numbers

def write_file(seq, output_file_path, tag):
    f = open(output_file_path, "w")
    f.write('>'+tag +'\n')
    f.write(seq)
    f.close()
    
def metadata(output_file_path, note):
    f = open(output_file_path, "a+")
    f.write(note + '\n')
    f.close()
    return
    
def write_args(output_file_path, R, run_number, seed_pop, ref, snv, rec, rebound_size):

    note = 'Settings used: '
    metadata(output_file_path, note)
    
    note = 'Basic R : ' + str(R)
    metadata(output_file_path, note)    
    
    note = 'Repeat time: ' + str(run_number)
    metadata(output_file_path, note)
    
    note = 'Starting viral population: ' + str(seed_pop)
    metadata(output_file_path, note)
    
    note = 'Reference: ' + str(ref)
    metadata(output_file_path, note)
    
    note = 'Mutation rate: ' + str(snv) + '; Recombination rate: ' + str(rec)
    metadata(output_file_path, note)
    
    note = 'Rebound size: ' + str(rebound_size)
    metadata(output_file_path, note)
    
    note = '---*---*---*---*---*---*---'
    metadata(output_file_path, note)
    
    return

def codon_lib_generator():
    Base1_list = ['T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G']
    Base2_list = ['T', 'T', 'T', 'T', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'T', 'T', 'T', 'T', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'T', 'T', 'T', 'T', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'T', 'T', 'T', 'T', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G']
    Base3_list = ['T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G']
    AA_list = ['F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', '*', '*', 'C', 'C', '*', 'W', 'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R', 'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R', 'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G']
    codon_lib = {}
    
    for i in range(len(Base1_list)):
        codon = Base1_list[i] + Base2_list[i] + Base3_list[i]
        codon_lib[codon] = AA_list[i]
    return codon_lib

    
## -------------------------------------------------------------------------
## STEP1 mutation
## -------------------------------------------------------------------------

def mutation(input_file, ori_seed, run_number, simulation_time, mutation_rate):
    command = '/home/swang6/Simulation/pipeline_test/scripts/mutation-simulator_V5.py ' + input_file +' -ori_seed ' + str(ori_seed) + ' -run_number ' + str(run_number) + ' -simulation_time ' + str(simulation_time) + ' args -sn ' + str(mutation_rate) + ' > shell.txt'
    os.system(command)
    return

def split_concatamer(input_file, n):
    seq = ''
    input_file_path = input_file + '.fa'
    with open(input_file_path, 'r') as f:
        for line in f:
            if '>' not in line:
                seq += line.rstrip()
    seq_list = seq.split('N'*n)
    return seq_list

def write_ms_sp_strain_file(seq_list, output_file):
    size = len(seq_list)
    output_file_path = output_file + '.fa'
    f = open(output_file_path, "w")
    for i in range(size):
        tagline = '>' + output_file.split('/')[-1] + '_' + str(i) + '\n'
        f.write(tagline)
        f.write(seq_list[i] + '\n')
    f.close()

## -------------------------------------------------------------------------
## STEP2 recombination
## -------------------------------------------------------------------------
def recombination(input_file, ori_seed, run_number, simulation_time, recombination_rate):
    command = '/home/swang6/Simulation/pipeline_test/scripts/mutation-simulator_V5.py ' + input_file +' -ori_seed ' + str(ori_seed) + ' -run_number ' + str(run_number) + ' -simulation_time ' + str(simulation_time) + ' it ' + str(recombination_rate) + ' > shell.txt'
    os.system(command)
    return

## -------------------------------------------------------------------------
## STEP3 fitness calculation and replication
## -------------------------------------------------------------------------

def find_different_nuc(codon_ref, codon_seq):
    
    compare_list = []
    
    for i, j in zip(codon_ref, codon_seq):
        compare_list.append(i==j)
        
    return compare_list

def get_fasta_tag(file):
    with open(file, 'r') as f:
        for line in f:
            if '>' in line:
                seq_tag = line.rstrip().rsplit('>')[1]
                    
    return seq_tag

def translator(nuc_seq, ref_sequence):
    # checked. Correct. 0522.
    import pysam
    from pysam import FastaFile
    
    # get HXB2_PR
    sequences_object = pysam.FastaFile(ref_sequence)
    seq_tag = get_fasta_tag(ref_sequence)
    
    aa_mutation_list = []
    nuc_mutation_list = []
    
    for i in range(0, 297, 3): # PR 297 nucleotides
        codon_ref = ''
        codon_seq = ''
        aa_ref = ''
        aa_seq = ''
        aa_mutation = ''
        nuc_mutation = ''
        
        pos_list = [i+1, i+2, i+3]
        start_position = 0
        nuc_pos = (pos_list[0]-start_position)//3 + 1 # get aa position
        region = 'PR'

        codon_ref = sequences_object.fetch(seq_tag, i, i+3)
        codon_seq = nuc_seq[i:i+3]
        aa_ref = codon_lib[codon_ref]
        aa_seq = codon_lib[codon_seq]

        if aa_seq != aa_ref:
        #if codon_seq != codon_ref:
            aa_mutation = region + '_' + aa_ref + str(nuc_pos) + aa_seq
            
            for i, j in zip(list(codon_ref), find_different_nuc(codon_ref, codon_seq)):
                if j == False:
                    nuc_mutation += str(i)
            for i, j in zip(pos_list, find_different_nuc(codon_ref, codon_seq)):
                if j == False:
                    nuc_mutation += '_'
                    nuc_mutation += str(i)
                    nuc_mutation += '_'
            for i, j in zip(list(codon_seq), find_different_nuc(codon_ref, codon_seq)):
                if j == False:
                    nuc_mutation += str(i)
            
            aa_mutation_list.append(aa_mutation)
            nuc_mutation_list.append(nuc_mutation)
        
    return aa_mutation_list

def get_fasta_seq(input_file):
    with open(input_file, 'r') as f:
        for line in f:
            if '>' not in line:
                seq = line.rstrip()
    return seq

def R_cal(fitness_score, p, MB_pasg):
    
    R = R0/(1 + e**(p-fitness_score)) * MB_pasg
    
    return R

def R_calculation(drug, seq, ref_sequence, df_scores, p, r, c, MB_DRM):
    
    mutations = translator(seq, ref_sequence)
    pDRMs = df_scores['DRM'].tolist() # all possible DRMs
                            
    DRM_CM_lib = {} # store the DRM-CM found in that virus
    neDRM = 0 # number of effective DRMs in the seq
    
    DRMs = list(set(mutations).intersection(set(pDRMs))) # DRMs in target sequence

    for DRM in DRMs:
        pCMs = df_scores['CM'].tolist()[pDRMs.index(DRM)][2:-2].split("', '") # possible CMs
        CMs_number = len(set(mutations).intersection(set(pCMs))) # CMs in target sequence
        if CMs_number > 1:
            DRM_CM_lib[DRM] = 1
        else:
            DRM_CM_lib[DRM] = CMs_number
        if df_scores[drug].tolist()[pDRMs.index(DRM)] == 1:
            neDRM = 1
        else: pass

    if p == 0:
        DR_score = 0
    else:
        DR_score = neDRM * r
        
    if neDRM == 1:
        if len(mutations) <= len(df_MB) - 1:
            MB_pasg = df_MB.iloc[len(mutations)]['K_MB']
        else: 
            MB_pasg = df_MB.iloc[-1]['K_MB']
    else:
        if len(mutations) <= len(df_MB) - 1:
            MB_pasg = df_MB.iloc[len(mutations)]['K_MB_DN']
        else:
            MB_pasg = df_MB.iloc[-1]['K_MB_DN']

    RC_score = sum(list(DRM_CM_lib.values()))*c

    relative_fitness =  DR_score + RC_score - len(DRMs) * MB_DRM
    
    R = R_cal(relative_fitness, p, MB_pasg)

    return R

def replication(input_file_path, ref_sequence, drug, df_scores, p, r, c, MB_DRM):
    
    line_count = 0
    seq = ''
    
    #ReadID_list = [] # store the ReadID for each seq
    
    progeny_number = 0
    progeny_number_list = [] # store the progeny number for each seq
    
    progeny_pool_size = 0 # store the total progeny number in the newly generated population
    
    progeny_list = []
    
    with open(input_file_path, 'r') as f:
        for line in f:
            line_count += 1
            if line_count % 2 == 0:
                seq = line.rstrip()
                #ReadID_list.append(line_count//2)
                
                R = R_calculation(drug, seq, ref_sequence, df_scores, p, r, c, MB_DRM)
                progeny_number = rng.poisson(R, 1)[0]
                progeny_number_list.append(progeny_number)
                progeny_pool_size += progeny_number
                
    for i in range(len(progeny_number_list)):
        progeny_list += [i+1]*progeny_number_list[i]
        
    return progeny_list, progeny_pool_size
    
def random_select(data, pick_number):
    # random pick n number from N number (n<=N)
    # return is a list
    
    import random as rnd
    
    total_number = len(data)
    select_result = []

    for i in sorted(rnd.sample(range(0, total_number), pick_number)):
        select_result.append(data[i])
    
    return select_result

def progeny_concatemer_sepNs(input_file, progeny_list, n):
    # input_file stores the strains after mutation and recombination
    # progeny_list stores the information of strains after replication
    # example: 1,1,2,4 --> progeny from strain 1, 1, 2, 4 in the input_file
    seq_list = []
    with open(input_file, 'r') as f:
        seq_list = []
        for line in f:
            if '>' not in line:
                seq_list.append(line.rstrip())
        progeny_concatemer = ''
        for progeny in progeny_list:
            progeny_concatemer += seq_list[progeny-1]
            progeny_concatemer += 'N'* n
        
        return progeny_concatemer[:-n]
        
def progeny_concatemer_sepNs_dh(input_file, progeny_list, n, progeny_pool_size):
    # input_file stores the strains after mutation and recombination
    # progeny_list stores the information of strains after replication
    # example: 1,1,2,4 --> progeny from strain 1, 1, 2, 4 in the input_file
    seq_list = []
    new_progeny_pool_size = 0
    with open(input_file, 'r') as f:
        seq_list = []
        for line in f:
            if '>' not in line:
                seq_list.append(line.rstrip())
        progeny_concatemer = ''
        for progeny in progeny_list:
            repeat = int(rng.poisson(200000/progeny_pool_size, 1)[0])
            new_progeny_pool_size += repeat
            for r in range(repeat):
                progeny_concatemer += seq_list[progeny-1]
                progeny_concatemer += 'N'* n
        return progeny_concatemer[:-n], new_progeny_pool_size

## -------------------------------------------------------------------------
## Tracking file
## -------------------------------------------------------------------------
    # track the generations and store the population size in each generation
    
def tracking(tracking_file, treatment, initial_pop_size, tracking_treatment_generation_VL):
    
    generation_list = [treatment + str(0)]
    VL_list = [initial_pop_size]
    
    for key, value in tracking_treatment_generation_VL.items():

        for i in range(len(value)):
            generation_list.append(key + str(i+1))
            VL_list.append(value[i])
        
    data = {'Generaion':generation_list, 'Population_size':VL_list}
    df = pd.DataFrame(data)
    
    if Path(tracking_file).is_file():
        df_previous = pd.read_csv(tracking_file)
        df = pd.concat([df_previous, df])
        df.to_csv(tracking_file, index = False)
    else:
        df.to_csv(tracking_file, index = False)
    
    return
    

## -------------------------------------------------------------------------

if __name__ == "__main__":
    try:
        if args.seed:
            ori_seed = int(args.seed)
            #print(f'Seed given: {ori_seed}')
        else:
            ori_seed = int(time.time())
            #print(f'Seed generated: {ori_seed}')
        
        # store the original seed used
        # args.o should end with /
        ori_seed_file = args.o + 'SEED.txt'
        f = open(ori_seed_file, "w")
        f_line = 'Original seed set:' + str(ori_seed)
        f.write(f_line)
        f.close()

        # create the RNG that you want to pass around
        rng = np.random.default_rng(ori_seed)
        # get the SeedSequence of the passed RNG
        ss = rng.bit_generator._seed_seq
        # create n initial independent states
        child_states = ss.spawn(args.run_number * args.redo_number)
        
        if not args.disc_simu:
            simulation_time_list = [i for i in range(args.start_number - 1, args.start_number - 1 + args.run_number)]
        else:
            with open(args.disc_simu,'r') as f:
                simulation_time_list = [int(i)-1 for i in f.readline().rstrip().rsplit(', ')]

        if len(simulation_time_list) != args.run_number:
            print('ERROR: Incorrect simulation number!')
            sys.exit(1)

        Parallel(n_jobs = args.cores)\
        (delayed(main)(simulation_time) for simulation_time in simulation_time_list)
         
    except KeyboardInterrupt:
        print("\n KeyboardInterrupt")
        sys.exit(1)
