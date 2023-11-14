"""
Shiyi Wang
Last Modified: 11-13-2023 3:37 PM
"""

import argparse
import os
import numpy as np
import pandas as pd
import pysam
import math
from math import e
from pysam import FastaFile
from pathlib import Path
from joblib import Parallel, delayed
import time
import sys
from pyfaidx import Fasta, Faidx



EOL = "\n"              # end of line terminator


def simu_treated(args):   

    # determine simulation run numbers: e.g., from Simulation 1 to Simulation 10

    simulation_time_list = []

    # protect against negative start_number and run_number
    if (args.disc_simu is None
        and args.start_number >= 1
        and args.run_number > 0):

        for i in range(args.start_number - 1, args.start_number - 1 + args.run_number):
            simulation_time_list.append(i)

    elif os.path.exists(args.disc_simu):

        with open(args.disc_simu,'r') as f:

            for tline in f:
                tline = tline.strip()

                simulation_time_list.append(int(tline) - 1)

    else:
        raise ValueError("Invalid disc_simu, start_number or run_number")

    Path(args.o).mkdir(parents=True, exist_ok=True)                                      
                                        
    # generate a random independent seed for each simulation run
    child_states = random_seed_generator(args.o,
                                         args.run_number *  args.redo_number,
                                         seed=args.seed)

    # store all arguments
    # save space for child_states and simulation_time_list
    simu_var = Variables(args.seed, args.mode, args.input_dir, args.run_number, args.start_number,
                         args.disc_simu, args.ref, args.g, args.sample_time,
                         args.score_info, args.treatment, args.redo_number, args.settings,
                         args.snv, args.rec, args.R, args.rebound_size, args.o, args.tag,
                         args.cores, child_states, simulation_time_list)

    Parallel(n_jobs = simu_var.cores)(
    delayed(simu_var.job)(simulation_time) for simulation_time in simulation_time_list)
    

def concatemer_sepNs(input_file_path):
    # concatenate all seq in a multi-FASTA file
    con_seq = ''
    count_seq = 0
    with open(input_file_path, 'r') as f:
        for line in f:
            if '>' not in line:
                con_seq += line.rstrip()
                con_seq += 'N' * 5
                count_seq += 1
    return (con_seq[:-5]), count_seq

    
def codon_lib_generator():
    Base1_list = ['T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T',
                  'T', 'T', 'T', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
                  'C', 'C', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
                  'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G',
                  'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G']
    Base2_list = ['T', 'T', 'T', 'T', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A',
                  'G', 'G', 'G', 'G', 'T', 'T', 'T', 'T', 'C', 'C', 'C', 'C',
                  'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'T', 'T', 'T', 'T',
                  'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G',
                  'T', 'T', 'T', 'T', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A',
                  'G', 'G', 'G', 'G']
    Base3_list = ['T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G',
                  'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G',
                  'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G',
                  'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G',
                  'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T', 'C', 'A', 'G',
                  'T', 'C', 'A', 'G']
    AA_list = ['F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', '*', '*', 'C',
               'C', '*', 'W', 'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H',
               'Q', 'Q', 'R', 'R', 'R', 'R', 'I', 'I', 'I', 'M', 'T', 'T', 'T',
               'T', 'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R', 'V', 'V', 'V', 'V',
               'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G']
    codon_lib = {}
    
    for i in range(len(Base1_list)):
        codon = Base1_list[i] + Base2_list[i] + Base3_list[i]
        codon_lib[codon] = AA_list[i]
    return codon_lib
    

class Variables:
        
    def __init__(self, seed, mode, input_dir, run_number, start_number,
                 disc_simu, ref, g, sample_time, score_info,
                 treatment, redo_number, settings, snv, rec, R,
                 rebound_size, o, tag, cores, child_states, simulation_time_list):
        self.seed = seed
        self.mode = mode
        self.input_dir = input_dir
        self.run_number = run_number
        self.start_number = start_number
        self.disc_simu = disc_simu
        self.ref = ref
        self.g = g
        self.sample_time = sample_time
        self.score_info = score_info
        self.treatment = treatment
        self.redo_number = redo_number
        self.settings = settings
        self.snv = snv
        self.rec = rec
        self.R = R
        self.rebound_size = rebound_size
        self.o = o
        self.tag = tag
        self.cores = cores
        self.child_states = child_states
        self.simulation_time_list = simulation_time_list

    def job(self, simulation_time):   
        
        # create the result folder for this simulation run if not existed
        output_folder = os.path.join(self.o, f"Simulation_time_{simulation_time + 1}")
        Path(output_folder).mkdir(parents=True, exist_ok=True)
        
        # name the file with starting materials
        sequences_file = os.path.join(output_folder, f"{self.tag}_simu.fa")
        
        # create a Metadata file
        Metadata_file = os.path.join(output_folder, 'Metadata.txt')
        
        # create a viral loading tracking file for this simulation run
        tracking_file = os.path.join(output_folder, 'VL_tracking_file.csv')
        tracking_treatment_generation_VL = {}
        
        # generate the starting viral population
        # store it in sequences_file
        # take notes of metadata
        initial_pop_size, concatemer, p, r, c, MB_DRM = self.preparation(simulation_time,
                                                                         output_folder,
                                                                         sequences_file,
                                                                         Metadata_file)
        
        # store the starting materials to "concatemer_side"
        concatemer_side = concatemer
        
        # perform the simulation for redo_number times unless reaching rebound
        # define the pointer to decide whether to continuously repeat the simulation run
        switch = False
        
        # count the number of repeats 
        redo_count = 0
        
        while not switch:
            switch, progeny_pool_size_list = Variables.each_repeat(self, Metadata_file,
                                                                   simulation_time, output_folder,
                                                                   sequences_file, p, r, c, MB_DRM,
                                                                   redo_count, switch)
            redo_count += 1
            
            if not switch: 
                # repeat
                # rewrite starting file using the concatemer_side
                note = (f"In Round {redo_count} : Virus is completely suppressed"
                        f" in treatment {self.treatment} with"
                        f" {progeny_pool_size_list[-1]/progeny_pool_size_list[-2]}")

                metadata(Metadata_file, note)

                write_file_split(concatemer = concatemer_side,
                                 output_file_path = sequences_file, tag = self.tag)

                if redo_count >= self.redo_number:
                    note = f"Virus fails to rebound in all attempts under treatment {self.treatment}"
                    metadata(Metadata_file, note)
                    return
            
        tracking_treatment_generation_VL[self.treatment] = progeny_pool_size_list  
        tracking(tracking_file, self.treatment, initial_pop_size, tracking_treatment_generation_VL)
    
    def preparation(self, simulation_time, output_folder, sequences_file, Metadata_file):

        if self.mode == 'init':
            # write all args to the metadata file.
            self.write_args(Metadata_file)

            # create an initial viral population
            if self.input_dir:
                initial_pop_size = 30000
                start_materials_fas = write_starting_pop(simulation_time, 
                                                         get_starting_materials(self.input_dir),
                                                         initial_pop_size,
                                                         sequences_file,
                                                         self.tag)
                concatemer, initial_pop_size = concatemer_sepNs(sequences_file)

                note = EOL.join([f"Initial population size: {initial_pop_size}",
                                 f"Starting materials collected from file: {start_materials_fas}"])
                metadata(Metadata_file, note)
                
            else:
                print('Error: missing initial viral population.')
                sys.exit()

        elif self.mode == 'cont':

            note = EOL.join(['---*---*---*---*---*---*---',
                             f'Continue the simulation on {time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())}'])
            metadata(Metadata_file, note)

            concatemer, initial_pop_size = concatemer_sepNs(sequences_file)

        if self.settings:
            with open(self.settings,'r') as f:
                p, r, c, MB_DRM = [float(i) for i in f.readline().rstrip().rsplit(', ')]
  
            note = EOL.join(['Settings: p, r, c, MB_DRM',
                             f"Settings: {p}, {r}, {c}, {MB_DRM}"])
            metadata(Metadata_file, note)
        else:
            print('Error: missing setting information.')
            sys.exit()

        if self.treatment:
            note = EOL.join([
            f"Treatment: {self.treatment}",
            f"---*---*---*---*---*---*---{EOL}"
            ])
            metadata(Metadata_file, note)       
        else:
            print('Error: missing treatment information.')
            sys.exit()

        return initial_pop_size, concatemer, p, r, c, MB_DRM

    def each_repeat(self, Metadata_file, simulation_time, output_folder,
                    sequences_file, p, r, c, MB_DRM, redo_count, switch):

        # generate a random seed
        switch = False
        rng = stochastic_function(self.child_states 
                                  [self.simulation_time_list.index(simulation_time) * 
                                  self.redo_number + redo_count])
        
        progeny_pool_size_list = [] 

        # simulate in g generations until a viral rebound or suppression is achieved
        for generation in range(1, self.g + 1):
            progeny_list, progeny_pool_size, recombined_mutated_sequences_file = (
                each_generation(sequences_file, self.snv, self.rec, self.score_info, 
                                self.ref, self.treatment, p, r, c, MB_DRM, self.R, rng))
            progeny_pool_size_list.append(progeny_pool_size)
            metadata(Metadata_file, 
            		 f"{generation}: {progeny_pool_size}")

            # check
            if (len(progeny_pool_size_list) >= 3 and progeny_pool_size < 100 and 
                progeny_pool_size_list[-2] < 100 and progeny_pool_size_list[-3] < 100):
                break

            elif (len(progeny_pool_size_list) >= 3 and 
                  progeny_pool_size > self.rebound_size and 
                  progeny_pool_size_list[-1] > progeny_pool_size_list[-2] and 
                  progeny_pool_size_list[-2] > progeny_pool_size_list[-3]):
                metadata(Metadata_file, 
                		(f"In treatment {self.treatment}:"
                         f" Virus rebound at generation {generation}"
                         f" with {progeny_pool_size_list[-1]/progeny_pool_size_list[-2]}"))

                # save the rebound viral population
                # e.g., ./HXB2_simu_g100_rebound.fa
                rebound_pop_file = os.path.join(output_folder, 
                                                f"{self.tag}_simu_{self.treatment}_"
                                                f"g{generation}_rebound.fa")
                write_pop_file(input_file = recombined_mutated_sequences_file, 
                               output_file = rebound_pop_file, 
                               progeny_list = progeny_list, 
                               tag = self.tag)
                switch = True            
                break

            else:
                # set a pseudo drug-holiday at generation 14
                # amplify the population size to 30,000
                if generation == 14:
                    new_progeny_pool_size = write_dh_pop_file(
                        input_file = recombined_mutated_sequences_file, 
                        output_file = sequences_file, 
                        progeny_list = progeny_list, 
                        progeny_pool_size = progeny_pool_size, 
                        tag = self.tag, 
                        rng = rng)
                    progeny_pool_size_list[-1] = new_progeny_pool_size  
                else:
                    # generate the progeny population
                    # ATTENTION: output file will cover the previous concatemer_output file
                    write_pop_file(input_file = recombined_mutated_sequences_file, 
                                   output_file = sequences_file, 
                                   progeny_list = progeny_list, 
                                   tag = self.tag)

                # save intermediate timepoint data
                # e.g., ./HXB2_simu_g100.fa
                if generation % self.sample_time == 0:
                    put_aside = os.path.join(output_folder, 
                                             f"{self.tag}_simu_{self.treatment}_"
                                             f"g{generation}.fa")
                    write_pop_file(input_file = recombined_mutated_sequences_file,
                                   output_file = put_aside, 
                                   progeny_list = progeny_list,
                                   tag = self.tag)

        return switch, progeny_pool_size_list

    
    def write_args(self, output_file_path):
    
        note = EOL.join(['Settings used: ',
                         f"Basic R : {self.R}",
                         f"Repeat time: {self.run_number}",
                         f"Starting materials stored in: {self.input_dir}",
                         f"Reference: {self.ref}"
                         f"Mutation rate: {self.snv}; Recombination rate: {self.rec}",
                         f"Rebound size: {self.rebound_size}",
                         f"---*---*---*---*---*---*---"])
    
        metadata(output_file_path, note)

            
def each_generation(sequences_file, snv, rec, score_info, ref, treatment, p, r, c, MB_DRM, R, rng):
    
    # mutate
    # ./HXB2_simu_ms
    mutated_sequences_file = f"{sequences_file.split('.fa')[0]}_ms.fa"
    mutator(input_file = sequences_file,
            output_file = mutated_sequences_file,
            mutation_rate = float(snv),
            rng = rng)

    # recombine
    # ./HXB2_simu_ms_split_ms_it.fa
    if rec != 0: # perform recombine
        recombined_mutated_sequences_file = (f"{sequences_file.split('.fa')[0]}_"
                                             f"ms_it.fa")
        recombination(input_file = mutated_sequences_file, 
                      output_file = recombined_mutated_sequences_file,
                      recombination_rate = float(rec), 
                      rng = rng)
    elif rec == 0: # skep the recombination process
        recombined_mutated_sequences_file = mutated_sequences_file
        
    # fitness calculation and viral replication
    df_scores = pd.read_csv(score_info)
    progeny_list, progeny_pool_size = replication(recombined_mutated_sequences_file, 
                                                  ref, treatment,
                                                  df_scores, 
                                                  p, r, c, MB_DRM, 
                                                  R, rng)
    
    return progeny_list, progeny_pool_size, recombined_mutated_sequences_file

    
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

def get_starting_materials(starting_material_directory):
    
    fasta_file_list = []
    for root, dirs, files in os.walk(starting_material_directory):
        for file in files:
            if file.endswith(".fa"):
                fasta_file_list.append(os.path.join(root, file))
                
    return fasta_file_list
    
def get_fasta_seq(input_file):
    with open(input_file, 'r') as f:
        for line in f:
            if '>' not in line:
                seq = line.rstrip()
    return seq

def metadata(output_file_path, note):
    with open(output_file_path, "a+") as f:
        f.write(note + EOL)

def mutator(input_file, output_file, mutation_rate, rng):
    seq_lib = read_fasta_file(input_file)
    seq_names = list(seq_lib.keys())

    with open(output_file, "w") as f:
        for seq_name in seq_names:
            seq = seq_lib[seq_name][:].seq
            seq_mut = mutate(seq, mutation_rate, rng)
            f.write(EOL.join([
                f">{seq_name}",
                f"{seq_mut}{EOL}"
            ]))

def mutate(seq, mutation_rate, rng):
    
    # get positions to mutate
    mean = len(seq) * mutation_rate
    mut_count = rng.poisson(mean, 1)[0]
    positions = sorted(rng.choice(np.arange(0, len(seq)), mut_count, replace=False))
    
    if len(positions) == 0:
        return seq
    else:
        mutant_seq = mutant_sequence(seq, positions, rng)
        return mutant_seq
      
def mutant_sequence(seq, positions, rng):
    mut_seq = seq
    for pos in positions:
        nucleotide_list = ['A', 'T', 'G', 'C']
        nucleotide_list.remove(seq[pos])
        alt = rng.choice(nucleotide_list, 1)[0]
        mut_seq = mut_seq[:pos] + alt + mut_seq[pos+1:]
    return mut_seq
    
def random_seed_generator(output_dir, total_number, seed=None):

    seed = np.random.SeedSequence(entropy=seed)

    # store the seed used in each simulation run
    with open(os.path.join(output_dir, "SEED.txt"), "w") as f:
        f.write(f"Seed set to {seed.entropy}{EOL}")

    return seed.spawn(total_number)

    
def read_fasta_file(input_file):
    seq_lib = Fasta(input_file)
    return seq_lib
    
def recombination(input_file, output_file, recombination_rate, rng):
    seq_lib = read_fasta_file(input_file)
    seq_names = list(seq_lib.keys())
    
    # pick 20% of the sequences to recombine
    seq_to_recombine = list(rng.choice(seq_names, int(len(seq_names)/6), replace=False))
    seq_names_norecombine = sorted(list(set(seq_names) - set(seq_to_recombine)))
    
    recombine_pairs_lib = {}
    # pair
    for i in range(int(len(seq_to_recombine)/2)):
        paired = rng.choice(seq_to_recombine, 2, replace=False)
        recombine_pairs_lib[paired[0]] = paired[1]
        seq_to_recombine.remove(paired[0])
        seq_to_recombine.remove(paired[1])
    
    
    with open(output_file, "w") as f:
    
        for seq_name in seq_names_norecombine:
            seq = seq_lib[seq_name][:].seq
            f.write(EOL.join([
                f">{seq_name}",
                f"{seq}{EOL}"
            ]))
            
        for key, value in recombine_pairs_lib.items():
            seq1 = seq_lib[key][:].seq
            seq2 = seq_lib[value][:].seq
            seq1_rec, seq2_rec = recombine(seq1, seq2, recombination_rate, rng)
            f.write(EOL.join([
                f">{key}",
                f"{seq1_rec}",
                f">{value}",
                f"{seq2_rec}{EOL}",
            ]))

def recombine(seq1, seq2, recombination_rate, rng):
    mean = (len(seq1) -1)* recombination_rate
    rec_count = rng.poisson(mean, 1)[0]
    positions = sorted(rng.choice(np.arange(1, len(seq1)), rec_count, replace=False))
    
    if len(positions) == 0:
        return seq1, seq2
    else:
        seq1_rec, seq2_rec = switch(seq1, seq2, positions)
        return seq1_rec, seq2_rec
        
def R_cal(R0, x, p):
    
    if p == 0:
        R = R0 * 2/(1+e**(-x))
    else:
        R = R0 * 2/(1+e**(p-x))
        
    return R

def R_calculation(drug, seq, ref_sequence, df_scores, p, r, c, MB_DRM, R0):
    
    mutations = translator(seq, ref_sequence)
    pDRMs = df_scores['DRM'].tolist() # all possible DRMs
                            
    DRM_CM_lib = {} # store the DRM-CM found in that virus
    neDRM = 0 # number of effective DRMs in the seq
    
    DRMs = list(set(mutations).intersection(set(pDRMs))) # DRMs in target sequence

    for DRM in DRMs:
        pCMs = df_scores['CM'].tolist()[pDRMs.index(DRM)][2:-2].split("', '") # possible CMs
        CMs_number = len(set(mutations).intersection(set(pCMs))) # CMs in target sequence
        DRM_CM_lib[DRM] = CMs_number
        if df_scores[drug].tolist()[pDRMs.index(DRM)] == 1:
            neDRM = 1
        else: pass

    x = 0
    if p == 0:
        x = 0 - MB_DRM * len(DRMs) + sum(list(DRM_CM_lib.values()))*c - MB_DRM/10 * len(mutations)
    else:
        x = neDRM * r - MB_DRM * len(DRMs) + sum(list(DRM_CM_lib.values()))*c - MB_DRM/10 * len(mutations)

    R = R_cal(R0, x, p)

    return R

def replication(input_file_path, ref_sequence, drug, df_scores, p, r, c, MB_DRM, R0, rng):
    
    line_count = 0
    seq = ''
    
    #ReadID_list = [] # store the ReadID for each seq
    
    progeny_number = 0
    progeny_number_list = [] # store the progeny number for each seq
    
    progeny_pool_size = 0 # store total progeny number in the newly generated population
    
    progeny_list = []
    
    with open(input_file_path, 'r') as f:
        for line in f:
            line_count += 1
            if line_count % 2 == 0:
                seq = line.rstrip()
                #ReadID_list.append(line_count//2)
                
                R = R_calculation(drug, seq, ref_sequence, df_scores, p, r, c, MB_DRM, R0)
                progeny_number = rng.poisson(R, 1)[0]
                progeny_number_list.append(progeny_number)
                progeny_pool_size += progeny_number
                
    for i in range(len(progeny_number_list)):
        progeny_list += [i+1]*progeny_number_list[i]
        
    return progeny_list, progeny_pool_size
        

def stochastic_function(seed): # generate RNG
    rng = np.random.default_rng(seed)
    return rng
    
def switch(seq1, seq2, positions):
    for pos in positions:
        seq1_list = [seq1[:pos], seq1[pos:]]
        seq2_list = [seq2[:pos], seq2[pos:]]
        seq1 = seq1_list[0] + seq2_list[1]
        seq2 = seq2_list[0] + seq1_list[1]
    return seq1, seq2

def translator(nuc_seq, ref_sequence):
    # checked. Correct. 0522.
    
    codon_lib = codon_lib_generator()
    
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
    
def write_dh_pop_file(input_file, output_file, progeny_list, progeny_pool_size, tag, rng):
    seq_list = []
    with open(input_file, 'r') as f:
        seq_list = []
        for line in f:
            if '>' not in line:
                seq_list.append(line.rstrip())
                
    read_count = 0
    with open(output_file, "w") as fo:
        for progeny in progeny_list:
            repeat = int(rng.poisson(30000/progeny_pool_size, 1)[0])
            for r in range(repeat):
                fo.write('>'+ tag + '_' + str(read_count + 1) + EOL)
                fo.write(seq_list[progeny-1] + EOL)
                read_count += 1
    return read_count
    
def write_pop_file(input_file, output_file, progeny_list, tag):
    seq_list = []
    with open(input_file, 'r') as f:
        seq_list = []
        for line in f:
            if '>' not in line:
                seq_list.append(line.rstrip())
                
    read_count = 0
    with open(output_file, "w") as fo:
        for progeny in progeny_list:
            fo.write('>'+ tag + '_' + str(read_count + 1) + EOL)
            fo.write(seq_list[progeny-1] + EOL)
            read_count += 1
    return
        
    
def write_file_split(concatemer, output_file_path, tag):
    seq_list = concatemer.split('N' * 5)

    with open(output_file_path, "w") as f:
        for i in range(len(seq_list)):
            f.write('>'+ tag + '_' + str(i + 1) + EOL)
            f.write(seq_list[i] + EOL)

def write_starting_pop(simulation_time, start_materials, pop_size, output_file_path, tag):
    # generate starting viral population
    start_materials_fas = start_materials[simulation_time]
    
    read_list = []
    with open(start_materials_fas, 'r') as f:
        for line in f:
            if '>' not in line:
                read_list.append(line.strip())

    read_count = 0
    repeats = int(pop_size/len(read_list))
    with open(output_file_path, "w") as fo:
        for read in read_list:
            for i in range(repeats):
                fo.write(EOL.join([
                    f">{tag}_{read_count + 1}",
                    f"{read}{EOL}"
                ]))
                read_count += 1
    return start_materials_fas

