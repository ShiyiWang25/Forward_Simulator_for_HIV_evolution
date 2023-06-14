# Individual virus-based forward simulator

**This is an individual virus-based forward simulator designed to simulate the acquisition of co-variant mutations in the HIV-1 protease sequence during resistance development.**
---

## Description

This forward simulator models the genetic changes occurring in individual HIV-1 protease genomes during drug-mediated selection over time, encompassing both mutations and the recombination process. It quantifies the fitness of each genome and simulates their evolution under drug pressure.  In addition to capturing Drug-Resistance Mutations (DRMs) that directly contribute to resistance, the simulator incorporates compensatory mutations that mitigate the viral fitness loss caused by DRMs. Throughout the simulation, the effect of each DRM on viral fitness is quantified, and the presence of a compensatory mutation within the same genome earns a bonus point. In this way, this individual virus-based simulator recaptures the selection of co-variant mutations in the HIV protease sequence during drug-mediated selection, allowing for the generation of simulated drug-experienced samples for analysis.

This individual virus-based forward simulator is written in _Python_.

## Features

* The current version of this simulator focuses on the 297-nucleotide HIV-1 _protease_ sequence.

![Image](https://github.com/ShiyiWang25/202306_Simulator/blob/main/Figures/Simulator_WorkFlow.png)
* This simulator simulates the viral evolution by generations (Figure A)
  * In each generation, individual genomes go through 4 steps sequentially: mutation, recombination, fitness calculation, and replication.
  * Mutation and recombination steps are performed using the _ARGS_ and _IT_ modules in the [_mutation-simulator(ver2.0.03)_](https://github.com/mkpython3/Mutation-Simulator), with specific modifications made to fulfill the purposes of this work.
  
* The simulator stops in either of the 2 conditions (Figure B):
  1. Simulated viral rebound: the size of the simulated population exceed a given threshold, e.g., 150,000 genomes
  2. Simulated viral suppression: the size of the simulated population is continuously lower than a given threshold, e.g., 100 genomes, for at least 3 generations.
  
* We designed linkages by generating co-variant mutation pairs composed of synthetic DRM and synthetic compensatory mutations in a 1-to-1 relationship.
 * Each co-variant mutation pair designed has completed resistance and high fitness, which is necessary and sufficient to lead to a simulated viral rebound.

## Installation

Create the environment from the [simulator.yml](https://github.com/ShiyiWang25/202306_Simulator/blob/main/Installation/simulator.yml) file through [Conda](https://anaconda.org/anaconda/conda):
```
conda env create -f simu.yml
```
Activate the environment before running the simulator:
```
conda activate Simu
```

## Dependencies
* Python = 3.6.10
* pandas = 1.0.5
* numpy = 1.19.5
* pysam = 0.16.0.1
* joblib = 1.1.0
* tqdm = 4.62.3

## Methods

This simulator starts with individual HIV-1 genomes stored in a FASTA file. 
In the [**Test Dataset**](https://github.com/ShiyiWang25/202306_Simulator/tree/main/Data), we provided three files needed to start one test simulation:
 1. starting materials for the simulator: [Simu_starting_sequences_2.fa](https://github.com/ShiyiWang25/202306_Simulator/blob/main/Data/materials/Simu_starting_sequences_2.fa)
  * we provided the 83 groups of starting materials for 83 independent simulations. 
  * Each group has 5 _protease_ sequences.
 2. the quantified effects of mutational burden on viral fitness: [kmb_unbiased_0122_4.csv](https://github.com/ShiyiWang25/202306_Simulator/blob/main/Data/materials/kmb_unbiased_0122_4.csv)
 3. synthetic co-variant mutation: [SimuV9_scoring_system_0130](https://github.com/ShiyiWang25/202306_Simulator/blob/main/Data/materials/SimuV9_scoring_system_0130.csv)
  * 24 pairs in total
  * equally assigned to 8 synthetic drug classes (from A to H)

The simulator can be executed with several arguments with a combination of all possible options as bellow:

### Command Line ARGS:

| ARGS | Description |
| --- | --- |
| `-seed_pop` | Starting materials for the simulator stored in a FASTA file |
| `-ref` | Reference sequence stored in a FASTA file |

### Command Line Options:



```
conda activate Simu
python3 ../scripts/Simu_V9_3_hpc_dh.py -run_number 1 -mode cont --start_number 32 -seed_pop ../materials/Simu_starting_sequences.fa -ref ../materials/HXB2_PR.fa --tag 0201 -o ./ -g 800 -R 2.6 -kmb ../scripts/kmb_unbiased_0122_4.csv --sample_time 900 -treatment B --redo_number 20 -settings ./settings.txt -score_info ../materials/SimuV9_scoring_system_0130.csv -rebound_size 200000 --cores 1
```









