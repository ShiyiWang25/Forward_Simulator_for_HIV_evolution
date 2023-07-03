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
  
* The simulator stops in either of the 2 conditions (Figure B):
  1. Simulated viral rebound: the size of the simulated population exceed a given threshold, e.g., 150,000 genomes
  2. Simulated viral suppression: the size of the simulated population is continuously lower than a given threshold, e.g., 100 genomes, for at least 3 generations.
  
* We designed linkages by generating co-variant mutation pairs composed of synthetic DRM and synthetic compensatory mutations in a 1-to-1 relationship.
 * Each co-variant mutation pair designed has completed resistance and high fitness, which is necessary and sufficient to lead to a simulated viral rebound.

## Installation

To install this script with its dependencies, run the following command in your command line:
```
pip install simulator-sw-0701 --extra-index-url=https://test.pypi.org/simple/
```

The latest script in construction can be viewed on [TestPyPI](https://test.pypi.org/project/pysam/): [simulator-sw-0701 0.0.5](https://test.pypi.org/project/simulator-sw-0701/0.0.5/)

## python_requires = >= 3.9
## Dependencies
* numpy
* pandas
* pysam
* joblib
* pyfaidx

## Test
The test folder provides some basic materials to perform a test run using this script. The materials provided are:
1. The reference protease sequences (or the 'wild-type' protease sequence): HXB2_PR.fa
   - all synthetic sequences generated will be aligned to this reference sequence, based on the alignment results their fitness will be determined
2. Synthetic Drug-Resistance and compensatory mutation pairs stored in one CSV file: score_system.csv
   - 8 synthetic drug pressure are provided, named from A to H
   - for each synthetic drug pressure, 3 pairs of synthetic mutation pairs are assigned
3. The values for variables used in this simulator: settings.txt
4. Eleven different viral populations for performing 11 independent simulation runs:
   - each viral population has 30,000 synthetic Drug-Naive protease sequences

Besides the agreements shown in the above example command, the simulator can be executed with several other arguments to provide a more detailed control of the viral evolution. All possible arguments are listed below:

### Command Line ARGS:

| ARGS | Description |
| --- | --- |
| `-seed_pop` | Starting materials for the simulator stored in a FASTA file |
| `-ref` | Reference sequence stored in a FASTA file |
| `-kmb` | Quantified effects of mutational burden on viral fitness in the protease genome |
| `-settings` | Values for all variables used in the simulator stored in a TXT file |
| `-score_info` | Input synthetic co-variant mutation pairs stored in a CSV file|
| `--tag` | Name tag for output files|
| `-o` | Output directory to store the simulation results |

### Command Line Options:

| Options | Description | Default |
| --- | --- | --- |
| `-run_number` | The number of simulations to run in parallel | 1 |
| `-mode` | Request the script to run new simulations from the input starting materials `init`, or continue the previous runs from existing simulation outputs `cont` | `init` |
| `--start_number` | Tell the script to use the materials in the `seed_pop` starting from `start_number` | 1 |
| `--disc_simu` | Instead of starting from a specific material group, define the groups of materials to use in the `seed_pop` | None |
| `-g` | Maxmium generations to simulate | 5 |
| `-R` | Basic reproductive number | 2 |
| `-snv` | Mutation rate | 0.000036 |
| `-rec` | Recombination rate | 0.002 |
| `--sample_time` | Save the intermediate simulated population every `sample_time` generations | 50 |
| `--redo_number` | Allow each simulation to be repeated for `redo_number` times before acquiring a simulated viral rebound | 5 |
| `-rebound_size` | Population size threshold to define a simulated viral rebound | 5000 |
| `-treatment` | Define the synthetic drug class to be used in the simulator | 'A' |
| `--cores` | Number of cores to use | 1 |

### Here are a few example command lines:
1. Run 5 independent simulations using materials #1 to #5 with basic a reproductive ratio of 2.6:
 - Allow the viral populations to evolve for 800 generations
 - Allow 10 attempts for each simulation to reach simulated viral rebound (simulated viral population size > 150000)
 - No sampling during each simulation (`sample_time` > `g`)
 - Use synthetic drug class 'A'
 - Use 5 cores in parallel
```
python3 Simu_V9_3_hpc_dh_2.py -seed_pop ../materials/Simu_starting_sequences.fa -ref ../materials/HXB2_PR.fa -kmb ../materials/kmb_unbiased_0122_4.csv -settings ./materials/settings.txt -score_info ../materials/SimuV9_scoring_system_0130.csv --tag Test  -o ../Outputs -run_number 5 -g 800 -R 2.6 --sample_time 900 -treatment A --redo_number 10 -rebound_size 150000 --cores 5
```
2. Run 2 independent simulations using materials #23 to #44 with a basic reproductive ratio of 2.6:
```
python3 Simu_V9_3_hpc_dh_2.py -seed_pop ../materials/Simu_starting_sequences.fa -ref ../materials/HXB2_PR.fa -kmb ../materials/kmb_unbiased_0122_4.csv -settings ./materials/settings.txt -score_info ../materials/SimuV9_scoring_system_0130.csv --tag Test -o ../Outputs --disc_simu  ../materials/disc_simu.txt --run_number 2 -g 800 -R 2.6 -sample_time 900 -treatment A --redo_number 10  -rebound_size 150000 --cores 2
```
3. Continuous the 2 simulations from the simulated outputs of **Example 2**:
 - switching to treatment B
```
python3 Simu_V9_3_hpc_dh_2.py -seed_pop ../materials/Simu_starting_sequences.fa -ref ../materials/HXB2_PR.fa  -kmb ../materials/kmb_unbiased_0122_4.csv --settings ./materials/settings.txt -score_info ../materials/SimuV9_scoring_system_0130.csv -tag Test -o ../Outputs -mode cont -run_number 2 --disc_simu  ../materials/disc_simu.txt  -g 800 -R 2.6 --sample_time 900 -treatment B --redo_number 10 -rebound_size 150000 --cores 2
```

### How to read the outputs:

The results of the simulations are to be saved in the output folder defined by `-o` in the command line.

In each independent simulation:

  * Each simulation will have its unique subfolder named after the number of the starting materials:
    - The outputs of the simulation starting from material group #11 will be stored in the subfolder Outputs/_Simulation_time_11_.
  
  * The simulator collects the viral population at the simulated viral rebound, the name of which composes of:
    - Tag
    - treatment
    - number of generations took to rebound
    - For one simulated viral rebound at generation 100 under synthetic drug class 'A' with tag 'Test': `Test_simu_A_g100_rebound.fa`
    
  * The size of the simulated population over generations is recorded in _VL_tracking_file.csv_.

  * Other Metadata of this simulation is stored in _Metadata.csv_:
    - Basic R
    - Repeat time
    - Starting viral population
    - Reference 
    - Mutation rate & Recombination rate
    - Rebound size
    - Initial population size
    - Input material group number
    - Settings used
    - Treatment
    - population size over generations
    - Simulation Outcome: 
      - Suppression Notice: In Round 0: Virus is completely suppressed in treatment A
      - Rebound Notice: In treatment A: Virus rebound at generation 100
   
 









