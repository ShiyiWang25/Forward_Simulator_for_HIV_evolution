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
pip install simulator-sw-0701==0.0.8 --extra-index-url=https://test.pypi.org/simple/
```

The latest script in construction can be viewed on [TestPyPI](https://test.pypi.org/project/pysam/) using the following link: [simulator-sw-0701 0.0.5](https://test.pypi.org/project/simulator-sw-0701/0.0.5/)

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
  
In the Test folder, run the following command line to initiate 3 independent simulation runs.
This test will stop in minutes as only 5 generations are allowed in each simulation run. The simulation outputs will be saved in `output/Simulation_time_1`, `output/Simulation_time_2`, and `output/Simulation_time_3`.
```
python3 -m  simulator_sw_0701 -i start_materials/ -run_number 3 -ref HXB2_PR.fa -score_info score_system.csv -treatment A -settings settings.txt -o ./output -tag test -cores 3 -mode init -g 5 -rebound_size 100000
```

Besides the agreements shown in the above example command, the simulator can be executed with several other arguments to provide a more detailed control of the viral evolution. All possible arguments are listed below:

### Command Line ARGS:

| Options | Description | Default |
| --- | --- | --- |
| `--seed` | Provide a random seed. If not provided, the script will generate one automatically | None |
| `-m`, `-mode` | Tells the script to run new simulations from the input starting materials by `init`, or continue the previous runs from existing simulation outputs by `cont` | `init` |
| `-i` | Provides the directory where Fasta files with starting materials are stored` | None |
| `--run_number` | Pefines the number of simulations to run | 1 |
| `--start_number` | Provide the index of the first simulation to run | 1 |
| `--disc_simu` | Provides the index of seperate simulation to run | None |
| `-r`, `--ref` | Provides the reference sequence | None |
| `-g` | Maximum number of generations in each simulation | 5 |
| `--sample_time` | Saves the intermediate simulated population every `sample_time` generations | 5 |
| `--score_info` | Provides synthetic mutation pairs act against each synthetic drug pressure | None |
| `--treatment` | Select one synthetic treatment from A to H | None |
| `--redo_number` | Allow `redo_number` attempts for each simulation before acquiring a simulated viral rebound | 1 |
| `--settings` | Assigns values to variables used | None |
| `--snv` | HIV-1 Mutation rate | 0.000036 |
| `--rec` | HIV-1 Recombination rate | 0.002 |
| `-R` | Basic reproductive number | 1.25 |
| `--rebound_size` | Defines the population size threshold to define a simulated viral rebound | 500 |
| `-o` | Define the folder to store the output files | None |
| `--tag` | The tag to name output files | None |
| `--cores` | Number of cores to use | 1 |


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
   
 









