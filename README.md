# Individual virus-based forward simulator

**This is an individual virus-based forward simulator designed to simulate the acquisition of co-variant mutations in the HIV-1 protease sequence during resistance development.**
---

## Description

This forward simulator models the genetic changes occurring in individual HIV-1 protease genomes during drug-mediated selection over time, encompassing both mutations and the recombination process. It quantifies the fitness of each genome and simulates their evolution under drug pressure.  In addition to capturing Drug-Resistance Mutations (DRMs) that directly contribute to resistance, the simulator incorporates compensatory mutations that mitigate the viral fitness loss caused by DRMs. Throughout the simulation, the effect of each DRM on viral fitness is quantified, and the presence of a compensatory mutation within the same genome earns a bonus point. In this way, this individual virus-based simulator recaptures the selection of co-variant mutations in the HIV protease sequence during drug-mediated selection, allowing for the generation of simulated drug-experienced samples for analysis.

This individual virus-based forward simulator is written in _Python_.

## Features

![Image]([Figures/Simulator_WorkFlow.png](https://github.com/ShiyiWang25/202306_Simulator/blob/main/Figures/Simulator_WorkFlow.png))

* The current version of this simulator focuses on the 297-nucleotide HIV-1 _protease_ sequence.
* This simulator simulates the viral evolution by generations
  * In each generation, individual genomes go through 4 steps sequentially: mutation, recombination, fitness calculation, and replication.
  
* The simulator stops in either of the 2 conditions:
  1. Simulated viral rebound: the size of the simulated population exceed a given threshold, e.g., 150,000 genomes
  2. Simulated viral suppression: the size of the simulated population is continuously lower than a given threshold, e.g., 100 genomes, for at least 3 generations.
  
* We designed linkages by generating co-variant mutation pairs composed of synthetic DRM and synthetic compensatory mutations in a 1-to-1 relationship.
 * Each co-variant mutation pair designed has completed resistance and high fitness, which is necessary and sufficient to lead to a simulated viral rebound.




