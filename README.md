# Local Search Solver

A solver based on local search.

## Description

The goal of this repository is to provide a simple framework to quickly implement heuristic algorithms based on local search.

For complex and evolving problems, implementing and extending local search algorithms based on complex neighborhoods quickly becomes cumbersome and time consuming. This often makes them unsuitable for practical use.
The algorithms of this repository are designed to get the best performances out of the simplest neighborhoods. Thus making them easier to implement and extend.

## Examples

### Packing

[Knapsack Problem with Conflicts](examples/knapsackwithconflicts.hpp)
* Literature:
  * "Local branching-based algorithms for the disjunctively constrained knapsack problem" (Akeb et al., 2011) [DOI](https://doi.org/10.1016/j.cie.2011.01.019)
  * "Bin Packing with Conflicts: A Generic Branch-and-Price Algorithm" (Sadykov et Vanderbeck, 2012) [DOI](https://doi.org/10.1287/ijoc.1120.0499)
  * "A Branch-and-Bound Algorithm for the Knapsack Problem with Conflict Graph" (Bettinelli et al., 2017) [DOI](https://doi.org/10.1287/ijoc.2016.0742)
  * "A new combinatorial branch-and-bound algorithm for the Knapsack Problem with Conflicts" (Coniglio et al., 2020) [DOI](https://doi.org/10.1016/j.ejor.2020.07.023)
* Algorithm:
  * Perturbation: force item `j` in/out of the solution
  * Local search neighborhoods:
    * Move item `j` in/out of the solution
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../treesearchdata/knapsackwithconflicts/data.csv -l knsapsackwithconflicts -t 60`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../treesearchdata/knapsackwithconflicts/data.csv -l knsapsackwithconflicts -b heuristiclong -t 62`

[Multidimensional Multiple-choice Knapsack Problem](examples/multidimensionalmultiplechoiceknapsack.hpp)

<details><summary>Literature</summary>
<p>

* Literature:
  * "Heuristic algorithms for the multiple-choice multidimensional knapsack problem" (Hifi et al., 2004) [DOI](https://doi.org/10.1057/palgrave.jors.2601796)
  * "A Reactive Local Search-Based Algorithm for the Multiple-Choice Multi-Dimensional Knapsack Problem" (Hifi et al., 2005) [DOI](https://doi.org/10.1007/s10589-005-3057-0)
  * "Development of core to solve the multidimensional multiple-choice knapsack problem" (Ghasemia et Razzaz, 2011) [DOI](https://doi.org/10.1016/j.cie.2010.12.001)
  * "Iterative semi-continuous relaxation heuristics for the multiple-choice multidimensional knapsack problem" (Crévits et al., 2012) [DOI](https://doi.org/10.1016/j.cor.2010.12.016)
  * "A hybrid heuristic for the multiple choice multidimensional knapsack problem" (Mansi et al., 2012) [DOI](https://doi.org/10.1080/0305215X.2012.717072)
  * "A fast and scalable multidimensional multiple-choice knapsack heuristic" (Shojaei et al., 2013) [DOI](https://doi.org/10.1145/2541012.2541014)
  * "A “reduce and solve” approach for the multiple-choice multidimensional knapsack problem" (Chen et Hao, 2014) [DOI](https://doi.org/10.1016/j.ejor.2014.05.025)
  * "Lagrangian heuristic-based neighbourhood search for the multiple-choice multi-dimensional knapsack problem" (Hifi et al., 2015) [DOI](https://doi.org/10.1080/0305215X.2014.982631)
  * "A set partitioning reformulation for the multiple-choice multidimensional knapsack problem" (Voß et Lalla-Ruiz, 2015) [DOI](https://doi.org/10.1080/0305215X.2015.1062094)
  * "A Core-Based Exact Algorithm for the Multidimensional Multiple Choice Knapsack Problem" (Mansini et Zanotti, 2020) [DOI](https://doi.org/10.1287/ijoc.2019.0909) 
  * "A two-phase kernel search variant for the multidimensional multiple-choice knapsack problem" (Lamanna et al., 2021) [DOI](https://doi.org/10.1016/j.ejor.2021.05.007)

</p>
</details>

* Algorithm:
  * Perturbation: force item `j` in the knapsack
  * Local search neighborhoods:
    * Add item `j` in the knapsack
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../localsearchdata/multidimensionalmultiplechoiceknapsack/data.csv -l multidimensionalmultiplechoiceknapsack -t 60`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../localsearchdata/multidimensionalmultiplechoiceknapsack/data.csv -l multidimensionalmultiplechoiceknapsack -b heuristiclong -t 62`

[Quadratic Assignment Problem](examples/quadraticassignment.hpp)

<details><summary>Literature</summary>
<p>

* Literature:
  * "Iterated local search for the quadratic assignment problem" (Stützle, 2006) [DOI](https://doi.org/10.1016/j.ejor.2005.01.066)
  * "A survey for the quadratic assignment problem" (Loiola et al., 2007) [DOI](https://doi.org/10.1016/j.ejor.2005.09.032)
  * "A branch-and-cut algorithm for quadratic assignment problems based on linearizations" (Erdoğan et Tansel, 2007) [DOI](https://doi.org/10.1016/j.cor.2005.05.027)
  * "Extensive experiments with hybrid genetic algorithms for the solution of the quadratic assignment problem" (Drezner, 2008) [DOI](https://doi.org/10.1016/j.cor.2006.05.004)
  * "A cooperative parallel tabu search algorithm for the quadratic assignment problem" (James et al., 2009) [DOI](https://doi.org/10.1016/j.ejor.2007.06.061)
  * "Multistart Tabu Search and Diversification Strategies for the Quadratic Assignment Problem" (James et al., 2009) [DOI](https://doi.org/10.1109/TSMCA.2009.2014556)
  * "An ejection chain algorithm for the quadratic assignment problem" (Rego et al., 2009) [DOI](https://doi.org/10.1002/net.20360)
  * "Breakout local search for the quadratic assignment problem" (Benlic et Hao, 2013) [DOI](https://doi.org/10.1016/j.amc.2012.10.106)
  * "Memetic search for the quadratic assignment problem" (Benlic et Hao, 2015) [DOI](https://doi.org/10.1016/j.eswa.2014.08.011)

</p>
</details>

* Algorithm:
  * Perturbation: swap two assignments
  * Local search neighborhoods:
    * Swap two assignments
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../localsearchdata/quadraticassignment/data.csv -l quadraticassignment -t 60`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../localsearchdata/quadraticassignment/data.csv -l quadraticassignment -b heuristiclong -t 62`

Generalized Assignment Problem from [fontanf/generalizedassignmentsolver](https://github.com/fontanf/generalizedassignmentsolver/blob/master/generalizedassignmentsolver/algorithms/localsearch.cpp)
* Algorithm:
  * Perturbation: force item `j` to be assigned to agent `i`
  * Local search neighborhoods:
    * Shift item `j` to agent `i`

### Scheduling

#### Single machine scheduling

[Single machine scheduling problem with sequence-dependent setup times, Total weighted Tardiness](examples/schedulingwithsdsttwt.hpp)
* Three field classification: `1 | sᵢⱼ | ∑wⱼTⱼ`

<details><summary>Literature</summary>
<p>

* Literature:
  * "Real-time scheduling of an automated manufacturing center" (Raman, 1989) [DOI](https://doi.org/10.1016/0377-2217(89)90332-9)
  * "A heuristic to minimize the total weighted tardiness with sequence-dependent setups" (Lee et al., 1997) [DOI](https://doi.org/10.1080/07408179708966311)
  * "Enhancing Stochastic Search Performance by Value-Biased Randomization of Heuristics" (Cicirello et Smith, 2005) [DOI](https://doi.org/10.1007/s10732-005-6997-8)
  * "Non-wrapping order crossover: an order preserving crossover operator that respects absolute position" (Cicirello, 2006) [DOI](https://doi.org/10.1145/1143997.1144177)
  * "An ant colony optimization for single-machine tardiness scheduling with sequence-dependent setups" (Liao et Juan, 2007) [DOI](https://doi.org/10.1016/j.cor.2005.07.020)
  * "Solving single-machine total weighted tardiness problems with sequence-dependent setup times by meta-heuristics" (Lin et Ying, 2007) [DOI](https://doi.org/10.1007/s00170-006-0693-1)
  * "Beam search algorithms for the single machine total weighted tardiness scheduling problem with sequence-dependent setups" (Valente et Alves, 2008) [DOI](https://doi.org/10.1016/j.cor.2006.11.004)
  * "A new discrete particle swarm optimization approach for the single-machine total weighted tardiness scheduling problem with sequence-dependent setup times" (Anghinolfi et Paolucci, 2009) [DOI](https://doi.org/10.1016/j.ejor.2007.10.044)
  * "A discrete differential evolution algorithm for the single machine total weighted tardiness problem with sequence dependent setup times" (Tasgetiren et al., 2009) [DOI](https://doi.org/10.1016/j.cor.2008.06.007)
  * "Parallel path relinking method for the single machine total weighted tardiness problem with sequence-dependent setups" (Bożejko, 2010) [DOI](https://doi.org/10.1007/s10845-009-0253-2)
  * "A variable neighborhood search for minimizing total weighted tardiness with sequence dependent setup times on a single machine" (Kirlik et Oguz, 2012) [DOI](https://doi.org/10.1016/j.cor.2011.08.022)
  * "A discrete electromagnetism-like mechanism for single machine total weighted tardiness problem with sequence-dependent setup times" (Chao et Liao, 2012) [DOI](https://www.sciencedirect.com/science/article/abs/pii/S1568494612002530)
  * "Neighborhood search procedures for single machine tardiness scheduling with sequence-dependent setups" (Liao et al., 2012) [DOI](https://doi.org/10.1016/j.tcs.2012.01.043)
  * "A hybrid genetic algorithm for the single machine scheduling problem with sequence-dependent setup times" (Sioud et al., 2012) [DOI](https://doi.org/10.1016/j.cor.2011.12.017)
  * "An exact algorithm for the single-machine total weighted tardiness problem with sequence-dependent setup times" (Tanaka et Araki, 2013) [DOI](https://doi.org/10.1016/j.cor.2012.07.004)
  * "Iterated Local Search for single-machine scheduling with sequence-dependent setup times to minimize total weighted tardiness" (Xu et al., 2014) [DOI](https://doi.org/10.1007/s10951-013-0351-z)
  * "A study of hybrid evolutionary algorithms for single machine scheduling problem with sequence-dependent setup times" (Xu et al., 2014) [DOI](https://doi.org/10.1016/j.cor.2014.04.009)
  * "An Iterated Local Search heuristic for the single machine total weighted tardiness scheduling problem with sequence-dependent setup times" (Subramanian et al., 2014) [DOI](https://doi.org/10.1080/00207543.2014.883472)
  * "An iterated greedy algorithm for the single-machine total weighted tardiness problem with sequence-dependent setup times" (Deng et Gu, 2014) [DOI](https://doi.org/10.1080/00207721.2012.723054)
  * "An improved scatter search algorithm for the single machine total weighted tardiness scheduling problem with sequence-dependent setup times" (Guo et Tang, 2015) [DOI](https://doi.org/10.1016/j.asoc.2014.12.030)
  * "Efficient local search limitation strategy for single machine total weighted tardiness scheduling with sequence-dependent setup times" (Subramanian et Farias., 2017) [DOI](https://doi.org/10.1016/j.cor.2016.10.008)

</p>
</details>

* Algorithm:
  * The algorithm is mainly based on the Iterated Local Search from "Efficient local search limitation strategy for single machine total weighted tardiness scheduling with sequence-dependent setup times" (Subramanian et Farias., 2017)
  * Perturbation: swap two blocs (double-bridge)
  * Local search neighborhoods:
    * Shift a bloc of `k` consecutive jobs, `k = 1..13`
    * Swap job `j1` and job `j2`
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../treesearchdata/schedulingwithsdsttwt/data.csv -l schedulingwithsdsttwt -t 60`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../treesearchdata/schedulingwithsdsttwt/data.csv -l schedulingwithsdsttwt -b heuristiclong -t 31`

#### Flow shop scheduling

[Permutation flow shop scheduling problem, Makespan](examples/permutationflowshopschedulingmakespan.hpp)
* Three field classification: `Fm | prmu | Cₘₐₓ`
* Literature:
  * "An Effective Hybrid Heuristic for Flow Shop Scheduling" (Zheng et Wang, 2003) [DOI](https://doi.org/10.1007/s001700300005)
  * "A simple and effective iterated greedy algorithm for the permutation flowshop scheduling problem" (Ruiz et Stützle, 2007) [DOI](https://doi.org/10.1016/j.ejor.2005.12.009)
  * "Cooperative metaheuristics for the permutation flowshop scheduling problem" (Vallada et Ruiz, 2009) [DOI](https://doi.org/10.1016/j.ejor.2007.11.049)
  * "A Variable Block Insertion Heuristic for Solving Permutation Flow Shop Scheduling Problem with Makespan Criterion" (Kizilay et al., 2019) [DOI](https://doi.org/10.3390/a12050100)
  * "A best-of-breed iterated greedy for the permutation flowshop scheduling problem with makespan objective" (Fernandez-Viagas, Framinan, 2019) [DOI](https://doi.org/10.1016/j.cor.2019.104767)
  * "A memetic algorithm with novel semi-constructive evolution operators for permutation flowshop scheduling problem" (Kurdi, 2020) [DOI](https://doi.org/10.1016/j.asoc.2020.106458)
  * "Iterative beam search algorithms for the permutation flowshop" (Libralesso et al., 2020)
* Algorithm:
  * Perturbation: swap two blocs (double-bridge)
  * Local search neighborhoods:
    * Move `k` consecutive jobs, `k = 1..4`
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../treesearchdata/permutationflowshopscheduling/data_makespan.csv -f "row['Dataset'] == 'vallada2015_large'" -l permutationflowshopschedulingmakespan --timelimitfield "Time limit"`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../treesearchdata/permutationflowshopscheduling/data_makespan.csv -f "row['Dataset'] == 'vallada2015_large' and int(row['Job number']) <= 100" -l permutationflowshopschedulingmakespan -b heuristiclong -t 500`

[Permutation flow shop scheduling problem, Total tardiness](examples/permutationflowshopschedulingtt.hpp)
* Three field classification: `Fm | prmu | ∑Tⱼ`
* Literature:
  * "Minimising total tardiness in the m-machine flowshop problem: A review and evaluation of heuristics and metaheuristics"  (Vallada et al., 2008) [DOI](https://doi.org/10.1016/j.cor.2006.08.016)
  * "Cooperative metaheuristics for the permutation flowshop scheduling problem" (Vallada et Ruiz, 2009) [DOI](https://doi.org/10.1016/j.ejor.2007.11.049)
  * "Genetic algorithms with path relinking for the minimum tardiness permutation flowshop problem" (Vallada et Ruiz, 2010) [DOI](https://doi.org/10.1016/j.omega.2009.04.002)
  * "NEH-based heuristics for the permutation flowshop scheduling problem to minimise total tardiness" (Fernandez-Viagas et Framinan, 2015) [DOI](https://doi.org/10.1016/j.cor.2015.02.002)
  * "Matheuristic algorithms for minimizing total tardiness in the m-machine flow-shop scheduling problem" (Ta et al., 2015) [DOI](https://doi.org/10.1007/s10845-015-1046-4)
  * "Iterated-greedy-based algorithms with beam search initialization for the permutation flowshop to minimise total tardiness" (Fernandez-Viagas et al., 2018) [DOI](https://doi.org/10.1016/j.eswa.2017.10.050)
* Algorithm:
  * Perturbation: force job `j` to start first or right after job `j_prev`
  * Local search neighborhoods:
    * Move `k` consecutive jobs, `k = 1..4`
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../treesearchdata/permutationflowshopscheduling/data_totaltardiness.csv -l permutationflowshopschedulingtt --timelimitfield "Time limit"`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../treesearchdata/permutationflowshopscheduling/data_totaltardiness.csv -l permutationflowshopschedulingtt -b heuristiclong -t 500`

#### Resource constrained scheduling

[ROADEF/EURO Challenge 2020: Maintenance Planning Problem](examples/roadef2020.hpp)
* Website: https://www.roadef.org/challenge/2020/en/
* Algorithm:
  * Perturbation: force intervention `j` to start at time `t_start`
  * Local search neighborhoods:
    * Shift intervention `j` to time `t_start`
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --main "./bazel-bin/examples/roadef2020_main" --csv "../localsearchdata/roadef2020/data.csv" -f "'A' not in row['Dataset']" -l roadef2020 -t 900`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py -b heuristiclong --csv ../localsearchdata/roadef2020/data.csv -l roadef2020 -t 920 -f "'A' not in row['Dataset']"`

## Usage, running examples from command line

[Download data](https://github.com/fontanf/treesearchsolver/releases/download/data/treesearchdata.7z)

Compile:
```shell
bazel build -- //...
```

Then, examples can be executed as follows:
```shell
./bazel-bin/examples/main -v -p knapsackwithconflicts -i ../treesearchdata/knapsackwithconflicts/C1/BPPC_1_0_1.txt_0.1 -t 5
```

## Usage, C++ library

See examples.

