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
* Algorithm:
  * Perturbation: force item `j` in the knapsack
  * Local search neighborhoods:
    * Add item `j` in the knapsack
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../localsearchdata/multidimensionalmultiplechoiceknapsack/data.csv -l multidimensionalmultiplechoiceknapsack -t 60`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../localsearchdata/multidimensionalmultiplechoiceknapsack/data.csv -l multidimensionalmultiplechoiceknapsack -b heuristiclong -t 62`

[Quadratic Assignment Problem](examples/quadraticassignment.hpp)
* Literature:
  * "Iterated local search for the quadratic assignment problem" (Stützle, 2006) [DOI](https://doi.org/10.1016/j.ejor.2005.01.066)
  * "A cooperative parallel tabu search algorithm for the quadratic assignment problem" (James et al., 2009) [DOI](https://doi.org/10.1016/j.ejor.2007.06.061)
  * "Multistart Tabu Search and Diversification Strategies for the Quadratic Assignment Problem" (James et al., 2009) [DOI](https://doi.org/10.1109/TSMCA.2009.2014556)
  * "An ejection chain algorithm for the quadratic assignment problem" (Rego et al., 2009) [DOI](https://doi.org/10.1002/net.20360)
  * "Breakout local search for the quadratic assignment problem" (Benlic et Hao, 2013) [DOI](https://doi.org/10.1016/j.amc.2012.10.106)
  * "Memetic search for the quadratic assignment problem" (Benlic et Hao, 2015) [DOI](https://doi.org/10.1016/j.eswa.2014.08.011)
* Algorithm:
  * Perturbation: swap two assignments
  * Local search neighborhoods:
    * Swap two assignments
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../localsearchdata/quadraticassignment/data.csv -l quadraticassignment -t 60`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../localsearchdata/quadraticassignment/data.csv -l quadraticassignment -b heuristiclong -t 62`

Generalized Assignment Problem from [fontanf/generalizedassignmentsolver](https://github.com/fontanf/generalizedassignmentsolver/blob/master/generalizedassignmentsolver/algorithms/localsearch.cpp)
* Algorithm:
  * Perturbation: for item `j` to be assigned to agent `i`
  * Local search neighborhoods:
    * Shift item `j` to agent `i`

### Scheduling

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
  * Perturbation: force job `j` to start first or right after job `j_prev`
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

