# Local Search Solver

A solver based on local search.

## Description

The goal of this repository is to provide a simple framework to quickly implement heuristic algorithms based on local search.

For complex and evolving problems, implementing and extending local search algorithms based on complex neighborhoods quickly becomes cumbersome and time consuming. This often makes them unsuitable for practical use.
The algorithms of this repository are designed to get the best performances out of the simplest neighborhoods. Thus making them easier to implement and extend.

Still, more complex neighborhoods can be implemented if better performances are needed.

In addition, the algorithms don't require parameter tuning.

## Examples

Data can be downloaded from [fontanf/orproblems](https://github.com/fontanf/orproblems)

### Packing

[Knapsack Problem with Conflicts](examples/knapsackwithconflicts.hpp)
* Algorithm:
  * Perturbation: force item `j` in/out of the solution
  * Local search neighborhoods:
    * Move item `j` in/out of the solution (and remove conflicting items)
    * Swap two non-conflicting items
    * Remove one item and add two non-conflicting items from its neighbors
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'hifi2006'" -l knapsackwithconflicts -t 300`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'hifi2006'" -l knapsackwithconflicts -b heuristiclong -t 310`
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'bettinelli2017'" -l knapsackwithconflicts -t 5`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'bettinelli2017'" -l knapsackwithconflicts -b heuristiclong -t 6`

[Multidimensional Multiple-choice Knapsack Problem](examples/multidimensionalmultiplechoiceknapsack.hpp)
* Algorithm:
  * Perturbation: force item `j` in the knapsack
  * Local search neighborhoods:
    * Add item `j` in the knapsack
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/multidimensionalmultiplechoiceknapsack/data.csv -l multidimensionalmultiplechoiceknapsack -t 60`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/multidimensionalmultiplechoiceknapsack/data.csv -l multidimensionalmultiplechoiceknapsack -b heuristiclong -t 62`

[Quadratic Assignment Problem](examples/quadraticassignment.hpp)
* Algorithm:
  * Perturbation: ejection chain
  * Local search neighborhoods:
    * Swap two assignments
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/quadraticassignment/data.csv -l quadraticassignment --timelimitfield "Time limit" -a "astarlocalsearch -x 6"`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/quadraticassignment/data.csv -l quadraticassignment -b heuristiclong -t 185`

Generalized Assignment Problem from [fontanf/generalizedassignmentsolver](https://github.com/fontanf/generalizedassignmentsolver/blob/master/generalizedassignmentsolver/algorithms/localsearch.cpp)
* Algorithm:
  * Perturbation: shift 8 random jobs
  * Local search neighborhoods:
    * Shift item `j` to agent `i`
  * The local search implementation avoids recomputing moves which have not change since the last neighborhood exploration

### Routing

[Travelling Salesman Problem](examples/travellingsalesman.hpp)
* Three field classification: `1 | sᵢⱼ | Cₘₐₓ`
* Algorithm:
  * Perturbation: swap two blocs (double-bridge)
  * Local search neighborhoods:
    * Shift a bloc of `k` consecutive vertices, `k = 1..8` (or-opt)
    * Swap vertex `j1` with vertex `j2`
    * Replace edges `(j1, j2)` and `(j3, j4)` by edges `(j1, j3)` and `(j2, j4)` (2-opt)
  * The local search implementation avoids recomputing moves which have not change since the last neighborhood exploration
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/travellingsalesman/data.csv -l travellingsalesman -t 60`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/travellingsalesman/data.csv -l travellingsalesman -b heuristiclong -t 62`

### Scheduling

#### Single machine scheduling

[Single machine scheduling problem with sequence-dependent setup times, Total weighted Tardiness](examples/schedulingwithsdsttwt.hpp)
* Algorithm:
  * The algorithm is mainly based on the Iterated Local Search from "Efficient local search limitation strategy for single machine total weighted tardiness scheduling with sequence-dependent setup times" (Subramanian et Farias., 2017)
  * Perturbation: swap two blocs (double-bridge)
  * Local search neighborhoods:
    * Shift a bloc of `k` consecutive jobs, `k = 1..13`
    * Swap job `j1` and job `j2`
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/schedulingwithsdsttwt/data.csv -l schedulingwithsdsttwt -t 60`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/schedulingwithsdsttwt/data.csv -l schedulingwithsdsttwt -b heuristiclong -t 31`

#### Flow shop scheduling

[Permutation flow shop scheduling problem, Makespan](examples/permutationflowshopschedulingmakespan.hpp)
* Algorithm:
  * Perturbation: swap two blocs (double-bridge)
  * Local search neighborhoods:
    * Move `k` consecutive jobs, `k = 1..4`
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/permutationflowshopscheduling/data_makespan.csv -f "row['Dataset'] == 'vallada2015_large'" -l permutationflowshopschedulingmakespan --timelimitfield "Time limit"`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/permutationflowshopscheduling/data_makespan.csv -f "row['Dataset'] == 'vallada2015_large' and int(row['Job number']) <= 100" -l permutationflowshopschedulingmakespan -b heuristiclong -t 500`

[Permutation flow shop scheduling problem, Total tardiness](examples/permutationflowshopschedulingtt.hpp)
* Algorithm:
  * Perturbation: swap two blocs (double-bridge)
  * Local search neighborhoods:
    * Shift a bloc of `k` consecutive jobs, `k = 1..3`
    * Swap job `j1` and job `j2`
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/permutationflowshopscheduling/data_totaltardiness.csv -l permutationflowshopschedulingtt --timelimitfield "Time limit"`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/permutationflowshopscheduling/data_totaltardiness.csv -l permutationflowshopschedulingtt -b heuristiclong -t 500`

#### Resource constrained scheduling

[ROADEF/EURO Challenge 2020: Maintenance Planning Problem](examples/roadef2020.hpp)
* Algorithm:
  * Perturbation: force intervention `j` to start at time `t_start`
  * Local search neighborhoods:
    * Shift intervention `j` to time `t_start`
* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --main "./bazel-bin/examples/roadef2020_main -w 0 -y 1 " --csv ../ordata/roadef2020/data/data.csv -l roadef2020 -t 900 -f "'A' not in row['Dataset']"`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py -b heuristiclong --csv ../ordata/roadef2020/data.csv -l roadef2020 -t 920 -f "'A' not in row['Dataset']"`

### Graphs

Maximum-Weight Independent Set Problem from [fontanf/stablesolver](https://github.com/fontanf/stablesolver/blob/master/stablesolver/algorithms/localsearch.cpp)
* Algorithm:
  * Perturbation: force vertex `v` in/out of the solution
  * Local search neighborhoods:
    * Move vertex `v` in/out of the solution (and remove conflicting vertices)
    * Remove one vertex and add two non-conflicting vertices from its neighbors

## Usage, running examples from command line

Compile:
```shell
bazel build -- //...
```

Then, examples can be executed as follows:
```shell
./bazel-bin/examples/main -v -p knapsackwithconflicts -i ../ordata/knapsackwithconflicts/C1/BPPC_1_0_1.txt_0.1 -t 5
```

## Usage, C++ library

See examples.

