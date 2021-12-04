# Local Search Solver

A solver based on local search.

## Description

The goal of this repository is to provide a simple framework to quickly implement heuristic algorithms based on local search.

For complex and evolving problems, implementing and extending local search algorithms based on complex neighborhoods quickly becomes cumbersome and time consuming. This often makes them unsuitable for practical use.
The algorithms of this repository are designed to get the best performances out of the simplest neighborhoods. Thus making them easier to implement and extend.

Still, more complex neighborhoods can be implemented if better performances are needed.

All algorithm require the Local Scheme to implement: `GlobalCost`, `global_cost_worst`, `Solution`, `initial_solution`, `global_cost` and `local_search`.

Implemented algorithms:
* Restarting Local Search `-a restarting_local_search`
  * Most basic algorithm, can be used as a baseline
* Iterated Local Search `-a iterated_local_search`
  * Single solution algorithm
  * Low implementation cost
  * Good for very large problems
  * Additional requirements: `Move`, `move_null`, `perturbations`, `apply_move`
* Best First Local Search `-a best_first_local_search`
  * Multi-solution algorithm
  * Higher implementation cost
  * Higher quality solutions
  * Additional requirements: `CompactSolution`, `compact_solution_hasher`, `solution2compact`, `compact2solution`, `Move`, `move_null`, `move_hasher`, `perturbations`, `apply_move`
* Genetic Local Search `-a genetic_local_search`
  * Population-based algorithm
  * Good when the ruggedness of the landscape is high
  * Additional requirements: `crossover` and `distance`

A specific implementation is also available for sequencing problems. The neighborhoods, crossovers and perturbations are already implemented, it is only required to provide an `append(solution, j)` method to use them:
* Perturbation: swap two blocs (double-bridge)
* Crossover algorithms:
  * OX crossover
  * SJOX crossover
  * SBOX crossover
* Local search neighborhoods:
  * Shift a bloc of `k` consecutive jobs
  * Swap a bloc of `k1` consecutive jobs with another bloc of `k2` consecutive jobs
  * Reverse a bloc of consecutive jobs
  * Shift and reverse a bloc of `k` consecutive jobs

## Examples

Data can be downloaded from [fontanf/orproblems](https://github.com/fontanf/orproblems)

### Packing

[Multidimensional Multiple-choice Knapsack Problem](examples/multidimensionalmultiplechoiceknapsack.hpp)
* Straightforward example: single neighborhood, simple perturbation, basic operators
* Algorithm:
  * Perturbation: force item `j` in the knapsack
  * Local search neighborhoods:
    * Add item `j` in the knapsack

<details><summary>Benchmarks</summary>
<p>

* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/multidimensionalmultiplechoiceknapsack/data.csv -l multidimensionalmultiplechoiceknapsack -t 60`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/multidimensionalmultiplechoiceknapsack/data.csv -l multidimensionalmultiplechoiceknapsack -b heuristiclong -t 62`

</p>
</details>

[Quadratic Assignment Problem](examples/quadraticassignment.hpp)
* Example which implements a problem specific acceleration strategy to compute the move costs
* Algorithm:
  * Perturbation: ejection chain
  * Crossover algorithm: UX crossover
  * Local search neighborhoods:
    * Swap two assignments

<details><summary>Benchmarks</summary>
<p>

* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/quadraticassignment/data.csv -l quadraticassignment --timelimitfield "Time limit" -a "astarlocalsearch -x 6"`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/quadraticassignment/data.csv -l quadraticassignment -b heuristiclong -t 185`

</p>
</details>

[Knapsack Problem with Conflicts](examples/knapsackwithconflicts.hpp)
* Example with multiple neigborhoods
* Algorithm:
  * Perturbation: force item `j` in/out of the solution
  * Crossover algorithm: double backbone-based crossover
  * Local search neighborhoods:
    * Move item `j` in/out of the solution (and remove conflicting items)
    * Swap two non-conflicting items
    * Remove one item and add two non-conflicting items from its neighbors

<details><summary>Benchmarks</summary>
<p>

* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'hifi2006'" -l knapsackwithconflicts -t 300`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'hifi2006'" -l knapsackwithconflicts -b heuristiclong -t 310`
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'bettinelli2017'" -l knapsackwithconflicts -t 5`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'bettinelli2017'" -l knapsackwithconflicts -b heuristiclong -t 6`

</p>
</details>

Generalized Assignment Problem from [fontanf/generalizedassignmentsolver](https://github.com/fontanf/generalizedassignmentsolver/blob/master/generalizedassignmentsolver/algorithms/localsearch.cpp)
* Example which implements a generic strategy for very large problems to avoid recomputing moves which have not change since the last neighborhood exploration
* Algorithm:
  * Perturbation: shift 8 random jobs
  * Local search neighborhoods:
    * Shift item `j` to agent `i`

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

<details><summary>Benchmarks</summary>
<p>

* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/travellingsalesman/data.csv -l travellingsalesman -t 60`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/travellingsalesman/data.csv -l travellingsalesman -b heuristiclong -t 62`

</p>
</details>

### Scheduling

#### Single machine scheduling

[Single machine scheduling problem with sequence-dependent setup times, Total weighted tardiness](examples/schedulingwithsdsttwt.hpp)
* Example of a sequencing problem

<details><summary>Benchmarks</summary>
<p>

* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/schedulingwithsdsttwt/data.csv -l schedulingwithsdsttwt -t 60`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/schedulingwithsdsttwt/data.csv -l schedulingwithsdsttwt -b heuristiclong -t 31`

</p>
</details>

#### Flow shop scheduling

[Permutation flow shop scheduling problem, Makespan](examples/permutationflowshopschedulingmakespan.hpp)
* This one is not considered as a sequencing problem since the dedicated acceleration strategy makes it possible to explore the neighborhoods more efficiently
* Algorithm:
  * Perturbation: swap two blocs (double-bridge)
  * Local search neighborhoods:
    * Move `k` consecutive jobs, `k = 1..4`

<details><summary>Benchmarks</summary>
<p>

* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/permutationflowshopscheduling/data_makespan.csv -f "row['Dataset'] == 'vallada2015_large'" -l permutationflowshopschedulingmakespan --timelimitfield "Time limit"`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/permutationflowshopscheduling/data_makespan.csv -f "row['Dataset'] == 'vallada2015_large' and int(row['Job number']) <= 100" -l permutationflowshopschedulingmakespan -b heuristiclong -t 500`

</p>
</details>

[Permutation flow shop scheduling problem, Total tardiness](examples/permutationflowshopschedulingtt.hpp)
* Example of a sequencing problem

<details><summary>Benchmarks</summary>
<p>

* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --csv ../ordata/permutationflowshopscheduling/data_totaltardiness.csv -l permutationflowshopschedulingtt --timelimitfield "Time limit"`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/permutationflowshopscheduling/data_totaltardiness.csv -l permutationflowshopschedulingtt -b heuristiclong -t 500`

</p>
</details>

#### Resource constrained scheduling

[ROADEF/EURO Challenge 2020: Maintenance Planning Problem](examples/roadef2020.hpp)
* Example of a problem with expensive move evaluations
* Algorithm:
  * Perturbation: force intervention `j` to start at time `t_start`
  * Local search neighborhoods:
    * Shift intervention `j` to time `t_start`

<details><summary>Benchmarks</summary>
<p>

* Benchmarks:
  * `python3 ../optimizationtools/optimizationtools/bench_run.py --main "./bazel-bin/examples/roadef2020_main -w 0 -y 1 " --csv ../ordata/roadef2020/data/data.csv -l roadef2020 -t 900 -f "'A' not in row['Dataset']"`
  * `python3 ../optimizationtools/optimizationtools/bench_process.py -b heuristiclong --csv ../ordata/roadef2020/data.csv -l roadef2020 -t 920 -f "'A' not in row['Dataset']"`

</p>
</details>

### Graphs

Maximum-Weight Independent Set Problem from [fontanf/stablesolver](https://github.com/fontanf/stablesolver/blob/master/stablesolver/algorithms/localsearch.cpp)
* Algorithm:
  * Perturbation: force vertex `v` in/out of the solution
  * Local search neighborhoods:
    * Move vertex `v` in/out of the solution (and remove conflicting vertices)
    * Remove one vertex and add two non-conflicting vertices from its neighbors

Graph Coloring Problem from [fontanf/coloringsolver](https://github.com/fontanf/coloringsolver/blob/master/coloringsolver/algorithms/localsearch.cpp)
* Example of a problem with a granular objective
* Algorithm:
  * Perturbation:
    * If the solution is feasible, merge two colors
    * If the solution is infeasible, change the color of a conflicting vertex
  * Local search neighborhoods:
    * Change the color of a conflicting vertex

## Usage, running examples from command line

Compile:
```shell
bazel build -- //...
```

Then, examples can be executed as follows:
```shell
./bazel-bin/examples/main -v -p knapsackwithconflicts -i ../ordata/knapsackwithconflicts/bettinelli2017/C1/BPPC_1_0_1.txt_0.1 -f bettinelli2017 -t 5
```

## Usage, C++ library

See examples.

