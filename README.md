# Local Search Solver

A solver based on local search.

## Description

The goal of this repository is to provide a simple framework to quickly implement heuristic algorithms based on local search.

For complex and evolving problems, implementing and extending local search algorithms based on complex neighborhoods quickly becomes cumbersome and time consuming. This often makes them unsuitable for practical use.
The algorithms of this repository are designed to get the best performances out of the simplest neighborhoods. Thus making them easier to implement and extend.

Still, more complex neighborhoods can be implemented if better performances are needed.

All algorithm require the Local Scheme to implement: `GlobalCost`, `Solution`, `initial_solution`, `global_cost` and `local_search`.

Implemented algorithms:
* Restarting Local Search `-a restarting_local_search`
  * Most basic algorithm, can be used as a baseline
* Iterated Local Search `-a iterated_local_search`
  * Single solution algorithm
  * Low implementation cost
  * Good for very large problems
  * Additional requirements: `Move`, `perturbations`, `apply_move`
* Best First Local Search `-a best_first_local_search`
  * Multi-solution algorithm
  * Higher implementation cost
  * Higher quality solutions
  * Additional requirements: `CompactSolution`, `compact_solution_hasher`, `solution2compact`, `compact2solution`, `Move`, `move_hasher`, `perturbations`, `apply_move`
* Genetic Local Search `-a genetic_local_search`
  * Population-based algorithm
  * Good when the ruggedness of the landscape is high
  * Additional requirements: `crossover` and `distance`

## Sequencing module

A specific implementation is also available for sequencing problems. The neighborhoods, crossovers and perturbations are already implemented, it is only required to provide an `append(solution, j)` method to use them:
* Local search neighborhoods:
  * Shift a block of `k` consecutive jobs
  * Swap a block of `k1` consecutive jobs with another block of `k2` consecutive jobs
  * Reverse a block of consecutive jobs
  * Shift and reverse a block of `k` consecutive jobs
  * Remove a job from the solution / Add a job into the solution
* Perturbations:
  * Swap two blocks of consecutive jobs (double-bridge)
  * Remove `k` jobs and re-insert them (ruin-and-recreate)
  * Force a job into the solution
* Crossover algorithms:
  * OX crossover
  * SJOX crossover
  * SBOX crossover

### Examples

[Sequential Ordering Problem](examples/sequentialordering.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/optimizationtools/bench_run.py --main ./bazel-bin/examples/sequentialordering_main --csv ../ordata/sequentialordering/data.csv -f "row['Dataset'] == 'tsplib'" -l "${DATE}_sequentialordering" -t 60
python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/sequentialordering/data.csv -f "row['Dataset'] == 'tsplib'" -l "${DATE}_sequentialordering" -b heuristiclong -t 62
```

![sequentialordering](img/sequentialordering_tsplib.png?raw=true "sequentialordering_tsplib")

</p>
</details>

[Single machine scheduling problem with sequence-dependent setup times, Total weighted tardiness](examples/schedulingwithsdsttwt.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/optimizationtools/bench_run.py --main ./bazel-bin/examples/schedulingwithsdsttwt_main --csv ../ordata/schedulingwithsdsttwt/data.csv -l "${DATE}_schedulingwithsdsttwt" -t 60
python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/schedulingwithsdsttwt/data.csv -l "${DATE}_schedulingwithsdsttwt" -b heuristiclong -t 62
```

![schedulingwithsdsttwt](img/schedulingwithsdsttwt.png?raw=true "schedulingwithsdsttwt")

</p>
</details>

[Permutation flow shop scheduling problem, Total completion time](examples/permutationflowshopschedulingtct.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/optimizationtools/bench_run.py --main ./bazel-bin/examples/permutationflowshopschedulingtct_main --csv ../ordata/flowshopscheduling/data_permutationflowshopschedulingtct.csv -l "${DATE}_permutationflowshopschedulingtct" --timelimitfield "Time limit"
python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/flowshopscheduling/data_permutationflowshopschedulingtct.csv -l "${DATE}_permutationflowshopschedulingtct" -b heuristiclong -t 500
```

![permutationflowshopschedulingtct](img/permutationflowshopschedulingtct.png?raw=true "permutationflowshopschedulingtct")

</p>
</details>

[Permutation flow shop scheduling problem, Total tardiness](examples/permutationflowshopschedulingtt.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/optimizationtools/bench_run.py --main ./bazel-bin/examples/permutationflowshopschedulingtt_main --csv ../ordata/flowshopscheduling/data_permutationflowshopschedulingtt.csv -l "${DATE}_permutationflowshopschedulingtt" --timelimitfield "Time limit"
python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/flowshopscheduling/data_permutationflowshopschedulingtt.csv -l "${DATE}_permutationflowshopschedulingtt" -b heuristiclong -t 2200
```

![permutationflowshopschedulingtt](img/permutationflowshopschedulingtt.png?raw=true "permutationflowshopschedulingtt")

</p>
</details>

[Single machine order acceptance and scheduling problem with time windows and sequence-dependent setup times, Total weighted tardiness](examples/orderacceptanceandscheduling.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
python3 ../optimizationtools/optimizationtools/bench_run.py --main ./bazel-bin/examples/orderacceptanceandscheduling_main --csv ../ordata/orderacceptanceandscheduling/data.csv -l orderacceptanceandscheduling -t 60
python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/orderacceptanceandscheduling/data.csv -l orderacceptanceandscheduling -b heuristiclong -t 61
```

</p>
</details>

[Time-dependent orienteering problem](examples/timedependentorienteering.hpp)

## Other examples

Data can be downloaded from [fontanf/orproblems](https://github.com/fontanf/orproblems)

### Packing

[Multidimensional Multiple-choice Knapsack Problem](examples/multidimensionalmultiplechoiceknapsack.hpp)
* Straightforward example: single neighborhood, simple perturbation, basic operators
* Algorithm:
  * Local search neighborhoods:
    * Add item `j` in the knapsack
  * Perturbation: force item `j` in the knapsack
  * Crossover algorithm

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/optimizationtools/bench_run.py --main ./bazel-bin/examples/multidimensionalmultiplechoiceknapsack_main --csv ../ordata/multidimensionalmultiplechoiceknapsack/data.csv -l "${DATE}_multidimensionalmultiplechoiceknapsack" -t 60
python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/multidimensionalmultiplechoiceknapsack/data.csv -l "${DATE}_multidimensionalmultiplechoiceknapsack" -b heuristiclong -t 62
```

![multidimensionalmultiplechoiceknapsack](img/multidimensionalmultiplechoiceknapsack.png?raw=true "multidimensionalmultiplechoiceknapsack")

</p>
</details>

[Quadratic Assignment Problem](examples/quadraticassignment.hpp)
* Example which implements a problem specific acceleration strategy to compute the move costs
* Algorithm:
  * Local search neighborhood: swap two assignments
  * Perturbation: ejection chain
  * Crossover algorithm: UX crossover

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/optimizationtools/bench_run.py --main ./bazel-bin/examples/quadraticassignment_main --csv ../ordata/quadraticassignment/data.csv -l "${DATE}_quadraticassignment" -t 60
python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/quadraticassignment/data.csv -l "${DATE}_quadraticassignment" -b heuristiclong -t 62
```

![quadraticassignment](img/quadraticassignment.png?raw=true "quadraticassignment")

</p>
</details>

[Knapsack Problem with Conflicts](examples/knapsackwithconflicts.hpp)
* Example with multiple neigborhoods
* Algorithm:
  * Local search neighborhoods:
    * Move item `j` in/out of the solution (and remove conflicting items)
    * Swap two non-conflicting items
    * Remove one item and add two non-conflicting items from its neighbors
  * Perturbation: force item `j` in/out of the solution
  * Crossover algorithm: double backbone-based crossover

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/optimizationtools/bench_run.py --main ./bazel-bin/examples/knapsackwithconflicts_main --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'hifi2006'" -l "${DATE}_knapsackwithconflicts_hifi2006" -t 60
python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'hifi2006'" -l "${DATE}_knapsackwithconflicts_hifi2006" -b heuristiclong -t 62
```

![knapsackwithconflicts](img/knapsackwithconflicts_hifi2006.png?raw=true "knapsackwithconflicts_hifi2006")

</p>
</details>

Generalized Assignment Problem from [fontanf/generalizedassignmentsolver](https://github.com/fontanf/generalizedassignmentsolver/blob/master/generalizedassignmentsolver/algorithms/localsearch.cpp)
* Example which implements a generic strategy for very large problems to avoid recomputing moves which have not change since the last neighborhood exploration
* Algorithm:
  * Local search neighborhoods: shift item `j` to agent `i`
  * Perturbation: shift 8 random jobs

### Routing

[Travelling Salesman Problem](examples/travellingsalesman.hpp)
* Three field classification: `1 | sᵢⱼ | Cₘₐₓ`
* Algorithm:
  * Local search neighborhoods:
    * Shift a block of `k` consecutive vertices, `k = 1..8` (or-opt)
    * Swap vertex `j1` with vertex `j2`
    * Replace edges `(j1, j2)` and `(j3, j4)` by edges `(j1, j3)` and `(j2, j4)` (2-opt)
  * Perturbation: swap two blocks (double-bridge)
  * The local search implementation avoids recomputing moves which have not change since the last neighborhood exploration

<details><summary>Benchmarks</summary>
<p>

```shell
python3 ../optimizationtools/optimizationtools/bench_run.py --main ./bazel-bin/examples/travellingsalesman_main --csv ../ordata/travellingsalesman/data.csv -l travellingsalesman -t 60
python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/travellingsalesman/data.csv -l travellingsalesman -b heuristiclong -t 62
```

</p>
</details>

### Scheduling

#### Flow shop scheduling

[Permutation flow shop scheduling problem, Makespan](examples/permutationflowshopschedulingmakespan.hpp)
* This one is not considered as a sequencing problem since the dedicated acceleration strategy makes it possible to explore the neighborhoods more efficiently
* Algorithm:
  * Local search neighborhood: move a block of `k` consecutive jobs, `k = 1..4`
  * Perturbation: swap two blocks (double-bridge)

<details><summary>Benchmarks</summary>
<p>

```shell
python3 ../optimizationtools/optimizationtools/bench_run.py --main ./bazel-bin/examples/permutationflowshopschedulingmakespan_main --csv ../ordata/permutationflowshopscheduling/data_makespan.csv -f "row['Dataset'] == 'vallada2015_large'" -l permutationflowshopschedulingmakespan --timelimitfield "Time limit"
python3 ../optimizationtools/optimizationtools/bench_process.py --csv ../ordata/permutationflowshopscheduling/data_makespan.csv -f "row['Dataset'] == 'vallada2015_large' and int(row['Job number']) <= 100" -l permutationflowshopschedulingmakespan -b heuristiclong -t 500
```

</p>
</details>

#### Resource constrained scheduling

[ROADEF/EURO Challenge 2020: Maintenance Planning Problem](examples/roadef2020.hpp)
* Example of a problem with expensive move evaluations
* Algorithm:
  * Local search neighborhood: shift intervention `j` to time `t_start`
  * Perturbation: force intervention `j` to start at time `t_start`

<details><summary>Benchmarks</summary>
<p>

```shell
python3 ../optimizationtools/optimizationtools/bench_run.py --main "./bazel-bin/examples/roadef2020_main -w 0 -y 1 " --csv ../ordata/roadef2020/data/data.csv -l roadef2020 -t 900 -f "'A' not in row['Dataset']"
python3 ../optimizationtools/optimizationtools/bench_process.py -b heuristiclong --csv ../ordata/roadef2020/data.csv -l roadef2020 -t 920 -f "'A' not in row['Dataset']"
```

</p>
</details>

### Graphs

Maximum-Weight Independent Set Problem from [fontanf/stablesolver](https://github.com/fontanf/stablesolver/blob/master/stablesolver/algorithms/localsearch.cpp)
* Algorithm:
  * Local search neighborhoods:
    * Move vertex `v` in/out of the solution (and remove conflicting vertices)
    * Remove one vertex and add two non-conflicting vertices from its neighbors
  * Perturbation: force vertex `v` in/out of the solution

Graph Coloring Problem from [fontanf/coloringsolver](https://github.com/fontanf/coloringsolver/blob/master/coloringsolver/algorithms/localsearch.cpp)
* Example of a problem with a granular objective
* Algorithm:
  * Local search neighborhood: change the color of a conflicting vertex
  * Perturbation:
    * If the solution is feasible, merge two colors
    * If the solution is infeasible, change the color of a conflicting vertex

## Usage, running examples from command line

Compile:
```shell
bazel build -- //...
```

Then, examples can be executed as follows:
```shell
./bazel-bin/examples/knapsackwithconflicts_main -v -i ../ordata/knapsackwithconflicts/bettinelli2017/sparse_corr/test_1000_1000_r0.001-0.dat -f bettinelli2017 -t 10 -c sol.txt
```

## Usage, C++ library

See examples.

