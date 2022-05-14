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
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/sequentialordering_main --csv ../ordata/sequentialordering/data.csv -f "row['Dataset'] == 'tsplib'" -l "${DATE}_sequentialordering" -t 60
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/sequentialordering/data.csv -f "row['Dataset'] == 'tsplib'" -l "${DATE}_sequentialordering" -b heuristiclong -t 62
```

![sequentialordering](img/sequentialordering_tsplib.png?raw=true "sequentialordering_tsplib")

</p>
</details>

[Single machine scheduling problem with sequence-dependent setup times, Total weighted tardiness](examples/schedulingwithsdsttwt.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/schedulingwithsdsttwt_main --csv ../ordata/schedulingwithsdsttwt/data.csv -l "${DATE}_schedulingwithsdsttwt" -t 60
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/schedulingwithsdsttwt/data.csv -l "${DATE}_schedulingwithsdsttwt" -b heuristiclong -t 62
```

![schedulingwithsdsttwt](img/schedulingwithsdsttwt.png?raw=true "schedulingwithsdsttwt")

</p>
</details>

[Permutation flow shop scheduling problem, Total completion time](examples/permutationflowshopschedulingtct.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/permutationflowshopschedulingtct_main --csv ../ordata/flowshopscheduling/data_permutationflowshopschedulingtct.csv -l "${DATE}_permutationflowshopschedulingtct" --timelimitfield "Time limit"
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/flowshopscheduling/data_permutationflowshopschedulingtct.csv -l "${DATE}_permutationflowshopschedulingtct" -b heuristiclong -t 500
```

![permutationflowshopschedulingtct](img/permutationflowshopschedulingtct.png?raw=true "permutationflowshopschedulingtct")

</p>
</details>

[Permutation flow shop scheduling problem, Total tardiness](examples/permutationflowshopschedulingtt.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/permutationflowshopschedulingtt_main --csv ../ordata/flowshopscheduling/data_permutationflowshopschedulingtt.csv -l "${DATE}_permutationflowshopschedulingtt" --timelimitfield "Time limit"
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/flowshopscheduling/data_permutationflowshopschedulingtt.csv -l "${DATE}_permutationflowshopschedulingtt" -b heuristiclong -t 2200
```

![permutationflowshopschedulingtt](img/permutationflowshopschedulingtt.png?raw=true "permutationflowshopschedulingtt")

</p>
</details>

[Single machine order acceptance and scheduling problem with time windows and sequence-dependent setup times, Total weighted tardiness](examples/orderacceptanceandscheduling.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/orderacceptanceandscheduling_main --csv ../ordata/orderacceptanceandscheduling/data.csv -l orderacceptanceandscheduling -t 60
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/orderacceptanceandscheduling/data.csv -l orderacceptanceandscheduling -b heuristiclong -t 61
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
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/multidimensionalmultiplechoiceknapsack_main --csv ../ordata/multidimensionalmultiplechoiceknapsack/data.csv -l "${DATE}_multidimensionalmultiplechoiceknapsack" -t 60
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/multidimensionalmultiplechoiceknapsack/data.csv -l "${DATE}_multidimensionalmultiplechoiceknapsack" -b heuristiclong -t 62
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
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/quadraticassignment_main --csv ../ordata/quadraticassignment/data.csv -l "${DATE}_quadraticassignment" -t 60
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/quadraticassignment/data.csv -l "${DATE}_quadraticassignment" -b heuristiclong -t 62
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
DATE=$(date '+%Y-%m-%d--%H-%M') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/knapsackwithconflicts_main --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'hifi2006'" -l "${DATE}_knapsackwithconflicts_hifi2006" -t 60
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'hifi2006'" -l "${DATE}_knapsackwithconflicts_hifi2006" -b heuristiclong -t 62
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
python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/travellingsalesman_main --csv ../ordata/travellingsalesman/data.csv -l travellingsalesman -t 60
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/travellingsalesman/data.csv -l travellingsalesman -b heuristiclong -t 62
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
python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/permutationflowshopschedulingmakespan_main --csv ../ordata/permutationflowshopscheduling/data_makespan.csv -f "row['Dataset'] == 'vallada2015_large'" -l permutationflowshopschedulingmakespan --timelimitfield "Time limit"
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/permutationflowshopscheduling/data_makespan.csv -f "row['Dataset'] == 'vallada2015_large' and int(row['Job number']) <= 100" -l permutationflowshopschedulingmakespan -b heuristiclong -t 500
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
python3 ../optimizationtools/scripts/bench_run.py --main "./bazel-bin/examples/roadef2020_main -w 0 -y 1 " --csv ../ordata/roadef2020/data/data.csv -l roadef2020 -t 900 -f "'A' not in row['Dataset']"
python3 ../optimizationtools/scripts/bench_process.py -b heuristiclong --csv ../ordata/roadef2020/data.csv -l roadef2020 -t 920 -f "'A' not in row['Dataset']"
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
./bazel-bin/examples/knapsackwithconflicts_main -v -i ../ordata/knapsackwithconflicts/bettinelli2017/sparse_corr/test_1000_1000_r0.001-0.dat -f bettinelli2017 -t 5 -c sol.txt --print-solution
```
```
=======================================
          Local Search Solver          
=======================================

Algorithm
---------
Best First Local Search

Parameters
----------
Maximum number of nodes:     -1
Seed:                        0
Maximum size of the pool:    1
Time limit:                  5

Local scheme parameters
-----------------------
Swap:                        1
(2,1)-swap:                  1

      Time                                   Value                                 Comment
      ----                                   -----                                 -------
     0.001                                -1310, 0                                 s0 (t0)
     0.002                                -1320, 0                           n1 d1 c1 (t0)
     0.005                                -1330, 0                           n4 d2 c3 (t0)
     0.006                                -1340, 0                           n5 d3 c3 (t0)
     0.007                                -1350, 0                           n7 d4 c4 (t0)
     0.008                                -1360, 0                           n8 d5 c5 (t0)
     0.009                                -1370, 0                           n9 d6 c2 (t0)
     0.011                                -1380, 0                          n11 d7 c8 (t0)
     0.013                                -1390, 0                         n15 d8 c12 (t0)
     0.014                                -1410, 0                         n16 d9 c11 (t0)
     0.015                                -1420, 0                        n17 d10 c12 (t0)
     0.016                                -1430, 0                        n18 d11 c13 (t0)
     0.017                                -1440, 0                        n19 d12 c13 (t0)
     0.019                                -1450, 0                        n23 d13 c17 (t0)
     0.021                                -1460, 0                        n24 d14 c18 (t0)
     0.022                                -1470, 0                        n25 d15 c18 (t0)
     0.028                                -1480, 0                        n31 d16 c26 (t0)
     0.029                                -1500, 0                        n32 d17 c26 (t0)
     0.029                                -1510, 0                        n33 d18 c26 (t0)
     0.030                                -1530, 0                        n35 d19 c28 (t0)
     0.032                                -1540, 0                        n36 d20 c27 (t0)
     0.033                                -1550, 0                        n37 d21 c29 (t0)
     0.034                                -1560, 0                        n38 d22 c26 (t0)
     0.034                                -1570, 0                        n39 d23 c33 (t0)
     0.053                                -1580, 0                        n54 d24 c48 (t0)
     0.059                                -1600, 0                        n60 d25 c53 (t0)
     0.068                                -1610, 0                        n68 d26 c62 (t0)
     0.099                                -1620, 0                        n96 d27 c89 (t0)
     0.113                                -1640, 0                       n110 d28 c-2 (t0)
     0.116                                -1650, 0                       n113 d29 c12 (t0)
     0.120                                -1660, 0                       n116 d30 c-2 (t0)
     0.139                                -1670, 0                       n130 d31 c26 (t0)
     0.140                                -1680, 0                       n131 d32 c28 (t0)
     0.152                                -1690, 0                       n143 d33 c-2 (t0)
     0.161                                -1700, 0                       n152 d34 c-2 (t0)
     0.165                                -1710, 0                       n155 d35 c51 (t0)
     0.169                                -1720, 0                       n159 d36 c54 (t0)
     0.220                                -1730, 0                       n191 d37 c84 (t0)
     0.232                                -1750, 0                       n202 d38 c-2 (t0)
     0.235                                -1760, 0                       n204 d39 c93 (t0)
     0.260                                -1770, 0                       n219 d40 c-2 (t0)
     0.268                                -1780, 0                       n224 d41 c-2 (t0)
     0.291                                -1790, 0                       n238 d42 c-2 (t0)
     0.492                                -1800, 0                       n349 d44 c-2 (t0)
     0.550                                -1810, 0                       n380 d45 c-2 (t0)
     0.554                                -1820, 0                       n384 d46 c78 (t0)
     0.560                                -1830, 0                       n389 d47 c83 (t0)
     0.586                                -1840, 0                       n408 d48 c-2 (t0)
     0.593                                -1850, 0                       n414 d49 c14 (t0)
     0.677                                -1860, 0                       n464 d50 c-2 (t0)
     0.688                                -1870, 0                       n470 d51 c67 (t0)
     0.691                                -1880, 0                       n472 d52 c66 (t0)
     0.696                                -1890, 0                       n477 d53 c-2 (t0)
     0.702                                -1900, 0                       n482 d54 c-2 (t0)
     0.707                                -1910, 0                       n486 d55 c-2 (t0)
     0.742                                -1920, 0                       n512 d56 c-2 (t0)
     0.761                                -1930, 0                       n526 d57 c-2 (t0)
     0.764                                -1940, 0                       n530 d58 c91 (t0)
     0.765                                -1950, 0                       n531 d59 c28 (t0)
     0.777                                -1960, 0                       n540 d60 c37 (t0)
     0.778                                -1970, 0                       n541 d61 c37 (t0)
     0.862                                -1980, 0                       n591 d62 c-2 (t0)
     0.882                                -1990, 0                       n601 d63 c96 (t0)
     1.225                                -2000, 0                       n756 d65 c-2 (t0)
     1.237                                -2010, 0                       n763 d66 c-2 (t0)
     1.310                                -2020, 0                       n797 d67 c91 (t0)
     1.516                                -2030, 0                       n883 d68 c-2 (t0)
     1.629                                -2040, 0                       n936 d69 c-2 (t0)
     1.679                                -2050, 0                       n963 d70 c-2 (t0)
     1.689                                -2060, 0                       n968 d71 c66 (t0)
     1.710                                -2070, 0                       n977 d72 c-2 (t0)
     1.712                                -2080, 0                       n978 d73 c69 (t0)
     1.719                                -2090, 0                       n982 d74 c-2 (t0)
     1.726                                -2100, 0                       n984 d75 c74 (t0)
     1.731                                -2110, 0                       n988 d76 c75 (t0)
     1.741                                -2120, 0                       n996 d77 c85 (t0)
     1.786                                -2130, 0                     n1017 d78 c115 (t0)
     1.787                                -2140, 0                     n1018 d79 c114 (t0)
     1.802                                -2150, 0                     n1026 d80 c124 (t0)
     1.818                                -2160, 0                     n1035 d81 c133 (t0)
     1.881                                -2170, 0                     n1065 d82 c165 (t0)
     2.036                                -2180, 0                     n1149 d83 c247 (t0)
     2.041                                -2190, 0                     n1152 d84 c248 (t0)
     2.042                                -2200, 0                     n1153 d85 c246 (t0)
     2.072                                -2210, 0                     n1171 d86 c268 (t0)
     2.079                                -2220, 0                     n1175 d87 c265 (t0)
     2.082                                -2230, 0                     n1177 d88 c271 (t0)
     2.099                                -2240, 0                     n1187 d89 c280 (t0)
     2.100                                -2250, 0                     n1188 d90 c279 (t0)
     2.108                                -2260, 0                     n1192 d91 c284 (t0)
     2.391                                -2270, 0                     n1328 d92 c424 (t0)
     2.403                                -2280, 0                     n1336 d93 c429 (t0)
     2.423                                -2290, 0                     n1348 d94 c444 (t0)
     2.454                                -2300, 0                     n1367 d95 c462 (t0)
     2.456                                -2310, 0                     n1368 d96 c457 (t0)
     3.512                                -2320, 0                      n1811 d97 c-2 (t0)
     3.584                                -2330, 0                     n1849 d98 c445 (t0)

Final statistics
----------------
Value:                      -2330, 0
Time:                       5.00216

Local scheme statistics
-----------------------
Toggle:                     0 / 0 / -nan%
Swap:                       0 / 0 / -nan%
(2-1)-swap:                 0 / 0 / -nan%

Solution
--------
Items: 0 1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 57 58 59 60 62 64 65 66 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 110 111 112 113 114 115 117 118 119 121 122 123 124 125 126 127 128 129 132 133 134 136 137 138 140 141 142 143 157 518
Weight: 1000 / 1000
Profit: 2330

Checker
-------
Number of Items:                133 / 1000
Number of duplicates:           0
Number of conflict violations:  0
Weight:                         1000 / 1000
Feasible:                       1
Profit:                         2330
```

## Usage, C++ library

See examples.

