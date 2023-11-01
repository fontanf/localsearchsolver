# LocalSearchSolver

A solver based on local search.

![localsearch](img/localsearch.jpg?raw=true "localsearch")

[image source](https://commons.wikimedia.org/wiki/File:Chocolate_Hills_overview.JPG)

## Description

The goal of this repository is to provide a simple framework to quickly implement heuristic algorithms based on local search.

For complex and evolving problems, implementing and extending local search algorithms based on complex neighborhoods quickly becomes cumbersome and time consuming. This often makes them unsuitable for practical use.
The algorithms of this repository are designed to get the best performances out of the simplest neighborhoods. Thus making them easier to implement and extend.

Still, more complex neighborhoods can be implemented if better performances are needed.

All algorithms require the Local Scheme to implement: `GlobalCost`, `Solution`, `initial_solution`, `global_cost` and `local_search`.

Implemented algorithms:
* Multi-start local search `-a restarting-local-search`
  * Most basic algorithm, can be used as a baseline
* Iterated local search `-a iterated-local-search`
  * Single solution algorithm
  * Low implementation cost
  * Good for very large problems
  * Additional requirements: `Perturbation`, `perturbations`, `apply_perturbation`
* Best first local search `-a best-first-local-search`
  * Multi-solution algorithm
  * Higher implementation cost
  * Higher quality solutions
  * Additional requirements: `CompactSolution`, `compact_solution_hasher`, `solution2compact`, `compact2solution`, `Perturbation`, `perturbation_hasher`, `perturbations`, `apply_perturbation`
* Genetic local search `-a genetic-local-search`
  * Population-based algorithm
  * Good when the ruggedness of the landscape is high
  * Additional requirements: `crossover` and `distance`

## Sequencing module

A specific implementation is also available for sequencing problems. The neighborhoods, perturbations and crossovers are already implemented. It is only required to provide a `void append(SequenceData&, ElementId) const` method to use them.

Following the framework described in "A unified solution framework for multi-attribute vehicle routing problems" (Vidal et al., 2014) [DOI](https://doi.org/10.1016/j.ejor.2013.09.045) , more efficient neighborhood explorations can be obtained by providing either a `void concatenate(SequenceData&, const SequenceData&) const` method or a `GlobalCost global_cost_concatenate(SequenceData&, const SequenceData&) const` method (when it is possible).

In case the `GlobalCost global_cost_concatenate(SequenceData&, const SequenceData&) const` method is provided, the `Swap` and `Reverse` neighborhoods won't take much avantage of it. Therefore, it might be better to disable these neighborhoods, in particular when optimizing a single sequence.

* Local search neighborhoods:
  * Intra:
    * Shift a block of `k` consecutive elements
    * Shift and reverse a block of `k` consecutive elements
    * Swap a block of `k1` consecutive elements with another block of `k2` consecutive elements
    * Reverse a block of consecutive elements (2-opt)
  * Inter:
    * Swap the tails of two sequences (2-opt\*)
    * Replace the tail of sequence `i1` with the reversed head of sequence `i2` and replace the head of sequence `i2` with the reversed tail of sequence `i1` (reversed 2-opt\*)
    * Shift a block of `k` consecutive elements from one sequence to another
    * Shift and reverse a block of `k` consecutive elements from one sequence to another
    * Swap a block of `k1` consecutive elements from a sequence with another block of `k2` consecutive elements from another sequence
    * Swap the sequences of two elements `j1` and `j2` and try to place them at promising positions in their new sequences (swap\*)
  * Sub-sequences:
    * Add an element into the solution
    * Remove an element from the solution
    * Replace an element from the solution by an element outside of the solution
  * Modes:
    * Shift a block of `k` consecutive elements and change the mode of the first one.
    * Swap the mode of element `j1` with the mode of another element `j2` from the same sequence
    * Swap an elemnt `j1` with another element `j2` from the same sequence and swap their modes
    * Increment the mode of element `j1` and decrement the mode of element `j2` from the same sequence
* Perturbations:
  * Swap two blocks of consecutive elements (double-bridge) (single sequence only)
  * Remove `k` elements and re-insert them (ruin-and-recreate)
  * Force an element into the solution
* Crossover algorithms:
  * Order crossover (OX): single or multiple sequences
  * Similar job order crossover (SJOX): single sequence
  * Similar block order crossover (SBOX): single sequence
  * Maximal preservative crossover (MPX): single sequence
  * Cycle crossover (CX): single sequence
  * Selective route exchange crossover 1 (SREX1): multiple sequences
  * Selective route exchange crossover 2 (SREX2): multiple sequences

### Examples

#### Single sequence

[Sequential ordering problem](examples/sequentialordering.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M-%S') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/sequentialordering_main --csv ../ordata/sequentialordering/data.csv -f "row['Dataset'] == 'tsplib'" -l "${DATE}_sequentialordering" -t 60
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/sequentialordering/data.csv -f "row['Dataset'] == 'tsplib'" -l "${DATE}_sequentialordering" -b heuristiclong -t 62
```

![sequentialordering](img/sequentialordering_tsplib.png?raw=true "sequentialordering_tsplib")

</p>
</details>

[Single machine scheduling problem with sequence-dependent setup times, total weighted tardiness](examples/schedulingwithsdsttwt.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M-%S') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/schedulingwithsdsttwt_main --csv ../ordata/schedulingwithsdsttwt/data.csv -l "${DATE}_schedulingwithsdsttwt" -t 60
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/schedulingwithsdsttwt/data.csv -l "${DATE}_schedulingwithsdsttwt" -b heuristiclong -t 62
```

![schedulingwithsdsttwt](img/schedulingwithsdsttwt.png?raw=true "schedulingwithsdsttwt")

</p>
</details>

[Permutation flow shop scheduling problem, total completion time](examples/permutationflowshopschedulingtct.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M-%S') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/permutationflowshopschedulingtct_main --csv ../ordata/flowshopscheduling/data_permutationflowshopschedulingtct.csv -l "${DATE}_permutationflowshopschedulingtct" --timelimitfield "Time limit"
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/flowshopscheduling/data_permutationflowshopschedulingtct.csv -l "${DATE}_permutationflowshopschedulingtct" -b heuristiclong -t 500
```

![permutationflowshopschedulingtct](img/permutationflowshopschedulingtct.png?raw=true "permutationflowshopschedulingtct")

</p>
</details>

[Permutation flow shop scheduling problem, total tardiness](examples/permutationflowshopschedulingtt.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M-%S') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/permutationflowshopschedulingtt_main --csv ../ordata/flowshopscheduling/data_permutationflowshopschedulingtt.csv -l "${DATE}_permutationflowshopschedulingtt" --timelimitfield "Time limit"
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/flowshopscheduling/data_permutationflowshopschedulingtt.csv -l "${DATE}_permutationflowshopschedulingtt" -b heuristiclong -t 2200
```

![permutationflowshopschedulingtt](img/permutationflowshopschedulingtt.png?raw=true "permutationflowshopschedulingtt")

</p>
</details>

#### Single sub-sequence

[Single machine order acceptance and scheduling problem with time windows and sequence-dependent setup times, total weighted tardiness](examples/orderacceptanceandscheduling.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/orderacceptanceandscheduling_main --csv ../ordata/orderacceptanceandscheduling/data.csv -l orderacceptanceandscheduling -t 60
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/orderacceptanceandscheduling/data.csv -l orderacceptanceandscheduling -b heuristiclong -t 61
```

</p>
</details>

[Time-dependent orienteering problem](examples/timedependentorienteering.hpp)

#### Multiple sequences

[Distributed permutation flow shop scheduling problem, total completion time](examples/distributedpfsstct.hpp)

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M-%S') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/permutationflowshopschedulingtct_main --csv ../ordata/flowshopscheduling/data_distributedpfsstct.csv -l "${DATE}_distributedpfsstct" --timelimitfield "Time limit"
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/flowshopscheduling/data_distributedpfsstct.csv -l "${DATE}_distributedpfsstct" -b heuristiclong -t 500
```

</p>
</details>

#### Single sequence, modes

[Traveling salesman problem with release dates](examples/travelingsalesmanwithreleasedates.hpp)

[Single machine batch scheduling problem, total weighted tardiness](examples/batchschedulingtotalweightedtardiness.hpp)

#### Single sequence, `sequence_data_init/concatenate`

[Traveling repairman problem](examples/travelingrepairman.hpp)

#### Multiple sequences, `sequence_data_init/concatenate`

[Capacitated vehicle routing problem](examples/capacitatedvehiclerouting.hpp)

[Vehicle routing problem with time windows](examples/vehicleroutingwithtimewindows.hpp)


## Other examples

Data can be downloaded from [fontanf/orproblems](https://github.com/fontanf/orproblems)

[Multidimensional multiple-choice knapsack problem](examples/multidimensionalmultiplechoiceknapsack.hpp)
* Straightforward example for a genetic local search: single neighborhood, basic operators
* Algorithm:
  * Local search neighborhoods:
    * Add item `j` in the knapsack
  * Crossover algorithm

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M-%S') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/multidimensionalmultiplechoiceknapsack_main --csv ../ordata/multidimensionalmultiplechoiceknapsack/data.csv -l "${DATE}_multidimensionalmultiplechoiceknapsack" -t 60
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/multidimensionalmultiplechoiceknapsack/data.csv -l "${DATE}_multidimensionalmultiplechoiceknapsack" -b heuristiclong -t 62
```

![multidimensionalmultiplechoiceknapsack](img/multidimensionalmultiplechoiceknapsack.png?raw=true "multidimensionalmultiplechoiceknapsack")

</p>
</details>

[Quadratic assignment problem](examples/quadraticassignment.hpp)
* Example which implements a problem specific acceleration strategy to compute the move costs
* Algorithm:
  * Local search neighborhood: swap two assignments
  * Perturbation: ejection chain
  * Crossover algorithm: UX crossover

<details><summary>Benchmarks</summary>
<p>

```shell
DATE=$(date '+%Y-%m-%d--%H-%M-%S') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/quadraticassignment_main --csv ../ordata/quadraticassignment/data.csv -l "${DATE}_quadraticassignment" -t 60
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/quadraticassignment/data.csv -l "${DATE}_quadraticassignment" -b heuristiclong -t 62
```

![quadraticassignment](img/quadraticassignment.png?raw=true "quadraticassignment")

</p>
</details>

[Knapsack problem with conflicts](examples/knapsackwithconflicts.hpp)
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
DATE=$(date '+%Y-%m-%d--%H-%M-%S') && python3 ../optimizationtools/scripts/bench_run.py --main ./bazel-bin/examples/knapsackwithconflicts_main --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'hifi2006'" -l "${DATE}_knapsackwithconflicts_hifi2006" -t 60
python3 ../optimizationtools/scripts/bench_process.py --csv ../ordata/knapsackwithconflicts/data.csv -f "row['Dataset'] == 'hifi2006'" -l "${DATE}_knapsackwithconflicts_hifi2006" -b heuristiclong -t 62
```

![knapsackwithconflicts](img/knapsackwithconflicts_hifi2006.png?raw=true "knapsackwithconflicts_hifi2006")

</p>
</details>

[Generalized assignment problem](https://github.com/fontanf/generalizedassignmentsolver/blob/master/generalizedassignmentsolver/algorithms/local_search.cpp) from [fontanf/generalizedassignmentsolver](https://github.com/fontanf/generalizedassignmentsolver)
* Example which implements a generic strategy for very large problems to avoid recomputing moves which have not change since the last neighborhood exploration
* Algorithm:
  * Local search neighborhoods: shift item `j` to agent `i`
  * Perturbation: shift 8 random jobs

[Permutation flow shop scheduling problem, makespan](examples/permutationflowshopschedulingmakespan.hpp)
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

[ROADEF/EURO challenge 2020: maintenance planning problem](examples/roadef2020.hpp)
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

[Maximum-weight independent set problem](https://github.com/fontanf/stablesolver/blob/master/stablesolver/stable/algorithms/localsearch.cpp) and [maximum-weight clique problem](https://github.com/fontanf/stablesolver/blob/master/stablesolver/clique/algorithms/localsearch.cpp) from [fontanf/stablesolver](https://github.com/fontanf/stablesolver)
* Algorithm:
  * Local search neighborhoods:
    * Move vertex `v` in/out of the solution (and remove conflicting vertices)
    * Remove one vertex and add two non-conflicting vertices from its neighbors
  * Perturbation: force vertex `v` in/out of the solution

## Usage, running examples from command line

Compile:
```shell
bazel build -- //...
```

Then, examples can be executed as follows:
```shell
./bazel-bin/examples/knapsackwithconflicts_main -v 1 -i ../ordata/knapsackwithconflicts/bettinelli2017/sparse_corr/test_1000_1000_r0.001-0.dat -f bettinelli2017 -t 5 -c sol.txt --print-solution 1
```
```
Instance
--------
Number of items:         1000
Capacity:                1000
Number of conflicts:     499
Weight ratio:            50.703
Average # of conflicts:  0.499

=======================================
           LocalSearchSolver           
=======================================

Algorithm
---------
Best first local search

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
     0.012                             1310.000000                                 s0 (t0)
     0.018                             1320.000000                           n1 d1 c1 (t0)
     0.030                             1330.000000                           n4 d2 c3 (t0)
     0.033                             1340.000000                           n5 d3 c3 (t0)
     0.039                             1350.000000                           n7 d4 c4 (t0)
     0.042                             1360.000000                           n8 d5 c5 (t0)
     0.045                             1370.000000                           n9 d6 c2 (t0)
     0.048                             1380.000000                          n11 d7 c8 (t0)
     0.053                             1390.000000                         n15 d8 c12 (t0)
     0.053                             1410.000000                         n16 d9 c11 (t0)
     0.055                             1420.000000                        n17 d10 c12 (t0)
     0.056                             1430.000000                        n18 d11 c13 (t0)
     0.057                             1440.000000                        n19 d12 c13 (t0)
     0.061                             1450.000000                        n23 d13 c17 (t0)
     0.063                             1460.000000                        n24 d14 c18 (t0)
     0.064                             1470.000000                        n25 d15 c18 (t0)
     0.070                             1480.000000                        n31 d16 c26 (t0)
     0.071                             1500.000000                        n32 d17 c26 (t0)
     0.072                             1510.000000                        n33 d18 c26 (t0)
     0.073                             1530.000000                        n35 d19 c28 (t0)
     0.074                             1540.000000                        n36 d20 c27 (t0)
     0.076                             1550.000000                        n37 d21 c29 (t0)
     0.076                             1560.000000                        n38 d22 c26 (t0)
     0.077                             1570.000000                        n39 d23 c33 (t0)
     0.095                             1580.000000                        n54 d24 c48 (t0)
     0.101                             1600.000000                        n60 d25 c53 (t0)
     0.110                             1610.000000                        n68 d26 c62 (t0)
     0.140                             1620.000000                        n96 d27 c89 (t0)
     0.154                             1640.000000                       n110 d28 c-2 (t0)
     0.157                             1650.000000                       n113 d29 c12 (t0)
     0.161                             1660.000000                       n116 d30 c-2 (t0)
     0.179                             1670.000000                       n130 d31 c26 (t0)
     0.180                             1680.000000                       n131 d32 c28 (t0)
     0.192                             1690.000000                       n143 d33 c-2 (t0)
     0.202                             1700.000000                       n152 d34 c-2 (t0)
     0.205                             1710.000000                       n155 d35 c51 (t0)
     0.209                             1720.000000                       n159 d36 c54 (t0)
     0.259                             1730.000000                       n191 d37 c84 (t0)
     0.270                             1750.000000                       n202 d38 c-2 (t0)
     0.274                             1760.000000                       n204 d39 c93 (t0)
     0.298                             1770.000000                       n219 d40 c-2 (t0)
     0.306                             1780.000000                       n224 d41 c-2 (t0)
     0.329                             1790.000000                       n238 d42 c-2 (t0)
     0.529                             1800.000000                       n349 d44 c-2 (t0)
     0.586                             1810.000000                       n380 d45 c-2 (t0)
     0.591                             1820.000000                       n384 d46 c78 (t0)
     0.597                             1830.000000                       n389 d47 c83 (t0)
     0.622                             1840.000000                       n408 d48 c-2 (t0)
     0.629                             1850.000000                       n414 d49 c14 (t0)
     0.714                             1860.000000                       n464 d50 c-2 (t0)
     0.724                             1870.000000                       n470 d51 c67 (t0)
     0.727                             1880.000000                       n472 d52 c66 (t0)
     0.732                             1890.000000                       n477 d53 c-2 (t0)
     0.739                             1900.000000                       n482 d54 c-2 (t0)
     0.744                             1910.000000                       n486 d55 c-2 (t0)
     0.778                             1920.000000                       n512 d56 c-2 (t0)
     0.797                             1930.000000                       n526 d57 c-2 (t0)
     0.800                             1940.000000                       n530 d58 c91 (t0)
     0.801                             1950.000000                       n531 d59 c28 (t0)
     0.814                             1960.000000                       n540 d60 c37 (t0)
     0.815                             1970.000000                       n541 d61 c37 (t0)
     0.898                             1980.000000                       n591 d62 c-2 (t0)
     0.919                             1990.000000                       n601 d63 c96 (t0)
     1.262                             2000.000000                       n756 d65 c-2 (t0)
     1.273                             2010.000000                       n763 d66 c-2 (t0)
     1.346                             2020.000000                       n797 d67 c91 (t0)
     1.552                             2030.000000                       n883 d68 c-2 (t0)
     1.666                             2040.000000                       n936 d69 c-2 (t0)
     1.716                             2050.000000                       n963 d70 c-2 (t0)
     1.726                             2060.000000                       n968 d71 c66 (t0)
     1.748                             2070.000000                       n977 d72 c-2 (t0)
     1.750                             2080.000000                       n978 d73 c69 (t0)
     1.757                             2090.000000                       n982 d74 c-2 (t0)
     1.763                             2100.000000                       n984 d75 c74 (t0)
     1.769                             2110.000000                       n988 d76 c75 (t0)
     1.779                             2120.000000                       n996 d77 c85 (t0)
     1.824                             2130.000000                     n1017 d78 c115 (t0)
     1.825                             2140.000000                     n1018 d79 c114 (t0)
     1.840                             2150.000000                     n1026 d80 c124 (t0)
     1.857                             2160.000000                     n1035 d81 c133 (t0)
     1.918                             2170.000000                     n1065 d82 c165 (t0)
     2.078                             2180.000000                     n1149 d83 c247 (t0)
     2.082                             2190.000000                     n1152 d84 c248 (t0)
     2.084                             2200.000000                     n1153 d85 c246 (t0)
     2.113                             2210.000000                     n1171 d86 c268 (t0)
     2.121                             2220.000000                     n1175 d87 c265 (t0)
     2.123                             2230.000000                     n1177 d88 c271 (t0)
     2.140                             2240.000000                     n1187 d89 c280 (t0)
     2.142                             2250.000000                     n1188 d90 c279 (t0)
     2.149                             2260.000000                     n1192 d91 c284 (t0)
     2.435                             2270.000000                     n1328 d92 c424 (t0)
     2.447                             2280.000000                     n1336 d93 c429 (t0)
     2.468                             2290.000000                     n1348 d94 c444 (t0)
     2.498                             2300.000000                     n1367 d95 c462 (t0)
     2.501                             2310.000000                     n1368 d96 c457 (t0)
     3.567                             2320.000000                      n1811 d97 c-2 (t0)
     3.640                             2330.000000                     n1849 d98 c445 (t0)

Final statistics
----------------
Value:                      2330.000000
Time:                       5.00283

Local scheme statistics
-----------------------
Toggle:            14361 / 11902 / 82.8772%
Swap:              14152 / 11721 / 82.8222%
(2-1)-swap:        10440 / 13 / 0.124521%

Solution
--------
Profit:            2330
Weight:            1000 / 1000
Number of items:   133 / 1000

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

