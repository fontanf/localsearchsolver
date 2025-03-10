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

[Sequential ordering problem](examples/sequencing/sequential_ordering_main.cpp)

[Single machine scheduling problem with sequence-dependent setup times, total weighted tardiness](examples/sequencing/scheduling_with_sdst_twt_main.cpp)

[Permutation flow shop scheduling problem, total completion time](examples/sequencing/permutation_flowshop_scheduling_tct_main.cpp)

[Permutation flow shop scheduling problem, total tardiness](examples/sequencing/permutation_flowshop_scheduling_tt_main.cpp)

#### Single sub-sequence

[Single machine order acceptance and scheduling problem with time windows and sequence-dependent setup times, total weighted tardiness](examples/sequencing/order_acceptance_and_scheduling_main.cpp)

[Time-dependent orienteering problem](examples/sequencing/time_dependent_orienteering_main.cpp)

#### Multiple sequences

[Distributed permutation flow shop scheduling problem, total completion time](examples/sequencing/distributed_pfss_tct_main.cpp)

#### Single sequence, modes

[Traveling salesman problem with release dates](examples/sequencing/traveling_salesman_with_release_dates_main.cpp)

[Single machine batch scheduling problem, total weighted tardiness](examples/sequencing/batch_scheduling_total_weighted_tardiness_main.cpp)

#### Single sequence, `sequence_data_init/concatenate`

[Traveling repairman problem](examples/sequencing/traveling_repairman_main.cpp)

#### Multiple sequences, `sequence_data_init/concatenate`

[Capacitated vehicle routing problem](examples/sequencing/capacitated_vehicle_routing_main.cpp)

[Vehicle routing problem with time windows](examples/sequencing/vehicle_routing_with_time_windows_main.cpp)

#### Multiple sub-sequences, `sequence_data_init/concatenate`

[Team orienteering problem](examples/sequencing/team_orienteering_main.cpp)


## Other examples

Data can be downloaded from [fontanf/orproblems](https://github.com/fontanf/orproblems)

[Multidimensional multiple-choice knapsack problem](examples/multidimensional_multiple_choice_knapsack_main.cpp)
* Straightforward example for a genetic local search: single neighborhood, basic operators
* Algorithm:
  * Local search neighborhoods:
    * Add item `j` in the knapsack
  * Crossover algorithm

[Quadratic assignment problem](examples/quadratic_assignment_main.cpp)
* Example which implements a problem specific acceleration strategy to compute the move costs
* Algorithm:
  * Local search neighborhood: swap two assignments
  * Perturbation: ejection chain
  * Crossover algorithm: UX crossover

[Knapsack problem with conflicts](examples/knapsack_with_conflicts_main.cpp)
* Example with multiple neigborhoods
* Algorithm:
  * Local search neighborhoods:
    * Move item `j` in/out of the solution (and remove conflicting items)
    * Swap two non-conflicting items
    * Remove one item and add two non-conflicting items from its neighbors
  * Perturbation: force item `j` in/out of the solution
  * Crossover algorithm: double backbone-based crossover

[Generalized assignment problem](https://github.com/fontanf/generalizedassignmentsolver/blob/master/generalizedassignmentsolver/algorithms/local_search.cpp) from [fontanf/generalizedassignmentsolver](https://github.com/fontanf/generalizedassignmentsolver)
* Example which implements a generic strategy for very large problems to avoid recomputing moves which have not change since the last neighborhood exploration
* Algorithm:
  * Local search neighborhoods: shift item `j` to agent `i`
  * Perturbation: shift 8 random jobs

[Permutation flow shop scheduling problem, makespan](examples/permutation_flowshop_scheduling_makespan_main.cpp)
* This one is not considered as a sequencing problem since the dedicated acceleration strategy makes it possible to explore the neighborhoods more efficiently
* Algorithm:
  * Local search neighborhood: move a block of `k` consecutive jobs, `k = 1..4`
  * Perturbation: swap two blocks (double-bridge)

[Maximum-weight independent set problem](https://github.com/fontanf/stablesolver/blob/master/src/stable/algorithms/local_search.cpp) and [maximum-weight clique problem](https://github.com/fontanf/stablesolver/blob/master/src/clique/algorithms/local_search.cpp) from [fontanf/stablesolver](https://github.com/fontanf/stablesolver)
* Algorithm:
  * Local search neighborhoods:
    * Move vertex `v` in/out of the solution (and remove conflicting vertices)
    * Remove one vertex and add two non-conflicting vertices from its neighbors
  * Perturbation: force vertex `v` in/out of the solution

## Usage, running examples from command line

Compile:
```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
cmake --install build --config Release --prefix install
```

Download data:
```shell
python3 scripts/download_data.py
```

Then, examples can be executed as follows:
```shell
./install/bin/localsearchsolver_sequential_ordering  --verbosity-level 1  --input ./data/sequential_ordering/soplib/R.200.100.1.sop  --format soplib  --algorithm best-first-local-search --maximum-number-of-nodes 100  --certificate solution.txt
```
```
=======================================
           LocalSearchSolver           
=======================================

Problem
-------
Sequential ordering problem
Number of locations:  200

Local scheme parameters
-----------------------
Neighborhoods
    Shift
        Block maximum length:                  8
    Swap
        Block maximum length:                  2
    Reverse:                                   0
    Shift-reverse
        Block maximum length:                  3
    Add/Remove:                                0
    Replace:                                   0
Perturbations
    Double-bridge
        Number of perturbations:               0
    Ruin-and-recreate
        Number of perturbations:               10
        Number of elements removed:            4
        Ruin random weight:                    1
        Ruin nearest weight:                   0
        Ruin adjacent string removal weight:   0
            Maximum string cardinality:        10
            Split rate:                        0.5
            Beta:                              0.01
        Recreate random:                       0
        Recreate best:                         0
    Force-add:                                 0
Crossovers
    Order crossover:                           0
    Similar job order crossover:               0
    Similar block order crossover:             0
    Selective route exchange crossover 1:      0
    Selective route exchange crossover 2:      0
    Maximal preservative crossover:            0

Algorithm
---------
Best first local search

Parameters
----------
Time limit:               inf
Messages
    Verbosity level:      1
    Standard output:      1
    File path:            
    # streams:            0
Logger
    Has logger:           0
    Standard error:       0
    File path:            
Solution pool size:       1
Seed:                     0
Goal:                     
Maximum number of nodes:  100

      Time                                   Value                                 Comment
      ----                                   -----                                 -------
     4.486                                  290, 0                                 s0 (t0)
     7.116                                  272, 0                           n2 d1 c2 (t0)
    10.021                                  236, 0                           n5 d2 c2 (t0)
    31.258                                  224, 0                          n34 d6 c5 (t0)
    38.488                                  209, 0                          n44 d7 c9 (t0)
    39.214                                  206, 0                          n45 d8 c0 (t0)
    51.860                                  204, 0                         n64 d11 c5 (t0)
    53.334                                  202, 0                         n66 d12 c1 (t0)
    54.232                                  199, 0                         n67 d13 c0 (t0)
    70.541                                  194, 0                         n93 d13 c6 (t0)

Algorithm statistics
--------------------
Value:            194, 0
Time (s):         74.0243
Number of nodes:  100

Local scheme statistics
-----------------------
General
    Initial solution:           1 / 1.5102e-05s / 1.5102e-05s
    Crossover time:             0 / 0s / -nans
    Local search time:          101 / 74.0166s / 0.732838s
    Local search iterations:    2839 / 28.1089
    Temporary structures time:  0.0108234s
Neighborhoods
    1-shift                     634 / 290 / 45.7413% / 6.81606s
    2-shift                     652 / 302 / 46.319% / 6.83057s
    3-shift                     657 / 284 / 43.2268% / 6.93018s
    4-shift                     670 / 317 / 47.3134% / 6.9931s
    5-shift                     618 / 305 / 49.3528% / 6.38262s
    6-shift                     645 / 319 / 49.4574% / 6.69405s
    7-shift                     638 / 309 / 48.4326% / 6.50337s
    8-shift                     659 / 311 / 47.1927% / 6.78887s
    (1,1)-swap                  571 / 76 / 13.31% / 2.53698s
    (2,1)-swap                  601 / 108 / 17.97% / 5.02293s
    (2,2)-swap                  572 / 76 / 13.2867% / 2.43386s
    2-shift-reverse             599 / 96 / 16.0267% / 5.4184s
    3-shift-reverse             597 / 46 / 7.70519% / 4.66367s

Local scheme solution
---------------------
Sequence 0: 0 16 105 172 177 107 91 119 80 158 164 104 151 42 74 187 53 175 99 147 70 27 136 153 56 87 118 123 125 196 114 26 28 6 30 36 143 195 82 20 127 92 135 148 109 55 72 193 24 160 25 197 113 5 157 63 8 76 149 21 38 191 51 22 45 88 59 52 3 188 111 134 106 9 37 40 130 103 128 10 54 81 186 146 166 145 180 131 120 95 62 61 184 100 190 154 44 183 169 178 129 17 78 165 18 116 139 23 173 126 85 75 50 66 67 58 174 69 39 60 159 64 98 93 141 79 140 110 144 14 65 189 133 115 101 4 156 43 117 34 102 12 83 108 29 94 132 170 89 121 46 182 19 168 73 97 171 162 194 122 90 33 84 152 96 77 86 198 161 112 11 68 142 13 35 179 163 32 181 176 150 192 137 124 71 2 49 47 57 41 15 185 48 31 155 7 167 138 1 199
    Cost: 194, 0
Total cost: 194, 0

Checker
-------
Number of Vertices:               200 / 200
Number of duplicates:             0
Number of precedence violations:  0
Feasible:                         1
Total distance:                   194
```

```shell
./install/bin/localsearchsolver_knapsack_with_conflicts  --verbosity-level 1  --input ./data/knapsack_with_conflicts/bettinelli2017/sparse_corr/test_1000_1000_r0.001-0.dat --format bettinelli2017  --time-limit 5  --certificate solution.txt
```
```
=======================================
           LocalSearchSolver           
=======================================

Problem
-------
Knapsack problem with conflicts
Number of items:         1000
Capacity:                1000
Number of conflicts:     499
Weight ratio:            50.703
Average # of conflicts:  0.499

Algorithm
---------
Best first local search

Parameters
----------
Time limit:               5
Messages
    Verbosity level:      1
    Standard output:      1
    File path:            
    # streams:            0
Logger
    Has logger:           0
    Standard error:       0
    File path:            
Solution pool size:       1
Seed:                     0
Goal:                     
Maximum number of nodes:  -1

      Time                                   Value                                 Comment
      ----                                   -----                                 -------
     0.001                             1310.000000                                 s0 (t0)
     0.002                             1320.000000                           n1 d1 c1 (t0)
     0.010                             1330.000000                           n4 d2 c3 (t0)
     0.010                             1340.000000                           n5 d3 c3 (t0)
     0.012                             1350.000000                           n7 d4 c4 (t0)
     0.013                             1360.000000                           n8 d5 c5 (t0)
     0.014                             1370.000000                           n9 d6 c2 (t0)
     0.016                             1380.000000                          n11 d7 c8 (t0)
     0.018                             1390.000000                         n15 d8 c12 (t0)
     0.018                             1410.000000                         n16 d9 c11 (t0)
     0.020                             1420.000000                        n17 d10 c12 (t0)
     0.020                             1430.000000                        n18 d11 c13 (t0)
     0.021                             1440.000000                        n19 d12 c13 (t0)
     0.024                             1450.000000                        n23 d13 c17 (t0)
     0.025                             1460.000000                        n24 d14 c18 (t0)
     0.026                             1470.000000                        n25 d15 c18 (t0)
     0.032                             1480.000000                        n31 d16 c26 (t0)
     0.033                             1500.000000                        n32 d17 c26 (t0)
     0.034                             1510.000000                        n33 d18 c26 (t0)
     0.035                             1530.000000                        n35 d19 c28 (t0)
     0.036                             1540.000000                        n36 d20 c27 (t0)
     0.038                             1550.000000                        n37 d21 c29 (t0)
     0.038                             1560.000000                        n38 d22 c26 (t0)
     0.039                             1570.000000                        n39 d23 c33 (t0)
     0.057                             1580.000000                        n54 d24 c48 (t0)
     0.063                             1600.000000                        n60 d25 c53 (t0)
     0.071                             1610.000000                        n68 d26 c62 (t0)
     0.101                             1620.000000                        n96 d27 c89 (t0)
     0.115                             1640.000000                       n110 d28 c-2 (t0)
     0.118                             1650.000000                       n113 d29 c12 (t0)
     0.122                             1660.000000                       n116 d30 c-2 (t0)
     0.140                             1670.000000                       n130 d31 c26 (t0)
     0.141                             1680.000000                       n131 d32 c28 (t0)
     0.153                             1690.000000                       n143 d33 c-2 (t0)
     0.162                             1700.000000                       n152 d34 c-2 (t0)
     0.166                             1710.000000                       n155 d35 c51 (t0)
     0.170                             1720.000000                       n159 d36 c54 (t0)
     0.219                             1730.000000                       n191 d37 c84 (t0)
     0.230                             1750.000000                       n202 d38 c-2 (t0)
     0.233                             1760.000000                       n204 d39 c93 (t0)
     0.258                             1770.000000                       n219 d40 c-2 (t0)
     0.266                             1780.000000                       n224 d41 c-2 (t0)
     0.289                             1790.000000                       n238 d42 c-2 (t0)
     0.485                             1800.000000                       n349 d44 c-2 (t0)
     0.542                             1810.000000                       n380 d45 c-2 (t0)
     0.547                             1820.000000                       n384 d46 c78 (t0)
     0.553                             1830.000000                       n389 d47 c83 (t0)
     0.578                             1840.000000                       n408 d48 c-2 (t0)
     0.585                             1850.000000                       n414 d49 c14 (t0)
     0.668                             1860.000000                       n464 d50 c-2 (t0)
     0.678                             1870.000000                       n470 d51 c67 (t0)
     0.681                             1880.000000                       n472 d52 c66 (t0)
     0.686                             1890.000000                       n477 d53 c-2 (t0)
     0.692                             1900.000000                       n482 d54 c-2 (t0)
     0.697                             1910.000000                       n486 d55 c-2 (t0)
     0.731                             1920.000000                       n512 d56 c-2 (t0)
     0.749                             1930.000000                       n526 d57 c-2 (t0)
     0.752                             1940.000000                       n530 d58 c91 (t0)
     0.754                             1950.000000                       n531 d59 c28 (t0)
     0.766                             1960.000000                       n540 d60 c37 (t0)
     0.767                             1970.000000                       n541 d61 c37 (t0)
     0.849                             1980.000000                       n591 d62 c-2 (t0)
     0.869                             1990.000000                       n601 d63 c96 (t0)
     1.204                             2000.000000                       n756 d65 c-2 (t0)
     1.215                             2010.000000                       n763 d66 c-2 (t0)
     1.287                             2020.000000                       n797 d67 c91 (t0)
     1.488                             2030.000000                       n883 d68 c-2 (t0)
     1.599                             2040.000000                       n936 d69 c-2 (t0)
     1.649                             2050.000000                       n963 d70 c-2 (t0)
     1.658                             2060.000000                       n968 d71 c66 (t0)
     1.680                             2070.000000                       n977 d72 c-2 (t0)
     1.682                             2080.000000                       n978 d73 c69 (t0)
     1.689                             2090.000000                       n982 d74 c-2 (t0)
     1.695                             2100.000000                       n984 d75 c74 (t0)
     1.701                             2110.000000                       n988 d76 c75 (t0)
     1.711                             2120.000000                       n996 d77 c85 (t0)
     1.755                             2130.000000                     n1017 d78 c115 (t0)
     1.756                             2140.000000                     n1018 d79 c114 (t0)
     1.771                             2150.000000                     n1026 d80 c124 (t0)
     1.787                             2160.000000                     n1035 d81 c133 (t0)
     1.848                             2170.000000                     n1065 d82 c165 (t0)
     2.001                             2180.000000                     n1149 d83 c247 (t0)
     2.006                             2190.000000                     n1152 d84 c248 (t0)
     2.007                             2200.000000                     n1153 d85 c246 (t0)
     2.036                             2210.000000                     n1171 d86 c268 (t0)
     2.044                             2220.000000                     n1175 d87 c265 (t0)
     2.046                             2230.000000                     n1177 d88 c271 (t0)
     2.063                             2240.000000                     n1187 d89 c280 (t0)
     2.064                             2250.000000                     n1188 d90 c279 (t0)
     2.072                             2260.000000                     n1192 d91 c284 (t0)
     2.349                             2270.000000                     n1328 d92 c424 (t0)
     2.361                             2280.000000                     n1336 d93 c429 (t0)
     2.381                             2290.000000                     n1348 d94 c444 (t0)
     2.411                             2300.000000                     n1367 d95 c462 (t0)
     2.414                             2310.000000                     n1368 d96 c457 (t0)
     3.451                             2320.000000                      n1811 d97 c-2 (t0)
     3.521                             2330.000000                     n1849 d98 c445 (t0)

Algorithm statistics
--------------------
Value:            2330.000000
Time (s):         5.00157
Number of nodes:  2484

Local scheme statistics
-----------------------
Toggle:            14746 / 12232 / 82.9513%
Swap:              14556 / 12070 / 82.9211%
(2-1)-swap:        10738 / 13 / 0.121065%

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
