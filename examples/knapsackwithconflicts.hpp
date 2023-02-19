/**
 * Knapsack Problem with Conflicts.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/knapsackwithconflicts.hpp
 *
 * Local Search:
 *
 * Three neighborhoods are implemented:
 * - Toggle: remove an item from the knapsack or add an item into the knapsack
 *   and remove its neighbors from the knapsack
 *   Complexity: O(n)
 * - Swap: swap an item from the knapsack with an item out of the knapsack
 *   which is not one of its neighbors (this is already considered by the
 *   Toggle neighborhood).
 *   Complexity: O(number of non-conflicts)
 * - (2-1)-swap: remove an item from the knapsack and add two of its neighbors
 *   which are not neighbors.
 *   This neighborhood has originaly been proposed for the Maximum Independent
 *   Set Problem and for the Maximum-Weight Independent Set Problem.
 *   Complexity: O(number of conflicts)
 *   "Fast local search for the maximum independent set problem" (Andrade et al., 2012)
 *   https://doi.org/10.1007/s10732-012-9196-4
 *   "A hybrid iterated local search heuristic for the maximum weight independent set problem" (Nogueira et al., 2018)
 *   https://doi.org/10.1007/s11590-017-1128-7
 *
 */

#pragma once

#include "localsearchsolver/common.hpp"

#include "orproblems/knapsackwithconflicts.hpp"

namespace localsearchsolver
{

namespace knapsackwithconflicts
{

using namespace orproblems::knapsackwithconflicts;

class LocalScheme
{

public:

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        /** Enable Swap neighborhood. */
        bool swap = true;
        /** Enable (2-1)-swap neighborhood. */
        bool swap_2_1 = true;
        bool shuffle_neighborhood_order = true;
    };

    LocalScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters),
        tabu_(instance_.number_of_items()),
        items_(instance.number_of_items()),
        neighbors_(instance_.number_of_items()),
        free_items_(instance_.number_of_items()),
        free_items_2_(instance_.number_of_items())
    {
        // Initialize items_.
        std::iota(items_.begin(), items_.end(), 0);
    }

    /*
     * Global cost.
     */

    /** Global cost: <Overweight, Profit, Weight>; */
    using GlobalCost = std::tuple<Weight, Profit>;

    inline Weight&       overweight(GlobalCost& global_cost) const { return std::get<0>(global_cost); }
    inline Profit&           profit(GlobalCost& global_cost) const { return std::get<1>(global_cost); }
    inline Weight  overweight(const GlobalCost& global_cost) const { return std::get<0>(global_cost); }
    inline Profit      profit(const GlobalCost& global_cost) const { return std::get<1>(global_cost); }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {
            0,
            -value,
        };
    }

    inline std::string to_string(const GlobalCost& global_cost) const
    {
        if (overweight(global_cost) > 0)
            return "-inf";
        return std::to_string(-profit(global_cost));
    }

    /*
     * Solutions.
     */

    struct SolutionItem
    {
        /**
         * in == true iff item j is in the solution.
         */
        bool in = false;

        /**
         * neighbor_weight = w iff the sum of the weights of the neighbors of j
         * which are in the solution is equal to w.
         */
        Weight neighbor_weight = 0;

        /**
         * neighbor_profit = p iff the sum of the profits of the neighbors of j
         * which are in the solution is equal to p.
         */
        Profit neighbor_profit = 0;
    };

    struct Solution
    {
        std::vector<SolutionItem> items;
        Profit profit = 0;
        Weight weight = 0;
    };

    inline Solution empty_solution() const
    {
        Solution solution;
        solution.items.resize(instance_.number_of_items());
        return solution;
    }

    inline Solution initial_solution(
            Counter,
            std::mt19937_64& generator)
    {
        Solution solution = empty_solution();
        std::shuffle(items_.begin(), items_.end(), generator);
        for (ItemId j: items_)
            if (cost_add(solution, j) < global_cost(solution))
                add(solution, j);
        return solution;
    }

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            overweight(solution),
            -solution.profit,
        };
    }

    /*
     * Local search.
     */

    struct Perturbation;

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Perturbation& tabu = Perturbation())
    {
        //std::cout << "> local search" << std::endl;
        //print(std::cout, solution);
        //std::cout << "tabu.j " << tabu.j << " " << tabu.add << std::endl;

        // Get neighborhoods.
        std::vector<Counter> neighborhoods = {0};
        if (parameters_.swap)
            neighborhoods.push_back(1);
        if (parameters_.swap_2_1)
            neighborhoods.push_back(2);

        // We forbid changing the status of tabu.j. This means that, if tabu.j
        // has been forced into the solution, we also need to forbid to add a
        // neighbor of tabu.j.
        std::fill(tabu_.begin(), tabu_.end(), 0);
        if (tabu.j != -1) {
            tabu_[tabu.j] = 1;
            if (contains(solution, tabu.j))
                for (ItemId j_neighbor: instance_.item(tabu.j).neighbors)
                    tabu_[j_neighbor] = 1;
        }

        Counter it = 0;
        (void)it;
        for (;; ++it) {
            //std::cout << "it " << it
            //    << " c " << to_string(global_cost(solution))
            //    << std::endl;
            //print(std::cout, solution);

            if (parameters_.shuffle_neighborhood_order)
                std::shuffle(neighborhoods.begin(), neighborhoods.end(), generator);

            bool improved = false;
            // Loop through neighborhoods.
            for (Counter neighborhood: neighborhoods) {
                switch (neighborhood) {
                case 0: { // Toggle neighborhood.
                    std::shuffle(items_.begin(), items_.end(), generator);
                    ItemId j_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (ItemId j: items_) {
                        if (tabu_[j] == 1)
                            continue;
                        GlobalCost c = (contains(solution, j))?
                            cost_remove(solution, j):
                            cost_add(solution, j);
                        if (c >= c_best)
                            continue;
                        if (j_best != -1 && !dominates(c, c_best))
                            continue;
                        j_best = j;
                        c_best = c;
                    }
                    if (j_best != -1) {
                        improved = true;
                        // Apply perturbation.
                        if (contains(solution, j_best)) {
                            remove(solution, j_best);
                        } else {
                            add(solution, j_best);
                        }
                        toggle_number_of_sucesses_++;
                    }
                    toggle_number_of_explorations_++;
                    break;
                } case 1: { // Swap neighborhood.
                    std::shuffle(items_.begin(), items_.end(), generator);
                    items_in_.clear();
                    items_out_.clear();
                    for (ItemId j: items_) {
                        if (contains(solution, j)) {
                            items_in_.push_back(j);
                        } else {
                            items_out_.push_back(j);
                        }
                    }

                    ItemId j_in_best = -1;
                    ItemId j_out_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (ItemId j_in: items_in_) {
                        if (tabu_[j_in] == 1)
                            continue;
                        neighbors_.clear();
                        for (ItemId j_neighbor: instance_.item(j_in).neighbors)
                            neighbors_.add(j_neighbor);
                        for (ItemId j_out: items_out_) {
                            if (tabu_[j_out] == 1)
                                continue;
                            if (neighbors_.contains(j_out))
                                continue;
                            GlobalCost c = cost_swap(solution, j_in, j_out, c_best);
                            if (c >= c_best)
                                continue;
                            if (j_in_best != -1 && !dominates(c, c_best))
                                continue;
                            j_in_best = j_in;
                            j_out_best = j_out;
                            c_best = c;
                        }
                    }
                    if (j_in_best != -1) {
                        improved = true;
                        // Apply perturbation.
                        remove(solution, j_in_best);
                        add(solution, j_out_best);
                        swap_number_of_sucesses_++;
                    }
                    swap_number_of_explorations_++;
                    break;
                } case 2: { // (2-1)-swap neighborhood.
                    std::shuffle(items_.begin(), items_.end(), generator);
                    // Get items inside the knapsack.
                    items_in_.clear();
                    for (ItemId j: items_)
                        if (contains(solution, j))
                            items_in_.push_back(j);

                    ItemId j_in_best = -1;
                    ItemId j_out_1_best = -1;
                    ItemId j_out_2_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (ItemId j_in: items_in_) {
                        if (tabu_[j_in] == 1)
                            continue;
                        // Update free_items_
                        free_items_.clear();
                        for (ItemId j_neighbor: instance_.item(j_in).neighbors)
                            if (tabu_[j_neighbor] == 0
                                    && solution.items[j_neighbor].neighbor_profit
                                    == instance_.item(j_in).profit)
                                free_items_.add(j_neighbor);
                        if (free_items_.size() <= 2)
                            continue;
                        free_items_.shuffle_in(generator);
                        remove(solution, j_in);
                        for (ItemId j_out_1: free_items_) {
                            free_items_2_.clear();
                            for (ItemId j: free_items_)
                                free_items_2_.add(j);
                            free_items_2_.remove(j_out_1);
                            for (ItemId j_neighbor: instance_.item(j_out_1).neighbors)
                                if (free_items_2_.contains(j_neighbor))
                                    free_items_2_.remove(j_neighbor);
                            if (free_items_2_.empty())
                                continue;
                            free_items_2_.shuffle_in(generator);
                            add(solution, j_out_1);
                            for (ItemId j_out_2: free_items_2_) {
                                GlobalCost c = cost_add(solution, j_out_2);
                                if (c >= c_best)
                                    continue;
                                if (j_in_best != -1 && !dominates(c, c_best))
                                    continue;
                                j_in_best = j_in;
                                j_out_1_best = j_out_1;
                                j_out_2_best = j_out_2;
                                c_best = c;
                            }
                            remove(solution, j_out_1);
                        }
                        add(solution, j_in);
                    }
                    if (j_in_best != -1) {
                        improved = true;
                        // Apply perturbation.
                        remove(solution, j_in_best);
                        add(solution, j_out_1_best);
                        add(solution, j_out_2_best);
                        swap_2_1_number_of_sucesses_++;
                    }
                    swap_2_1_number_of_explorations_++;
                    break;
                }
                }
                if (improved)
                    break;
            }
            if (!improved)
                break;
        }
        //print(std::cout, solution);
    }

    /*
     * Genetic local search.
     */

    inline Solution crossover(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        Solution solution = empty_solution();
        std::vector<ItemId> items;
        for (ItemId j = 0; j < instance_.number_of_items(); ++j) {
            if (contains(solution_parent_1, j)
                    && contains(solution_parent_2, j)) {
                // Add items which are in both parents.
                add(solution, j);
            } else if (contains(solution_parent_1, j)
                    || contains(solution_parent_2, j)) {
                // Store items which are in one parent.
                items.push_back(j);
            }
        }
        // Add some of the items which are in one parent.
        std::shuffle(items.begin(), items.end(), generator);
        for (ItemId j: items)
            if (cost_add(solution, j) < global_cost(solution))
                add(solution, j);
        return solution;
    }

    inline ItemId distance(
            const Solution& solution_1,
            const Solution& solution_2) const
    {
        ItemId d = 0;
        for (ItemId j = 0; j < instance_.number_of_items(); ++j)
            if (contains(solution_1, j) != contains(solution_2, j))
                d++;
        return d;
    }

    /*
     * Iterated local search.
     */

    struct Perturbation
    {
        Perturbation(): j(-1), global_cost(worst<GlobalCost>()) { }

        ItemId j;
        bool add;
        GlobalCost global_cost;
    };

    inline std::vector<Perturbation> perturbations(
            const Solution& solution,
            std::mt19937_64&)
    {
        std::vector<Perturbation> perturbations;
        for (ItemId j = 0; j < instance_.number_of_items(); ++j) {
            Perturbation perturbation;
            perturbation.j = j;
            if (contains(solution, j)) {
                perturbation.global_cost = global_cost(solution);
                perturbation.add = false;
            } else {
                perturbation.global_cost = cost_add(solution, j);
                overweight(perturbation.global_cost) = overweight(solution);
                perturbation.add = true;
            }
            perturbations.push_back(perturbation);
        }
        return perturbations;
    }

    inline void apply_perturbation(
            Solution& solution,
            const Perturbation& perturbation,
            std::mt19937_64&) const
    {
        if (perturbation.add) {
            add(solution, perturbation.j);
        } else {
            remove(solution, perturbation.j);
        }
    }

    /*
     * Best first local search.
     */

    using CompactSolution = std::vector<bool>;

    struct CompactSolutionHasher
    {
        std::hash<CompactSolution> hasher;

        inline bool operator()(
                const std::shared_ptr<CompactSolution>& compact_solution_1,
                const std::shared_ptr<CompactSolution>& compact_solution_2) const
        {
            return *compact_solution_1 == *compact_solution_2;
        }

        inline std::size_t operator()(
                const std::shared_ptr<CompactSolution>& compact_solution) const
        {
            return hasher(*compact_solution);
        }
    };

    inline CompactSolutionHasher compact_solution_hasher() const { return CompactSolutionHasher(); }

    CompactSolution solution2compact(const Solution& solution)
    {
        std::vector<bool> items(instance_.number_of_items(), false);
        for (ItemId j = 0; j < instance_.number_of_items(); ++j)
            if (solution.items[j].in)
                items[j] = true;
        return items;
    }

    Solution compact2solution(const CompactSolution& compact_solution)
    {
        auto solution = empty_solution();
        for (ItemId j = 0; j < instance_.number_of_items(); ++j)
            if (compact_solution[j])
                add(solution, j);
        return solution;
    }

    struct PerturbationHasher
    {
        std::hash<ItemId> hasher;
        std::hash<bool> hasher_2;

        inline bool hashable(const Perturbation&) const { return true; }

        inline bool operator()(
                const Perturbation& perturbation_1,
                const Perturbation& perturbation_2) const
        {
            return perturbation_1.j == perturbation_2.j && perturbation_1.add == perturbation_2.add;
        }

        inline std::size_t operator()(
                const Perturbation& perturbation) const
        {
            size_t hash = hasher(perturbation.j);
            optimizationtools::hash_combine(hash, hasher_2(perturbation.add));
            return hash;
        }
    };

    inline PerturbationHasher perturbation_hasher() const { return PerturbationHasher(); }

    /*
     * Outputs.
     */

    std::ostream& print(
            std::ostream &os,
            const Solution& solution,
            int verbosity_level)
    {
        if (verbosity_level >= 1) {
            ItemId number_of_items = 0;
            for (ItemId j = 0; j < instance_.number_of_items(); ++j)
                if (contains(solution, j))
                    number_of_items++;
            os << "Profit:            " << solution.profit << std::endl;
            os << "Weight:            " << solution.weight << " / " << instance_.capacity() << std::endl;
            os << "Number of items:   " << number_of_items << " / " << instance_.number_of_items() << std::endl;
        }
        if (verbosity_level >= 2) {
            os << std::endl
                << std::setw(12) << "Item"
                << std::setw(12) << "Profit"
                << std::setw(12) << "Weight"
                << std::setw(12) << "Efficiency"
                << std::setw(12) << "# conflicts"
                << std::endl
                << std::setw(12) << "----"
                << std::setw(12) << "------"
                << std::setw(12) << "------"
                << std::setw(12) << "----------"
                << std::setw(12) << "-----------"
                << std::endl;
            for (ItemId j = 0; j < instance_.number_of_items(); ++j) {
                if (!contains(solution, j))
                    continue;
                os << std::setw(12) << j
                    << std::setw(12) << instance_.item(j).profit
                    << std::setw(12) << instance_.item(j).weight
                    << std::setw(12) << (double)instance_.item(j).profit / instance_.item(j).weight
                    << std::setw(12) << instance_.item(j).neighbors.size()
                    << std::endl;
            }
        }
        return os;
    }

    inline void write(
            const Solution& solution,
            std::string certificate_path) const
    {
        if (certificate_path.empty())
            return;
        std::ofstream cert(certificate_path);
        if (!cert.good()) {
            throw std::runtime_error(
                    "Unable to open file \"" + certificate_path + "\".");
        }

        for (ItemId j = 0; j < instance_.number_of_items(); ++j)
            if (contains(solution, j))
                cert << j << " ";
    }

    void print_parameters(
            optimizationtools::Info& info) const
    {
        info.os() << ""
                << "Swap:                        " << parameters_.swap << std::endl
                << "(2,1)-swap:                  " << parameters_.swap_2_1 << std::endl;
    }

    void print_statistics(
            optimizationtools::Info& info) const
    {
        info.os()
            << "Toggle:            "
            << toggle_number_of_explorations_
            << " / " << toggle_number_of_sucesses_
            << " / " << (double)toggle_number_of_sucesses_ / toggle_number_of_explorations_ * 100 << "%"
            << std::endl;
        info.output->json["Algorithm"]["ToggleNumberOfExplorations"]
            = toggle_number_of_explorations_;
        info.output->json["Algorithm"]["ToggleNumberOfSuccesses"]
            = toggle_number_of_explorations_;
        if (parameters_.swap) {
            info.os()
                << "Swap:              "
                << swap_number_of_explorations_
                << " / " << swap_number_of_sucesses_
                << " / " << (double)swap_number_of_sucesses_ / swap_number_of_explorations_ * 100 << "%"
                << std::endl;
            info.output->json["Algorithm"]["SwapNumberOfExplorations"]
                = swap_number_of_explorations_;
            info.output->json["Algorithm"]["SwapNumberOfSuccesses"]
                = swap_number_of_explorations_;
        }
        if (parameters_.swap_2_1) {
            info.os()
                << "(2-1)-swap:        "
                << swap_2_1_number_of_explorations_
                << " / " << swap_2_1_number_of_sucesses_
                << " / " << (double)swap_2_1_number_of_sucesses_ / swap_2_1_number_of_explorations_ * 100 << "%"
                << std::endl;
            info.output->json["Algorithm"]["Swap2,1NumberOfExplorations"]
                = swap_2_1_number_of_explorations_;
            info.output->json["Algorithm"]["Swap2,1NumberOfSuccesses"]
                = swap_2_1_number_of_explorations_;
        }
    }

private:

    /*
     * Manipulate solutions.
     */

    inline bool contains(const Solution& solution, ItemId j) const
    {
        return solution.items[j].in;
    }

    inline Weight overweight(const Solution& solution) const
    {
        return std::max((Weight)0, solution.weight - instance_.capacity());
    }

    inline void add(Solution& solution, ItemId j) const
    {
        solution.items[j].in = true;
        Weight w = instance_.item(j).weight;
        Profit p = instance_.item(j).profit;
        solution.weight += w;
        solution.profit += p;
        for (ItemId j_neighbor: instance_.item(j).neighbors) {
            if (contains(solution, j_neighbor))
                remove(solution, j_neighbor);
            solution.items[j_neighbor].neighbor_weight += w;
            solution.items[j_neighbor].neighbor_profit += p;
        }
    }

    inline void remove(Solution& solution, ItemId j) const
    {
        solution.items[j].in = false;
        Weight w = instance_.item(j).weight;
        Profit p = instance_.item(j).profit;
        solution.weight -= w;
        solution.profit -= p;
        for (ItemId j_neighbor: instance_.item(j).neighbors) {
            solution.items[j_neighbor].neighbor_weight -= w;
            solution.items[j_neighbor].neighbor_profit -= p;
        }
    }

    /*
     * Evaluate moves.
     */

    inline GlobalCost cost_remove(const Solution& solution, ItemId j) const
    {
        Weight w = solution.weight - instance_.item(j).weight;
        Profit p = solution.profit - instance_.item(j).profit;
        return {std::max((Weight)0, w - instance_.capacity()), -p};
    }

    inline GlobalCost cost_add(const Solution& solution, ItemId j) const
    {
        Weight w = solution.weight
            + instance_.item(j).weight
            - solution.items[j].neighbor_weight;
        Profit p = solution.profit
            + instance_.item(j).profit
            - solution.items[j].neighbor_profit;
        return {std::max((Weight)0, w - instance_.capacity()), -p};
    }

    inline GlobalCost cost_swap(
            const Solution& solution,
            ItemId j_in,
            ItemId j_out,
            GlobalCost) const
    {
        Weight w = solution.weight
                    - instance_.item(j_in).weight
                    + instance_.item(j_out).weight
                    - solution.items[j_out].neighbor_weight;
        Profit p = solution.profit
                    - instance_.item(j_in).profit
                    + instance_.item(j_out).profit
                    - solution.items[j_out].neighbor_profit;
        return {std::max((Weight)0, w - instance_.capacity()), -p};
    }

    /*
     * Private attributes.
     */

    /** Instance. */
    const Instance& instance_;

    /** Parameters. */
    Parameters parameters_;

    std::vector<int8_t> tabu_;

    std::vector<ItemId> items_;

    std::vector<ItemId> items_in_;

    std::vector<ItemId> items_out_;

    optimizationtools::IndexedSet neighbors_;

    optimizationtools::IndexedSet free_items_;

    optimizationtools::IndexedSet free_items_2_;

    /*
     * Statistics.
     */

    Counter toggle_number_of_explorations_ = 0;

    Counter toggle_number_of_sucesses_ = 0;

    Counter swap_number_of_explorations_ = 0;

    Counter swap_number_of_sucesses_ = 0;

    Counter swap_2_1_number_of_explorations_ = 0;

    Counter swap_2_1_number_of_sucesses_ = 0;

};

}

}

