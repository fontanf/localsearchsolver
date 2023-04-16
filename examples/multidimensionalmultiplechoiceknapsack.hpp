/**
 * Multidimansional Multiple-Choice Knapsack Problem.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/multidimensionalmultiplechoiceknapsack.hpp
 */

#pragma once

#include "localsearchsolver/common.hpp"

#include "orproblems/multidimensionalmultiplechoiceknapsack.hpp"

#include "optimizationtools/containers/indexed_set.hpp"

namespace localsearchsolver
{

namespace multidimensionalmultiplechoiceknapsack
{

using namespace orproblems::multidimensionalmultiplechoiceknapsack;

class LocalScheme
{

public:

    /*
     * Constructors and destructor.
     */

    LocalScheme(const Instance& instance):
        instance_(instance),
        groups_(instance.number_of_groups()),
        items_(instance.largest_group_size())
    {
        std::iota(groups_.begin(), groups_.end(), 0);
        std::iota(items_.begin(), items_.end(), 0);
    }

    /*
     * Global cost.
     */

    /** Global cost: <Overweight, Profit>; */
    using GlobalCost = std::tuple<Weight, Profit>;

    inline Weight&       overweight(GlobalCost& global_cost) const { return std::get<0>(global_cost); }
    inline Profit&           profit(GlobalCost& global_cost) const { return std::get<1>(global_cost); }
    inline Weight  overweight(const GlobalCost& global_cost) const { return std::get<0>(global_cost); }
    inline Profit      profit(const GlobalCost& global_cost) const { return std::get<1>(global_cost); }

    /*
     * Solution.
     */

    struct Solution
    {
        std::vector<ItemId> items;
        std::vector<Weight> weights;
        Weight overweight = 0;
        Profit profit = 0;
    };

    inline Solution empty_solution() const
    {
        Solution solution;
        solution.items.resize(instance_.number_of_groups(), -1);
        solution.weights.resize(instance_.number_of_resources(), 0);
        return solution;
    }

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            solution.overweight,
            -solution.profit,
        };
    }

    inline Solution initial_solution(
            Counter,
            std::mt19937_64& generator) const
    {
        Solution solution = empty_solution();
        for (GroupId group_id = 0; group_id < instance_.number_of_groups(); ++group_id) {
            std::uniform_int_distribution<ItemId> d(0, instance_.number_of_items(group_id) - 1);
            ItemId item_id = d(generator);
            add(solution, group_id, item_id);
        }
        return solution;
    }

    /*
     * Local search.
     */

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator)
    {
        Counter it = 0;
        (void)it;
        for (;; ++it) {
            //std::cout << "it " << it << " cost " << to_string(global_cost(solution)) << std::endl;
            std::shuffle(items_.begin(), items_.end(), generator);
            std::shuffle(groups_.begin(), groups_.end(), generator);
            GroupId group_id_best = -1;
            ItemId item_id_best = -1;
            GlobalCost c_best = global_cost(solution);
            for (GroupId group_id: groups_) {
                ItemId item_id_old = solution.items[group_id];
                remove(solution, group_id, item_id_old);
                for (ItemId item_id: items_) {
                    if (item_id >= instance_.number_of_items(group_id))
                        continue;
                    if (item_id == item_id_old)
                        continue;
                    GlobalCost c = cost_add(solution, group_id, item_id);
                    if (c >= c_best)
                        continue;
                    if (group_id_best != -1 && !dominates(c, c_best))
                        continue;
                    group_id_best = group_id;
                    item_id_best = item_id;
                    c_best = c;
                }
                add(solution, group_id, item_id_old);
            }
            if (group_id_best == -1)
                break;
            remove(solution, group_id_best, solution.items[group_id_best]);
            add(solution, group_id_best, item_id_best);
        }
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
        std::uniform_int_distribution<int> distribution(0, 1);
        for (GroupId group_id = 0; group_id < instance_.number_of_groups(); ++group_id) {
            if (distribution(generator) == 0) {
                add(solution, group_id, solution_parent_1.items[group_id]);
            } else {
                add(solution, group_id, solution_parent_2.items[group_id]);
            }
        }
        return solution;
    }

    inline ItemId distance(
            const Solution& solution_1,
            const Solution& solution_2) const
    {
        ItemId d = 0;
        for (GroupId group_id = 0; group_id < instance_.number_of_groups(); ++group_id)
            if (solution_1.items[group_id] != solution_2.items[group_id])
                d++;
        return d;
    }

    /*
     * Outputs.
     */

    std::ostream& print(
            std::ostream &os,
            const Solution& solution,
            int verbosity_level)
    {
        if (verbosity_level >= 1) {
            os << "Overweight:        " << solution.overweight << std::endl;
            os << "Profit:            " << solution.profit << std::endl;
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

        for (ItemId j: solution.items)
            cert << j << " ";
    }

private:

    /*
     * Manipulate solutions.
     */

    inline void add(
            Solution& solution,
            GroupId group_id,
            ItemId item_id) const
    {
        assert(solution.items[group_id] == -1);
        // Update weights.
        for (ResourceId r = 0; r < instance_.number_of_resources(); ++r) {
            Weight w_max = instance_.capacity(r);
            Weight w = instance_.item(group_id, item_id).weights[r];
            if (solution.weights[r] >= w_max) {
                solution.overweight += w;
            } else if (solution.weights[r] + w <= w_max) {
            } else {
                solution.overweight += (solution.weights[r] + w - w_max);
            }
            solution.weights[r] += w;
        }
        // Update profit.
        solution.profit += instance_.item(group_id, item_id).profit;
        // Update items.
        solution.items[group_id] = item_id;
    }

    inline void remove(
            Solution& solution,
            GroupId group_id,
            ItemId item_id) const
    {
        assert(solution.items[group_id] == item_id);
        // Update weights.
        for (ResourceId r = 0; r < instance_.number_of_resources(); ++r) {
            Weight w_max = instance_.capacity(r);
            Weight w = instance_.item(group_id, item_id).weights[r];
            if (solution.weights[r] - w >= w_max) {
                solution.overweight -= w;
            } else if (solution.weights[r] <= w_max) {
            } else {
                solution.overweight -= (solution.weights[r] - w_max);
            }
            solution.weights[r] -= w;
        }
        // Update profit.
        solution.profit -= instance_.item(group_id, item_id).profit;
        // Update items.
        solution.items[group_id] = -1;
    }

    /*
     * Evaluate moves.
     */

    inline GlobalCost cost_add(
            const Solution& solution,
            GroupId group_id,
            ItemId item_id) const
    {
        GlobalCost gc = global_cost(solution);
        assert(solution.items[group_id] == -1);
        // Update overweigt.
        for (ResourceId r = 0; r < instance_.number_of_resources(); ++r) {
            Weight w_max = instance_.capacity(r);
            Weight w = instance_.item(group_id, item_id).weights[r];
            if (solution.weights[r] >= w_max) {
                overweight(gc) += w;
            } else if (solution.weights[r] + w <= w_max) {
            } else {
                overweight(gc) += (solution.weights[r] + w - w_max);
            }
        }
        // Update profit.
        profit(gc) -= instance_.item(group_id, item_id).profit;
        return gc;
    }

    /*
     * Private attributes.
     */

    /** Instance. */
    const Instance& instance_;

    std::vector<GroupId> groups_;

    std::vector<ItemId> items_;

};

}

}

