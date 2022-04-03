#pragma once

/**
 * Multidimansional Multiple-Choice Knapsack Problem..
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/multidimensionalmultiplechoiceknapsack.hpp
 *
 */

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

    /** Global cost: <Number of groups, Overweight, Profit>; */
    using GlobalCost = std::tuple<GroupId, Weight, Profit>;

    inline GroupId&       number_of_groups(GlobalCost& global_cost) const { return std::get<0>(global_cost); }
    inline Weight&              overweight(GlobalCost& global_cost) const { return std::get<1>(global_cost); }
    inline Profit&                  profit(GlobalCost& global_cost) const { return std::get<2>(global_cost); }
    inline GroupId  number_of_groups(const GlobalCost& global_cost) const { return std::get<0>(global_cost); }
    inline Weight         overweight(const GlobalCost& global_cost) const { return std::get<1>(global_cost); }
    inline Profit             profit(const GlobalCost& global_cost) const { return std::get<2>(global_cost); }

    /*
     * Solutions.
     */

    using CompactSolution = std::vector<ItemId>;

    struct CompactSolutionHasher
    {
        std::hash<ItemId> hasher;

        inline bool operator()(
                const std::shared_ptr<CompactSolution>& compact_solution_1,
                const std::shared_ptr<CompactSolution>& compact_solution_2) const
        {
            return *compact_solution_1 == *compact_solution_2;
        }

        inline std::size_t operator()(
                const std::shared_ptr<CompactSolution>& compact_solution) const
        {
            size_t hash = 0;
            for (ItemId j: *compact_solution)
                optimizationtools::hash_combine(hash, hasher(j));
            return hash;
        }
    };

    inline CompactSolutionHasher compact_solution_hasher() const { return CompactSolutionHasher(); }

    struct Solution
    {
        std::vector<ItemId> items;
        std::vector<Weight> weights;
        GroupId number_of_groups = 0;
        Weight overweight = 0;
        Profit profit = 0;
    };

    CompactSolution solution2compact(const Solution& solution)
    {
        return solution.items;
    }

    Solution compact2solution(const CompactSolution& compact_solution)
    {
        Solution solution = empty_solution();
        for (GroupId group_id = 0; group_id < instance_.number_of_groups(); ++group_id) {
            ItemId j = compact_solution[group_id];
            if (j != -1)
                add(solution, group_id, j);
        }
        return solution;
    }

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
    };

    LocalScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters),
        groups_(instance.number_of_groups()),
        items_(instance.largest_group_size())
    {
        std::iota(groups_.begin(), groups_.end(), 0);
        std::iota(items_.begin(), items_.end(), 0);
    }

    LocalScheme(const LocalScheme& local_scheme):
        LocalScheme(local_scheme.instance_, local_scheme.parameters_) { }

    virtual ~LocalScheme() { }

    /*
     * Initial solutions.
     */

    inline Solution empty_solution() const
    {
        Solution solution;
        solution.items.resize(instance_.number_of_groups(), -1);
        solution.weights.resize(instance_.number_of_resources(), 0);
        return solution;
    }

    inline Solution initial_solution(
            Counter,
            std::mt19937_64& generator) const
    {
        Solution solution = empty_solution();
        for (GroupId group_id = 0; group_id < instance_.number_of_groups(); ++group_id) {
            std::uniform_int_distribution<ItemId> d(0, instance_.number_of_items(group_id) - 1);
            ItemId j = d(generator);
            add(solution, group_id, j);
        }
        return solution;
    }

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

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            -solution.number_of_groups,
            solution.overweight,
            -solution.profit,
        };
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
     * Local search.
     */

    struct Move
    {
        Move(): group_id(-1), global_cost(worst<GlobalCost>()) { }

        GroupId group_id;
        ItemId j;
        GlobalCost global_cost;
    };

    struct MoveHasher
    {
        std::hash<ItemId> hasher;

        inline bool hashable(const Move&) const { return true; }

        inline bool operator()(
                const Move& move_1,
                const Move& move_2) const
        {
            return move_1.group_id == move_2.group_id
                && move_1.j == move_2.j;
        }

        inline std::size_t operator()(
                const Move& move) const
        {
            std::size_t hash = hasher(move.group_id);
            optimizationtools::hash_combine(hash, move.j);
            return hash;
        }
    };

    inline MoveHasher move_hasher() const { return MoveHasher(); }

    inline std::vector<Move> perturbations(
            Solution& solution,
            std::mt19937_64&)
    {
        std::vector<Move> moves;
        for (GroupId group_id = 0; group_id < instance_.number_of_groups(); ++group_id) {
            ItemId j_old = solution.items[group_id];
            remove(solution, group_id, j_old);
            for (ItemId j = 0; j < instance_.number_of_items(group_id); ++j) {
                if (j == j_old)
                    continue;
                Move move;
                move.group_id = group_id;
                move.j = j;
                move.global_cost = cost_add(solution, group_id, j);
                moves.push_back(move);
            }
            add(solution, group_id, j_old);
        }
        return moves;
    }

    inline void apply_move(Solution& solution, const Move& move) const
    {
        remove(solution, move.group_id, solution.items[move.group_id]);
        add(solution, move.group_id, move.j);
    }

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Move& tabu = Move())
    {
        Counter it = 0;
        for (;; ++it) {
            //std::cout << "it " << it << " cost " << to_string(global_cost(solution)) << std::endl;
            std::shuffle(items_.begin(), items_.end(), generator);
            std::shuffle(groups_.begin(), groups_.end(), generator);
            GroupId group_id_best = -1;
            ItemId j_best = -1;
            GlobalCost c_best = global_cost(solution);
            for (GroupId group_id: groups_) {
                if (tabu.group_id != -1
                        && group_id == tabu.group_id)
                    continue;
                ItemId j_old = solution.items[group_id];
                remove(solution, group_id, j_old);
                for (ItemId j: items_) {
                    if (j >= instance_.number_of_items(group_id))
                        continue;
                    if (j == j_old)
                        continue;
                    GlobalCost c = cost_add(solution, group_id, j);
                    if (c >= c_best)
                        continue;
                    if (group_id_best != -1 && !dominates(c, c_best))
                        continue;
                    group_id_best = group_id;
                    j_best = j;
                    c_best = c;
                }
                add(solution, group_id, j_old);
            }
            if (group_id_best == -1)
                break;
            remove(solution, group_id_best, solution.items[group_id_best]);
            add(solution, group_id_best, j_best);
        }
    }

    /*
     * Outputs.
     */

    std::ostream& print(
            std::ostream &os,
            const Solution& solution) const
    {
        os << "number of groups: " << solution.number_of_groups << std::endl;
        os << "items:";
        for (ItemId j: solution.items)
            os << " " << j;
        os << std::endl;
        os << "weights:";
        for (Weight w: solution.weights)
            os << " " << w;
        os << std::endl;
        os << "overweight: " << solution.overweight << std::endl;
        os << "profit: " << solution.profit << std::endl;
        return os;
    }

    inline void write(
            const Solution& solution,
            std::string filepath) const
    {
        if (filepath.empty())
            return;
        std::ofstream cert(filepath);
        if (!cert.good()) {
            std::cerr << "\033[31m" << "ERROR, unable to open file \"" << filepath << "\"" << "\033[0m" << std::endl;
            return;
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
            ItemId j) const
    {
        assert(solution.items[group_id] == -1);
        // Update number_of_groups.
        solution.number_of_groups++;
        // Update weights.
        for (ResourceId r = 0; r < instance_.number_of_resources(); ++r) {
            Weight w_max = instance_.capacity(r);
            Weight w = instance_.item(group_id, j).weights[r];
            if (solution.weights[r] >= w_max) {
                solution.overweight += w;
            } else if (solution.weights[r] + w <= w_max) {
            } else {
                solution.overweight += (solution.weights[r] + w - w_max);
            }
            solution.weights[r] += w;
        }
        // Update profit.
        solution.profit += instance_.item(group_id, j).profit;
        // Update items.
        solution.items[group_id] = j;
    }

    inline void remove(
            Solution& solution,
            GroupId group_id,
            ItemId j) const
    {
        assert(solution.items[group_id] == j);
        // Update number_of_groups.
        solution.number_of_groups--;
        // Update weights.
        for (ResourceId r = 0; r < instance_.number_of_resources(); ++r) {
            Weight w_max = instance_.capacity(r);
            Weight w = instance_.item(group_id, j).weights[r];
            if (solution.weights[r] - w >= w_max) {
                solution.overweight -= w;
            } else if (solution.weights[r] <= w_max) {
            } else {
                solution.overweight -= (solution.weights[r] - w_max);
            }
            solution.weights[r] -= w;
        }
        // Update profit.
        solution.profit -= instance_.item(group_id, j).profit;
        // Update items.
        solution.items[group_id] = -1;
    }

    /*
     * Evaluate moves.
     */

    inline GlobalCost cost_add(
            const Solution& solution,
            GroupId group_id,
            ItemId j) const
    {
        GlobalCost gc = global_cost(solution);
        assert(solution.items[group_id] == -1);
        // Update number_of_groups.
        number_of_groups(gc)--;
        // Update overweigt.
        for (ResourceId r = 0; r < instance_.number_of_resources(); ++r) {
            Weight w_max = instance_.capacity(r);
            Weight w = instance_.item(group_id, j).weights[r];
            if (solution.weights[r] >= w_max) {
                overweight(gc) += w;
            } else if (solution.weights[r] + w <= w_max) {
            } else {
                overweight(gc) += (solution.weights[r] + w - w_max);
            }
        }
        // Update profit.
        profit(gc) -= instance_.item(group_id, j).profit;
        return gc;
    }

    /*
     * Private attributes.
     */

    const Instance& instance_;
    Parameters parameters_;

    std::vector<GroupId> groups_;
    std::vector<ItemId> items_;

};

}

}

