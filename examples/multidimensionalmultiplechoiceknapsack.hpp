#pragma once

/**
 * Multidimansional Multiple-Choice Knapsack Problem..
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/multidimensionalmultiplechoiceknapsack.hpp
 *
 * TODO
 *
 */

#include "localsearchsolver/common.hpp"

#include "orproblems/multidimensionalmultiplechoiceknapsack.hpp"

#include "optimizationtools/indexed_set.hpp"

namespace localsearchsolver
{

namespace multidimensionalmultiplechoiceknapsack
{

using namespace orproblems::multidimensionalmultiplechoiceknapsack;

class LocalScheme
{

public:

    /** Global cost: <Group number, Overweight, Profit>; */
    using GlobalCost = std::tuple<GroupId, Weight, Profit>;

    inline GroupId&       group_number(GlobalCost& global_cost) const { return std::get<0>(global_cost); }
    inline Weight&          overweight(GlobalCost& global_cost) const { return std::get<1>(global_cost); }
    inline Profit&              profit(GlobalCost& global_cost) const { return std::get<2>(global_cost); }
    inline GroupId  group_number(const GlobalCost& global_cost) const { return std::get<0>(global_cost); }
    inline Weight     overweight(const GlobalCost& global_cost) const { return std::get<1>(global_cost); }
    inline Profit         profit(const GlobalCost& global_cost) const { return std::get<2>(global_cost); }

    static GlobalCost global_cost_worst()
    {
        return {
            std::numeric_limits<GroupId>::max(),
            std::numeric_limits<Weight>::max(),
            std::numeric_limits<Profit>::max(),
        };
    }

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
        GroupId group_number = 0;
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
        for (ItemId j: compact_solution)
            if (j != -1)
                add(solution, j);
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
        items_(instance.group_number())
    {
        for (GroupId group_id = 0; group_id < instance.group_number(); ++group_id)
            items_[group_id] = instance.items(group_id);
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
        solution.items.resize(instance_.group_number(), -1);
        solution.weights.resize(instance_.resource_number(), 0);
        return solution;
    }

    inline Solution initial_solution(
            Counter,
            std::mt19937_64& generator) const
    {
        Solution solution = empty_solution();
        for (GroupId group_id = 0; group_id < instance_.group_number(); ++group_id) {
            std::uniform_int_distribution<ItemId> d(0, instance_.items(group_id).size() - 1);
            ItemId j = instance_.items(group_id)[d(generator)];
            add(solution, j);
        }
        return solution;
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            -solution.group_number,
            solution.overweight,
            -solution.profit,
        };
    }

    /*
     * Local search.
     */

    struct Move
    {
        ItemId j;
        GlobalCost global_cost;
    };

    static Move move_null() { return {-1, global_cost_worst()}; }

    struct MoveHasher
    {
        std::hash<ItemId> hasher;

        inline bool hashable(const Move&) const { return true; }

        inline bool operator()(
                const Move& move_1,
                const Move& move_2) const
        {
            return move_1.j == move_2.j;
        }

        inline std::size_t operator()(
                const Move& move) const
        {
            return hasher(move.j);
        }
    };

    inline MoveHasher move_hasher() const { return MoveHasher(); }

    inline std::vector<Move> perturbations(
            Solution& solution,
            std::mt19937_64&)
    {
        std::vector<Move> moves;
        for (GroupId group_id = 0; group_id < instance_.group_number(); ++group_id) {
            ItemId j_old = solution.items[group_id];
            remove(solution, j_old);
            for (ItemId j: instance_.items(group_id)) {
                if (j == j_old)
                    continue;
                Move move;
                move.j = j;
                move.global_cost = cost_add(solution, j);
                moves.push_back(move);
            }
            add(solution, j_old);
        }
        return moves;
    }

    inline void apply_move(Solution& solution, const Move& move) const
    {
        GroupId group_id = instance_.item(move.j).group_id;
        remove(solution, solution.items[group_id]);
        add(solution, move.j);
    }

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Move& tabu = move_null())
    {
        Counter it = 0;
        for (;; ++it) {
            //std::cout << "it " << it << " cost " << to_string(global_cost(solution)) << std::endl;
            std::shuffle(items_.begin(), items_.end(), generator);
            ItemId j_best = -1;
            GlobalCost c_best = global_cost(solution);
            for (auto& items: items_) {
                std::shuffle(items.begin(), items.end(), generator);
                GroupId group_id = instance_.item(items[0]).group_id;
                if (tabu.j != -1 && group_id == instance_.item(tabu.j).group_id)
                    continue;
                ItemId j_old = solution.items[group_id];
                remove(solution, j_old);
                for (ItemId j: items) {
                    if (j == j_old)
                        continue;
                    GlobalCost c = cost_add(solution, j);
                    if (c >= c_best)
                        continue;
                    if (j_best != -1 && !dominates(c, c_best))
                        continue;
                    j_best = j;
                    c_best = c;
                }
                add(solution, j_old);
            }
            if (j_best == -1)
                break;
            GroupId group_id = instance_.item(j_best).group_id;
            remove(solution, solution.items[group_id]);
            add(solution, j_best);
        }
    }

    /*
     * Outputs.
     */

    std::ostream& print(
            std::ostream &os,
            const Solution& solution) const
    {
        os << "group number: " << solution.group_number << std::endl;
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

    inline void add(Solution& solution, ItemId j) const
    {
        assert(j >= 0);
        assert(j < instance_.item_number());
        GroupId group_id = instance_.item(j).group_id;
        assert(solution.items[group_id] == -1);
        // Update group_number.
        solution.group_number++;
        // Update weights.
        for (ResourceId r = 0; r < instance_.resource_number(); ++r) {
            Weight w_max = instance_.capacity(r);
            Weight w = instance_.item(j).weights[r];
            if (solution.weights[r] >= w_max) {
                solution.overweight += w;
            } else if (solution.weights[r] + w <= w_max) {
            } else {
                solution.overweight += (solution.weights[r] + w - w_max);
            }
            solution.weights[r] += w;
        }
        // Update profit.
        solution.profit += instance_.item(j).profit;
        // Update items.
        solution.items[group_id] = j;
    }

    inline void remove(Solution& solution, ItemId j) const
    {
        GroupId group_id = instance_.item(j).group_id;
        assert(solution.items[group_id] == j);
        // Update group_number.
        solution.group_number--;
        // Update weights.
        for (ResourceId r = 0; r < instance_.resource_number(); ++r) {
            Weight w_max = instance_.capacity(r);
            Weight w = instance_.item(j).weights[r];
            if (solution.weights[r] - w >= w_max) {
                solution.overweight -= w;
            } else if (solution.weights[r] <= w_max) {
            } else {
                solution.overweight -= (solution.weights[r] - w_max);
            }
            solution.weights[r] -= w;
        }
        // Update profit.
        solution.profit -= instance_.item(j).profit;
        // Update items.
        solution.items[group_id] = -1;
    }

    /*
     * Evaluate moves.
     */

    inline GlobalCost cost_add(
            const Solution& solution,
            ItemId j) const
    {
        GlobalCost gc = global_cost(solution);
        assert(solution.items[instance_.item(j).group_id] == -1);
        // Update group_number.
        group_number(gc)--;
        // Update overweigt.
        for (ResourceId r = 0; r < instance_.resource_number(); ++r) {
            Weight w_max = instance_.capacity(r);
            Weight w = instance_.item(j).weights[r];
            if (solution.weights[r] >= w_max) {
                overweight(gc) += w;
            } else if (solution.weights[r] + w <= w_max) {
            } else {
                overweight(gc) += (solution.weights[r] + w - w_max);
            }
        }
        // Update profit.
        profit(gc) -= instance_.item(j).profit;
        return gc;
    }

    /*
     * Private attributes.
     */

    const Instance& instance_;
    Parameters parameters_;

    std::vector<std::vector<ItemId>> items_;

};

}

}

