#pragma once

#include "localsearchsolver/common.hpp"

#include "external/treesearchsolver/examples/knapsackwithconflicts.hpp"

namespace localsearchsolver
{

namespace knapsackwithconflicts
{

using ItemId = treesearchsolver::knapsackwithconflicts::ItemId;
using ItemPos = treesearchsolver::knapsackwithconflicts::ItemPos;
using Profit = treesearchsolver::knapsackwithconflicts::Profit;
using Weight = treesearchsolver::knapsackwithconflicts::Weight;
using Instance = treesearchsolver::knapsackwithconflicts::Instance;

class LocalScheme
{

public:

    /** Global cost: <Overweight, Profit>; */
    using GlobalCost = std::tuple<Weight, Profit>;

    inline Weight&       overweight(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Profit&           profit(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline Weight  overweight(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Profit      profit(const GlobalCost& global_cost) { return std::get<1>(global_cost); }

    static GlobalCost global_cost_worst()
    {
        return {
            std::numeric_limits<Weight>::max(),
            std::numeric_limits<Profit>::max(),
        };
    }

    /*
     * Solutions.
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

    struct Solution
    {
        std::vector<bool> items;
        Profit profit = 0;
        Weight weight = 0;
    };

    CompactSolution solution2compact(const Solution& solution)
    {
        return solution.items;
    }

    Solution compact2solution(const CompactSolution& compact_solution)
    {
        auto solution = empty_solution();
        for (ItemId j = 0; j < instance_.item_number(); ++j)
            if (compact_solution[j])
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
        items_(instance.item_number())
    {
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
        solution.items.resize(instance_.item_number(), false);
        return solution;
    }

    inline Solution initial_solution(
            Counter,
            std::mt19937_64& generator) const
    {
        Solution solution = empty_solution();
        std::uniform_int_distribution<int> d(0, 1);
        for (ItemId j = 0; j < instance_.item_number(); ++j)
            if (d(generator) == 0)
                add(solution, j);
        return solution;
    }

    /*
     * Solution properties.
     */

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
            size_t hash = hasher(move.j);
            return hash;
        }
    };

    inline MoveHasher move_hasher() const { return MoveHasher(); }

    inline std::vector<Move> perturbations(
            const Solution& solution,
            std::mt19937_64&)
    {
        std::vector<Move> moves;
        for (ItemId j: items_) {
            GlobalCost c = (contains(solution, j))?
                cost_remove(solution, j, global_cost_worst()):
                cost_add(solution, j, global_cost_worst());
            Move move;
            move.j = j;
            move.global_cost = c;
            moves.push_back(move);
        }
        return moves;
    }

    inline void apply_move(Solution& solution, const Move& move) const
    {
        if (contains(solution, move.j)) {
            remove(solution, move.j);
        } else {
            add(solution, move.j);
        }
    }

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Move& tabu = move_null())
    {
        for (;;) {
            std::shuffle(items_.begin(), items_.end(), generator);
            ItemId j_best = -1;
            GlobalCost c_best = global_cost(solution);
            for (ItemId j: items_) {
                if (j == tabu.j)
                    continue;
                GlobalCost c = (contains(solution, j))?
                    cost_remove(solution, j, c_best):
                    cost_add(solution, j, c_best);
                if (c >= c_best)
                    continue;
                if (j_best != -1 && !dominates(c, c_best))
                    continue;
                j_best = j;
                c_best = c;
            }
            if (j_best == -1)
                break;
            if (contains(solution, j_best)) {
                remove(solution, j_best);
            } else {
                add(solution, j_best);
            }
        }
    }

    /*
     * Outputs.
     */

    std::ostream& print(
            std::ostream &os,
            const Solution& solution)
    {
        os << "items:";
        for (ItemId j = 0; j < instance_.item_number(); ++j)
            if (contains(solution, j))
                os << " " << j;
        os << std::endl;
        os << "weight: " << solution.weight << " / " << instance_.capacity() << std::endl;
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

        for (ItemId j = 0; j < instance_.item_number(); ++j)
            if (contains(solution, j))
                cert << j << " ";
    }

private:

    /*
     * Manipulate solutions.
     */

    inline bool contains(const Solution& solution, ItemId j) const
    {
        return solution.items[j];
    }

    inline Weight overweight(const Solution& solution) const
    {
        return std::max((Weight)0, solution.weight - instance_.capacity());
    }

    inline void add(Solution& solution, ItemId j) const
    {
        solution.items[j] = true;
        solution.weight += instance_.item(j).weight;
        solution.profit += instance_.item(j).profit;
        for (ItemId j_neighbor: instance_.item(j).neighbors)
            if (contains(solution, j_neighbor))
                remove(solution, j_neighbor);
    }

    inline void remove(Solution& solution, ItemId j) const
    {
        solution.items[j] = false;
        solution.weight -= instance_.item(j).weight;
        solution.profit -= instance_.item(j).profit;
    }

    /*
     * Evaluate moves.
     */

    inline GlobalCost cost_remove(const Solution& solution, ItemId j, GlobalCost) const
    {
        return {
            std::max((Weight)0, (solution.weight - instance_.item(j).weight) - instance_.capacity()),
            - (solution.profit - instance_.item(j).profit),
        };
    }

    inline GlobalCost cost_add(const Solution& solution, ItemId j, GlobalCost) const
    {
        Profit profit = solution.profit + instance_.item(j).profit;
        Weight weight = solution.weight + instance_.item(j).weight;
        for (ItemId j_neighbor: instance_.item(j).neighbors) {
            if (contains(solution, j_neighbor)) {
                profit -= instance_.item(j_neighbor).profit;
                weight -= instance_.item(j_neighbor).weight;
            }
        }
        return {
            std::max((Weight)0, weight - instance_.capacity()),
            -profit,
        };
    }

    /*
     * Private attributes.
     */

    const Instance& instance_;
    Parameters parameters_;

    std::vector<ItemId> items_;

};

}

}

