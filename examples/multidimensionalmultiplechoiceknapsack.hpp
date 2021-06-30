#pragma once

#include "localsearchsolver/common.hpp"

#include "optimizationtools/indexed_set.hpp"

namespace localsearchsolver
{

namespace multidimensionalmultiplechoiceknapsack
{

using ItemId = int64_t;
using GroupId = int64_t;
using Profit = int64_t;
using ResourceId = int64_t;
using Weight = int64_t;

struct Item
{
    ItemId id;
    GroupId group_id;
    Profit profit;
    std::vector<Weight> weights;
};

class Instance
{

public:

    Instance() { }
    void add_resource(Weight capacity) { capacities_.push_back(capacity); }
    void add_item(GroupId group_id, Profit profit)
    {
        Item item;
        item.id = items_.size();
        item.group_id = group_id;
        item.profit = profit;
        item.weights.resize(resource_number(), 0);
        items_.push_back(item);
        while ((GroupId)groups_.size() <= group_id)
            groups_.push_back({});
        groups_[group_id].push_back(item.id);
    }
    void set_weight(ItemId j, ResourceId r, Weight weight) { items_[j].weights[r] = weight; }

    Instance(std::string instance_path, std::string format = "")
    {
        std::ifstream file(instance_path);
        if (!file.good()) {
            std::cerr << "\033[31m" << "ERROR, unable to open file \"" << instance_path << "\"" << "\033[0m" << std::endl;
            assert(false);
            return;
        }
        if (format == "" || format == "khan2002") {
            read_khan2002(file);
        } else if (format == "shojaei2013") {
            read_shojaei2013(file);
        } else if (format == "mansi2013") {
            read_mansi2013(file);
        } else {
            std::cerr << "\033[31m" << "ERROR, unknown instance format \"" << format << "\"" << "\033[0m" << std::endl;
        }
        file.close();
    }

    virtual ~Instance() { }

    inline ItemId item_number() const { return items_.size(); }
    inline GroupId group_number() const { return groups_.size(); }
    inline ResourceId resource_number() const { return capacities_.size(); }
    inline const Item& item(ItemId j) const { return items_[j]; }
    inline const std::vector<ItemId>& items(GroupId group_id) const { return groups_[group_id]; }
    inline Weight capacity(ResourceId r) const { return capacities_[r]; }

    std::pair<bool, Profit> check(std::string certificate_path)
    {
        std::ifstream file(certificate_path);
        if (!file.good()) {
            std::cerr << "\033[31m" << "ERROR, unable to open file \"" << certificate_path << "\"" << "\033[0m" << std::endl;
            assert(false);
            return {false, 0};
        }

        optimizationtools::IndexedSet groups(group_number());
        std::vector<Weight> weights(resource_number(), 0);
        Profit profit = 0;
        GroupId duplicates = 0;
        ItemId j = 0;
        while (file >> j) {
            if (groups.contains(item(j).group_id)) {
                duplicates++;
                std::cout << "GroupId " << item(j).group_id << " already in the knapsack." << std::endl;
            }
            groups.add(item(j).group_id);
            for (ResourceId r = 0; r < resource_number(); ++r)
                weights[r] += item(j).weights[r];
            profit += item(j).profit;
            std::cout << "Item: " << j
                << "; Group: " << item(j).group_id
                << "; Profit: " << item(j).profit
                << std::endl;
        }
        Weight overweight = 0;
        for (ResourceId r = 0; r < resource_number(); ++r)
            if (weights[r] > capacity(r))
                overweight += (weights[r] - capacity(r));
        bool feasible
            = (duplicates == 0)
            && (groups.size() == group_number())
            && (overweight == 0);
        std::cout << "---" << std::endl;
        std::cout << "Groups:                     " << groups.size() << " / " << group_number() << std::endl;
        std::cout << "Duplicates:                 " << duplicates << std::endl;
        std::cout << "Overweight:                 " << overweight << std::endl;
        std::cout << "Feasible:                   " << feasible << std::endl;
        std::cout << "Profit:                     " << profit << std::endl;
        return {feasible, profit};
    }

private:

    void read_khan2002(std::ifstream& file)
    {
        ItemId group_number = -1;
        ItemId group_size = -1;
        ResourceId resource_number = -1;
        file >> group_number >> group_size >> resource_number;

        ResourceId capacity = -1;
        for (ResourceId r = 0; r < resource_number; ++r) {
            file >> capacity;
            add_resource(capacity);
        }

        std::string tmp;
        Profit profit = -1;
        Weight weight = -1;
        ItemId j = 0;
        for (GroupId group_id = 0; group_id < group_number; group_id++) {
            file >> tmp;
            for (ItemId j_pos = 0; j_pos < group_size; ++j_pos) {
                file >> profit;
                add_item(group_id, profit);
                for (ResourceId r = 0; r < resource_number; ++r) {
                    file >> weight;
                    set_weight(j, r, weight);
                }
                ++j;
            }
        }
    }

    void read_shojaei2013(std::ifstream& file)
    {
        ItemId group_number = -1;
        ResourceId resource_number = -1;
        file >> group_number >> resource_number;

        ResourceId capacity = -1;
        for (ResourceId r = 0; r < resource_number; ++r) {
            file >> capacity;
            add_resource(capacity);
        }

        ItemId group_size = -1;
        Profit profit = -1;
        Weight weight = -1;
        ItemId j = 0;
        for (GroupId group_id = 0; group_id < group_number; group_id++) {
            file >> group_size;
            for (ItemId j_pos = 0; j_pos < group_size; ++j_pos) {
                file >> profit;
                add_item(group_id, profit);
                for (ResourceId r = 0; r < resource_number; ++r) {
                    file >> weight;
                    set_weight(j, r, weight);
                }
                ++j;
            }
        }
    }

    void read_mansi2013(std::ifstream& file)
    {
        ItemId group_number = -1;
        ItemId group_size = -1;
        ResourceId resource_number = -1;
        file >> group_number >> group_size >> resource_number;

        ResourceId capacity = -1;
        for (ResourceId r = 0; r < resource_number; ++r) {
            file >> capacity;
            add_resource(capacity);
        }

        Profit profit = -1;
        Weight weight = -1;
        ItemId j = 0;
        for (GroupId group_id = 0; group_id < group_number; group_id++) {
            for (ItemId j_pos = 0; j_pos < group_size; ++j_pos) {
                file >> profit;
                add_item(group_id, profit);
                for (ResourceId r = 0; r < resource_number; ++r) {
                    file >> weight;
                    set_weight(j, r, weight);
                }
                ++j;
            }
        }
    }

    std::vector<Weight> capacities_;
    std::vector<Item> items_;
    std::vector<std::vector<ItemId>> groups_;

};

std::ostream& operator<<(
        std::ostream &os, const Instance& instance)
{
    os << "item number " << instance.item_number() << std::endl;
    os << "group number " << instance.group_number() << std::endl;
    os << "resource number " << instance.resource_number() << std::endl;
    os << "capacities";
    for (ResourceId r = 0; r < instance.resource_number(); ++r)
        os << " " << instance.capacity(r);
    os << std::endl;
    os << "items" << std::endl;
    for (ItemId j = 0; j < instance.item_number(); ++j) {
        os << "item " << instance.item(j).id
            << " group " << instance.item(j).group_id
            << " profit " << instance.item(j).profit
            << " weights";
        for (ResourceId r = 0; r < instance.resource_number(); ++r)
            os << " " << instance.item(j).weights[r];
        os << std::endl;
    }
    return os;
}

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

    inline std::vector<Move> perturbations(Solution& solution) const
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

