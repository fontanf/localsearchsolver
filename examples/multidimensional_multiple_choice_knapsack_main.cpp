/**
 * Multidimansional multiple-choice knapsack problem
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/multidimensional_multiple_choice_knapsack.hpp
 */

#include "read_args.hpp"

#include "localsearchsolver/common.hpp"

#include "orproblems/packing/multidimensional_multiple_choice_knapsack.hpp"

using namespace localsearchsolver;
using namespace orproblems::multidimensional_multiple_choice_knapsack;

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
     * Local search
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
     * Outputs
     */

    void instance_format(
            std::ostream& os,
            int verbosity_level) const
    {
        os << "Multidimansional multiple-choice knapsack problem" << std::endl;
        instance_.format(os, verbosity_level);
    }

    void solution_format(
            const Solution& solution,
            std::ostream& os,
            int verbosity_level) const
    {
        if (verbosity_level >= 1) {
            os
                << "Overweight:  " << solution.overweight << std::endl
                << "Profit:      " << solution.profit << std::endl
                ;
        }
    }

    void solution_write(
            const Solution& solution,
            const std::string& certificate_path) const
    {
        if (certificate_path.empty())
            return;
        std::ofstream file(certificate_path);
        if (!file.good()) {
            throw std::runtime_error(
                    "Unable to open file \"" + certificate_path + "\".");
        }

        for (ItemId item_id: solution.items)
            file << item_id << " ";
    }

private:

    /*
     * Manipulate solutions
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
     * Evaluate moves
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
     * Private attributes
     */

    /** Instance. */
    const Instance& instance_;

    std::vector<GroupId> groups_;

    std::vector<ItemId> items_;

};

int main(int argc, char *argv[])
{
    // Create command line options.
    boost::program_options::options_description desc = setup_args();
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        std::cout << desc << std::endl;;
        throw "";
    }
    try {
        boost::program_options::notify(vm);
    } catch (const boost::program_options::required_option& e) {
        std::cout << desc << std::endl;;
        throw "";
    }

    // Create instance.
    InstanceBuilder instance_builder;
    instance_builder.read(
            vm["input"].as<std::string>(),
            vm["format"].as<std::string>());
    const Instance instance = instance_builder.build();

    // Create local scheme.
    LocalScheme local_scheme(instance);

    // Run algorithm.
    std::string algorithm = vm["algorithm"].as<std::string>();
    auto output =
        (algorithm == "multi-start-local-search")?
        run_multi_start_local_search(local_scheme, vm):
        //(algorithm == "iterated-local-search")?
        //run_iterated_local_search(local_scheme, vm):
        //(algorithm == "best-first-local-search")?
        //run_best_first_local_search(local_scheme, vm):
        run_genetic_local_search(local_scheme, vm);

    // Run checker.
    if (vm["print-checker"].as<int>() > 0
            && vm["certificate"].as<std::string>() != "") {
        std::cout << std::endl
            << "Checker" << std::endl
            << "-------" << std::endl;
        instance.check(
                vm["certificate"].as<std::string>(),
                std::cout,
                vm["print-checker"].as<int>());
    }

    return 0;
}
