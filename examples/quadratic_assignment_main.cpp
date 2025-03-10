/**
 * Quadratic assignment problem
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/quadratic_assignment.hpp
 *
 * TODO
 *
 */

#include "read_args.hpp"

#include "localsearchsolver/common.hpp"

#include "orproblems/assignment/quadratic_assignment.hpp"

#include "optimizationtools/containers/indexed_set.hpp"
#include "optimizationtools//utils//utils.hpp"

using namespace localsearchsolver;
using namespace orproblems::quadratic_assignment;

class LocalScheme
{

public:

    struct Parameters
    {
        double crossover_ux_weight = 1.0;
        double crossover_swap_weight = 0.0;
        double crossover_insert_weight = 0.0;
    };

    /*
     * Constructors and destructor
     */

    LocalScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters),
        facilities_(instance.number_of_facilities()),
        facility_pairs_(instance.number_of_facilities())
    {
        // Initialize facilities_.
        std::iota(facilities_.begin(), facilities_.end(), 0);
        // Initialize facility_pairs_.
        for (FacilityId facility_id_1 = 0; facility_id_1 < instance.number_of_facilities(); ++facility_id_1)
            for (FacilityId facility_id_2 = facility_id_1 + 1; facility_id_2 < instance.number_of_facilities(); ++facility_id_2)
                facility_pairs_.push_back({facility_id_1, facility_id_2});
    }

    /*
     * Global cost
     */

    /** Global cost: <Cost>; */
    using GlobalCost = std::tuple<Cost>;

    inline Cost&       cost(GlobalCost& global_cost) const { return std::get<0>(global_cost); }
    inline Cost  cost(const GlobalCost& global_cost) const { return std::get<0>(global_cost); }

    inline GlobalCost global_cost_goal(double value) const
    {
        return {
            value,
        };
    }

    /*
     * Solutions
     */

    struct Solution
    {
        std::vector<LocationId> locations;
        Cost cost = 0;
    };

    inline Solution empty_solution() const
    {
        Solution solution;
        solution.locations.resize(instance_.number_of_facilities(), -1);
        return solution;
    }

    inline Solution initial_solution(
            Counter,
            std::mt19937_64& generator) const
    {
        Solution solution = empty_solution();
        std::vector<LocationId> locations(instance_.number_of_facilities(), -1);
        std::iota(locations.begin(), locations.end(), 0);
        std::shuffle(locations.begin(), locations.end(), generator);
        for (FacilityId facility_id = 0; facility_id < instance_.number_of_facilities(); ++facility_id)
            add(solution, facility_id, locations[facility_id]);
        return solution;
    }

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            solution.cost,
        };
    }

    /*
     * Local search
     */

    struct Perturbation;

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Perturbation& = Perturbation())
    {
        Counter it = 0;
        (void)it;
        for (;; ++it) {
            //std::cout << "it " << it
            //    << " cost " << to_string(global_cost(solution))
            //    << std::endl;
            std::shuffle(facility_pairs_.begin(), facility_pairs_.end(), generator);
            FacilityId facility_id_1_best = -1;
            FacilityId facility_id_2_best = -1;
            GlobalCost c_best = global_cost(solution);
            for (auto p: facility_pairs_) {
                FacilityId facility_id_1 = p.first;
                FacilityId facility_id_2 = p.second;
                GlobalCost c = cost_swap(
                        solution,
                        facility_id_1,
                        facility_id_2);
                if (c >= c_best)
                    continue;
                if (facility_id_1 != -1 && !dominates(c, c_best))
                    continue;
                facility_id_1_best = facility_id_1;
                facility_id_2_best = facility_id_2;
                c_best = c;
            }
            if (facility_id_1_best == -1)
                break;
            LocationId location_id_1 = solution.locations[facility_id_1_best];
            LocationId location_id_2 = solution.locations[facility_id_2_best];
            remove(solution, facility_id_1_best);
            remove(solution, facility_id_2_best);
            add(solution, facility_id_1_best, location_id_2);
            add(solution, facility_id_2_best, location_id_1);
        }
    }

    /*
     * Genetic local search
     */

    inline Solution crossover_ux(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        Solution solution = empty_solution();

        std::shuffle(facilities_.begin(), facilities_.end(), generator);
        optimizationtools::IndexedSet free_locations(instance_.number_of_facilities());
        for (LocationId location_id = 0;
                location_id < instance_.number_of_facilities();
                ++location_id)
            free_locations.add(location_id);

        // If a facility is assigned to the same location in both parent
        // solutions, keep it assigned to the same location in the child
        // solution.
        for (FacilityId facility_id: facilities_) {
            LocationId location_id_1 = solution_parent_1.locations[facility_id];
            LocationId location_id_2 = solution_parent_2.locations[facility_id];
            if (location_id_1 == location_id_2) {
                add(solution, facility_id, location_id_1);
                free_locations.remove(location_id_1);
            }
        }

        // If a location to which a facility is assigned in one of the parent
        // solution is still free, assign the facility to it in the child
        // solution.
        for (FacilityId facility_id: facilities_) {
            if (solution.locations[facility_id] != -1)
                continue;
            LocationId location_id_1 = solution_parent_1.locations[facility_id];
            LocationId location_id_2 = solution_parent_2.locations[facility_id];
            if (!free_locations.contains(location_id_1)
                    && !free_locations.contains(location_id_2))
                continue;
            LocationId location_id = -1;
            if (free_locations.contains(location_id_1)
                    && free_locations.contains(location_id_2)) {
                std::uniform_int_distribution<LocationId> d(0, 1);
                LocationId pos = d(generator);
                location_id = (pos == 0)? location_id_1: location_id_2;
            } else if (free_locations.contains(location_id_1)) {
                location_id = location_id_1;
            } else {
                location_id = location_id_2;
            }
            add(solution, facility_id, location_id);
            free_locations.remove(location_id);
        }

        // Assign a random location to the remaining facilities.
        for (FacilityId facility_id: facilities_) {
            if (solution.locations[facility_id] != -1)
                continue;
            free_locations.shuffle_in(generator);
            LocationId location_id = *free_locations.begin();
            add(solution, facility_id, location_id);
            free_locations.remove(location_id);
        }

        return solution;
    }

    /**
     * References:
     * - "A greedy genetic algorithm for the quadratic assignment problem"
     *   (Ahuja et al., 2000)
     *   https://doi.org/10.1016/S0305-0548(99)00067-2
     */
    inline Solution crossover_swap(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        Solution solution = empty_solution();
        (void)solution_parent_1;
        (void)solution_parent_2;
        (void)generator;
        return solution;
    }

    /**
     * References:
     * - "A greedy genetic algorithm for the quadratic assignment problem"
     *   (Ahuja et al., 2000)
     *   https://doi.org/10.1016/S0305-0548(99)00067-2
     */
    inline Solution crossover_insert(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        Solution solution = empty_solution();
        (void)solution_parent_1;
        (void)solution_parent_2;
        (void)generator;
        return solution;
    }

    inline Solution crossover(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        std::discrete_distribution<Counter> d_crossover({
                parameters_.crossover_ux_weight,
                parameters_.crossover_swap_weight,
                parameters_.crossover_insert_weight,
                });
        Counter x = d_crossover(generator);
        switch (x) {
        case 1: {
            return crossover_swap(solution_parent_1, solution_parent_2, generator);
        } case 2: {
            return crossover_insert(solution_parent_1, solution_parent_2, generator);
        } default: {
            return crossover_ux(solution_parent_1, solution_parent_2, generator);
        }
        }
    }

    inline FacilityId distance(
            const Solution& solution_1,
            const Solution& solution_2) const
    {
        FacilityId d = 0;
        for (FacilityId facility_id = 0;
                facility_id < instance_.number_of_facilities();
                ++facility_id) {
            if (solution_1.locations[facility_id]
                    != solution_2.locations[facility_id])
                d++;
        }
        return d;
    }

    /*
     * Iterated local search
     */

    struct Perturbation
    {
        Perturbation(): global_cost(worst<GlobalCost>()) { }

        std::vector<FacilityId> facilities;
        GlobalCost global_cost;
    };

    inline std::vector<Perturbation> perturbations(
            const Solution& solution,
            std::mt19937_64& generator)
    {
        std::vector<Perturbation> perturbations;
        for (FacilityId facility_id = 0;
                facility_id < instance_.number_of_facilities();
                ++facility_id) {
            Solution solution_tmp = solution;
            std::vector<FacilityId> facility_ids = {facility_id};
            FacilityId facility_id_1 = facility_id;
            for (Counter c = 0; c < 16; ++c) {
                FacilityId facility_id_2_best = -1;
                GlobalCost c_best = worst<GlobalCost>();
                std::shuffle(facilities_.begin(), facilities_.end(), generator);
                for (FacilityId facility_id_2: facilities_) {
                    if (std::find(facility_ids.begin(), facility_ids.end(), facility_id_2)
                            != facility_ids.end())
                        continue;
                    GlobalCost c = cost_swap(
                            solution_tmp,
                            facility_id_1,
                            facility_id_2);
                    if (c >= c_best)
                        continue;
                    facility_id_2_best = facility_id_2;
                    c_best = c;
                }
                if (facility_id_2_best == -1)
                    break;
                LocationId location_id_1 = solution_tmp.locations[facility_id_1];
                LocationId location_id_2 = solution_tmp.locations[facility_id_2_best];
                remove(solution_tmp, facility_id_1);
                remove(solution_tmp, facility_id_2_best);
                add(solution_tmp, facility_id_1, location_id_2);
                add(solution_tmp, facility_id_2_best, location_id_1);
                facility_ids.push_back(facility_id_2_best);
                facility_id_1 = facility_id_2_best;
            }

            Perturbation perturbation;
            perturbation.facilities = facility_ids;
            perturbation.global_cost = global_cost(solution_tmp);
            perturbations.push_back(perturbation);
        }
        return perturbations;
    }

    inline void apply_perturbation(
            Solution& solution,
            const Perturbation& perturbation,
            std::mt19937_64&) const
    {
        LocationId location_id_1 = solution.locations[perturbation.facilities[0]];
        remove(solution, perturbation.facilities[0]);
        for (FacilityId facility_pos = 0;
                facility_pos < (Counter)perturbation.facilities.size() - 1;
                ++facility_pos) {
            LocationId location_id = solution.locations[perturbation.facilities[facility_pos + 1]];
            remove(solution, perturbation.facilities[facility_pos + 1]);
            add(solution, perturbation.facilities[facility_pos], location_id);
        }
        add(solution, perturbation.facilities.back(), location_id_1);
        assert(global_cost(solution) == perturbation.global_cost);
    }

    /*
     * Best first local search
     */

    using CompactSolution = std::vector<LocationId>;

    struct CompactSolutionHasher
    {
        std::hash<LocationId> hasher;

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
            for (LocationId location_id: *compact_solution)
                optimizationtools::hash_combine(hash, hasher(location_id));
            return hash;
        }
    };

    inline CompactSolutionHasher compact_solution_hasher() const { return CompactSolutionHasher(); }

    CompactSolution solution2compact(const Solution& solution)
    {
        return solution.locations;
    }

    Solution compact2solution(const CompactSolution& compact_solution)
    {
        Solution solution;
        solution.locations = compact_solution;
        for (FacilityId facility_id_1 = 0; facility_id_1 < instance_.number_of_facilities(); ++facility_id_1) {
            LocationId location_id_1 = solution.locations[facility_id_1];
            if (location_id_1 == -1)
                continue;
            solution.cost += instance_.flow(facility_id_1, facility_id_1)
                * instance_.distance(location_id_1, location_id_1);
            for (FacilityId facility_id_2 = facility_id_1 + 1;
                    facility_id_2 < instance_.number_of_facilities(); ++facility_id_2) {
                LocationId location_id_2 = solution.locations[facility_id_2];
                if (location_id_2 == -1)
                    continue;
                solution.cost += instance_.flow(facility_id_1, facility_id_2)
                    * instance_.distance(location_id_1, location_id_2);
                solution.cost += instance_.flow(facility_id_2, facility_id_1)
                    * instance_.distance(location_id_2, location_id_1);
            }
        }
        return solution;
    }

    struct PerturbationHasher
    {
        std::hash<FacilityId> hasher;

        inline bool hashable(const Perturbation&) const { return true; }

        inline bool operator()(
                const Perturbation& perturbation_1,
                const Perturbation& perturbation_2) const
        {
            return perturbation_1.facilities[0] == perturbation_2.facilities[0];
        }

        inline std::size_t operator()(
                const Perturbation& perturbation) const
        {
            size_t hash = hasher(perturbation.facilities[0]);
            return hash;
        }
    };

    inline PerturbationHasher perturbation_hasher() const { return PerturbationHasher(); }

    /*
     * Outputs
     */

    void instance_format(
            std::ostream& os,
            int verbosity_level) const
    {
        os << "Quadratic assignment problem" << std::endl;
        instance_.format(os, verbosity_level);
    }

    void solution_format(
            std::ostream& os,
            const Solution& solution,
            int verbosity_level) const
    {
        (void)verbosity_level;
        os << "locations:";
        for (LocationId location_id: solution.locations)
            os << " " << location_id;
        os << std::endl;
        os << "cost: " << solution.cost << std::endl;
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

        for (LocationId location_id: solution.locations)
            file << location_id << " ";
    }

private:

    /*
     * Manipulate solutions
     */

    inline void add(Solution& solution, FacilityId facility_id, LocationId location_id) const
    {
        assert(solution.locations[facility_id] == -1);
        assert(std::find(solution.locations.begin(), solution.locations.end(), location_id) == solution.locations.end());
        solution.locations[facility_id] = location_id;
        solution.cost += instance_.flow(facility_id, facility_id)
            * instance_.distance(location_id, location_id);
        for (FacilityId facility_id_2 = 0; facility_id_2 < instance_.number_of_facilities(); facility_id_2++) {
            LocationId location_id_2 = solution.locations[facility_id_2];
            if (location_id_2 == -1 || location_id_2 == location_id)
                continue;
            solution.cost += instance_.flow(facility_id, facility_id_2)
                * instance_.distance(location_id, location_id_2);
            solution.cost += instance_.flow(facility_id_2, facility_id)
                * instance_.distance(location_id_2, location_id);
        }
    }

    inline void remove(Solution& solution, FacilityId facility_id) const
    {
        LocationId location_id = solution.locations[facility_id];
        assert(location_id != -1);
        solution.cost -= instance_.flow(facility_id, facility_id)
            * instance_.distance(location_id, location_id);
        for (FacilityId facility_id_2 = 0; facility_id_2 < instance_.number_of_facilities(); facility_id_2++) {
            LocationId location_id_2 = solution.locations[facility_id_2];
            if (location_id_2 == -1 || location_id_2 == location_id)
                continue;
            solution.cost -= instance_.flow(facility_id, facility_id_2)
                * instance_.distance(location_id, location_id_2);
            solution.cost -= instance_.flow(facility_id_2, facility_id)
                * instance_.distance(location_id_2, location_id);
        }
        solution.locations[facility_id] = -1;
    }

    /*
     * Evaluate moves
     */

    inline GlobalCost cost_swap(
            const Solution& solution,
            FacilityId facility_id_1,
            FacilityId facility_id_2) const
    {
        LocationId location_id_1 = solution.locations[facility_id_1];
        LocationId location_id_2 = solution.locations[facility_id_2];
        assert(location_id_1 != -1);
        assert(location_id_2 != -1);
        GlobalCost gc = global_cost(solution);
        cost(gc)
            += instance_.flow(facility_id_1, facility_id_1)
            * (instance_.distance(location_id_2, location_id_2)
                    - instance_.distance(location_id_1, location_id_1))
            + instance_.flow(facility_id_2, facility_id_2)
            * (instance_.distance(location_id_1, location_id_1)
                    - instance_.distance(location_id_2, location_id_2))
            + instance_.flow(facility_id_1, facility_id_2)
            * (instance_.distance(location_id_2, location_id_1)
                    - instance_.distance(location_id_1, location_id_2))
            + instance_.flow(facility_id_2, facility_id_1)
            * (instance_.distance(location_id_1, location_id_2)
                    - instance_.distance(location_id_2, location_id_1));
        for (FacilityId facility_id_3 = 0;
                facility_id_3 < instance_.number_of_facilities();
                ++facility_id_3) {
            if (facility_id_3 == facility_id_1)
                continue;
            if (facility_id_3 == facility_id_2)
                continue;
            LocationId location_id_3 = solution.locations[facility_id_3];
            cost(gc)
                += instance_.flow(facility_id_1, facility_id_3)
                * (instance_.distance(location_id_2, location_id_3)
                        - instance_.distance(location_id_1, location_id_3))
                + instance_.flow(facility_id_3, facility_id_1)
                * (instance_.distance(location_id_3, location_id_2)
                        - instance_.distance(location_id_3, location_id_1))
                + instance_.flow(facility_id_2, facility_id_3)
                * (instance_.distance(location_id_1, location_id_3)
                        - instance_.distance(location_id_2, location_id_3))
                + instance_.flow(facility_id_3, facility_id_2)
                * (instance_.distance(location_id_3, location_id_1)
                        - instance_.distance(location_id_3, location_id_2));
        }
        return gc;
    }

    /*
     * Private attributes
     */

    /** Instance. */
    const Instance& instance_;

    /** Parameters. */
    Parameters parameters_;

    std::vector<FacilityId> facilities_;

    std::vector<std::pair<FacilityId, FacilityId>> facility_pairs_;

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
    LocalScheme::Parameters parameters;
    LocalScheme local_scheme(instance, parameters);

    // Run algorithm.
    std::string algorithm = vm["algorithm"].as<std::string>();
    auto output =
        (algorithm == "multi-start-local-search")?
        run_multi_start_local_search(local_scheme, vm):
        (algorithm == "iterated-local-search")?
        run_iterated_local_search(local_scheme, vm):
        (algorithm == "best-first-local-search")?
        run_best_first_local_search(local_scheme, vm):
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
