/**
 * ROADEF/EURO Challenge 2020: Maintenance Planning Problem.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/roadef2020.hpp
 *
 * TODO
 *
 */

#pragma once

#include "localsearchsolver/common.hpp"

#include "orproblems/roadef2020.hpp"

namespace localsearchsolver
{

namespace roadef2020
{

using namespace orproblems::roadef2020;

class LocalScheme
{

public:

    /*
     * Constructors and destructor.
     */

    LocalScheme(
            const Instance& instance):
        instance_(instance),
        interventions_(instance.number_of_interventions()),
        times_(instance.horizon())
    {
        // Initialize interventions_.
        std::iota(interventions_.begin(), interventions_.end(), 0);
        // Initialize times_.
        std::iota(times_.begin(), times_.end(), 0);
        // Initialize risks_ and sorted_scenarios_.
        ScenarioId s_max = 0;
        for (Time time_cur = 0; time_cur < instance.horizon(); ++time_cur)
            if (s_max < instance_.number_of_scenarios(time_cur))
                s_max = instance_.number_of_scenarios(time_cur);
        risks_.resize(s_max, 0);
        sorted_scenarios_.resize(s_max, 0);
    }

    /*
     * Glocal cost.
     */

    /**
     * Global cost:
     * - Number of conflicts
     * - Overwork
     * - Underwork
     * - Cost
     */
    using GlobalCost = std::tuple<ExclusionId, Workload, Workload, Cost>;

    inline ExclusionId&       number_of_conflicts(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Workload&                     overwork(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline Workload&                    underwork(GlobalCost& global_cost) { return std::get<2>(global_cost); }
    inline Cost&                             cost(GlobalCost& global_cost) { return std::get<3>(global_cost); }
    inline ExclusionId  number_of_conflicts(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Workload                overwork(const GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline Workload               underwork(const GlobalCost& global_cost) { return std::get<2>(global_cost); }
    inline Cost                        cost(const GlobalCost& global_cost) { return std::get<3>(global_cost); }

    /*
     * Solutions.
     */

    struct SolutionTimeStep
    {
        /** Workload for each resource, workloads[r]. */
        std::vector<Workload> workloads;

        /** Risk for each scenario, risks[scenario]. */
        std::vector<Risk> risks;

        /** Scenarios sorted by risk. */
        std::vector<ScenarioId> sorted_scenarios;

        /** Sum of risks of all scenarios. */
        Risk risk_sum = 0;
    };

    struct SolutionConflict
    {
        /** First intervention of the conflict. */
        InterventionId intervention_id_1 = -1;

        /** Second intervention of the conflict. */
        InterventionId intervention_id_2 = -1;

        /** Time of the conflict. */
        Time time_cur = -1;
    };

    struct Solution
    {
        /** Start date for each intervention, -1 if not in solution. */
        std::vector<Time> intervention_starts;

        /** Informations for each time step. */
        std::vector<SolutionTimeStep> time_steps;

        /** List of pair of interventions in conflict. */
        std::vector<SolutionConflict> conflicts;

        /** Mean cost. */
        Cost mean_cost = 0;

        /** Expected excess. */
        Cost expected_excess = 0;

        /** Underwork. */
        Workload underwork = 0;

        /** Overwork. */
        Workload overwork = 0;
    };

    Solution empty_solution() const;

    Solution initial_solution(
            Counter initial_solution_id,
            std::mt19937_64& generator);

    inline double real_cost(const Solution& solution) const
    {
        return (double)cost(solution)
            / instance_.horizon()
            / instance_.least_common_multiple()
            / instance_.alpha_multiplier()
            / instance_.risk_multiplier();
    }

    inline bool feasible(const Solution& solution) const
    {
        return (solution.conflicts.size() == 0
                && solution.overwork == 0
                && solution.underwork == 0);
    }

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            solution.conflicts.size(),
            solution.overwork,
            solution.underwork,
            cost(solution),
        };
    }

    /*
     * Local search.
     */

    struct Perturbation;

    void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Perturbation& tabu = Perturbation());

    /*
     * Iterated local search.
     */

    struct Perturbation
    {
        Perturbation(): intervention_id(-1), time_start(-1), global_cost(worst<GlobalCost>()) { }

        InterventionId intervention_id;
        Time time_start;
        GlobalCost global_cost;
    };

    std::vector<Perturbation> perturbations(
            Solution& solution,
            std::mt19937_64&);

    inline void apply_perturbation(
            Solution& solution,
            const Perturbation& perturbation,
            std::mt19937_64&) const
    {
        remove(solution, perturbation.intervention_id);
        add(solution, perturbation.intervention_id, perturbation.time_start);
    }

    /*
     * Best first local search.
     */

    using CompactSolution = std::vector<Time>;

    struct CompactSolutionHasher
    {
        std::hash<Time> hasher;

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
            for (Time t: *compact_solution)
                optimizationtools::hash_combine(hash, hasher(t));
            return hash;
        }
    };

    inline CompactSolutionHasher compact_solution_hasher() const { return CompactSolutionHasher(); }

    CompactSolution solution2compact(const Solution& solution)
    {
        return solution.intervention_starts;
    }

    Solution compact2solution(const CompactSolution& compact_solution)
    {
        auto solution = empty_solution();
        for (InterventionId j = 0; j < instance_.number_of_interventions(); ++j)
            if (compact_solution[j] != -1)
                add(solution, j, compact_solution[j]);
        return solution;
    }

    struct PerturbationHasher
    {
        std::hash<InterventionId> hasher_1;
        std::hash<Time> hasher_2;

        inline bool hashable(const Perturbation&) const { return true; }

        inline bool operator()(
                const Perturbation& perturbation_1,
                const Perturbation& perturbation_2) const
        {
            return perturbation_1.intervention_id == perturbation_2.intervention_id
                && perturbation_1.time_start == perturbation_2.time_start;
        }

        inline std::size_t operator()(
                const Perturbation& perturbation) const
        {
            size_t hash = hasher_1(perturbation.intervention_id);
            optimizationtools::hash_combine(hash, hasher_2(perturbation.time_start));
            return hash;
        }
    };

    inline PerturbationHasher perturbation_hasher() const { return PerturbationHasher(); }

    /*
     * Outputs.
     */

    std::ostream& print(
            std::ostream &os,
            const Solution& solution)
    {
        for (InterventionId intervention_id = 0;
                intervention_id < instance_.number_of_interventions();
                ++intervention_id) {
            os << "intervention_id " << intervention_id
                << " time_start " << solution.intervention_starts[intervention_id]
                << std::endl;
        }
        for (Time time_cur = 0; time_cur < instance_.horizon(); ++time_cur) {
            os << "t " << time_cur;
            for (ResourceId r = 0; r < instance_.number_of_resources(); ++r)
                os << " " << solution.time_steps[time_cur].workloads[r]
                    << "/" << instance_.workload_max(r, time_cur);
            os << std::endl;
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

        for (InterventionId j = 0; j < instance_.number_of_interventions(); ++j)
            cert << instance_.intervention_name(j)
                << " " << solution.intervention_starts[j] + 1 << std::endl;
    }

private:

    /*
     * Manipulate solutions.
     */

    inline double cost(const Solution& solution) const
    {
        return (instance_.alpha_1() * solution.mean_cost
                + instance_.alpha_2() * solution.expected_excess);
    }

    /** Add an intervention to a solution. */
    void add(
            Solution& solution,
            InterventionId intervention_id,
            Time time_start) const;

    /** Remove an intervention from a solution. */
    void remove(
            Solution& solution,
            InterventionId intervention_id) const;

    /*
     * Evaluate moves.
     */

    /** Compute the cost of adding an intervention. */
    GlobalCost cost_add(
            const Solution& solution,
            InterventionId intervention_id,
            Time time_start,
            const GlobalCost& cutoff);

    /*
     * Private attributes.
     */

    /** Instance. */
    const Instance& instance_;

    /*
     * Temporary structures.
     */

    std::vector<InterventionId> interventions_;

    std::vector<Time> times_;

    std::vector<Risk> risks_;

    std::vector<ScenarioId> sorted_scenarios_;

};

}

}
