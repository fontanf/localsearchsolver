#pragma once

#include "localsearchsolver/common.hpp"

namespace localsearchsolver
{

namespace roadef2020
{

using Cost = int64_t;
using SeasonId = int64_t;
using ScenarioId = int64_t;
using ResourceId = int64_t;
using ResourcePos = int64_t;
using InterventionId = int16_t;
using InterventionPos = int16_t;
using ExclusionId = int64_t;
using Time = int16_t;
using Workload = int64_t;
using Risk = int64_t;

struct Exclusion
{
    InterventionId j1;
    InterventionId j2;
    SeasonId season;

    InterventionId j(InterventionId j_cur) const { return (j_cur == j1)? j2: j1; }
};

struct InterventionRisk
{
    Time t_first; // default: t_cur + 1
    std::vector<std::vector<Risk>> risks;
    std::vector<std::vector<double>> risks_double;

    inline Risk risk(Time t_start, ScenarioId s) const
    {
        if (t_start < t_first
                || t_start - t_first >= (Time)risks.size())
            return 0;
        return risks[t_start - t_first][s];
    }
};

struct InterventionResourceWorkload
{
    Time t_first = 0; // default: t_cur + 1
    std::vector<Workload> workloads;
    std::vector<double> workloads_double;

    inline Workload workload(Time t_start) const
    {
        if (t_start < t_first
                || t_start - t_first >= (Time)workloads.size())
            return 0;
        return workloads[t_start - t_first];
    }
};

struct InterventionResource
{
    ResourceId r;
    std::vector<InterventionResourceWorkload> workloads;

    inline Workload workload(Time t_cur, Time t_start) const
    {
        return workloads[t_cur].workload(t_start);
    }
};

struct Intervention
{
    std::string name;
    std::vector<Time> deltas;
    /** risk[t_cur]. */
    std::vector<InterventionRisk> risks;
    Time t_start_max;
    /** workloads[r_pos]. */
    std::vector<InterventionResource> resources;
    /** exclusions[season] = {j1, j2, ..., j3}. */
    std::vector<std::vector<ExclusionId>> exclusions;
    InterventionPos exclusion_number = 0;
    std::vector<int8_t> fixed_assignments;
    double duration_mean = 0;
    std::vector<double> workload_mean;
};

struct Resource
{
    std::string name;
    /* Minimum usage of the resource at time t. */
    std::vector<Workload> min;
    /* Maximum usage of the resource at time t. */
    std::vector<Workload> max;
    std::vector<double> min_double;
    std::vector<double> max_double;
    Workload multiplier = 1;
    std::vector<bool> unconstrained;
};

class Instance
{

public:

    Instance() { }

    /** Create instance from file. */
    Instance(std::string filepath, std::string format);

    /**
     * Return a new instance containing only the subset of interventions of
     * 'interventions'.
     */
    Instance reduced_instance(const std::vector<bool>& interventions) const;

    /*
     * getters
     */

    inline double                        alpha() const { return alpha_; }
    inline double                     quantile() const { return quantile_; }
    inline Time                        horizon() const { return horizon_; }
    inline InterventionId  intervention_number() const { return interventions_.size(); }
    inline ResourceId          resource_number() const { return resources_.size(); }
    inline SeasonId              season_number() const { return season_strings_.size(); }
    inline ExclusionId        exclusion_number() const { return exclusions_.size(); }

    inline ScenarioId    least_common_multiple() const { return least_common_multiple_; }
    inline Risk                        alpha_1() const { return alpha_1_; }
    inline Risk                        alpha_2() const { return alpha_2_; }
    inline Risk               alpha_multiplier() const { return alpha_multiplier_; }
    inline Risk                risk_multiplier() const { return risk_multiplier_; }

    inline ScenarioId  scenario_number(Time t) const { return scenario_numbers_[t]; }
    inline SeasonId             season(Time t) const { return (t < (Time)seasons_.size())? seasons_[t]: -1; }

    inline Time                start_max(InterventionId j) const { return interventions_[j].t_start_max; }
    inline Time                 duration(InterventionId j, Time t) const { return interventions_[j].deltas[t]; }
    inline double          duration_mean(InterventionId j) const { return interventions_[j].duration_mean; }
    inline double          workload_mean(InterventionId j, ResourceId r) const { return interventions_[j].workload_mean[r]; }
    inline ExclusionId  exclusion_number(InterventionId j) const { return interventions_[j].exclusion_number; }
    inline Risk                     risk(InterventionId j, Time t_cur, Time t_start, ScenarioId s) const { return interventions_[j].risks[t_cur].risk(t_start, s); }
    inline ResourcePos   resource_number(InterventionId j) const { return interventions_[j].resources.size(); }
    inline ResourceId           resource(InterventionId j, ResourcePos r_pos) const { return interventions_[j].resources[r_pos].r; }
    inline Workload             workload(InterventionId j, ResourcePos r_pos, Time t_cur, Time t_start) const { return interventions_[j].resources[r_pos].workload(t_cur, t_start); }
    inline const std::vector<ExclusionId>& exclusions(InterventionId j, SeasonId season) const { return interventions_[j].exclusions[season]; }

    inline Workload         workload_min(ResourceId r, Time t) const { return resources_[r].min[t]; }
    inline Workload         workload_max(ResourceId r, Time t) const { return resources_[r].max[t]; }
    inline Workload  resource_multiplier(ResourceId r) const { return resources_[r].multiplier; }

    inline const Exclusion& exclusion(ExclusionId e) const { return exclusions_[e]; }

    std::string intervention_name(InterventionId j) const { return interventions_[j].name; }
    std::string     resource_name(ResourceId r)     const { return resources_[r].name; }
    std::string       season_name(SeasonId season)  const { return season_strings_[season]; }

    InterventionId intervention(std::string intervention_name) const { return string_to_interv_.at(intervention_name); }

    int fixed(InterventionId j, Time t_start) const { return interventions_[j].fixed_assignments[t_start]; }
    bool unconstrained(ResourceId r, Time t_cur) const { return resources_[r].unconstrained[t_cur]; }

    const Instance* reduced_instance() const { return reduced_instance_.get(); }

private:

    /*
     * Attributes
     */

    double alpha_;
    double quantile_;
    Time horizon_;
    std::vector<Intervention> interventions_;
    std::vector<Resource> resources_;
    /** scenario_numbers_[t]. */
    std::vector<ScenarioId> scenario_numbers_;
    /** seasons_[t]. */
    std::vector<SeasonId> seasons_;
    std::vector<Exclusion> exclusions_;

    std::map<std::string, InterventionId> string_to_interv_;
    std::map<std::string, ResourceId>     string_to_resource_;
    std::map<std::string, SeasonId>       string_to_season_;
    std::vector<std::string> season_strings_;

    ScenarioId least_common_multiple_ = 1;
    Risk alpha_1_;
    Risk alpha_2_;
    Risk alpha_multiplier_ = 1;
    Risk risk_multiplier_ = 1;

    std::unique_ptr<Instance> reduced_instance_ = nullptr;

    void fix_assignments();

};

std::ostream& operator<<(std::ostream& os, const Instance& instance);

class LocalScheme
{

public:

    /** Global cost: <Interventions, Conflicts, Overwork, Underwork, Cost>; */
    using GlobalCost = std::tuple<InterventionId, ExclusionId, Workload, Workload, Cost>;

    inline InterventionId&      intervention_number(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline ExclusionId&             conflict_number(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline Workload&                       overwork(GlobalCost& global_cost) { return std::get<2>(global_cost); }
    inline Workload&                      underwork(GlobalCost& global_cost) { return std::get<3>(global_cost); }
    inline Cost&                               cost(GlobalCost& global_cost) { return std::get<4>(global_cost); }
    inline InterventionId intervention_number(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline ExclusionId        conflict_number(const GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline Workload                  overwork(const GlobalCost& global_cost) { return std::get<2>(global_cost); }
    inline Workload                 underwork(const GlobalCost& global_cost) { return std::get<3>(global_cost); }
    inline Cost                          cost(const GlobalCost& global_cost) { return std::get<4>(global_cost); }

    static GlobalCost global_cost_worst()
    {
        return {
            std::numeric_limits<InterventionId>::max(),
            std::numeric_limits<ExclusionId>::max(),
            std::numeric_limits<Workload>::max(),
            std::numeric_limits<Workload>::max(),
            std::numeric_limits<Cost>::max(),
        };
    }

    /*
     * Solutions.
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
        InterventionId j1 = -1;
        InterventionId j2 = -1;
        Time t_cur = -1;
    };

    struct Solution
    {
        /** Number of interventions. */
        InterventionId intervention_number = 0;
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

    CompactSolution solution2compact(const Solution& solution)
    {
        return solution.intervention_starts;
    }

    Solution compact2solution(const CompactSolution& compact_solution)
    {
        auto solution = empty_solution();
        for (InterventionId j = 0; j < instance_.intervention_number(); ++j)
            if (compact_solution[j] != -1)
                add(solution, j, compact_solution[j]);
        return solution;
    }

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        /** Enable reduced instance. */
        double reduced_instance_time = -1;

        /** Enable repair procedure. */
        int repair = -1;
    };

    LocalScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters),
        interventions_(instance.intervention_number()),
        times_(instance.horizon())
    {
        if (parameters_.reduced_instance_time == -1)
            parameters_.reduced_instance_time = 60;
        if (parameters_.repair == -1)
            parameters_.repair = 1;

        std::iota(interventions_.begin(), interventions_.end(), 0);
        std::iota(times_.begin(), times_.end(), 0);
        ScenarioId s_max = 0;
        for (Time t_cur = 0; t_cur < instance.horizon(); ++t_cur)
            if (s_max < instance_.scenario_number(t_cur))
                s_max = instance_.scenario_number(t_cur);
        risks_.resize(s_max, 0);
        sorted_scenarios_.resize(s_max, 0);
    }

    LocalScheme(const LocalScheme& local_scheme):
        LocalScheme(local_scheme.instance_, local_scheme.parameters_) { }

    virtual ~LocalScheme() { }

    /*
     * Initial solutions.
     */

    Solution empty_solution() const;

    Solution initial_solution(
            Counter initial_solution_id,
            std::mt19937_64& generator);

    /*
     * Solution properties.
     */

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
        return (solution.intervention_number == instance_.intervention_number()
                && solution.conflicts.size() == 0
                && solution.overwork == 0
                && solution.underwork == 0);
    }

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            -solution.intervention_number,
            solution.conflicts.size(),
            solution.overwork,
            solution.underwork,
            cost(solution),
        };
    }

    /*
     * Local search.
     */

    struct Move
    {
        InterventionId j;
        Time t_start;
        GlobalCost global_cost;
    };

    static Move move_null() { return {-1, -1, global_cost_worst()}; };

    struct MoveHasher
    {
        std::hash<InterventionId> hasher_1;
        std::hash<Time> hasher_2;

        inline bool operator()(
                const Move& move_1,
                const Move& move_2) const
        {
            return move_1.j == move_2.j && move_1.t_start == move_2.t_start;
        }

        inline std::size_t operator()(
                const Move& move) const
        {
            size_t hash = hasher_1(move.j);
            optimizationtools::hash_combine(hash, hasher_2(move.t_start));
            return hash;
        }
    };

    inline MoveHasher move_hasher() const { return MoveHasher(); }

    std::vector<Move> perturbations(Solution& solution);

    inline void apply_move(Solution& solution, const Move& move) const
    {
        remove(solution, move.j);
        add(solution, move.j, move.t_start);
    }

    void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Move& tabu = move_null());

    /*
     * Outputs.
     */

    std::ostream& print(
            std::ostream &os,
            const Solution& solution)
    {
        for (InterventionId j = 0; j < instance_.intervention_number(); ++j) {
            os << "j " << j << " t_start " << solution.intervention_starts[j] << std::endl;
        }
        for (Time t_cur = 0; t_cur < instance_.horizon(); ++t_cur) {
            os << "t " << t_cur;
            for (ResourceId r = 0; r < instance_.resource_number(); ++r)
                os << " " << solution.time_steps[t_cur].workloads[r] << "/" << instance_.workload_max(r, t_cur);
            os << std::endl;
        }
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

        for (InterventionId j = 0; j < instance_.intervention_number(); ++j)
            cert << instance_.intervention_name(j)
                << " " << solution.intervention_starts[j] + 1 << std::endl;
    }

private:

    /*
     * Initial solutions.
     */

    std::pair<int, Solution> repair(
            const Solution& infeasible_solution);

    /*
     * Manipulate solutions.
     */

    inline double cost(const Solution& solution) const
    {
        return (instance_.alpha_1() * solution.mean_cost
                + instance_.alpha_2() * solution.expected_excess);
    }

    void add(Solution& solution, InterventionId j, Time t_start) const;

    void remove(Solution& solution, InterventionId j) const;

    /*
     * Evaluate moves.
     */

    GlobalCost cost_add(
            const Solution& solution,
            InterventionId j,
            Time t_start,
            const GlobalCost& cutoff);

    /*
     * Private attributes.
     */

    const Instance& instance_;
    Parameters parameters_;

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
