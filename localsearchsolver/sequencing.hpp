#pragma once

/**
 * Sequencing problems.
 */

#include "localsearchsolver/common.hpp"

namespace localsearchsolver
{

namespace sequencing
{

using JobId = int64_t;
using JobPos = int64_t;

struct Parameters
{
    JobPos shift_bloc_maximum_length = 3;
    JobPos shift_maximum_distance = 1024;

    bool swap = true;
    JobPos swap_maximum_distance = 1024;

    bool reverse = false;
    JobPos reverse_maximum_length = 1024;

    bool shuffle_neighborhood_order = true;
    Counter number_of_perturbations = 10;

    double crossover_ox_weight = 1;
    double crossover_sjox_weight = 0;
    double crossover_sbox_weight = 0;
};

template <typename LocalScheme0>
class LocalScheme
{

public:

    using GlobalCost = typename LocalScheme0::GlobalCost;

    static GlobalCost global_cost_worst() { return LocalScheme0::global_cost_worst(); }

    using Solution = typename LocalScheme0::Solution;

    using CompactSolution = std::vector<JobId>;

    struct CompactSolutionHasher
    {
        std::hash<JobId> hasher;

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
            for (JobId j: *compact_solution)
                optimizationtools::hash_combine(hash, hasher(j));
            return hash;
        }
    };

    inline CompactSolutionHasher compact_solution_hasher() const { return CompactSolutionHasher(); }

    Solution compact2solution(const CompactSolution& compact_solution)
    {
        auto solution = local_scheme_0_.empty_solution();
        for (JobId j: compact_solution)
            local_scheme_0_.append(solution, j);
        return solution;
    }

    CompactSolution solution2compact(const Solution& solution)
    {
        return local_scheme_0_.jobs(solution);
    }

    /*
     * Constructors and destructor.
     */

    LocalScheme(
            LocalScheme0& local_scheme_0,
            Parameters parameters):
        local_scheme_0_(local_scheme_0),
        parameters_(parameters),
        positions1_(local_scheme_0_.number_of_jobs()),
        positions2_(local_scheme_0_.number_of_jobs()),
        global_costs_shift_(local_scheme_0_.number_of_jobs() + 1, local_scheme_0_.global_cost_worst()),
        global_costs_swap_(local_scheme_0_.number_of_jobs())
    {
        std::iota(positions1_.begin(), positions1_.end(), 0);
        std::iota(positions2_.begin(), positions2_.end(), 0);
        for (JobPos pos_1 = 0; pos_1 < local_scheme_0_.number_of_jobs(); ++pos_1)
            for (JobPos pos_2 = pos_1 + 1; pos_2 < local_scheme_0_.number_of_jobs(); ++pos_2)
                pairs_.push_back({pos_1, pos_2});
        for (JobPos pos_1 = 0; pos_1 < local_scheme_0_.number_of_jobs(); ++pos_1)
            global_costs_swap_[pos_1].resize(
                    local_scheme_0_.number_of_jobs() - pos_1,
                    local_scheme_0_.global_cost_worst());
    }

    LocalScheme(const LocalScheme& sequencing_scheme):
        LocalScheme(sequencing_scheme.local_scheme_0_, sequencing_scheme.parameters_) { }

    virtual ~LocalScheme() { }

    inline Solution empty_solution() const { return local_scheme_0_.empty_solution(); }

    inline Solution initial_solution(
            Counter initial_solution_id,
            std::mt19937_64& generator)
    {
        return local_scheme_0_.initial_solution(initial_solution_id, generator);
    }

    inline Solution crossover(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        std::discrete_distribution<Counter> d_crossover({
                parameters_.crossover_ox_weight,
                parameters_.crossover_sjox_weight,
                parameters_.crossover_sbox_weight,
                });
        Counter x = d_crossover(generator);
        switch (x) {
        case 0: {
            return crossover_ox(solution_parent_1, solution_parent_2, generator);
        } case 1: {
            return crossover_sjox(solution_parent_1, solution_parent_2, generator);
        } case 2: {
            return crossover_sbox(solution_parent_1, solution_parent_2, generator);
        } default: {
            return crossover_ox(solution_parent_1, solution_parent_2, generator);
        }
        }
    }

    /**
     * Generate a new solution from two parent solutions using the OX crossover
     * operator.
     *
     * The OX crossover operator consists in selecting a random substring from
     * the first parent, copying this substring into the child while leaving
     * the rest of the positions empty, and finally completing the child’s
     * missing positions circularly with the visits from the second parent,
     * starting from the end cutting point.
     *
     * References:
     * - "A simple and effective hybrid genetic search for the job sequencing
     *   and tool switching problem" (Mecler et al., 2021)
     *   https://doi.org/10.1016/j.cor.2020.105153
     */
    inline Solution crossover_ox(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        JobPos n = local_scheme_0_.jobs(solution_parent_1).size();

        std::vector<JobPos> edges = optimizationtools::bob_floyd<JobPos>(
                2, n + 1, generator);
        std::sort(edges.begin(), edges.end());

        JobPos pos_1 = edges[0];
        JobPos pos_2 = edges[1];

        std::vector<uint8_t> in_substring(n, false);
        for (JobPos pos = pos_1; pos < pos_2; ++pos) {
            JobId j = local_scheme_0_.jobs(solution_parent_1)[pos];
            in_substring[j] = true;
        }

        Solution solution = empty_solution();
        for (JobPos pos = 0; pos < n; ++pos) {
            if ((JobPos)local_scheme_0_.jobs(solution).size() == pos_1) {
                for (JobPos p = pos_1; p < pos_2; ++p) {
                    JobId j = local_scheme_0_.jobs(solution_parent_1)[p];
                    local_scheme_0_.append(solution, j);
                }
            }
            JobId j = local_scheme_0_.jobs(solution_parent_2)[pos];
            if (in_substring[j])
                continue;
            local_scheme_0_.append(solution, j);
        }
        if ((JobPos)local_scheme_0_.jobs(solution).size() == pos_1) {
            for (JobPos p = pos_1; p < pos_2; ++p) {
                JobId j = local_scheme_0_.jobs(solution_parent_1)[p];
                local_scheme_0_.append(solution, j);
            }
        }

        return solution;
    }

    /**
     * Generate a new solution from two parent solutions using the SJOX
     * crossover operator.
     *
     * References:
     * - "Two new robust genetic algorithms for the flowshop scheduling
     *   problem" (Ruiz et al., 2006)
     *   https://doi.org/10.1016/j.omega.2004.12.006
     */
    inline Solution crossover_sjox(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        JobPos n = local_scheme_0_.jobs(solution_parent_1).size();

        Solution solution = empty_solution();
        std::vector<JobPos> positions(n, -1);

        // Add jobs from parent_1 up to a given cut point.
        std::uniform_int_distribution<JobPos> d_point(1, n);
        JobPos pos_0 = d_point(generator);
        for (JobPos pos = 0; pos < pos_0; ++pos) {
            JobId j = local_scheme_0_.jobs(solution_parent_1)[pos];
            positions[j] = pos;
            local_scheme_0_.append(solution, j);
        }

        // Add jobs from parent_2 keeping the relative order.
        for (JobPos pos = 0; pos < n; ++pos) {
            // Add jobs which have the same positions in both parents.
            for (;;) {
                JobPos p = local_scheme_0_.jobs(solution).size();
                if (p == n)
                    break;
                JobId j1 = local_scheme_0_.jobs(solution_parent_1)[p];
                JobId j2 = local_scheme_0_.jobs(solution_parent_2)[p];
                if (j1 == j2) {
                    positions[j1] = p;
                    local_scheme_0_.append(solution, j1);
                    continue;
                }
                break;
            }
            JobId j = local_scheme_0_.jobs(solution_parent_2)[pos];
            if (positions[j] != -1)
                continue;
            positions[j] = local_scheme_0_.jobs(solution).size();
            local_scheme_0_.append(solution, j);
        }

        return solution;
    }

    /**
     * Generate a new solution from two parent solutions using the SBOX
     * crossover operator.
     *
     * References:
     * - "Two new robust genetic algorithms for the flowshop scheduling
     *   problem" (Ruiz et al., 2006)
     *   https://doi.org/10.1016/j.omega.2004.12.006
     */
    inline Solution crossover_sbox(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        JobPos n = local_scheme_0_.jobs(solution_parent_1).size();

        Solution solution = empty_solution();
        std::vector<JobPos> positions(n, -1);

        // Add jobs from parent_1 up to a given cut point.
        std::uniform_int_distribution<JobPos> d_point(1, n);
        JobPos pos_0 = d_point(generator);
        for (JobPos pos = 0; pos < pos_0; ++pos) {
            JobId j = local_scheme_0_.jobs(solution_parent_1)[pos];
            positions[j] = pos;
            local_scheme_0_.append(solution, j);
        }

        // Add jobs from parent_2 keeping the relative order.
        for (JobPos pos = 0; pos < n; ++pos) {
            // Add jobs which have the same positions in both parents.
            for (;;) {
                JobPos p = local_scheme_0_.jobs(solution).size();
                if (p == n)
                    break;
                JobId j1 = local_scheme_0_.jobs(solution_parent_1)[p];
                JobId j2 = local_scheme_0_.jobs(solution_parent_2)[p];
                if (p <= n - 1) {
                    JobId j1_next = local_scheme_0_.jobs(solution_parent_1)[p + 1];
                    JobId j2_next = local_scheme_0_.jobs(solution_parent_2)[p + 1];
                    if (j1 == j2 && j1_next == j2_next) {
                        positions[j1] = p;
                        local_scheme_0_.append(solution, j1);
                        continue;
                    }
                }
                if (p >= 1) {
                    JobId j1_prev = local_scheme_0_.jobs(solution_parent_1)[p - 1];
                    JobId j2_prev = local_scheme_0_.jobs(solution_parent_2)[p - 1];
                    if (j1 == j2 && j1_prev == j2_prev) {
                        positions[j1] = p;
                        local_scheme_0_.append(solution, j1);
                        continue;
                    }
                }
                break;
            }
            JobId j = local_scheme_0_.jobs(solution_parent_2)[pos];
            if (positions[j] != -1)
                continue;
            positions[j] = local_scheme_0_.jobs(solution).size();
            local_scheme_0_.append(solution, j);
        }

        return solution;
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return local_scheme_0_.global_cost(solution);
    }

    /**
     * Return the distance between two solutions.
     *
     * The distance between two solutions (represented as job permutations) is
     * defined as the number of different edges between them (broken-pairs
     * distance).
     *
     * References:
     * - "A simple and effective hybrid genetic search for the job sequencing
     *   and tool switching problem" (Mecler et al., 2021)
     *   https://doi.org/10.1016/j.cor.2020.105153
     */
    inline JobPos distance(
            const Solution& solution_1,
            const Solution& solution_2) const
    {
        JobPos n = local_scheme_0_.jobs(solution_1).size();
        std::vector<JobId> next_1(n, -1);
        for (JobPos pos = 0; pos < n - 1; ++pos) {
            JobId j = local_scheme_0_.jobs(solution_1)[pos];
            JobId j_next = local_scheme_0_.jobs(solution_1)[pos + 1];
            next_1[j] = j_next;
        }

        JobPos d = 0;
        for (JobPos pos = 0; pos < n - 1; ++pos) {
            JobId j = local_scheme_0_.jobs(solution_2)[pos];
            JobId j_next = local_scheme_0_.jobs(solution_2)[pos + 1];
            if (j_next != next_1[j])
                d++;
        }
        return d;
    }

    /*
     * Local search.
     */

    struct Move
    {
        JobPos pos_1;
        JobPos pos_2;
        JobPos pos_3;
        JobPos pos_4;
        GlobalCost global_cost;
    };

    static Move move_null() { return {-1, -1, -1, -1, LocalScheme0::global_cost_worst()}; }

    struct MoveHasher
    {
        inline bool hashable(const Move&) const { return false; }
        inline bool operator()(const Move&, const Move&) const { return false; }
        inline std::size_t operator()(const Move&) const { return 0; }
    };

    inline MoveHasher move_hasher() const { return MoveHasher(); }

    inline std::vector<Move> perturbations(
            const Solution& solution,
            std::mt19937_64& generator)
    {
        std::vector<Move> moves;
        for (Counter perturbation = 0; perturbation < parameters_.number_of_perturbations; ++perturbation) {
            std::vector<JobPos> edges = optimizationtools::bob_floyd<JobPos>(
                    4, local_scheme_0_.jobs(solution).size() + 1, generator);
            std::sort(edges.begin(), edges.end());
            Move move;
            move.pos_1 = edges[0];
            move.pos_2 = edges[1];
            move.pos_3 = edges[2];
            move.pos_4 = edges[3];
            assert(move.pos_1 >= 0);
            assert(move.pos_4 <= (JobPos)local_scheme_0_.jobs(solution).size());
            move.global_cost = global_cost(solution);
            moves.push_back(move);
        }
        return moves;
    }

    inline void apply_move(Solution& solution, const Move& move)
    {
        solution_tmp_ = local_scheme_0_.empty_solution();
        for (JobPos pos = 0; pos < move.pos_1; ++pos)
            local_scheme_0_.append(solution_tmp_, local_scheme_0_.jobs(solution)[pos]);
        for (JobPos pos = move.pos_3; pos < move.pos_4; ++pos)
            local_scheme_0_.append(solution_tmp_, local_scheme_0_.jobs(solution)[pos]);
        for (JobPos pos = move.pos_2; pos < move.pos_3; ++pos)
            local_scheme_0_.append(solution_tmp_, local_scheme_0_.jobs(solution)[pos]);
        for (JobPos pos = move.pos_1; pos < move.pos_2; ++pos)
            local_scheme_0_.append(solution_tmp_, local_scheme_0_.jobs(solution)[pos]);
        for (JobPos pos = move.pos_4; pos < (JobPos)local_scheme_0_.jobs(solution).size(); ++pos)
            local_scheme_0_.append(solution_tmp_, local_scheme_0_.jobs(solution)[pos]);
        solution = solution_tmp_;
        assert((JobPos)local_scheme_0_.jobs(solution).size() <= local_scheme_0_.number_of_jobs());
    }

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Move& = move_null())
    {
        //if (tabu.j != -1)
        //    std::cout << "j " << tabu.j << " j_prev " << tabu.j_prev
        //        << std::endl;
        //print(std::cout, solution);
        //std::cout << to_string(global_cost(solution)) << std::endl;

        // Get neighborhoods.
        std::vector<Counter> neighborhoods;
        if (parameters_.swap)
            neighborhoods.push_back(0);
        if (parameters_.reverse)
            neighborhoods.push_back(-1);
        for (JobPos bloc_size = 1; bloc_size <= parameters_.shift_bloc_maximum_length; ++bloc_size)
            neighborhoods.push_back(bloc_size);

        Counter it = 0;
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
                case 0: { // Swap neighborhood.
                    std::shuffle(pairs_.begin(), pairs_.end(), generator);
                    compute_cost_swap(solution);
                    JobPos pos_1_best = -1;
                    JobPos pos_2_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (auto pair: pairs_) {
                        GlobalCost c = global_costs_swap_[pair.first][pair.second - pair.first - 1];
                        if (c >= c_best)
                            continue;
                        if (pos_1_best != -1 && !dominates(c, c_best))
                            continue;
                        pos_1_best = pair.first;
                        pos_2_best = pair.second;
                        c_best = c;
                    }
                    if (pos_1_best != -1) {
                        //std::cout << "swap"
                        //    << " pos_1_best " << pos_1_best
                        //    << " pos_2_best " << pos_2_best
                        //    << std::endl;
                        improved = true;
                        // Apply best move.
                        solution_tmp_ = local_scheme_0_.empty_solution();
                        for (JobPos pos = 0; pos < local_scheme_0_.number_of_jobs(); ++pos) {
                            JobId j = local_scheme_0_.jobs(solution)[pos];
                            if (pos == pos_1_best)
                                j = local_scheme_0_.jobs(solution)[pos_2_best];
                            if (pos == pos_2_best)
                                j = local_scheme_0_.jobs(solution)[pos_1_best];
                            local_scheme_0_.append(solution_tmp_, j);
                        }
                        solution = solution_tmp_;
                        assert(local_scheme_0_.global_cost(solution) == c_best);
                    }
                    break;
                } case -1: { // Reverse neighborhood.
                    std::shuffle(pairs_.begin(), pairs_.end(), generator);
                    compute_cost_reverse(solution);
                    JobPos pos_1_best = -1;
                    JobPos pos_2_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (auto pair: pairs_) {
                        GlobalCost c = global_costs_swap_[pair.first][pair.second - pair.first - 1];
                        if (c >= c_best)
                            continue;
                        if (pos_1_best != -1 && !dominates(c, c_best))
                            continue;
                        pos_1_best = pair.first;
                        pos_2_best = pair.second;
                        c_best = c;
                    }
                    if (pos_1_best != -1) {
                        //std::cout << "reverse"
                        //    << " pos_1_best " << pos_1_best
                        //    << " pos_2_best " << pos_2_best
                        //    << std::endl;
                        improved = true;
                        // Apply best move.
                        solution_tmp_ = local_scheme_0_.empty_solution();
                        for (JobPos pos = 0; pos < pos_1_best; ++pos) {
                            JobId j = local_scheme_0_.jobs(solution)[pos];
                            local_scheme_0_.append(solution_tmp_, j);
                        }
                        for (JobPos pos = pos_2_best; pos >= pos_1_best; --pos) {
                            JobId j = local_scheme_0_.jobs(solution)[pos];
                            local_scheme_0_.append(solution_tmp_, j);
                        }
                        for (JobPos pos = pos_2_best + 1; pos < local_scheme_0_.number_of_jobs(); ++pos) {
                            JobId j = local_scheme_0_.jobs(solution)[pos];
                            local_scheme_0_.append(solution_tmp_, j);
                        }
                        solution = solution_tmp_;
                        assert(local_scheme_0_.global_cost(solution) == c_best);
                    }
                    break;
                } default: { // Shift neighborhood.
                    JobPos bloc_size = neighborhood;
                    std::shuffle(positions1_.begin(), positions1_.end(), generator);
                    std::shuffle(positions2_.begin(), positions2_.end(), generator);
                    JobPos pos_best = -1;
                    JobPos pos_new_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (JobPos pos: positions1_) {
                        if (pos > (JobPos)local_scheme_0_.jobs(solution).size() - bloc_size)
                            continue;
                        compute_cost_shift(solution, pos, bloc_size);
                        for (JobPos pos_new: positions2_) {
                            if (pos == pos_new || pos_new > (JobPos)local_scheme_0_.jobs(solution).size() - bloc_size)
                                continue;
                            GlobalCost c = global_costs_shift_[pos_new];
                            if (c >= c_best)
                                continue;
                            if (pos_best != -1 && !dominates(c, c_best))
                                continue;
                            pos_best = pos;
                            pos_new_best = pos_new;
                            c_best = c;
                        }
                    }
                    if (pos_best != -1) {
                        //std::cout << "shift " << neighborhood
                        //    << " pos_best " << pos_best
                        //    << " pos_new_best " << pos_new_best
                        //    << std::endl;
                        improved = true;
                        assert(c_best < local_scheme_0_.global_cost(solution));
                        // Apply best move.
                        //std::cout << bloc_size << std::endl;
                        solution_tmp_ = local_scheme_0_.empty_solution();
                        if (pos_best > pos_new_best) {
                            for (JobPos p = 0; p < pos_new_best; ++p)
                                local_scheme_0_.append(solution_tmp_, local_scheme_0_.jobs(solution)[p]);
                            for (JobPos p = pos_best; p < pos_best + bloc_size; ++p)
                                local_scheme_0_.append(solution_tmp_, local_scheme_0_.jobs(solution)[p]);
                            for (JobPos p = pos_new_best; p < pos_best; ++p)
                                local_scheme_0_.append(solution_tmp_, local_scheme_0_.jobs(solution)[p]);
                            for (JobPos p = pos_best + bloc_size; p < (JobPos)local_scheme_0_.jobs(solution).size(); ++p)
                                local_scheme_0_.append(solution_tmp_, local_scheme_0_.jobs(solution)[p]);
                        } else {
                            for (JobPos p = 0; p < pos_best; ++p)
                                local_scheme_0_.append(solution_tmp_, local_scheme_0_.jobs(solution)[p]);
                            for (JobPos p = pos_best + bloc_size; p < pos_new_best + bloc_size; ++p)
                                local_scheme_0_.append(solution_tmp_, local_scheme_0_.jobs(solution)[p]);
                            for (JobPos p = pos_best; p < pos_best + bloc_size; ++p)
                                local_scheme_0_.append(solution_tmp_, local_scheme_0_.jobs(solution)[p]);
                            for (JobPos p = pos_new_best + bloc_size; p < (JobPos)local_scheme_0_.jobs(solution).size(); ++p)
                                local_scheme_0_.append(solution_tmp_, local_scheme_0_.jobs(solution)[p]);
                        }
                        solution = solution_tmp_;
                        assert((JobPos)local_scheme_0_.jobs(solution).size() <= local_scheme_0_.number_of_jobs());
                        assert(local_scheme_0_.global_cost(solution) == c_best);
                    }
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
     * Outputs.
     */

    std::ostream& print(
            std::ostream &os,
            const Solution& solution) const
    {
        return local_scheme_0_.print(os, solution);
    }

    inline void write(
            const Solution& solution,
            std::string certificate_path) const
    {
        local_scheme_0_.write(solution, certificate_path);
    }

private:

    /*
     * Evaluate moves.
     */

    inline void compute_cost_shift(
            const Solution& solution,
            JobPos pos,
            JobPos size)
    {
        // Initialize solution_cur_.
        Solution solution_cur_ = local_scheme_0_.empty_solution();
        // Reset global_costs_shift_.
        std::fill(
                global_costs_shift_.begin(),
                global_costs_shift_.end(),
                local_scheme_0_.global_cost_worst());

        // Loop through all new positions.
        JobPos pos_min = std::max(
                (JobPos)0,
                pos - parameters_.shift_maximum_distance);
        JobPos pos_max = std::min(
                (JobPos)local_scheme_0_.jobs(solution).size() - size,
                pos + parameters_.shift_maximum_distance);
        for (JobPos pos_new = pos_min; pos_new <= pos_max; ++pos_new) {
            // Initialize solution_tmp_.
            Solution solution_tmp_ = solution_cur_;

            // Add bloc to times_.
            for (JobPos p = pos; p < pos + size; ++p) {
                // Add job to solution_tmp_.
                JobId j = local_scheme_0_.jobs(solution)[p];
                local_scheme_0_.append(solution_tmp_, j);
            }

            // Add the remaining jobs to solution_tmp.
            JobPos p0 = (pos_new < pos)? pos_new: pos_new + size;
            for (JobPos p = p0; p < (JobPos)local_scheme_0_.jobs(solution).size(); ++p) {
                // Skip jobs from the previously added bloc.
                if (pos <= p && p < pos + size)
                    continue;
                // Add job to solution_tmp_.
                JobId j = local_scheme_0_.jobs(solution)[p];
                local_scheme_0_.append(solution_tmp_, j);
            }
            global_costs_shift_[pos_new] = global_cost(solution_tmp_);

            // Stop condition.
            if (pos_new == (JobPos)local_scheme_0_.jobs(solution).size() - size)
                break;

            // Add j1 to solution_cur_.
            assert(p0 < (JobPos)local_scheme_0_.jobs(solution).size());
            JobId j1 = local_scheme_0_.jobs(solution)[p0];
            local_scheme_0_.append(solution_cur_, j1);
        }
    }

    inline void compute_cost_swap(const Solution& solution)
    {
        // Initialize solution_cur_.
        Solution solution_cur_ = local_scheme_0_.empty_solution();
        // Reset global_costs_swap_.
        for (JobPos pos_1 = 0; pos_1 < (JobPos)local_scheme_0_.jobs(solution).size(); ++pos_1)
            std::fill(
                    global_costs_swap_[pos_1].begin(),
                    global_costs_swap_[pos_1].end(),
                    local_scheme_0_.global_cost_worst());

        // Loop through all pairs.
        for (JobPos pos_1 = 0; pos_1 < (JobPos)local_scheme_0_.jobs(solution).size(); ++pos_1) {
            JobPos pos_max = std::min(
                    (JobPos)local_scheme_0_.jobs(solution).size(),
                    pos_1 + parameters_.swap_maximum_distance);
            for (JobPos pos_2 = pos_1 + 1; pos_2 < pos_max; ++pos_2) {
                // Initialize solution_tmp_.
                Solution solution_tmp_ = solution_cur_;
                // Add remaining jobs.
                for (JobPos pos = pos_1; pos < (JobPos)local_scheme_0_.jobs(solution).size(); ++pos) {
                    JobId j = local_scheme_0_.jobs(solution)[pos];
                    // If j1 or j2, swap.
                    if (pos == pos_1)
                        j = local_scheme_0_.jobs(solution)[pos_2];
                    if (pos == pos_2)
                        j = local_scheme_0_.jobs(solution)[pos_1];
                    // Add job to solution_tmp_.
                    local_scheme_0_.append(solution_tmp_, j);
                }
                global_costs_swap_[pos_1][pos_2 - pos_1 - 1] = local_scheme_0_.global_cost(solution_tmp_);
            }

            // Add j1 to solution_cur_.
            JobId j1 = local_scheme_0_.jobs(solution)[pos_1];
            local_scheme_0_.append(solution_cur_, j1);
        }
    }

    inline void compute_cost_reverse(const Solution& solution)
    {
        // Initialize solution_cur_.
        Solution solution_cur_ = local_scheme_0_.empty_solution();
        // Reset global_costs_swap_.
        for (JobPos pos_1 = 0; pos_1 < (JobPos)local_scheme_0_.jobs(solution).size(); ++pos_1)
            std::fill(
                    global_costs_swap_[pos_1].begin(),
                    global_costs_swap_[pos_1].end(),
                    local_scheme_0_.global_cost_worst());

        // Loop through all pairs.
        for (JobPos pos_1 = 0; pos_1 < (JobPos)local_scheme_0_.jobs(solution).size(); ++pos_1) {
            JobPos pos_max = std::min(
                    (JobPos)local_scheme_0_.jobs(solution).size(),
                    pos_1 + parameters_.reverse_maximum_length);
            for (JobPos pos_2 = pos_1 + 2; pos_2 < pos_max; ++pos_2) {
                // Initialize solution_tmp_.
                Solution solution_tmp_ = solution_cur_;
                // Add reverse sequence.
                for (JobPos pos = pos_2; pos >= pos_1; --pos) {
                    JobId j = local_scheme_0_.jobs(solution)[pos];
                    // Add job to solution_tmp_.
                    local_scheme_0_.append(solution_tmp_, j);
                }
                // Add remaining jobs.
                for (JobPos pos = pos_2 + 1; pos < (JobPos)local_scheme_0_.jobs(solution).size(); ++pos) {
                    JobId j = local_scheme_0_.jobs(solution)[pos];
                    // Add job to solution_tmp_.
                    local_scheme_0_.append(solution_tmp_, j);
                }
                global_costs_swap_[pos_1][pos_2 - pos_1 - 1] = local_scheme_0_.global_cost(solution_tmp_);
                //std::cout << to_string(global_costs_swap_[pos_1][pos_2 - pos_1 - 1]) << std::endl;
            }

            // Add j1 to solution_cur_.
            JobId j1 = local_scheme_0_.jobs(solution)[pos_1];
            local_scheme_0_.append(solution_cur_, j1);
        }
    }

    /*
     * Private attributes.
     */

    LocalScheme0& local_scheme_0_;
    Parameters parameters_;

    std::vector<JobPos> positions1_;
    std::vector<JobPos> positions2_;
    std::vector<std::pair<JobPos, JobPos>> pairs_;

    std::vector<GlobalCost> global_costs_shift_;
    std::vector<std::vector<GlobalCost>> global_costs_swap_;

    Solution solution_cur_;
    Solution solution_tmp_;

};

}

}

