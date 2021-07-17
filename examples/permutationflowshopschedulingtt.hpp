#pragma once

#include "localsearchsolver/common.hpp"

#include "external/treesearchsolver/examples/permutationflowshopschedulingtt.hpp"

namespace localsearchsolver
{

namespace permutationflowshopschedulingtt
{

using JobId = treesearchsolver::permutationflowshopschedulingtt::JobId;
using JobPos = treesearchsolver::permutationflowshopschedulingtt::JobPos;
using MachineId = treesearchsolver::permutationflowshopschedulingtt::MachineId;
using Time = treesearchsolver::permutationflowshopschedulingtt::Time;
using Instance = treesearchsolver::permutationflowshopschedulingtt::Instance;

class LocalScheme
{

public:

    /** Global cost: <Job number, Total tardiness>; */
    using GlobalCost = std::tuple<JobId, Time>;

    inline JobId&           job_number(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time&       total_tardiness(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline JobId      job_number(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time  total_tardiness(const GlobalCost& global_cost) { return std::get<1>(global_cost); }

    static GlobalCost global_cost_worst()
    {
        return {
            std::numeric_limits<JobId>::max(),
            std::numeric_limits<Time>::max(),
        };
    }

    /*
     * Solutions.
     */

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

    struct Solution
    {
        std::vector<JobId> jobs;
        Time total_tardiness = 0;
    };

    CompactSolution solution2compact(const Solution& solution)
    {
        return solution.jobs;
    }

    Solution compact2solution(const CompactSolution& compact_solution)
    {
        auto solution = empty_solution();
        compute(solution, compact_solution);
        return solution;
    }

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        JobPos bloc_size_max = 3;
        bool swap = true;
        bool shuffle_neighborhood_order = true;
        Counter perturbation_number = 10;
    };

    LocalScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters),
        positions1_(instance.job_number()),
        positions2_(instance.job_number()),
        times_(instance_.machine_number(), 0),
        heads_(instance.machine_number(), 0),
        total_tardinesses_shift_(instance.job_number() + 1, 0),
        total_tardinesses_swap_(instance.job_number())
    {
        std::iota(positions1_.begin(), positions1_.end(), 0);
        std::iota(positions2_.begin(), positions2_.end(), 0);
        for (JobPos pos_1 = 0; pos_1 < instance_.job_number(); ++pos_1)
            for (JobPos pos_2 = pos_1 + 1; pos_2 < instance_.job_number(); ++pos_2)
                pairs_.push_back({pos_1, pos_2});
        for (JobPos pos_1 = 0; pos_1 < instance_.job_number(); ++pos_1)
            total_tardinesses_swap_[pos_1].resize(instance_.job_number() - pos_1, 0);
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
        return solution;
    }

    inline Solution initial_solution(
            Counter,
            std::mt19937_64& generator)
    {
        std::vector<JobId> jobs(instance_.job_number());
        std::iota(jobs.begin(), jobs.end(), 0);
        std::shuffle(jobs.begin(), jobs.end(), generator);
        return compact2solution(jobs);
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            -solution.jobs.size(),
            solution.total_tardiness,
        };
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

    static Move move_null() { return {-1, -1, -1, -1, global_cost_worst()}; }

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
        for (Counter perturbation = 0; perturbation < parameters_.perturbation_number; ++perturbation) {
            std::vector<JobPos> edges = optimizationtools::bob_floyd<JobPos>(
                    4, solution.jobs.size() + 1, generator);
            std::sort(edges.begin(), edges.end());
            Move move;
            move.pos_1 = edges[0];
            move.pos_2 = edges[1];
            move.pos_3 = edges[2];
            move.pos_4 = edges[3];
            assert(move.pos_1 >= 0);
            assert(move.pos_4 <= (JobPos)solution.jobs.size());
            move.global_cost = global_cost(solution);
            moves.push_back(move);
        }
        return moves;
    }

    inline void apply_move(Solution& solution, const Move& move)
    {
        std::vector<JobId> jobs;
        for (JobPos pos = 0; pos < move.pos_1; ++pos)
            jobs.push_back(solution.jobs[pos]);
        for (JobPos pos = move.pos_3; pos < move.pos_4; ++pos)
            jobs.push_back(solution.jobs[pos]);
        for (JobPos pos = move.pos_2; pos < move.pos_3; ++pos)
            jobs.push_back(solution.jobs[pos]);
        for (JobPos pos = move.pos_1; pos < move.pos_2; ++pos)
            jobs.push_back(solution.jobs[pos]);
        for (JobPos pos = move.pos_4; pos < (JobPos)solution.jobs.size(); ++pos)
            jobs.push_back(solution.jobs[pos]);
        assert((JobPos)jobs.size() <= instance_.job_number());
        compute(solution, jobs);
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
        for (JobPos bloc_size = 1; bloc_size <= parameters_.bloc_size_max; ++bloc_size)
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
                        GlobalCost c = {
                            -solution.jobs.size(),
                            total_tardinesses_swap_[pair.first][pair.second - pair.first - 1]};
                        if (c >= c_best)
                            continue;
                        if (pos_1_best != -1 && !dominates(c, c_best))
                            continue;
                        pos_1_best = pair.first;
                        pos_2_best = pair.second;
                        c_best = c;
                    }
                    if (pos_1_best != -1) {
                        improved = true;
                        // Apply best move.
                        std::vector<JobId> jobs;
                        for (JobPos pos = 0; pos < instance_.job_number(); ++pos) {
                            if (pos == pos_1_best) {
                                jobs.push_back(solution.jobs[pos_2_best]);
                            } else if (pos == pos_2_best) {
                                jobs.push_back(solution.jobs[pos_1_best]);
                            } else {
                                jobs.push_back(solution.jobs[pos]);
                            }
                        }
                        compute(solution, jobs);
                        assert(solution.total_tardiness == total_tardiness(c_best));
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
                        if (pos > (JobPos)solution.jobs.size() - bloc_size)
                            continue;
                        compute_cost_shift(solution, pos, bloc_size);
                        for (JobPos pos_new: positions2_) {
                            if (pos == pos_new || pos_new > (JobPos)solution.jobs.size() - bloc_size)
                                continue;
                            GlobalCost c = {-solution.jobs.size(), total_tardinesses_shift_[pos_new]};
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
                        improved = true;
                        assert(total_tardiness(c_best) < solution.total_tardiness);
                        // Apply best move.
                        //std::cout << bloc_size << std::endl;
                        std::vector<JobId> jobs;
                        if (pos_best > pos_new_best) {
                            for (JobPos p = 0; p < pos_new_best; ++p)
                                jobs.push_back(solution.jobs[p]);
                            for (JobPos p = pos_best; p < pos_best + bloc_size; ++p)
                                jobs.push_back(solution.jobs[p]);
                            for (JobPos p = pos_new_best; p < pos_best; ++p)
                                jobs.push_back(solution.jobs[p]);
                            for (JobPos p = pos_best + bloc_size; p < (JobPos)solution.jobs.size(); ++p)
                                jobs.push_back(solution.jobs[p]);
                        } else {
                            for (JobPos p = 0; p < pos_best; ++p)
                                jobs.push_back(solution.jobs[p]);
                            for (JobPos p = pos_best + bloc_size; p < pos_new_best + bloc_size; ++p)
                                jobs.push_back(solution.jobs[p]);
                            for (JobPos p = pos_best; p < pos_best + bloc_size; ++p)
                                jobs.push_back(solution.jobs[p]);
                            for (JobPos p = pos_new_best + bloc_size; p < (JobPos)solution.jobs.size(); ++p)
                                jobs.push_back(solution.jobs[p]);
                        }
                        assert((JobPos)jobs.size() <= instance_.job_number());
                        compute(solution, jobs);
                        if (solution.total_tardiness != total_tardiness(c_best)) {
                            std::cout << "pos_best " << pos_best
                                << " pos_new_best " << pos_new_best
                                << " size " << bloc_size
                                << std::endl;
                            std::cout << total_tardiness(c_best) << std::endl;
                            std::cout << solution.total_tardiness << std::endl;
                            print(std::cout, solution);
                        }
                        assert(solution.total_tardiness == total_tardiness(c_best));
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
            const Solution& solution)
    {
        os << "jobs:";
        for (JobId j: solution.jobs)
            os << " " << j;
        os << std::endl;
        os << "total tardiness: " << solution.total_tardiness << std::endl;
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

        for (JobId j: solution.jobs)
            cert << j << " ";
    }

private:

    /*
     * Manipulate solutions.
     */

    inline void compute(
            Solution& solution,
            const std::vector<JobId>& jobs)
    {
        MachineId m = instance_.machine_number();
        // Update jobs.
        solution.jobs = jobs;
        // Initialize times_ and total tardiness.
        solution.total_tardiness = 0;
        std::fill(times_.begin(), times_.end(), 0);
        // Loop through all jobs.
        for (JobId j: solution.jobs) {
            // Update times_.
            times_[0] = times_[0] + instance_.job(j).processing_times[0];
            for (MachineId i = 1; i < m; ++i) {
                if (times_[i - 1] > times_[i]) {
                    times_[i] = times_[i - 1] + instance_.job(j).processing_times[i];
                } else {
                    times_[i] = times_[i] + instance_.job(j).processing_times[i];
                }
            }
            // Update total tardiness.
            if (times_[m - 1] > instance_.job(j).due_date)
                solution.total_tardiness += (times_[m - 1] - instance_.job(j).due_date);
        }
    }

    /*
     * Evaluate moves.
     */

    inline void compute_cost_shift(
            const Solution& solution,
            JobPos pos,
            JobPos size)
    {
        MachineId m = instance_.machine_number();

        // Initialize heads_ and head total tardiness.
        Time total_tardiness_cur = 0;
        std::fill(heads_.begin(), heads_.end(), 0);
        // Initialize total_tardinesses_shift_ high enough in case of early
        // termination.
        std::fill(
                total_tardinesses_shift_.begin(),
                total_tardinesses_shift_.end(),
                solution.total_tardiness + 1);

        // Loop through all new positions.
        for (JobPos pos_new = 0; pos_new <= (JobPos)solution.jobs.size() - size; ++pos_new) {
            // Initialize times_ and total tardiness.
            total_tardinesses_shift_[pos_new] = total_tardiness_cur;
            times_ = heads_;

            // Add bloc to times_.
            for (JobPos p = pos; p < pos + size; ++p) {
                // Check for early termination.
                if (total_tardinesses_shift_[pos_new] >= solution.total_tardiness)
                    break;
                // Update times_.
                JobId j = solution.jobs[p];
                times_[0] = times_[0] + instance_.job(j).processing_times[0];
                for (MachineId i = 1; i < m; ++i) {
                    if (times_[i - 1] > times_[i]) {
                        times_[i] = times_[i - 1] + instance_.job(j).processing_times[i];
                    } else {
                        times_[i] = times_[i] + instance_.job(j).processing_times[i];
                    }
                }
                // Update move total tardiness.
                if (times_[m - 1] > instance_.job(j).due_date)
                    total_tardinesses_shift_[pos_new] += (times_[m - 1] - instance_.job(j).due_date);
            }

            // Add the remaining jobs to times_.
            JobPos p0 = (pos_new < pos)? pos_new: pos_new + size;
            for (JobPos p = p0; p < (JobPos)solution.jobs.size(); ++p) {
                // Check for early termination.
                if (total_tardinesses_shift_[pos_new] >= solution.total_tardiness)
                    break;
                // Skip jobs from the previously added bloc.
                if (pos <= p && p < pos + size)
                    continue;
                // Update times_.
                JobId j = solution.jobs[p];
                times_[0] = times_[0] + instance_.job(j).processing_times[0];
                for (MachineId i = 1; i < m; ++i) {
                    if (times_[i - 1] > times_[i]) {
                        times_[i] = times_[i - 1] + instance_.job(j).processing_times[i];
                    } else {
                        times_[i] = times_[i] + instance_.job(j).processing_times[i];
                    }
                }
                // Update move total tardiness.
                if (times_[m - 1] > instance_.job(j).due_date)
                    total_tardinesses_shift_[pos_new] += (times_[m - 1] - instance_.job(j).due_date);
            }

            // Stop condition.
            if (pos_new == (JobPos)solution.jobs.size() - size)
                break;

            // Add j1 to heads_.
            assert(p0 < (JobPos)solution.jobs.size());
            JobId j1 = solution.jobs[p0];
            heads_[0] = heads_[0] + instance_.job(j1).processing_times[0];
            for (MachineId i = 1; i < m; ++i) {
                if (heads_[i - 1] > heads_[i]) {
                    heads_[i] = heads_[i - 1] + instance_.job(j1).processing_times[i];
                } else {
                    heads_[i] = heads_[i] + instance_.job(j1).processing_times[i];
                }
            }
            // Update head total tardiness.
            if (heads_[m - 1] > instance_.job(j1).due_date) {
                total_tardiness_cur += (heads_[m - 1] - instance_.job(j1).due_date);
                // Check for early termination.
                if (total_tardiness_cur >= solution.total_tardiness)
                    break;
            }

        }
    }

    inline void compute_cost_swap(const Solution& solution)
    {
        MachineId m = instance_.machine_number();
        // Initialize heads_ and head total tardiness.
        Time total_tardiness_cur = 0;
        std::fill(heads_.begin(), heads_.end(), 0);
        // Initialize total_tardinesses_swap_ high enough in case of early
        // termination.
        for (JobPos pos_1 = 0; pos_1 < (JobPos)solution.jobs.size(); ++pos_1)
            std::fill(
                    total_tardinesses_swap_[pos_1].begin(),
                    total_tardinesses_swap_[pos_1].end(),
                    solution.total_tardiness + 1);
        // Loop through all pairs.
        for (JobPos pos_1 = 0; pos_1 < (JobPos)solution.jobs.size(); ++pos_1) {
            for (JobPos pos_2 = pos_1 + 1; pos_2 < (JobPos)solution.jobs.size(); ++pos_2) {
                // Initialize times_ and total tardiness.
                times_ = heads_;
                total_tardinesses_swap_[pos_1][pos_2 - pos_1 - 1] = total_tardiness_cur;
                // Add remaining jobs.
                for (JobPos pos = pos_1; pos < (JobPos)solution.jobs.size(); ++pos) {
                    // Check for early termination.
                    if (total_tardinesses_swap_[pos_1][pos_2 - pos_1 - 1]
                            >= solution.total_tardiness)
                        break;
                    JobId j = solution.jobs[pos];
                    // If j1 or j2, swap.
                    if (pos == pos_1)
                        j = solution.jobs[pos_2];
                    if (pos == pos_2)
                        j = solution.jobs[pos_1];
                    // Update times_.
                    times_[0] = times_[0] + instance_.job(j).processing_times[0];
                    for (MachineId i = 1; i < m; ++i) {
                        if (times_[i - 1] > times_[i]) {
                            times_[i] = times_[i - 1] + instance_.job(j).processing_times[i];
                        } else {
                            times_[i] = times_[i] + instance_.job(j).processing_times[i];
                        }
                    }
                    // Update move total tardiness.
                    if (times_[m - 1] > instance_.job(j).due_date)
                        total_tardinesses_swap_[pos_1][pos_2 - pos_1 - 1]
                            += (times_[m - 1] - instance_.job(j).due_date);
                }
            }
            // Add job j1 to heads_.
            JobId j1 = solution.jobs[pos_1];
            heads_[0] = heads_[0] + instance_.job(j1).processing_times[0];
            for (MachineId i = 1; i < m; ++i) {
                if (heads_[i - 1] > heads_[i]) {
                    heads_[i] = heads_[i - 1] + instance_.job(j1).processing_times[i];
                } else {
                    heads_[i] = heads_[i] + instance_.job(j1).processing_times[i];
                }
            }
            // Update head total tardiness.
            if (heads_[m - 1] > instance_.job(j1).due_date) {
                total_tardiness_cur += (heads_[m - 1] - instance_.job(j1).due_date);
                // Check for early termination.
                if (total_tardiness_cur >= solution.total_tardiness)
                    break;
            }
        }
    }

    /*
     * Private attributes.
     */

    const Instance& instance_;
    Parameters parameters_;

    std::vector<JobPos> positions1_;
    std::vector<JobPos> positions2_;
    std::vector<std::pair<JobPos, JobPos>> pairs_;
    std::vector<Time> times_;

    std::vector<Time> heads_;
    std::vector<Time> total_tardinesses_shift_;
    std::vector<std::vector<Time>> total_tardinesses_swap_;

};

}

}

