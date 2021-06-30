#pragma once

#include "localsearchsolver/common.hpp"

#include "external/treesearchsolver/examples/permutationflowshopschedulingmakespan.hpp"

namespace localsearchsolver
{

namespace permutationflowshopschedulingmakespan
{

using JobId = treesearchsolver::permutationflowshopschedulingmakespan::JobId;
using JobPos = treesearchsolver::permutationflowshopschedulingmakespan::JobPos;
using MachineId = treesearchsolver::permutationflowshopschedulingmakespan::MachineId;
using Time = treesearchsolver::permutationflowshopschedulingmakespan::Time;
using Instance = treesearchsolver::permutationflowshopschedulingmakespan::Instance;

class LocalScheme
{

public:

    /** Global cost: <Job number, Makespan>; */
    using GlobalCost = std::tuple<JobId, Time>;

    inline JobId&       job_number(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time&          makespan(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline JobId  job_number(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time     makespan(const GlobalCost& global_cost) { return std::get<1>(global_cost); }

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

    /**
     * solution[n] is the first job scheduled.
     * solution[j] is the index of the job scheduled after job j, n if job
     * j is the last job, or -1 if job j is not in the solution.
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
        Time makespan = 0;
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
        JobPos bloc_size_max = 8;
        bool shuffle_neighborhood_order = true;
    };

    LocalScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters),
        positions1_(instance.job_number()),
        positions2_(instance.job_number()),
        times_(instance_.machine_number(), 0),
        heads_(instance.job_number() + 1),
        tails_(instance.job_number() + 1),
        completion_times_(instance.job_number() + 1)
    {
        std::iota(positions1_.begin(), positions1_.end(), 0);
        std::iota(positions2_.begin(), positions2_.end(), 0);
        for (JobId j = 0; j < instance_.job_number() + 1; ++j) {
            heads_[j] = std::vector<Time>(instance_.machine_number(), 0);
            tails_[j] = std::vector<Time>(instance_.machine_number(), 0);
            completion_times_[j] = std::vector<Time>(instance_.machine_number(), 0);
        }
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
        Solution solution = empty_solution();
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
            solution.makespan,
        };
    }

    /*
     * Local search.
     */

    struct Move
    {
        JobPos pos;
        JobPos pos_new;
        JobId j;
        JobId j_prev;
        GlobalCost global_cost;
    };

    static Move move_null() { return {-1, -1, -1, -1, global_cost_worst()}; }

    struct MoveHasher
    {
        std::hash<JobId> hasher;

        inline bool operator()(
                const Move& move_1,
                const Move& move_2) const
        {
            return move_1.j == move_2.j
                && move_1.j_prev == move_2.j_prev;
        }

        inline std::size_t operator()(
                const Move& move) const
        {
            size_t hash = hasher(move.j);
            optimizationtools::hash_combine(hash, hasher(move.j_prev));
            return hash;
        }
    };

    inline MoveHasher move_hasher() const { return MoveHasher(); }

    inline std::vector<Move> perturbations(const Solution& solution)
    {
        MachineId m = instance_.machine_number();
        std::vector<Move> moves;
        for (JobPos pos = 0; pos < (JobPos)solution.jobs.size(); ++pos) {
            JobId j = solution.jobs[pos];
            compute_structures(solution, pos, 1);
            for (JobPos pos_new = 0; pos_new < (JobPos)solution.jobs.size(); ++pos_new) {
                if (pos_new == pos)
                    continue;
                JobId j_prev = -1;
                if (pos_new > 0) {
                    JobPos p = (pos_new < pos)? pos_new - 1: pos_new;
                    j_prev = solution.jobs[p];
                }
                Time makespan = 0;
                for (MachineId i = 0; i < m; ++i)
                    makespan = std::max(makespan,
                            completion_times_[pos_new][i]
                                     + tails_[pos_new][i]);
                Move move;
                move.pos = pos;
                move.pos_new = pos_new;
                move.j = j;
                move.j_prev = j_prev;
                move.global_cost = {-solution.jobs.size(), makespan};
                moves.push_back(move);
            }
        }
        return moves;
    }

    inline void apply_move(Solution& solution, const Move& move)
    {
        insert(solution, move.pos, 1, move.pos_new);
        assert(solution.makespan == makespan(move.global_cost));
    }

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Move& tabu = move_null())
    {
        MachineId m = instance_.machine_number();
        Counter it = 0;
        std::vector<Counter> neighborhoods;
        for (JobPos bloc_size = 1; bloc_size <= parameters_.bloc_size_max; ++bloc_size)
            neighborhoods.push_back(bloc_size);
        for (;; ++it) {
            //std::cout << "it " << it
            //    << " c " << to_string(global_cost(solution))
            //    << std::endl;
            //print(std::cout, solution);

            if (parameters_.shuffle_neighborhood_order)
                std::shuffle(neighborhoods.begin(), neighborhoods.end(), generator);
            bool improved = false;
            // Loop through neighborhoods.
            for (Counter bloc_size: neighborhoods) {
                std::shuffle(positions1_.begin(), positions1_.end(), generator);
                std::shuffle(positions2_.begin(), positions2_.end(), generator);
                JobPos pos_best = -1;
                JobPos pos_new_best = -1;
                GlobalCost c_best = global_cost(solution);
                for (JobPos pos: positions1_) {
                    if (pos > (JobPos)solution.jobs.size() - bloc_size)
                        continue;
                    // We insert neither tabu.j nor tabu.j_prev.
                    JobId j_cur = solution.jobs[pos];
                    if (j_cur == tabu.j)
                        continue;
                    if (bloc_size == 1 && j_cur == tabu.j_prev)
                        continue;

                    compute_structures(solution, pos, bloc_size);
                    for (JobPos pos_new: positions2_) {
                        if (pos == pos_new || pos_new > (JobPos)solution.jobs.size() - bloc_size)
                            continue;
                        JobId j_prev = -1;
                        if (pos_new > 0) {
                            JobPos p = (pos_new < pos)? pos_new - 1: pos_new + bloc_size - 1;
                            //std::cout
                            //    << "pos " << pos
                            //    << " bloc_size " << bloc_size
                            //    << " pos_new " << pos_new
                            //    << " p " << p
                            //    << std::endl;
                            assert(p >= 0);
                            assert(p < (JobPos)solution.jobs.size());
                            j_prev = solution.jobs[p];
                        }
                        // We don't put anything between tabu.j_prev and tabu.j.
                        if (tabu.j_prev != -1 && j_prev == tabu.j_prev)
                            continue;
                        Time makespan = 0;
                        for (MachineId i = 0; i < m; ++i)
                            makespan = std::max(makespan,
                                    completion_times_[pos_new][i]
                                             + tails_[pos_new][i]);
                        GlobalCost c = {-solution.jobs.size(), makespan};
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
                    assert(makespan(c_best) < solution.makespan);
                    // Apply best move.
                    //std::cout << "pos_best " << pos_best
                    //    << " pos_new_best " << pos_new_best
                    //    << " size " << bloc_size
                    //    << std::endl;
                    insert(solution, pos_best, bloc_size, pos_new_best);
                    if (solution.makespan != makespan(c_best)) {
                        std::cout << "pos_best " << pos_best
                            << " pos_new_best " << pos_new_best
                            << " size " << bloc_size
                            << std::endl;
                        std::cout << makespan(c_best) << std::endl;
                        std::cout << solution.makespan << std::endl;
                        for (MachineId i = 0; i < m; ++i) {
                            std::cout << "i " << i
                                << " " << heads_[((pos_new_best <= pos_best)? pos_new_best: pos_new_best - bloc_size)][i]
                                << " " << completion_times_[((pos_new_best <= pos_best)? pos_new_best: pos_new_best - bloc_size)][i]
                                << " " << tails_[((pos_new_best <= pos_best)? pos_new_best: pos_new_best - bloc_size)][i]
                                << std::endl;
                        }
                        print(std::cout, solution);
                    }
                    assert(solution.makespan == makespan(c_best));
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
        os << "makespan: " << solution.makespan << std::endl;
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
        solution.jobs = jobs;
        std::fill(times_.begin(), times_.end(), 0);
        for (JobId j: solution.jobs) {
            times_[0] = times_[0] + instance_.processing_time(j, 0);
            for (MachineId i = 1; i < m; ++i) {
                if (times_[i - 1] > times_[i]) {
                    times_[i] = times_[i - 1] + instance_.processing_time(j, i);
                } else {
                    times_[i] = times_[i] + instance_.processing_time(j, i);
                }
            }
        }
        solution.makespan = times_[m - 1];
    }

    inline void insert(
            Solution& solution,
            JobPos pos,
            JobPos size,
            JobPos pos_new)
    {
        assert(pos != pos_new);
        assert(pos >= 0);
        assert(pos <= (JobPos)instance_.job_number() - size);
        assert(pos_new <= (JobPos)instance_.job_number() - size);
        assert(pos_new <= (JobPos)instance_.job_number() - size);

        std::vector<JobId> jobs;
        if (pos > pos_new) {
            for (JobPos p = 0; p < pos_new; ++p)
                jobs.push_back(solution.jobs[p]);
            for (JobPos p = pos; p < pos + size; ++p)
                jobs.push_back(solution.jobs[p]);
            for (JobPos p = pos_new; p < pos; ++p)
                jobs.push_back(solution.jobs[p]);
            for (JobPos p = pos + size; p < (JobPos)solution.jobs.size(); ++p)
                jobs.push_back(solution.jobs[p]);
        } else {
            for (JobPos p = 0; p < pos; ++p)
                jobs.push_back(solution.jobs[p]);
            for (JobPos p = pos + size; p < pos_new + size; ++p)
                jobs.push_back(solution.jobs[p]);
            for (JobPos p = pos; p < pos + size; ++p)
                jobs.push_back(solution.jobs[p]);
            for (JobPos p = pos_new + size; p < (JobPos)solution.jobs.size(); ++p)
                jobs.push_back(solution.jobs[p]);
        }
        assert((JobPos)jobs.size() <= instance_.job_number());
        compute(solution, jobs);
    }

    /*
     * Evaluate moves.
     */

    inline void compute_structures(
            const Solution& solution,
            JobPos pos,
            JobPos size)
    {
        MachineId m = instance_.machine_number();

        // Compute heads_.
        for (JobPos pos_new = 0; pos_new < (JobPos)solution.jobs.size() - size; ++pos_new) {
            JobId j = solution.jobs[((pos_new < pos)? pos_new: pos_new + size)];
            heads_[pos_new + 1][0] = heads_[pos_new][0]
                + instance_.processing_time(j, 0);
            for (MachineId i = 1; i < m; ++i) {
                if (heads_[pos_new + 1][i - 1] > heads_[pos_new][i]) {
                    heads_[pos_new + 1][i] = heads_[pos_new + 1][i - 1]
                        + instance_.processing_time(j, i);
                } else {
                    heads_[pos_new + 1][i] = heads_[pos_new][i]
                        + instance_.processing_time(j, i);
                }
            }
        }

        // Compute completion_times_.
        for (JobPos pos_new = 0; pos_new <= (JobPos)solution.jobs.size() - size; ++pos_new) {
            for (MachineId i = 0; i < m; ++i)
                completion_times_[pos_new][i] = heads_[pos_new][i];
            for (JobPos p0 = pos; p0 < pos + size; ++p0) {
                JobId j0 = solution.jobs[p0];
                completion_times_[pos_new][0] = completion_times_[pos_new][0]
                    + instance_.processing_time(j0, 0);
                for (MachineId i = 1; i < m; ++i) {
                    if (completion_times_[pos_new][i] > completion_times_[pos_new][i - 1]) {
                        completion_times_[pos_new][i] = completion_times_[pos_new][i]
                            + instance_.processing_time(j0, i);
                    } else {
                        completion_times_[pos_new][i] = completion_times_[pos_new][i - 1]
                            + instance_.processing_time(j0, i);
                    }
                }
            }
        }

        // Update tails_.
        for (MachineId i = m - 1; i >= 0; --i)
            tails_[solution.jobs.size() - size][i] = 0;
        for (JobPos pos_new = solution.jobs.size() - size - 1; pos_new >= 0; --pos_new) {
            JobId j = solution.jobs[((pos_new < pos)? pos_new: pos_new + size)];
            assert(j >= 0);
            assert(j < instance_.job_number());
            tails_[pos_new][m - 1] = tails_[pos_new + 1][m - 1]
                + instance_.processing_time(j, m - 1);
            for (MachineId i = m - 2; i >= 0; --i) {
                if (tails_[pos_new][i + 1] > tails_[pos_new + 1][i]) {
                    tails_[pos_new][i] = tails_[pos_new][i + 1]
                        + instance_.processing_time(j, i);
                } else {
                    tails_[pos_new][i] = tails_[pos_new + 1][i]
                        + instance_.processing_time(j, i);
                }
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
    std::vector<Time> times_;

    std::vector<std::vector<Time>> heads_;
    std::vector<std::vector<Time>> tails_;
    std::vector<std::vector<Time>> completion_times_;

};

}

}

