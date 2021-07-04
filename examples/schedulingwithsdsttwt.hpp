#pragma once

#include "localsearchsolver/common.hpp"

#include "external/treesearchsolver/examples/schedulingwithsdsttwt.hpp"

namespace localsearchsolver
{

namespace schedulingwithsdsttwt
{

using JobId = treesearchsolver::schedulingwithsdsttwt::JobId;
using JobPos = treesearchsolver::schedulingwithsdsttwt::JobPos;
using Time = treesearchsolver::schedulingwithsdsttwt::Time;
using Weight = treesearchsolver::schedulingwithsdsttwt::Weight;
using Instance = treesearchsolver::schedulingwithsdsttwt::Instance;

class LocalScheme
{

public:

    /** Global cost: <Job number, Total weighted tardiness>; */
    using GlobalCost = std::tuple<JobId, Weight>;

    inline JobId&                      job_number(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Weight&       total_weighted_tardiness(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline JobId                 job_number(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Weight  total_weighted_tardiness(const GlobalCost& global_cost) { return std::get<1>(global_cost); }

    static GlobalCost global_cost_worst()
    {
        return {
            std::numeric_limits<JobId>::max(),
            std::numeric_limits<Weight>::max(),
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
        Time total_weighted_tardiness = 0;
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
        JobPos bloc_size_max = 4;
        bool shuffle_neighborhood_order = true;
    };

    LocalScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters),
        positions1_(instance.job_number()),
        positions2_(instance.job_number()),
        total_tardinesses_(instance.job_number() + 1, 0)
    {
        std::iota(positions1_.begin(), positions1_.end(), 0);
        std::iota(positions2_.begin(), positions2_.end(), 0);
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
            solution.total_weighted_tardiness,
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
        std::vector<Move> moves;
        for (JobPos pos = 0; pos < instance_.job_number(); ++pos) {
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
                //if (total_tardinesses_[pos] == solution.total_tardiness)
                //    continue;
                Move move;
                move.pos = pos;
                move.pos_new = pos_new;
                move.j = j;
                move.j_prev = j_prev;
                move.global_cost = {-solution.jobs.size(), total_tardinesses_[pos_new]};
                moves.push_back(move);
            }
        }
        return moves;
    }

    inline void apply_move(Solution& solution, const Move& move)
    {
        insert(solution, move.pos, 1, move.pos_new);
        if (solution.total_weighted_tardiness != total_weighted_tardiness(move.global_cost)) {
            std::cout << "j " << move.j
                << " j_prev " << move.j_prev
                << std::endl;
            std::cout << total_weighted_tardiness(move.global_cost) << std::endl;
            std::cout << solution.total_weighted_tardiness << std::endl;
        }
        assert(solution.total_weighted_tardiness == total_weighted_tardiness(move.global_cost));
    }

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Move& tabu = move_null())
    {
        //if (tabu.j != -1)
        //    std::cout << "j " << tabu.j << " j_prev " << tabu.j_prev
        //        << std::endl;
        //print(std::cout, solution);
        //std::cout << to_string(global_cost(solution)) << std::endl;

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
                            j_prev = solution.jobs[p];
                        }
                        // We don't put anything between tabu.j_prev and tabu.j.
                        if (tabu.j_prev != -1 && j_prev == tabu.j_prev)
                            continue;
                        GlobalCost c = {-solution.jobs.size(), total_tardinesses_[pos_new]};
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
                    assert(total_weighted_tardiness(c_best) < solution.total_weighted_tardiness);
                    // Apply best move.
                    //std::cout << bloc_size << std::endl;
                    insert(solution, pos_best, bloc_size, pos_new_best);
                    if (solution.total_weighted_tardiness != total_weighted_tardiness(c_best)) {
                        std::cout << "pos_best " << pos_best
                            << " pos_new_best " << pos_new_best
                            << " size " << bloc_size
                            << std::endl;
                        std::cout << total_weighted_tardiness(c_best) << std::endl;
                        std::cout << solution.total_weighted_tardiness << std::endl;
                        print(std::cout, solution);
                    }
                    assert(solution.total_weighted_tardiness == total_weighted_tardiness(c_best));
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
        os << "total weighted tardiness: " << solution.total_weighted_tardiness << std::endl;
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
        JobId n = instance_.job_number();
        solution.jobs = jobs;
        solution.total_weighted_tardiness = 0;
        Time t = 0;
        JobId j_prev = n;
        for (JobId j: solution.jobs) {
            t += instance_.setup_time(j_prev, j);
            t += instance_.job(j).processing_time;
            if (t > instance_.job(j).due_date)
                solution.total_weighted_tardiness
                    += instance_.job(j).weight
                    * (t - instance_.job(j).due_date);
            j_prev = j;
        }
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
        JobId n = instance_.job_number();
        Time total_tardiness_cur = 0;
        Time t_cur = 0;
        JobId j_prev_cur = n;
        for (JobPos pos_new = 0; pos_new <= (JobPos)solution.jobs.size() - size; ++pos_new) {
            // Add bloc.
            total_tardinesses_[pos_new] = total_tardiness_cur;
            Time t = t_cur;
            JobId j_prev = j_prev_cur;
            for (JobPos p = pos; p < pos + size; ++p) {
                JobId j = solution.jobs[p];
                t += instance_.setup_time(j_prev, j);
                t += instance_.job(j).processing_time;
                if (t > instance_.job(j).due_date)
                    total_tardinesses_[pos_new]
                        += instance_.job(j).weight
                        * (t - instance_.job(j).due_date);
                j_prev = j;
            }

            // Add remaining jobs to times_.
            JobPos p0 = (pos_new < pos)? pos_new: pos_new + size;
            for (JobPos p = p0; p < (JobPos)solution.jobs.size(); ++p) {
                if (pos <= p && p < pos + size)
                    continue;
                JobId j = solution.jobs[p];
                t += instance_.setup_time(j_prev, j);
                t += instance_.job(j).processing_time;
                if (t > instance_.job(j).due_date)
                    total_tardinesses_[pos_new]
                        += instance_.job(j).weight
                        * (t - instance_.job(j).due_date);
                j_prev = j;
            }

            // Add j1 to heads_.
            if (pos_new == (JobPos)solution.jobs.size() - size)
                break;
            assert(p0 < (JobPos)solution.jobs.size());
            JobId j1 = solution.jobs[p0];
            t_cur += instance_.setup_time(j_prev_cur, j1);
            t_cur += instance_.job(j1).processing_time;
            if (t_cur > instance_.job(j1).due_date)
                total_tardiness_cur
                    += instance_.job(j1).weight
                    * (t_cur - instance_.job(j1).due_date);
            j_prev_cur = j1;

        }
    }

    /*
     * Private attributes.
     */

    const Instance& instance_;
    Parameters parameters_;

    std::vector<JobPos> positions1_;
    std::vector<JobPos> positions2_;

    std::vector<Time> total_tardinesses_;

};

}

}

