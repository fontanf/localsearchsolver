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
        JobPos bloc_size_max = 13;
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
        total_weighted_tardinesses_shift_(instance.job_number() + 1, 0),
        total_weighted_tardinesses_swap_(instance.job_number())
    {
        std::iota(positions1_.begin(), positions1_.end(), 0);
        std::iota(positions2_.begin(), positions2_.end(), 0);
        for (JobPos pos_1 = 0; pos_1 < instance_.job_number(); ++pos_1)
            for (JobPos pos_2 = pos_1 + 1; pos_2 < instance_.job_number(); ++pos_2)
                pairs_.push_back({pos_1, pos_2});
        for (JobPos pos_1 = 0; pos_1 < instance_.job_number(); ++pos_1)
            total_weighted_tardinesses_swap_[pos_1].resize(instance_.job_number() - pos_1, 0);
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
        JobPos pos_1;
        JobPos pos_2;
        JobPos pos_3;
        JobPos pos_4;
        GlobalCost global_cost;
    };

    static Move move_null() { return {-1, -1, -1, -1, global_cost_worst()}; }

    struct MoveHasher
    {
        std::hash<JobPos> hasher;

        inline bool operator()(const Move&, const Move&) const { return false; }

        inline std::size_t operator()(const Move& move) const
        {
            size_t hash = hasher(move.pos_1);
            optimizationtools::hash_combine(hash, hasher(move.pos_2));
            optimizationtools::hash_combine(hash, hasher(move.pos_3));
            optimizationtools::hash_combine(hash, hasher(move.pos_4));
            return hash;
        }
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
        //if (tabu.pos_1 != -1)
        //    std::cout << tabu.pos_1  << " " << tabu.pos_2 << " " << tabu.pos_3 << " " << tabu.pos_4 << std::endl;
        //print(std::cout, solution);
        //std::cout << to_string(global_cost(solution)) << std::endl;

        Counter it = 0;
        std::vector<Counter> neighborhoods;
        if (parameters_.swap)
            neighborhoods.push_back(0);
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
                            total_weighted_tardinesses_swap_[pair.first][pair.second - pair.first - 1]};
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
                            GlobalCost c = {
                                -solution.jobs.size(),
                                total_weighted_tardinesses_shift_[pos_new]};
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
                        compute(solution, jobs);

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

    /*
     * Evaluate moves.
     */

    inline void compute_cost_shift(
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
            total_weighted_tardinesses_shift_[pos_new] = total_tardiness_cur;
            Time t = t_cur;
            JobId j_prev = j_prev_cur;
            for (JobPos p = pos; p < pos + size; ++p) {
                JobId j = solution.jobs[p];
                t += instance_.setup_time(j_prev, j);
                t += instance_.job(j).processing_time;
                if (t > instance_.job(j).due_date)
                    total_weighted_tardinesses_shift_[pos_new]
                        += instance_.job(j).weight
                        * (t - instance_.job(j).due_date);
                j_prev = j;
            }

            // Add remaining jobs to times_.
            JobPos p0 = (pos_new < pos)? pos_new: pos_new + size;
            for (JobPos p = p0; p < (JobPos)solution.jobs.size(); ++p) {
                if (total_weighted_tardinesses_shift_[pos_new] > solution.total_weighted_tardiness)
                    break;
                if (pos <= p && p < pos + size)
                    continue;
                JobId j = solution.jobs[p];
                t += instance_.setup_time(j_prev, j);
                t += instance_.job(j).processing_time;
                if (t > instance_.job(j).due_date)
                    total_weighted_tardinesses_shift_[pos_new]
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

    inline void compute_cost_swap(const Solution& solution)
    {
        JobId n = instance_.job_number();
        Time total_tardiness_cur = 0;
        Time t_cur = 0;
        JobId j_prev_cur = n;
        for (JobPos pos_1 = 0; pos_1 < (JobPos)solution.jobs.size(); ++pos_1) {
            for (JobPos pos_2 = pos_1 + 1; pos_2 < (JobPos)solution.jobs.size(); ++pos_2) {
                // Add j2 ... j1 ...
                Time t = t_cur;
                JobId j_prev = j_prev_cur;
                total_weighted_tardinesses_swap_[pos_1][pos_2 - pos_1 - 1] = total_tardiness_cur;
                for (JobPos pos = pos_1; pos < (JobPos)solution.jobs.size(); ++pos) {
                    JobId j = solution.jobs[pos];
                    if (pos == pos_1)
                        j = solution.jobs[pos_2];
                    if (pos == pos_2)
                        j = solution.jobs[pos_1];
                    t += instance_.setup_time(j_prev, j);
                    t += instance_.job(j).processing_time;
                    if (t > instance_.job(j).due_date)
                        total_weighted_tardinesses_swap_[pos_1][pos_2 - pos_1 - 1]
                            += instance_.job(j).weight
                            * (t - instance_.job(j).due_date);
                    j_prev = j;
                }
            }
            // Add j1.
            JobId j1 = solution.jobs[pos_1];
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
    std::vector<std::pair<JobPos, JobPos>> pairs_;

    std::vector<Weight> total_weighted_tardinesses_shift_;
    std::vector<std::vector<Weight>> total_weighted_tardinesses_swap_;

};

}

}

