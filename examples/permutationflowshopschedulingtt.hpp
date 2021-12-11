#pragma once

/**
 * Permutation flow shop scheduling problem, Total tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/permutationflowshopschedulingtt.hpp
 *
 */

#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing2.hpp"

#include "orproblems/permutationflowshopschedulingtt.hpp"

namespace localsearchsolver
{

namespace permutationflowshopschedulingtt
{

using namespace orproblems::permutationflowshopschedulingtt;

class LocalScheme
{

public:

    /** Global cost: <Total tardiness>; */
    using GlobalCost = std::tuple<Time>;

    inline Time&       total_tardiness(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Time  total_tardiness(const GlobalCost& global_cost) { return std::get<0>(global_cost); }

    /*
     * Solutions.
     */

    struct Solution
    {
        std::vector<JobId> jobs;
        std::vector<Time> times;
        Time total_tardiness = 0;
    };

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shift_bloc_maximum_length = 7;
            sequencing_parameters.swap_bloc_maximum_length = 2;
            sequencing_parameters.reverse = true;
            sequencing_parameters.shift_reverse_bloc_maximum_length = 4;
        }

        sequencing2::Parameters sequencing_parameters;
    };

    LocalScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters) { }

    LocalScheme(const LocalScheme& local_scheme):
        LocalScheme(local_scheme.instance_, local_scheme.parameters_) { }

    virtual ~LocalScheme() { }

    /*
     * Initial solutions.
     */

    inline Solution empty_solution() const
    {
        Solution solution;
        solution.times = std::vector<Time>(instance_.number_of_machines(), 0);
        solution.total_tardiness = std::numeric_limits<Time>::max();
        return solution;
    }

    inline Solution initial_solution(
            Counter,
            std::mt19937_64& generator)
    {
        std::vector<JobId> jobs(instance_.number_of_jobs());
        std::iota(jobs.begin(), jobs.end(), 0);
        std::shuffle(jobs.begin(), jobs.end(), generator);
        Solution solution = empty_solution();
        for (JobId j: jobs)
            append(solution, j);
        return solution;
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            solution.total_tardiness,
        };
    }

    /*
     * Methods required by sequencing2::LocalScheme.
     */

    inline JobPos number_of_jobs() const { return instance_.number_of_jobs(); }

    inline const std::vector<JobId>& jobs(const Solution& solution) const { return solution.jobs; }

    inline void append(
            Solution& solution,
            JobId j) const
    {
        if (solution.jobs.size() == 0)
            solution.total_tardiness = 0;
        MachineId m = instance_.number_of_machines();
        // Update jobs.
        solution.jobs.push_back(j);
        // Update times.
        solution.times[0] = solution.times[0] + instance_.job(j).processing_times[0];
        for (MachineId i = 1; i < m; ++i) {
            if (solution.times[i - 1] > solution.times[i]) {
                solution.times[i] = solution.times[i - 1] + instance_.job(j).processing_times[i];
            } else {
                solution.times[i] = solution.times[i] + instance_.job(j).processing_times[i];
            }
        }
        // Update total tardiness.
        if (solution.times[m - 1] > instance_.job(j).due_date)
            solution.total_tardiness += (solution.times[m - 1] - instance_.job(j).due_date);
    }

    /*
     * Outputs.
     */

    std::ostream& print(
            std::ostream &os,
            const Solution& solution) const
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
     * Private attributes.
     */

    const Instance& instance_;
    Parameters parameters_;

};

}

}

