#pragma once

/**
 * Single machine scheduling problem with sequence-dependent setup times, Total
 * weighted Tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/schedulingwithsdsttwt.hpp
 *
 */


#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing.hpp"

#include "orproblems/schedulingwithsdsttwt.hpp"

namespace localsearchsolver
{

namespace schedulingwithsdsttwt
{

using namespace orproblems::schedulingwithsdsttwt;

class LocalScheme
{

public:

    /** Global cost: <Total weighted tardiness>; */
    using GlobalCost = std::tuple<Weight>;

    inline Weight&       total_weighted_tardiness(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Weight  total_weighted_tardiness(const GlobalCost& global_cost) { return std::get<0>(global_cost); }

    static GlobalCost global_cost_worst()
    {
        return {
            std::numeric_limits<Weight>::max(),
        };
    }

    /*
     * Solutions.
     */

    struct Solution
    {
        std::vector<JobId> jobs;
        Time time = 0;
        Weight total_weighted_tardiness = 0;
    };

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shift_bloc_maximum_length = 13;
            sequencing_parameters.swap = true;
            sequencing_parameters.shuffle_neighborhood_order = true;
            sequencing_parameters.number_of_perturbations = 10;
        }

        sequencing::Parameters sequencing_parameters;
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
        solution.total_weighted_tardiness = std::numeric_limits<Weight>::max();
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
            solution.total_weighted_tardiness,
        };
    }

    /*
     * Methods required by sequencing::LocalScheme.
     */

    inline JobPos number_of_jobs() const { return instance_.number_of_jobs(); }

    inline const std::vector<JobId>& jobs(const Solution& solution) const { return solution.jobs; }

    inline void append(
            Solution& solution,
            JobId j) const
    {
        if (solution.jobs.size() == 0)
            solution.total_weighted_tardiness = 0;
        // Update time.
        JobId j_prev = (solution.jobs.size() > 0)?
            solution.jobs.back():
            instance_.number_of_jobs();
        solution.time += instance_.setup_time(j_prev, j);
        solution.time += instance_.job(j).processing_time;
        // Update jobs.
        solution.jobs.push_back(j);
        // Update total weighted tardiness.
        if (solution.time > instance_.job(j).due_date)
            solution.total_weighted_tardiness
                += instance_.job(j).weight
                * (solution.time - instance_.job(j).due_date);
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
     * Private attributes.
     */

    const Instance& instance_;
    Parameters parameters_;

};

}

}

