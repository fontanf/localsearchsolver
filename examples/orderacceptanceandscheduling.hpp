#pragma once

/**
 * Single machine order acceptance and scheduling problem with time windows and
 * sequence-dependent setup times, Total weighted tardiness.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/orderacceptanceandscheduling.hpp
 *
 */


#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing2.hpp"

#include "orproblems/orderacceptanceandscheduling.hpp"

namespace localsearchsolver
{

namespace orderacceptanceandscheduling
{

using namespace orproblems::orderacceptanceandscheduling;

class LocalScheme
{

public:

    /** Global cost: <Reversed time, Total weighted tardiness - Profit>; */
    using GlobalCost = std::tuple<Time, Weight>;

    inline Time&        reverse_time(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Weight&         objective(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline Time  reversed_time(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Weight    objective(const GlobalCost& global_cost) { return std::get<1>(global_cost); }

    /*
     * Solutions.
     */

    struct Solution
    {
        std::vector<JobId> sequence;
        Time time = 0;
        Time reversed_time_curr = 0;
        Time reversed_time_full = 0;
        Weight total_weighted_tardiness_curr = 0;
        Weight total_weighted_tardiness_full = 0;
        Profit profit = 0;
    };

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shift_block_maximum_length = 7;
            sequencing_parameters.swap_block_maximum_length = 5;
            sequencing_parameters.reverse = true;
            sequencing_parameters.shift_reverse_block_maximum_length = 6;
            sequencing_parameters.add_remove = true;

            sequencing_parameters.double_bridge_number_of_perturbations = 0;
            sequencing_parameters.ruin_and_recreate_number_of_perturbations = 0;
            sequencing_parameters.force_add = true;
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
        solution.profit += instance_.job(0).profit;
        solution.profit += instance_.job(instance_.number_of_jobs() - 1).profit;
        return solution;
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            solution.reversed_time_full,
            solution.total_weighted_tardiness_full - solution.profit,
        };
    }

    /*
     * Methods required by sequencing::LocalScheme.
     */

    inline JobPos number_of_elements() const { return instance_.number_of_jobs() - 2; }

    inline GlobalCost bound(const Solution& solution) const
    {
        return {
            solution.reversed_time_curr,
            std::numeric_limits<Profit>::min(),
        };
    }

    inline void append(
            Solution& solution,
            JobId j) const
    {
        // Update time.
        Time rj = instance_.job(j + 1).release_date;
        if (solution.time < rj)
            solution.time = rj;
        JobId j_prev = (solution.sequence.size() > 0)?
            solution.sequence.back():
            -1;
        solution.time += instance_.setup_time(j_prev + 1, j + 1);
        solution.time += instance_.job(j + 1).processing_time;
        // Update reversed_time.
        Time dj = instance_.job(j + 1).deadline;
        if (solution.time > dj) {
            solution.reversed_time_curr += (solution.time - dj);
            solution.time = dj;
        }
        // Update jobs.
        solution.sequence.push_back(j);
        // Update total weighted tardiness.
        if (solution.time > instance_.job(j + 1).due_date)
            solution.total_weighted_tardiness_curr
                += instance_.job(j + 1).weight
                * (solution.time - instance_.job(j + 1).due_date);
        // Update profit.
        solution.profit += instance_.job(j + 1).profit;

        solution.reversed_time_full = solution.reversed_time_curr;
        solution.total_weighted_tardiness_full = solution.total_weighted_tardiness_curr;
        Time time_full = solution.time;
        JobId jn = instance_.number_of_jobs() - 1;
        Time rjn = instance_.job(jn).release_date;
        if (time_full < rjn)
            time_full = rjn;
        time_full += instance_.setup_time(j + 1, jn);
        time_full += instance_.job(jn).processing_time;
        Time djn = instance_.job(jn).deadline;
        if (time_full > djn) {
            solution.reversed_time_full += (time_full - djn);
            time_full = djn;
        }
        if (time_full > instance_.job(jn).due_date)
            solution.total_weighted_tardiness_full
                += instance_.job(jn).weight
                * (time_full - instance_.job(jn).due_date);
    }

    std::ostream& print(
            std::ostream &os,
            const Solution& solution) const
    {
        Solution solution_tmp_;
        JobId j_prev = -1;
        for (JobId j: solution.sequence) {
            append(solution_tmp_, j);
            std::cout << "j " << j
                << " rj " << instance_.job(j + 1).release_date
                << " st " << instance_.setup_time(j_prev + 1, j + 1)
                << " pj " << instance_.job(j + 1).processing_time
                << " time " << solution_tmp_.time
                << " dj " << instance_.job(j + 1).due_date
                << " dj " << instance_.job(j + 1).deadline
                << " rev_cur " << solution.reversed_time_curr
                << " rev_full " << solution.reversed_time_full
                << " profit " << solution_tmp_.profit
                << " twt_cur " << solution_tmp_.total_weighted_tardiness_curr
                << " twt_full " << solution_tmp_.total_weighted_tardiness_full
                << std::endl;
            j_prev = j;
        }
        os << "sequence:";
        for (JobId j: solution.sequence)
            os << " " << j;
        os << std::endl;
        os << "cost: " << to_string(global_cost(solution)) << std::endl;
        return os;
    }

    /*
     * Outputs.
     */

    inline void write(
            const Solution& solution,
            std::string certificate_path) const
    {
        if (certificate_path.empty())
            return;
        std::ofstream cert(certificate_path);
        if (!cert.good()) {
            std::cerr << "\033[31m" << "ERROR, unable to open file \"" << certificate_path << "\"" << "\033[0m" << std::endl;
            return;
        }

        for (JobId j: solution.sequence)
            cert << j + 1 << " ";
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

