#pragma once

/**
 * Time-dependent orienteering problem.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/timedependentorienteering.hpp
 *
 */


#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing2.hpp"

#include "orproblems/timedependentorienteering.hpp"

namespace localsearchsolver
{

namespace timedependentorienteering
{

using namespace orproblems::timedependentorienteering;

class LocalScheme
{

public:

    /** Global cost: <Overtime, Profit, Total time>; */
    using GlobalCost = std::tuple<Time, Profit, Time>;

    inline Time&         overtime(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Profit&         profit(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline Time&       total_time(GlobalCost& global_cost) { return std::get<2>(global_cost); }
    inline Time    overtime(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline Profit    profit(const GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline Time  total_time(const GlobalCost& global_cost) { return std::get<2>(global_cost); }

    /*
     * Solutions.
     */

    struct Solution
    {
        std::vector<LocationId> sequence;
        Time time_cur = 0;
        Time time_full = 0;
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
        return solution;
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            std::max((Time)0, solution.time_full - instance_.maximum_duration()),
            -solution.profit,
            solution.time_full,
        };
    }

    /*
     * Methods required by sequencing::LocalScheme.
     */

    inline LocationPos number_of_elements() const { return instance_.number_of_locations() - 2; }

    inline GlobalCost bound(const Solution& solution) const
    {
        return {
            std::max((Time)0, solution.time_full - instance_.maximum_duration()),
            std::numeric_limits<Profit>::lowest(),
            solution.time_full,
        };
    }

    inline void append(
            Solution& solution,
            LocationId j) const
    {
        // Update time_cur.
        LocationId j_prev = (solution.sequence.size() > 0)?
            solution.sequence.back():
            -1;
        //std::cout << solution.time_cur << std::endl;
        solution.time_cur = instance_.arrival_time(j_prev + 1, j + 1, solution.time_cur);
        //std::cout << solution.time_cur << std::endl;
        // Update sequence.
        solution.sequence.push_back(j);
        // Update profit.
        solution.profit += instance_.location(j + 1).profit;
        // Update time_full.
        LocationId jn = instance_.number_of_locations() - 1;
        solution.time_full = instance_.arrival_time(j + 1, jn, solution.time_cur);
    }

    std::ostream& print(
            std::ostream &os,
            const Solution& solution) const
    {
        Solution solution_tmp_;
        for (LocationId j: solution.sequence) {
            append(solution_tmp_, j);
            os << "j " << j
                << " time_curr " << solution_tmp_.time_cur
                << " time_full " << solution_tmp_.time_full
                << " pj " << instance_.location(j + 1).profit
                << " pS " << solution_tmp_.profit
                << std::endl;
        }
        os << "sequence:";
        for (LocationId j: solution.sequence)
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

        for (LocationId j: solution.sequence)
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

