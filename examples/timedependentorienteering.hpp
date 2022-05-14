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
     * Sequences.
     */

    struct Sequence
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
     * Sequence properties.
     */

    inline GlobalCost global_cost(const Sequence& sequence) const
    {
        return {
            std::max((Time)0, sequence.time_full - instance_.maximum_duration()),
            -sequence.profit,
            sequence.time_full,
        };
    }

    /*
     * Methods required by sequencing::LocalScheme.
     */

    inline LocationPos number_of_elements() const { return instance_.number_of_locations() - 2; }

    inline GlobalCost bound(const Sequence& sequence) const
    {
        return {
            std::max((Time)0, sequence.time_full - instance_.maximum_duration()),
            std::numeric_limits<Profit>::lowest(),
            sequence.time_full,
        };
    }

    inline void append(
            Sequence& sequence,
            LocationId j) const
    {
        // Update time_cur.
        LocationId j_prev = (sequence.sequence.size() > 0)?
            sequence.sequence.back():
            -1;
        sequence.time_cur = instance_.arrival_time(j_prev + 1, j + 1, sequence.time_cur);
        // Update sequence.
        sequence.sequence.push_back(j);
        // Update profit.
        sequence.profit += instance_.location(j + 1).profit;
        // Update time_full.
        LocationId jn = instance_.number_of_locations() - 1;
        sequence.time_full = instance_.arrival_time(j + 1, jn, sequence.time_cur);
    }

    std::ostream& print(
            std::ostream &os,
            const Sequence& sequence) const
    {
        Sequence sequence_tmp_;
        for (LocationId j: sequence.sequence) {
            append(sequence_tmp_, j);
            os << "j " << j
                << " time_curr " << sequence_tmp_.time_cur
                << " time_full " << sequence_tmp_.time_full
                << " pj " << instance_.location(j + 1).profit
                << " pS " << sequence_tmp_.profit
                << std::endl;
        }
        os << "sequence:";
        for (LocationId j: sequence.sequence)
            os << " " << j;
        os << std::endl;
        os << "cost: " << to_string(global_cost(sequence)) << std::endl;
        return os;
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

