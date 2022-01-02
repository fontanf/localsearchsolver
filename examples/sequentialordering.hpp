#pragma once

/**
 * Sequential Ordering Problem..
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/sequentialordering.hpp
 *
 */


#include "localsearchsolver/common.hpp"
#include "localsearchsolver/sequencing2.hpp"

#include "orproblems/sequentialordering.hpp"

namespace localsearchsolver
{

namespace sequentialordering
{

using namespace orproblems::sequentialordering;

class LocalScheme
{

public:

    /** Global cost: <Number of vertices, Number of precedence violations, Total distance>; */
    using GlobalCost = std::tuple<VertexPos, VertexPos, Distance>;

    inline VertexPos&                    number_of_vertices(GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline VertexPos&       number_of_precedence_violations(GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline Distance&                                 length(GlobalCost& global_cost) { return std::get<2>(global_cost); }
    inline VertexPos               number_of_vertices(const GlobalCost& global_cost) { return std::get<0>(global_cost); }
    inline VertexPos  number_of_precedence_violations(const GlobalCost& global_cost) { return std::get<1>(global_cost); }
    inline Distance                            length(const GlobalCost& global_cost) { return std::get<2>(global_cost); }

    /*
     * Solutions.
     */

    struct Solution
    {
        std::vector<VertexId> sequence;
        Distance length = 0;
        VertexPos number_of_precedence_violations = 0;
        std::vector<uint8_t> contains;
    };

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
        Parameters()
        {
            sequencing_parameters.shift_block_maximum_length = 8;
            sequencing_parameters.swap_block_maximum_length = 2;
            sequencing_parameters.reverse = true;
            sequencing_parameters.shift_reverse_block_maximum_length = 3;
            sequencing_parameters.shuffle_neighborhood_order = true;
            sequencing_parameters.double_bridge_number_of_perturbations = 0;
            sequencing_parameters.ruin_and_recreate_number_of_perturbations = 10;
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
        solution.contains = std::vector<uint8_t>(
                instance_.number_of_vertices(), 0);
        return solution;
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            -solution.sequence.size(),
            solution.number_of_precedence_violations,
            solution.length,
        };
    }

    inline GlobalCost global_cost_cutoff(double cutoff) const
    {
        return {
            -instance_.number_of_vertices(),
            0,
            cutoff,
        };
    }

    /*
     * Methods required by sequencing::LocalScheme.
     */

    inline VertexPos number_of_elements() const { return instance_.number_of_vertices(); }

    inline GlobalCost bound(const Solution& solution) const
    {
        return {
            -instance_.number_of_vertices(),
            solution.number_of_precedence_violations,
            solution.length,
        };
    }

    inline void append(
            Solution& solution,
            VertexPos j) const
    {
        // Update number_of_precedence_violations.
        for (VertexId j: instance_.predecessors(j))
            if (!solution.contains[j])
                solution.number_of_precedence_violations++;
        // Update time.
        if (solution.sequence.size() > 0) {
            VertexId j_prev = solution.sequence.back();
            Distance d = instance_.distance(j_prev, j);
            if (d == std::numeric_limits<Distance>::max())
                d = instance_.distance(j, j_prev);
            solution.length += d;
        }
        // Update contains.
        solution.contains[j] = 1;
        // Update sequence.
        solution.sequence.push_back(j);
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

