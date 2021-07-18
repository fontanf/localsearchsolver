#pragma once

/**
 * Quadratic Assignment Problem.
 *
 * Problem description:
 * See https://github.com/fontanf/orproblems/blob/main/orproblems/quadraticassignment.hpp
 *
 * TODO
 *
 */

#include "localsearchsolver/common.hpp"

#include "orproblems/quadraticassignment.hpp"

#include "optimizationtools/indexed_set.hpp"

namespace localsearchsolver
{

namespace quadraticassignment
{

using namespace orproblems::quadraticassignment;

class LocalScheme
{

public:

    /** Global cost: <Facility number, Cost>; */
    using GlobalCost = std::tuple<FacilityId, Cost>;

    inline FacilityId&       facility_number(GlobalCost& global_cost) const { return std::get<0>(global_cost); }
    inline Cost&                        cost(GlobalCost& global_cost) const { return std::get<1>(global_cost); }
    inline FacilityId  facility_number(const GlobalCost& global_cost) const { return std::get<0>(global_cost); }
    inline Cost                   cost(const GlobalCost& global_cost) const { return std::get<1>(global_cost); }

    static GlobalCost global_cost_worst()
    {
        return {
            std::numeric_limits<FacilityId>::max(),
            std::numeric_limits<Cost>::max(),
        };
    }

    /*
     * Solutions.
     */

    using CompactSolution = std::vector<LocationId>;

    struct CompactSolutionHasher
    {
        std::hash<LocationId> hasher;

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
            for (LocationId location_id: *compact_solution)
                optimizationtools::hash_combine(hash, hasher(location_id));
            return hash;
        }
    };

    inline CompactSolutionHasher compact_solution_hasher() const { return CompactSolutionHasher(); }

    struct Solution
    {
        std::vector<LocationId> locations;
        FacilityId facility_number = 0;
        Cost cost = 0;
    };

    CompactSolution solution2compact(const Solution& solution)
    {
        return solution.locations;
    }

    Solution compact2solution(const CompactSolution& compact_solution)
    {
        Solution solution;
        solution.locations = compact_solution;
        for (FacilityId facility_id_1 = 0; facility_id_1 < instance_.facility_number(); ++facility_id_1) {
            if (solution.locations[facility_id_1] == -1)
                continue;
            solution.facility_number++;
            for (FacilityId facility_id_2 = 0; facility_id_2 < instance_.facility_number(); ++facility_id_2) {
                if (solution.locations[facility_id_2] == -1)
                    continue;
                solution.cost += instance_.flow(facility_id_1, facility_id_2)
                    * instance_.distance(
                            solution.locations[facility_id_1],
                            solution.locations[facility_id_2]);
            }
        }
        return solution;
    }

    /*
     * Constructors and destructor.
     */

    struct Parameters
    {
    };

    LocalScheme(
            const Instance& instance,
            Parameters parameters):
        instance_(instance),
        parameters_(parameters),
        facilities_(instance.facility_number())
    {
        for (FacilityId facility_id_1 = 0; facility_id_1 < instance.facility_number(); ++facility_id_1)
            for (FacilityId facility_id_2 = facility_id_1 + 1; facility_id_2 < instance.facility_number(); ++facility_id_2)
                facilities_.push_back({facility_id_1, facility_id_2});
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
        solution.locations.resize(instance_.facility_number(), -1);
        return solution;
    }

    inline Solution initial_solution(
            Counter,
            std::mt19937_64& generator) const
    {
        Solution solution = empty_solution();
        std::vector<LocationId> locations(instance_.facility_number(), -1);
        std::iota(locations.begin(), locations.end(), 0);
        std::shuffle(locations.begin(), locations.end(), generator);
        for (FacilityId facility_id = 0; facility_id < instance_.facility_number(); ++facility_id)
            add(solution, facility_id, locations[facility_id]);
        return solution;
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return {
            -solution.facility_number,
            solution.cost,
        };
    }

    /*
     * Local search.
     */

    struct Move
    {
        FacilityId facility_id_1;
        FacilityId facility_id_2;
        GlobalCost global_cost;
    };

    static Move move_null() { return {-1, -1, global_cost_worst()}; }

    struct MoveHasher
    {
        std::hash<FacilityId> hasher;

        inline bool hashable(const Move&) const { return true; }

        inline bool operator()(
                const Move& move_1,
                const Move& move_2) const
        {
            return move_1.facility_id_1 == move_2.facility_id_1
                && move_1.facility_id_2 == move_2.facility_id_2;
        }

        inline std::size_t operator()(
                const Move& move) const
        {
            size_t hash = hasher(move.facility_id_1);
            optimizationtools::hash_combine(hash, hasher(move.facility_id_2));
            return hash;
        }
    };

    inline MoveHasher move_hasher() const { return MoveHasher(); }

    inline std::vector<Move> perturbations(
            const Solution& solution,
            std::mt19937_64&)
    {
        std::vector<Move> moves;
        for (FacilityId facility_id_1 = 0;
                facility_id_1 < instance_.facility_number();
                ++facility_id_1) {
            for (FacilityId facility_id_2 = facility_id_1 + 1;
                    facility_id_2 < instance_.facility_number();
                    ++facility_id_2) {
                Move move;
                move.facility_id_1 = facility_id_1;
                move.facility_id_2 = facility_id_2;
                move.global_cost = cost_swap(solution, facility_id_1, facility_id_2);
                moves.push_back(move);
            }
        }
        return moves;
    }

    inline void apply_move(Solution& solution, const Move& move) const
    {
        LocationId location_id_1 = solution.locations[move.facility_id_1];
        LocationId location_id_2 = solution.locations[move.facility_id_2];
        remove(solution, move.facility_id_1);
        remove(solution, move.facility_id_2);
        add(solution, move.facility_id_1, location_id_2);
        add(solution, move.facility_id_2, location_id_1);
    }

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Move& tabu = move_null())
    {
        Counter it = 0;
        for (;; ++it) {
            //std::cout << "it " << it << " cost " << to_string(global_cost(solution)) << std::endl;
            std::shuffle(facilities_.begin(), facilities_.end(), generator);
            FacilityId facility_id_1_best = -1;
            FacilityId facility_id_2_best = -1;
            GlobalCost c_best = global_cost(solution);
            for (auto p: facilities_) {
                FacilityId facility_id_1 = p.first;
                FacilityId facility_id_2 = p.second;
                if (facility_id_1 == tabu.facility_id_1
                        || facility_id_1 == tabu.facility_id_2
                        || facility_id_2 == tabu.facility_id_1
                        || facility_id_2 == tabu.facility_id_2)
                    continue;
                GlobalCost c = cost_swap(
                        solution,
                        facility_id_1,
                        facility_id_2);
                if (c >= c_best)
                    continue;
                if (facility_id_1 != -1 && !dominates(c, c_best))
                    continue;
                facility_id_1_best = facility_id_1;
                facility_id_2_best = facility_id_2;
                c_best = c;
            }
            if (facility_id_1_best == -1)
                break;
            LocationId location_id_1 = solution.locations[facility_id_1_best];
            LocationId location_id_2 = solution.locations[facility_id_2_best];
            remove(solution, facility_id_1_best);
            remove(solution, facility_id_2_best);
            add(solution, facility_id_1_best, location_id_2);
            add(solution, facility_id_2_best, location_id_1);
        }
    }

    /*
     * Outputs.
     */

    std::ostream& print(
            std::ostream &os,
            const Solution& solution) const
    {
        os << "locations:";
        for (LocationId location_id: solution.locations)
            os << " " << location_id;
        os << std::endl;
        os << "cost: " << solution.cost << std::endl;
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

        for (LocationId location_id: solution.locations)
            cert << location_id << " ";
    }

private:

    /*
     * Manipulate solutions.
     */

    inline void add(Solution& solution, FacilityId facility_id, LocationId location_id) const
    {
        assert(solution.locations[facility_id] == -1);
        solution.locations[facility_id] = location_id;
        for (FacilityId facility_id_2 = 0; facility_id_2 < instance_.facility_number(); facility_id_2++) {
            if (solution.locations[facility_id_2] == -1)
                continue;
            solution.cost += instance_.flow(facility_id, facility_id_2)
                * instance_.distance(
                        solution.locations[facility_id],
                        solution.locations[facility_id_2]);
            solution.cost += instance_.flow(facility_id_2, facility_id)
                * instance_.distance(
                        solution.locations[facility_id_2],
                        solution.locations[facility_id]);
        }
        solution.facility_number++;
    }

    inline void remove(Solution& solution, FacilityId facility_id) const
    {
        assert(solution.locations[facility_id] != -1);
        for (FacilityId facility_id_2 = 0; facility_id_2 < instance_.facility_number(); facility_id_2++) {
            if (solution.locations[facility_id_2] == -1)
                continue;
            solution.cost -= instance_.flow(facility_id, facility_id_2)
                * instance_.distance(
                        solution.locations[facility_id],
                        solution.locations[facility_id_2]);
            solution.cost -= instance_.flow(facility_id_2, facility_id)
                * instance_.distance(
                        solution.locations[facility_id_2],
                        solution.locations[facility_id]);
        }
        solution.locations[facility_id] = -1;
        solution.facility_number--;
    }

    /*
     * Evaluate moves.
     */

    inline GlobalCost cost_swap(
            const Solution& solution,
            FacilityId facility_id_1,
            FacilityId facility_id_2) const
    {
        LocationId location_id_1 = solution.locations[facility_id_1];
        LocationId location_id_2 = solution.locations[facility_id_2];
        assert(location_id_1 != -1);
        assert(location_id_2 != -1);
        GlobalCost gc = global_cost(solution);
        cost(gc)
            += instance_.flow(facility_id_1, facility_id_1)
            * (instance_.distance(location_id_2, location_id_2)
                    - instance_.distance(location_id_1, location_id_1))
            + instance_.flow(facility_id_2, facility_id_2)
            * (instance_.distance(location_id_1, location_id_1)
                    - instance_.distance(location_id_2, location_id_2))
            + instance_.flow(facility_id_1, facility_id_2)
            * (instance_.distance(location_id_2, location_id_1)
                    - instance_.distance(location_id_1, location_id_2))
            + instance_.flow(facility_id_2, facility_id_1)
            * (instance_.distance(location_id_1, location_id_2)
                    - instance_.distance(location_id_2, location_id_1));
        for (FacilityId facility_id_3 = 0;
                facility_id_3 < instance_.facility_number();
                ++facility_id_3) {
            if (facility_id_3 == facility_id_1)
                continue;
            if (facility_id_3 == facility_id_2)
                continue;
            LocationId location_id_3 = solution.locations[facility_id_3];
            cost(gc)
                += instance_.flow(facility_id_1, facility_id_3)
                * (instance_.distance(location_id_2, location_id_3)
                        - instance_.distance(location_id_1, location_id_3))
                + instance_.flow(facility_id_3, facility_id_1)
                * (instance_.distance(location_id_3, location_id_2)
                        - instance_.distance(location_id_3, location_id_1))
                + instance_.flow(facility_id_2, facility_id_3)
                * (instance_.distance(location_id_1, location_id_3)
                        - instance_.distance(location_id_2, location_id_3))
                + instance_.flow(facility_id_3, facility_id_2)
                * (instance_.distance(location_id_3, location_id_1)
                        - instance_.distance(location_id_3, location_id_2));
        }
        return gc;
    }

    /*
     * Private attributes.
     */

    const Instance& instance_;
    Parameters parameters_;

    std::vector<std::pair<FacilityId, FacilityId>> facilities_;

};

}

}

