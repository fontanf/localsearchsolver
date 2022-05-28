#pragma once

/**
 * Sequencing problems.
 */

#include "localsearchsolver/common.hpp"

#include "optimizationtools/containers/indexed_set.hpp"

namespace localsearchsolver
{

namespace sequencing2
{

using SequenceId = int64_t;
using SequencePos = int64_t;
using ElementId = int64_t;
using ElementPos = int64_t;

enum class Perturbations
{
    None,
    DoubleBridge,
    RuinAndRecreate,
    ForceAdd,
};

enum class Neighborhoods
{
    Shift,
    SwapK,
    SwapK1K2,
    Reverse,
    ShiftReverse,

    Add,
    Remove,

    InterTwoOpt,
    InterShift,
    InterSwap,
    InterShiftReverse,
};

struct Parameters
{

    bool shuffle_neighborhood_order = true;

    /*
     * Neighborhoods - Intra.
     */

    ElementPos shift_block_maximum_length = 3;
    ElementPos shift_maximum_distance = 1024;

    ElementPos swap_maximum_distance = 1024;
    ElementPos swap_block_maximum_length = 2;

    bool reverse = false;
    ElementPos reverse_maximum_length = 1024;

    ElementPos shift_reverse_block_maximum_length = 0;

    bool add_remove = false;

    /*
     * Neighborhoods - Inter.
     */

    bool inter_two_opt = false;

    ElementPos inter_shift_block_maximum_length = 0;

    ElementPos inter_swap_block_maximum_length = 0;

    ElementPos inter_shift_reverse_block_maximum_length = 0;

    /*
     * Perturbations.
     */

    Counter double_bridge_number_of_perturbations = 10;

    Counter ruin_and_recreate_number_of_perturbations = 0;
    ElementPos ruin_and_recreate_number_of_elements_removed = 4;

    bool force_add = false;

    /*
     * Crossovers.
     */

    double crossover_ox_weight = 1;
    double crossover_sjox_weight = 0;
    double crossover_sbox_weight = 0;
};

template <typename LocalScheme0>
class LocalScheme
{

public:

    using GlobalCost = typename LocalScheme0::GlobalCost;

    using Sequence = typename LocalScheme0::Sequence;

    struct Solution
    {
        std::vector<Sequence> sequences;
        GlobalCost global_cost;
    };

    using CompactSolution = std::vector<std::vector<ElementId>>;

    struct CompactSolutionHasher
    {
        std::hash<ElementId> hasher;

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
            for (const auto& sequence: *compact_solution) {
                size_t hash_tmp = 0;
                for (ElementId j: sequence)
                    optimizationtools::hash_combine(hash_tmp, hasher(j));
                optimizationtools::hash_combine(hash, hash_tmp);
            }
            return hash;
        }
    };

    inline CompactSolutionHasher compact_solution_hasher() const { return CompactSolutionHasher(); }

    Solution compact2solution(const CompactSolution& compact_solution)
    {
        auto solution = empty_solution();
        for (SequenceId i = 0; i < number_of_sequences(); ++i)
            for (ElementId j: compact_solution[i])
                local_scheme_0_.append(solution.sequences[i], j);
        compute_global_cost(solution);
        return solution;
    }

    CompactSolution solution2compact(const Solution& solution)
    {
        CompactSolution compact_solution(number_of_sequences());
        for (SequenceId i = 0; i < number_of_sequences(); ++i)
            compact_solution[i] = solution.sequences[i].sequence;
        return compact_solution;
    }

    /*
     * Constructors and destructor.
     */

    LocalScheme(
            LocalScheme0& local_scheme_0,
            Parameters parameters):
        local_scheme_0_(local_scheme_0),
        parameters_(parameters),
        sequences_1_(number_of_sequences()),
        sequences_2_(number_of_sequences()),
        elements_(local_scheme_0_.number_of_elements()),
        positions_1_(local_scheme_0_.number_of_elements()),
        positions_2_(local_scheme_0_.number_of_elements()),
        global_costs_1d_(local_scheme_0_.number_of_elements() + 1, worst<GlobalCost>()),
        global_costs_2d_1_(local_scheme_0_.number_of_elements()),
        global_costs_2d_2_(
                local_scheme_0_.number_of_elements(),
                std::vector<GlobalCost>(
                    local_scheme_0_.number_of_elements(),
                    worst<GlobalCost>())),
        solution_cur_(empty_solution()),
        solution_tmp_(empty_solution())
    {
        std::iota(sequences_1_.begin(), sequences_1_.end(), 0);
        std::iota(sequences_2_.begin(), sequences_2_.end(), 0);
        std::iota(positions_1_.begin(), positions_1_.end(), 0);
        std::iota(positions_2_.begin(), positions_2_.end(), 0);
        for (ElementPos pos_1 = 0; pos_1 < local_scheme_0_.number_of_elements(); ++pos_1)
            for (ElementPos pos_2 = pos_1 + 1; pos_2 < local_scheme_0_.number_of_elements(); ++pos_2)
                pairs_1_.push_back({pos_1, pos_2});
        for (ElementPos pos_1 = 0; pos_1 < local_scheme_0_.number_of_elements(); ++pos_1)
            for (ElementPos pos_2 = 0; pos_2 < local_scheme_0_.number_of_elements(); ++pos_2)
                if (pos_1 != pos_2)
                    pairs_2_.push_back({pos_1, pos_2});
        for (ElementPos pos_1 = 0; pos_1 < local_scheme_0_.number_of_elements(); ++pos_1)
            global_costs_2d_1_[pos_1].resize(
                    local_scheme_0_.number_of_elements() - pos_1,
                    worst<GlobalCost>());

        // Initialize statistics structures.
        shift_number_of_explorations_ = std::vector<Counter>(
                parameters_.shift_block_maximum_length + 1, 0);
        shift_number_of_sucesses_ = std::vector<Counter>(
                parameters_.shift_block_maximum_length + 1, 0);
        swap_number_of_explorations_ = std::vector<std::vector<Counter>>(
                parameters_.swap_block_maximum_length + 1);
        swap_number_of_sucesses_ = std::vector<std::vector<Counter>>(
                parameters_.swap_block_maximum_length + 1);
        for (ElementPos p = 0; p <= parameters_.swap_block_maximum_length; ++p) {
            swap_number_of_explorations_[p] = std::vector<Counter>(p + 1, 0);
            swap_number_of_sucesses_[p] = std::vector<Counter>(p + 1, 0);
        }
        shift_reverse_number_of_explorations_ = std::vector<Counter>(
                parameters_.shift_reverse_block_maximum_length + 1, 0);
        shift_reverse_number_of_sucesses_ = std::vector<Counter>(
                parameters_.shift_reverse_block_maximum_length + 1, 0);
        inter_shift_number_of_explorations_ = std::vector<Counter>(
                parameters_.inter_shift_block_maximum_length + 1, 0);
        inter_shift_number_of_sucesses_ = std::vector<Counter>(
                parameters_.inter_shift_block_maximum_length + 1, 0);
        inter_swap_number_of_explorations_ = std::vector<std::vector<Counter>>(
                parameters_.inter_swap_block_maximum_length + 1);
        inter_swap_number_of_sucesses_ = std::vector<std::vector<Counter>>(
                parameters_.inter_swap_block_maximum_length + 1);
        for (ElementPos p = 0; p <= parameters_.inter_swap_block_maximum_length; ++p) {
            inter_swap_number_of_explorations_[p] = std::vector<Counter>(p + 1, 0);
            inter_swap_number_of_sucesses_[p] = std::vector<Counter>(p + 1, 0);
        }
        inter_shift_reverse_number_of_explorations_ = std::vector<Counter>(
                parameters_.inter_shift_reverse_block_maximum_length + 1, 0);
        inter_shift_reverse_number_of_sucesses_ = std::vector<Counter>(
                parameters_.inter_shift_reverse_block_maximum_length + 1, 0);
    }

    LocalScheme(const LocalScheme& sequencing_scheme):
        LocalScheme(sequencing_scheme.local_scheme_0_, sequencing_scheme.parameters_) { }

    virtual ~LocalScheme() { }

    /*
     * Number of sequences.
     */

    template<typename, typename T>
    struct HasNumberOfSequencesMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasNumberOfSequencesMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().number_of_sequences(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    SequencePos number_of_sequences(
            std::false_type) const
    {
        return 1;
    }

    SequencePos number_of_sequences(
            std::true_type) const
    {
        return local_scheme_0_.number_of_sequences();;
    }

    /**
     * Return the number of sequences in the problem.
     *
     * If LocalScheme0 does not implement a method:
     * 'SequencePos number_of_sequences()',
     * then it returns 1.
     */
    SequencePos number_of_sequences() const
    {
        return number_of_sequences(
                std::integral_constant<
                    bool,
                    HasNumberOfSequencesMethod<LocalScheme0, SequencePos()>::value>());
    }

    /*
     * Empty sequence.
     */

    template<typename, typename T>
    struct HasEmptySequenceMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasEmptySequenceMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().empty_sequence(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    Sequence empty_sequence(
            SequenceId,
            std::false_type) const
    {
        return Sequence();
    }

    Sequence empty_sequence(
            SequenceId i,
            std::true_type) const
    {
        return local_scheme_0_.empty_sequence(i);;
    }

    Sequence empty_sequence(SequenceId i) const
    {
        return empty_sequence(
                i,
                std::integral_constant<
                    bool,
                    HasEmptySequenceMethod<LocalScheme0, Sequence(SequenceId)>::value>());
    }

    inline Solution empty_solution() const
    {
        Solution solution;
        for (SequenceId i = 0; i < number_of_sequences(); ++i)
            solution.sequences.push_back(empty_sequence(i));
        return solution;
    }

    /*
     * Initial solution.
     */

    template<typename, typename T>
    struct HasInitialSolutionMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasInitialSolutionMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().initial_solution(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    Solution initial_solution(
            Counter initial_solution_id,
            std::mt19937_64& generator,
            std::false_type)
    {
        SequenceId m = number_of_sequences();
        switch (initial_solution_id) {
        case 1: {
            // Fix a random order, find the best position.
            std::vector<ElementId> elements(local_scheme_0_.number_of_elements());
            std::iota(elements.begin(), elements.end(), 0);
            std::shuffle(elements.begin(), elements.end(), generator);
            Solution solution = empty_solution();
            for (ElementId j: elements) {
                GlobalCost c_best = global_cost(solution);
                SequenceId i_best = -1;
                ElementPos pos_best = -1;
                std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                std::shuffle(positions_2_.begin(), positions_2_.end(), generator);
                for (SequenceId i: sequences_1_) {
                    Sequence& sequence = solution.sequences[i];
                    compute_cost_add(solution, i, j);
                    for (ElementPos pos: positions_2_) {
                        if (pos > (ElementPos)sequence.sequence.size())
                            continue;
                        GlobalCost c = global_costs_1d_[pos];
                        if (c >= c_best)
                            continue;
                        if (pos_best != -1 && !dominates(c, c_best))
                            continue;
                        i_best = i;
                        pos_best = pos;
                        c_best = c;
                    }
                }
                if (i_best != -1) {
                    for (SequenceId i = 0; i < number_of_sequences(); ++i)
                        if (i != i_best)
                            solution_tmp_.sequences[i] = solution.sequences[i];
                    const Sequence& sequence = solution.sequences[i_best];
                    Sequence& sequence_tmp = solution_tmp_.sequences[i_best];
                    sequence_tmp = empty_sequence(i_best);
                    for (ElementPos p = 0; p < pos_best; ++p)
                        local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                    local_scheme_0_.append(sequence_tmp, j);
                    for (ElementPos p = pos_best; p < (ElementPos)sequence.sequence.size(); ++p)
                        local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                    compute_global_cost(solution_tmp_);
                    solution = solution_tmp_;
                }
            }
            compute_global_cost(solution);
            return solution;
        } case 2: {
            // Insert at the end, find the best element.
            std::vector<uint8_t> contains(local_scheme_0_.number_of_elements(), 0);
            Solution solution = empty_solution();
            for (;;) {
                GlobalCost c_best = global_cost(solution);
                //std::cout << to_string(c_best) << std::endl;
                ElementId j_best = -1;
                SequenceId i_best = -1;
                std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                std::shuffle(positions_1_.begin(), positions_1_.end(), generator);
                for (SequenceId i: sequences_1_) {
                    Sequence& sequence = solution.sequences[i];
                    for (ElementId j: positions_1_) {
                        if (contains[j])
                            continue;
                        sequence_tmp_1_ = sequence;
                        local_scheme_0_.append(sequence_tmp_1_, j);
                        GlobalCost c = global_cost(sequence_tmp_1_);
                        if (c >= c_best)
                            continue;
                        c_best = c;
                        j_best = j;
                        i_best = i;
                    }
                }
                if (j_best != -1) {
                    local_scheme_0_.append(solution.sequences[i_best], j_best);
                    contains[j_best] = 1;
                    continue;
                }
                break;
            }
            compute_global_cost(solution);
            return solution;
        } case 3: {
            // Find the best element x position.
            // Warning, this one can be expensive.
            std::vector<uint8_t> contains(local_scheme_0_.number_of_elements(), 0);
            Solution solution = empty_solution();
            for (;;) {
                GlobalCost c_best = global_cost(solution);
                //std::cout << to_string(c_best) << std::endl;
                SequenceId i_best = -1;
                ElementId j_best = -1;
                ElementPos pos_best = -1;
                std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                std::shuffle(positions_1_.begin(), positions_1_.end(), generator);
                std::shuffle(positions_2_.begin(), positions_2_.end(), generator);
                for (SequenceId i: sequences_1_) {
                    Sequence& sequence = solution.sequences[i];
                    for (ElementId j: positions_1_) {
                        if (contains[j])
                            continue;
                        compute_cost_add(solution, i, j);
                        for (ElementPos pos: positions_2_) {
                            if (pos > (ElementPos)sequence.sequence.size())
                                continue;
                            GlobalCost c = global_costs_1d_[pos];
                            if (c >= c_best)
                                continue;
                            pos_best = pos;
                            c_best = c;
                            j_best = j;
                            i_best = i;
                        }
                    }
                }
                if (i_best != -1) {
                    for (SequenceId i = 0; i < number_of_sequences(); ++i)
                        if (i != i_best)
                            solution_tmp_.sequences[i] = solution.sequences[i];
                    const Sequence& sequence = solution.sequences[i_best];
                    Sequence& sequence_tmp = solution_tmp_.sequences[i_best];
                    sequence_tmp = empty_sequence(i_best);
                    for (ElementPos p = 0; p < pos_best; ++p)
                        local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                    local_scheme_0_.append(sequence_tmp, j_best);
                    for (ElementPos p = pos_best; p < (ElementPos)sequence.sequence.size(); ++p)
                        local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                    compute_global_cost(solution_tmp_);
                    solution = solution_tmp_;
                    contains[j_best] = 1;
                    continue;
                }
                break;
            }
            compute_global_cost(solution);
            return solution;
        } default: {
            // Random permutation.
            elements_.fill();
            elements_.shuffle(generator);
            std::uniform_int_distribution<SequenceId> di(0, m - 1);
            Solution solution = empty_solution();
            for (ElementId j: elements_) {
                SequenceId i = di(generator);
                local_scheme_0_.append(solution.sequences[i], j);
            }
            compute_global_cost(solution);
            return solution;
        }
        }
    }

    Solution initial_solution(
            Counter initial_solution_id,
            std::mt19937_64& generator,
            std::true_type)
    {
        if (initial_solution_id >= 4)
            return local_scheme_0_.initial_solution(initial_solution_id, generator);
        return initial_solution(initial_solution_id, generator, false);
    }

    inline Solution initial_solution(
            Counter initial_solution_id,
            std::mt19937_64& generator)
    {
        return initial_solution(
                initial_solution_id,
                generator,
                std::integral_constant<
                    bool,
                    HasInitialSolutionMethod<LocalScheme0,
                    Solution(
                        Counter,
                        std::mt19937_64&)>::value>());
    }

    /*
     * Crossovers.
     */

    inline Solution crossover(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        std::discrete_distribution<Counter> d_crossover({
                parameters_.crossover_ox_weight,
                parameters_.crossover_sjox_weight,
                parameters_.crossover_sbox_weight,
                });
        Counter x = d_crossover(generator);
        switch (x) {
        case 0: {
            return crossover_ox(solution_parent_1, solution_parent_2, generator);
        } case 1: {
            return crossover_sjox(solution_parent_1, solution_parent_2, generator);
        } case 2: {
            return crossover_sbox(solution_parent_1, solution_parent_2, generator);
        } default: {
            return crossover_ox(solution_parent_1, solution_parent_2, generator);
        }
        }
    }

    /**
     * Generate a new solution from two parent solutions using the OX crossover
     * operator.
     *
     * The OX crossover operator consists in selecting a random substring from
     * the first parent, copying this substring into the child while leaving
     * the rest of the positions empty, and finally completing the child’s
     * missing positions circularly with the visits from the second parent,
     * starting from the end cutting point.
     *
     * References:
     * - "A simple and effective hybrid genetic search for the job sequencing
     *   and tool switching problem" (Mecler et al., 2021)
     *   https://doi.org/10.1016/j.cor.2020.105153
     */
    inline Solution crossover_ox(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        Solution solution = empty_solution();

        const Sequence& sequence_parent_1 = solution_parent_1.sequences[0];
        const Sequence& sequence_parent_2 = solution_parent_2.sequences[0];
        Sequence& sequence = solution.sequences[0];

        ElementPos n = sequence_parent_1.sequence.size();

        std::vector<ElementPos> edges = optimizationtools::bob_floyd<ElementPos>(
                2, n + 1, generator);
        std::sort(edges.begin(), edges.end());

        ElementPos pos_1 = edges[0];
        ElementPos pos_2 = edges[1];

        std::vector<uint8_t> in_substring(n, false);
        for (ElementPos pos = pos_1; pos < pos_2; ++pos) {
            ElementId j = sequence_parent_1.sequence[pos];
            in_substring[j] = true;
        }

        for (ElementPos pos = 0; pos < n; ++pos) {
            if ((ElementPos)sequence.sequence.size() == pos_1) {
                for (ElementPos p = pos_1; p < pos_2; ++p) {
                    ElementId j = sequence_parent_1.sequence[p];
                    local_scheme_0_.append(sequence, j);
                }
            }
            ElementId j = sequence_parent_2.sequence[pos];
            if (in_substring[j])
                continue;
            local_scheme_0_.append(sequence, j);
        }
        if ((ElementPos)sequence.sequence.size() == pos_1) {
            for (ElementPos p = pos_1; p < pos_2; ++p) {
                ElementId j = sequence_parent_1.sequence[p];
                local_scheme_0_.append(sequence, j);
            }
        }

        return solution;
    }

    /**
     * Generate a new solution from two parent solutions using the SJOX
     * crossover operator.
     *
     * References:
     * - "Two new robust genetic algorithms for the flowshop scheduling
     *   problem" (Ruiz et al., 2006)
     *   https://doi.org/10.1016/j.omega.2004.12.006
     */
    inline Solution crossover_sjox(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        Solution solution = empty_solution();

        const Sequence& sequence_parent_1 = solution_parent_1.sequences[0];
        const Sequence& sequence_parent_2 = solution_parent_2.sequences[0];
        Sequence& sequence = solution.sequences[0];

        ElementPos n = sequence_parent_1.sequence.size();
        std::vector<ElementPos> positions(n, -1);

        // Add elements from parent_1 up to a given cut point.
        std::uniform_int_distribution<ElementPos> d_point(1, n);
        ElementPos pos_0 = d_point(generator);
        for (ElementPos pos = 0; pos < pos_0; ++pos) {
            ElementId j = sequence_parent_1.sequence[pos];
            positions[j] = pos;
            local_scheme_0_.append(sequence, j);
        }

        // Add elements from parent_2 keeping the relative order.
        for (ElementPos pos = 0; pos < n; ++pos) {
            // Add elements which have the same positions in both parents.
            for (;;) {
                ElementPos p = sequence.sequence.size();
                if (p == n)
                    break;
                ElementId j1 = sequence_parent_1.sequence[p];
                ElementId j2 = sequence_parent_2.sequence[p];
                if (j1 == j2) {
                    positions[j1] = p;
                    local_scheme_0_.append(sequence, j1);
                    continue;
                }
                break;
            }
            ElementId j = sequence_parent_2.sequence[pos];
            if (positions[j] != -1)
                continue;
            positions[j] = sequence.sequence.size();
            local_scheme_0_.append(sequence, j);
        }

        return solution;
    }

    /**
     * Generate a new solution from two parent solutions using the SBOX
     * crossover operator.
     *
     * References:
     * - "Two new robust genetic algorithms for the flowshop scheduling
     *   problem" (Ruiz et al., 2006)
     *   https://doi.org/10.1016/j.omega.2004.12.006
     */
    inline Solution crossover_sbox(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        Solution solution = empty_solution();

        const Sequence& sequence_parent_1 = solution_parent_1.sequences[0];
        const Sequence& sequence_parent_2 = solution_parent_2.sequences[0];
        Sequence& sequence = solution.sequences[0];

        ElementPos n = sequence_parent_1.sequence.size();
        std::vector<ElementPos> positions(n, -1);

        // Add elements from parent_1 up to a given cut point.
        std::uniform_int_distribution<ElementPos> d_point(1, n);
        ElementPos pos_0 = d_point(generator);
        for (ElementPos pos = 0; pos < pos_0; ++pos) {
            ElementId j = sequence_parent_1.sequence[pos];
            positions[j] = pos;
            local_scheme_0_.append(sequence, j);
        }

        // Add elements from parent_2 keeping the relative order.
        for (ElementPos pos = 0; pos < n; ++pos) {
            // Add elements which have the same positions in both parents.
            for (;;) {
                ElementPos p = sequence.sequence.size();
                if (p == n)
                    break;
                ElementId j1 = sequence_parent_1.sequence[p];
                ElementId j2 = sequence_parent_2.sequence[p];
                if (p <= n - 1) {
                    ElementId j1_next = sequence_parent_1.sequence[p + 1];
                    ElementId j2_next = sequence_parent_2.sequence[p + 1];
                    if (j1 == j2 && j1_next == j2_next) {
                        positions[j1] = p;
                        local_scheme_0_.append(sequence, j1);
                        continue;
                    }
                }
                if (p >= 1) {
                    ElementId j1_prev = sequence_parent_1.sequence[p - 1];
                    ElementId j2_prev = sequence_parent_2.sequence[p - 1];
                    if (j1 == j2 && j1_prev == j2_prev) {
                        positions[j1] = p;
                        local_scheme_0_.append(sequence, j1);
                        continue;
                    }
                }
                break;
            }
            ElementId j = sequence_parent_2.sequence[pos];
            if (positions[j] != -1)
                continue;
            positions[j] = sequence.sequence.size();
            local_scheme_0_.append(sequence, j);
        }

        return solution;
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return solution.global_cost;
    }

    inline GlobalCost global_cost(const Sequence& sequence) const
    {
        return local_scheme_0_.global_cost(sequence);
    }

    inline GlobalCost bound(const Sequence& sequence) const
    {
        return local_scheme_0_.bound(sequence);
    }

    GlobalCost global_cost_goal(double value) const
    {
        return localsearchsolver::global_cost_goal(local_scheme_0_, value);
    }

    /**
     * Return the distance between two solutions.
     *
     * The distance between two solutions (represented as element permutations)
     * is defined as the number of different edges between them (broken-pairs
     * distance).
     *
     * References:
     * - "A simple and effective hybrid genetic search for the job sequencing
     *   and tool switching problem" (Mecler et al., 2021)
     *   https://doi.org/10.1016/j.cor.2020.105153
     */
    inline ElementPos distance(
            const Solution& solution_1,
            const Solution& solution_2) const
    {
        const Sequence& sequence_1 = solution_1.sequences[0];
        const Sequence& sequence_2 = solution_2.sequences[0];

        ElementPos n = sequence_1.sequence.size();
        std::vector<ElementId> next_1(n, -1);
        for (ElementPos pos = 0; pos < n - 1; ++pos) {
            ElementId j = sequence_1.sequence[pos];
            ElementId j_next = sequence_1.sequence[pos + 1];
            next_1[j] = j_next;
        }

        ElementPos d = 0;
        for (ElementPos pos = 0; pos < n - 1; ++pos) {
            ElementId j = sequence_2.sequence[pos];
            ElementId j_next = sequence_2.sequence[pos + 1];
            if (j_next != next_1[j])
                d++;
        }
        return d;
    }

    /*
     * Local search.
     */

    struct Move
    {
        Move(): type(Perturbations::None), global_cost(worst<GlobalCost>()) { }

        /** Type of perturbation. */
        Perturbations type;
        ElementId force_add_j = -1;
        GlobalCost global_cost;
    };

    struct MoveHasher
    {
        std::hash<ElementPos> hasher;

        inline bool hashable(const Move& move) const
        {
            if (move.type == Perturbations::ForceAdd)
                return true;
            return false;
        }

        inline bool operator()(
                const Move& move_1,
                const Move& move_2) const
        {
            if (move_1.type == Perturbations::ForceAdd
                    && move_2.type == Perturbations::ForceAdd
                    && move_1.force_add_j == move_2.force_add_j)
                return true;
            return false;
        }

        inline std::size_t operator()(
                const Move& move) const
        {
            return hasher(move.force_add_j);
        }
    };

    inline MoveHasher move_hasher() const { return MoveHasher(); }

    /*
     * Perturbations.
     */

    inline std::vector<Move> perturbations(
            const Solution& solution,
            std::mt19937_64& generator)
    {
        std::vector<Move> moves;

        // Double-bridge.
        for (Counter perturbation = 0;
                perturbation < parameters_.double_bridge_number_of_perturbations;
                ++perturbation) {
            Move move;
            move.type = Perturbations::DoubleBridge;
            move.global_cost = global_cost(solution);
            moves.push_back(move);
        }

        // Ruin-and-recreate.
        for (Counter perturbation = 0;
                perturbation < parameters_.ruin_and_recreate_number_of_perturbations;
                ++perturbation) {
            Move move;
            move.type = Perturbations::RuinAndRecreate;
            move.global_cost = global_cost(solution);
            moves.push_back(move);
        }

        // Force-add.
        if (parameters_.force_add) {
            std::vector<uint8_t> contains(local_scheme_0_.number_of_elements(), 0);
            for (SequenceId i = 0; i < number_of_sequences(); ++i)
                for (ElementId j: solution.sequences[i].sequence)
                    contains[j] = 1;
            for (ElementId j = 0; j < local_scheme_0_.number_of_elements(); ++j) {
                if (contains[j])
                    continue;
                Move move;
                move.type = Perturbations::ForceAdd;
                move.force_add_j = j;
                move.global_cost = global_cost(solution);
                moves.push_back(move);
            }
        }

        std::shuffle(moves.begin(), moves.end(), generator);
        return moves;
    }

    inline void apply_move(
            Solution& solution,
            const Move& move,
            std::mt19937_64& generator)
    {
        switch (move.type) {
        case Perturbations::None: {
            break;

        } case Perturbations::DoubleBridge: {
            std::uniform_int_distribution<SequenceId> distribution(
                    0, number_of_sequences() - 1);
            SequenceId i = distribution(generator);
            auto positions = optimizationtools::bob_floyd<ElementPos>(
                    4,
                    solution.sequences[i].sequence.size() + 1,
                    generator);
            std::sort(positions.begin(), positions.end());

            const Sequence& sequence = solution.sequences[i];
            Sequence& sequence_tmp = solution_tmp_.sequences[i];
            sequence_tmp = empty_sequence(i);
            for (ElementPos pos = 0; pos < positions[0]; ++pos)
                local_scheme_0_.append(sequence_tmp, sequence.sequence[pos]);
            for (ElementPos pos = positions[2]; pos < positions[3]; ++pos)
                local_scheme_0_.append(sequence_tmp, sequence.sequence[pos]);
            for (ElementPos pos = positions[1]; pos < positions[2]; ++pos)
                local_scheme_0_.append(sequence_tmp, sequence.sequence[pos]);
            for (ElementPos pos = positions[0]; pos < positions[1]; ++pos)
                local_scheme_0_.append(sequence_tmp, sequence.sequence[pos]);
            for (ElementPos pos = positions[3]; pos < (ElementPos)sequence.sequence.size(); ++pos)
                local_scheme_0_.append(sequence_tmp, sequence.sequence[pos]);
            compute_global_cost(solution_tmp_);
            solution = solution_tmp_;
            break;

        } case Perturbations::RuinAndRecreate: {
            ElementId n = local_scheme_0_.number_of_elements();
            std::vector<SequenceId> sequences(n, -1);
            std::vector<ElementPos> positions(n, -1);
            elements_.clear();
            for (SequenceId i = 0; i < number_of_sequences(); ++i) {
                const Sequence& sequence = solution.sequences[i];
                ElementPos seq_size = (ElementPos)sequence.sequence.size();
                for (ElementPos pos = 0; pos < seq_size; ++pos) {
                    ElementId j = sequence.sequence[pos];
                    elements_.add(j);
                    sequences[j] = i;
                    positions[j] = pos;
                }
            }

            ElementPos number_of_elements_removed = std::min(
                    parameters_.ruin_and_recreate_number_of_elements_removed,
                    elements_.size());
            elements_.shuffle_in(generator);
            while (elements_.size() > number_of_elements_removed) {
                ElementId j = *(elements_.begin() + (elements_.size() - 1));
                elements_.remove(j);
            }

            Solution solution_cur_ = empty_solution();
            for (SequenceId i = 0; i < number_of_sequences(); ++i) {
                const Sequence& sequence = solution.sequences[i];
                Sequence& sequence_cur = solution_cur_.sequences[i];
                for (ElementPos pos = 0; pos < (ElementPos)sequence.sequence.size(); ++pos) {
                    ElementId j = sequence.sequence[pos];
                    if (!elements_.contains(j))
                        local_scheme_0_.append(sequence_cur, j);
                }
            }
            // Add back removed elements.
            for (ElementId j: elements_) {
                std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                std::shuffle(positions_2_.begin(), positions_2_.end(), generator);
                GlobalCost c_best = global_cost(solution_cur_);
                SequenceId i_best = -1;
                ElementPos pos_best = -1;
                for (SequenceId i: sequences_1_) {
                    ElementPos seq_size = (ElementPos)solution_cur_.sequences[i].sequence.size();
                    compute_cost_add(solution_cur_, i, j);
                    for (ElementPos pos_new: positions_2_) {
                        if (pos_new > seq_size)
                            continue;
                        // Avoid inserting j at the same place as in the
                        // original solution.
                        if (i == sequences[j]) {
                            if (pos_new == 0 && positions[j] == 0)
                                continue;
                            if (pos_new != 0
                                    && positions[j] != 0
                                    && solution_cur_.sequences[i].sequence[pos_new - 1]
                                    == solution.sequences[i].sequence[positions[j] - 1])
                                continue;
                        }
                        GlobalCost c = global_costs_1d_[pos_new];
                        if (c >= c_best)
                            continue;
                        if (pos_best != -1 && !dominates(c, c_best))
                            continue;
                        i_best = i;
                        pos_best = pos_new;
                        c_best = c;
                    }
                }
                if (i_best != -1) {
                    for (SequenceId i = 0; i < number_of_sequences(); ++i)
                        if (i != i_best)
                            solution_tmp_.sequences[i] = solution_cur_.sequences[i];
                    const Sequence& sequence_cur = solution_cur_.sequences[i_best];
                    Sequence& sequence_tmp = solution_tmp_.sequences[i_best];
                    sequence_tmp = empty_sequence(i_best);
                    for (ElementPos p = 0; p < pos_best; ++p)
                        local_scheme_0_.append(sequence_tmp, sequence_cur.sequence[p]);
                    local_scheme_0_.append(sequence_tmp, j);
                    for (ElementPos p = pos_best; p < (ElementPos)sequence_cur.sequence.size(); ++p)
                        local_scheme_0_.append(sequence_tmp, sequence_cur.sequence[p]);
                    compute_global_cost(solution_tmp_);
                    solution_cur_ = solution_tmp_;
                }
            }
            solution = solution_cur_;
            break;

        } case Perturbations::ForceAdd: {
            ElementId j = move.force_add_j;
            GlobalCost c_best = global_cost(solution);
            std::mt19937_64 generator(0);
            std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
            std::shuffle(positions_2_.begin(), positions_2_.end(), generator);
            SequenceId i_best = -1;
            ElementPos pos_best = -1;
            for (SequenceId i: sequences_1_) {
                ElementPos seq_size = (ElementPos)solution.sequences[i].sequence.size();
                compute_cost_add(solution, i, j);
                for (ElementPos pos: positions_2_) {
                    if (pos > seq_size)
                        continue;
                    GlobalCost c = global_costs_1d_[pos];
                    if (pos_best != -1 && !dominates(c, c_best))
                        continue;
                    i_best = i;
                    pos_best = pos;
                    c_best = c;
                }
            }
            if (i_best != -1) {
                for (SequenceId i = 0; i < number_of_sequences(); ++i)
                    if (i != i_best)
                        solution_tmp_.sequences[i] = solution.sequences[i];
                const Sequence& sequence = solution.sequences[i_best];
                Sequence& sequence_tmp = solution_tmp_.sequences[i_best];
                sequence_tmp = empty_sequence(i_best);
                for (ElementPos p = 0; p < pos_best; ++p)
                    local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                local_scheme_0_.append(sequence_tmp, j);
                for (ElementPos p = pos_best; p < (ElementPos)sequence.sequence.size(); ++p)
                    local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                compute_global_cost(solution_tmp_);
                solution = solution_tmp_;
            }
            break;
        }
        }
    }

    /*
     * Local Search.
     */

    inline void local_search(
            Solution& solution,
            std::mt19937_64& generator,
            const Move& perturbation = Move())
    {
        //if (tabu.j != -1)
        //    std::cout << "j " << tabu.j << " j_prev " << tabu.j_prev
        //        << std::endl;
        //print(std::cout, solution);
        //std::cout << to_string(global_cost(solution)) << std::endl;

        SequenceId m = number_of_sequences();

        // Get neighborhoods.
        std::vector<std::tuple<Neighborhoods, ElementPos, ElementPos>> neighborhoods;
        for (ElementPos block_size = 1;
                block_size <= parameters_.shift_block_maximum_length;
                ++block_size) {
            neighborhoods.push_back({Neighborhoods::Shift, block_size, 0});
        }
        for (ElementPos block_size = 1;
                block_size <= parameters_.swap_block_maximum_length;
                ++block_size) {
            neighborhoods.push_back({Neighborhoods::SwapK, block_size, 0});
        }
        for (ElementPos block_size_1 = 1;
                block_size_1 <= parameters_.swap_block_maximum_length;
                ++block_size_1) {
            for (ElementPos block_size_2 = 1;
                    block_size_2 < block_size_1;
                    ++block_size_2) {
                neighborhoods.push_back({
                        Neighborhoods::SwapK1K2, block_size_1, block_size_2});
            }
        }
        if (parameters_.reverse)
            neighborhoods.push_back({Neighborhoods::Reverse, 0, 0});
        for (ElementPos block_size = 2;
                block_size <= parameters_.shift_reverse_block_maximum_length;
                ++block_size) {
            neighborhoods.push_back({
                    Neighborhoods::ShiftReverse, block_size, 0});
        }
        if (parameters_.add_remove) {
            neighborhoods.push_back({Neighborhoods::Add, 0, 0});
            neighborhoods.push_back({Neighborhoods::Remove, 0, 0});
        }
        if (m > 1) {
            if (parameters_.inter_two_opt)
                neighborhoods.push_back({Neighborhoods::InterTwoOpt, 0, 0});
            for (ElementPos block_size = 1;
                    block_size <= parameters_.inter_shift_block_maximum_length;
                    ++block_size) {
                neighborhoods.push_back({
                        Neighborhoods::InterShift, block_size, 0});
            }
            for (ElementPos block_size_1 = 1;
                    block_size_1 <= parameters_.inter_swap_block_maximum_length;
                    ++block_size_1) {
                for (ElementPos block_size_2 = 1;
                        block_size_2 <= block_size_1;
                        ++block_size_2) {
                    neighborhoods.push_back({
                            Neighborhoods::InterSwap,
                            block_size_1,
                            block_size_2});
                }
            }
        }

        Counter it = 0;
        for (;; ++it) {
            //std::cout << "it " << it
            //    << " c " << to_string(global_cost(solution))
            //    << std::endl;
            //print(std::cout, solution);

            if (parameters_.shuffle_neighborhood_order)
                std::shuffle(neighborhoods.begin(), neighborhoods.end(), generator);
            std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
            bool improved = false;
            // Loop through neighborhoods.
            for (auto neighborhood: neighborhoods) {
                switch (std::get<0>(neighborhood)) {

                case Neighborhoods::Shift: {
                    ElementPos block_size = std::get<1>(neighborhood);
                    std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                    std::shuffle(positions_1_.begin(), positions_1_.end(), generator);
                    std::shuffle(positions_2_.begin(), positions_2_.end(), generator);
                    SequenceId i_best = -1;
                    ElementPos pos_best = -1;
                    ElementPos pos_new_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (SequenceId i: sequences_1_) {
                        ElementPos seq_size = (ElementPos)solution.sequences[i].sequence.size();
                        if (seq_size <= block_size)
                            break;
                        for (ElementPos pos: positions_1_) {
                            if (pos > seq_size - block_size)
                                continue;
                            compute_cost_shift(solution, i, pos, block_size);
                            for (ElementPos pos_new: positions_2_) {
                                if (pos == pos_new || pos_new > seq_size - block_size)
                                    continue;
                                GlobalCost c = global_costs_1d_[pos_new];
                                if (c >= c_best)
                                    continue;
                                if (i_best != -1 && !dominates(c, c_best))
                                    continue;
                                i_best = i;
                                pos_best = pos;
                                pos_new_best = pos_new;
                                c_best = c;
                            }
                        }
                    }
                    // If an improving move has been found.
                    if (i_best != -1) {
                        //std::cout << "shift " << block_size
                        //    << " i_best " << i_best
                        //    << " pos_best " << pos_best
                        //    << " pos_new_best " << pos_new_best
                        //    << std::endl;

                        // Check that the move is improving.
                        if (c_best >= global_cost(solution)) {
                            throw std::logic_error(
                                    std::to_string(block_size) + "-shift."
                                    + " Best cost is worse than current cost:\n"
                                    + "* Best cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Current cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        improved = true;

                        // Apply best move.
                        for (SequenceId i = 0; i < m; ++i)
                            if (i != i_best)
                                solution_tmp_.sequences[i] = solution.sequences[i];
                        const Sequence& sequence = solution.sequences[i_best];
                        Sequence& sequence_tmp = solution_tmp_.sequences[i_best];
                        sequence_tmp = empty_sequence(i_best);
                        if (pos_best > pos_new_best) {
                            for (ElementPos p = 0; p < pos_new_best; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_best; p < pos_best + block_size; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_new_best; p < pos_best; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_best + block_size; p < (ElementPos)sequence.sequence.size(); ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                        } else {
                            for (ElementPos p = 0; p < pos_best; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_best + block_size; p < pos_new_best + block_size; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_best; p < pos_best + block_size; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_new_best + block_size; p < (ElementPos)sequence.sequence.size(); ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                        }
                        compute_global_cost(solution_tmp_);
                        solution = solution_tmp_;
                        // Check new current solution cost.
                        if (global_cost(solution) != c_best) {
                            throw std::logic_error(
                                    std::to_string(block_size)
                                    + "-shift. Costs do not match:\n"
                                    + "* Expected new cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Actual new cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        // Update statistics.
                        shift_number_of_sucesses_[block_size]++;
                    }
                    // Update statistics.
                    shift_number_of_explorations_[block_size]++;
                    break;

                } case Neighborhoods::SwapK: {
                    ElementPos block_size = std::get<1>(neighborhood);
                    std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                    std::shuffle(pairs_1_.begin(), pairs_1_.end(), generator);
                    SequenceId i_best = -1;
                    ElementPos pos_1_best = -1;
                    ElementPos pos_2_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (SequenceId i: sequences_1_) {
                        ElementPos seq_size = (ElementPos)solution.sequences[i].sequence.size();
                        if (seq_size <= block_size + block_size)
                            break;
                        compute_cost_swap(solution, i, block_size);
                        for (auto pair: pairs_1_) {
                            GlobalCost c = global_costs_2d_1_[pair.first][pair.second - pair.first - 1];
                            if (c >= c_best)
                                continue;
                            if (i_best != -1 && !dominates(c, c_best))
                                continue;
                            i_best = i;
                            pos_1_best = pair.first;
                            pos_2_best = pair.second;
                            c_best = c;
                        }
                    }
                    if (i_best != -1) {
                        //std::cout << "swap"
                        //    << " pos_1_best " << pos_1_best
                        //    << " pos_2_best " << pos_2_best
                        //    << std::endl;

                        // Check that the move is improving.
                        if (c_best >= global_cost(solution)) {
                            throw std::logic_error(
                                    std::to_string(block_size)
                                    + "-swap. Best cost is worse than current cost:\n"
                                    + "* Best cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Current cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        improved = true;

                        // Apply best move.
                        for (SequenceId i = 0; i < m; ++i)
                            if (i != i_best)
                                solution_tmp_.sequences[i] = solution.sequences[i];
                        const Sequence& sequence = solution.sequences[i_best];
                        Sequence& sequence_tmp = solution_tmp_.sequences[i_best];
                        sequence_tmp = empty_sequence(i_best);
                        for (ElementPos pos = 0; pos < (ElementPos)sequence.sequence.size(); ++pos) {
                            ElementId j = sequence.sequence[pos];
                            if (pos_1_best <= pos && pos < pos_1_best + block_size) {
                                ElementPos diff = pos - pos_1_best;
                                j = sequence.sequence[pos_2_best + diff];
                            }
                            if (pos_2_best <= pos && pos < pos_2_best + block_size) {
                                ElementPos diff = pos - pos_2_best;
                                j = sequence.sequence[pos_1_best + diff];
                            }
                            local_scheme_0_.append(sequence_tmp, j);
                        }
                        compute_global_cost(solution_tmp_);
                        solution = solution_tmp_;
                        // Check new current solution cost.
                        if (global_cost(solution) != c_best) {
                            throw std::logic_error(
                                    std::to_string(block_size)
                                    + "-swap. Costs do not match:\n"
                                    + "* Expected new cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Actual new cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        // Update statistics.
                        swap_number_of_sucesses_[block_size][block_size]++;
                    }
                    // Update statistics.
                    swap_number_of_explorations_[block_size][block_size]++;
                    break;

                } case Neighborhoods::SwapK1K2: {
                    ElementPos block_size_1 = std::get<1>(neighborhood);
                    ElementPos block_size_2 = std::get<2>(neighborhood);
                    std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                    std::shuffle(pairs_2_.begin(), pairs_2_.end(), generator);
                    SequenceId i_best = -1;
                    ElementPos pos_1_best = -1;
                    ElementPos pos_2_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (SequenceId i: sequences_1_) {
                        ElementPos seq_size = (ElementPos)solution.sequences[i].sequence.size();
                        if (seq_size <= block_size_1 + block_size_2)
                            break;
                        compute_cost_swap(solution, i, block_size_1, block_size_2);
                        for (auto pair: pairs_2_) {
                            GlobalCost c = global_costs_2d_2_[pair.first][pair.second];
                            if (c >= c_best)
                                continue;
                            if (i_best != -1 && !dominates(c, c_best))
                                continue;
                            i_best = i;
                            pos_1_best = pair.first;
                            pos_2_best = pair.second;
                            c_best = c;
                        }
                    }
                    if (i_best != -1) {
                        //std::cout << "swap "
                        //    << block_size_1 << " " << block_size_2
                        //    << " pos_1_best " << pos_1_best
                        //    << " pos_2_best " << pos_2_best
                        //    << " c_best " << to_string(c_best)
                        //    << std::endl;

                        // Check that the move is improving.
                        if (c_best >= global_cost(solution)) {
                            throw std::logic_error(
                                    "(" + std::to_string(block_size_1)
                                    + "," + std::to_string(block_size_2)
                                    + ")-swap. Best cost is worse than current cost:\n"
                                    + "* Best cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Current cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        improved = true;

                        // Apply best move.
                        for (SequenceId i = 0; i < m; ++i)
                            if (i != i_best)
                                solution_tmp_.sequences[i] = solution.sequences[i];
                        const Sequence& sequence = solution.sequences[i_best];
                        Sequence& sequence_tmp = solution_tmp_.sequences[i_best];
                        sequence_tmp = empty_sequence(i_best);
                        if (pos_1_best < pos_2_best) {
                            for (ElementPos p = 0; p < pos_1_best; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_2_best; p < pos_2_best + block_size_2; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_1_best + block_size_1; p < pos_2_best; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_1_best; p < pos_1_best + block_size_1; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_2_best + block_size_2; p < (ElementPos)sequence.sequence.size(); ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                        } else {
                            for (ElementPos p = 0; p < pos_2_best; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_1_best; p < pos_1_best + block_size_1; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_2_best + block_size_2; p < pos_1_best; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_2_best; p < pos_2_best + block_size_2; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_1_best + block_size_1; p < (ElementPos)sequence.sequence.size(); ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                        }
                        compute_global_cost(solution_tmp_);
                        solution = solution_tmp_;
                        // Check new current solution cost.
                        if (global_cost(solution) != c_best) {
                            throw std::logic_error(
                                    "(" + std::to_string(block_size_1)
                                    + "," + std::to_string(block_size_2)
                                    + ")-swap. Costs do not match:\n"
                                    + "* Expected new cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Actual new cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        // Update statistics.
                        swap_number_of_sucesses_[block_size_1][block_size_2]++;
                    }
                    // Update statistics.
                    swap_number_of_explorations_[block_size_1][block_size_2]++;
                    break;

                } case Neighborhoods::Reverse: {
                    std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                    std::shuffle(pairs_1_.begin(), pairs_1_.end(), generator);
                    SequenceId i_best = -1;
                    ElementPos pos_1_best = -1;
                    ElementPos pos_2_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (SequenceId i: sequences_1_) {
                        compute_cost_reverse(solution, i);
                        for (auto pair: pairs_1_) {
                            GlobalCost c = global_costs_2d_1_[pair.first][pair.second - pair.first - 1];
                            if (c >= c_best)
                                continue;
                            if (pos_1_best != -1 && !dominates(c, c_best))
                                continue;
                            i_best = i;
                            pos_1_best = pair.first;
                            pos_2_best = pair.second;
                            c_best = c;
                        }
                    }
                    if (i_best != -1) {
                        //std::cout << "reverse"
                        //    << " i_best " << i_best
                        //    << " pos_1_best " << pos_1_best
                        //    << " pos_2_best " << pos_2_best
                        //    << std::endl;

                        // Check that the move is improving.
                        if (c_best >= global_cost(solution)) {
                            throw std::logic_error(
                                    "Reverse. Best cost is worse than current cost:\n"
                                    "* Best cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Current cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        improved = true;

                        // Apply best move.
                        for (SequenceId i = 0; i < m; ++i)
                            if (i != i_best)
                                solution_tmp_.sequences[i] = solution.sequences[i];
                        const Sequence& sequence = solution.sequences[i_best];
                        Sequence& sequence_tmp = solution_tmp_.sequences[i_best];
                        sequence_tmp = empty_sequence(i_best);
                        for (ElementPos pos = 0; pos < pos_1_best; ++pos) {
                            ElementId j = sequence.sequence[pos];
                            local_scheme_0_.append(sequence_tmp, j);
                        }
                        for (ElementPos pos = pos_2_best; pos >= pos_1_best; --pos) {
                            ElementId j = sequence.sequence[pos];
                            local_scheme_0_.append(sequence_tmp, j);
                        }
                        for (ElementPos pos = pos_2_best + 1; pos < (ElementPos)sequence.sequence.size(); ++pos) {
                            ElementId j = sequence.sequence[pos];
                            local_scheme_0_.append(sequence_tmp, j);
                        }
                        compute_global_cost(solution_tmp_);
                        solution = solution_tmp_;
                        // Check new current solution cost.
                        if (global_cost(solution) != c_best) {
                            throw std::logic_error(
                                    "Reverse. Costs do not match:\n"
                                    "* Expected new cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Actual new cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        // Update statistics.
                        reverse_number_of_sucesses_++;
                    }
                    // Update statistics.
                    reverse_number_of_explorations_++;
                    break;

                } case Neighborhoods::ShiftReverse: {
                    ElementPos block_size = std::get<1>(neighborhood);
                    std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                    std::shuffle(positions_1_.begin(), positions_1_.end(), generator);
                    std::shuffle(positions_2_.begin(), positions_2_.end(), generator);
                    SequenceId i_best = -1;
                    ElementPos pos_best = -1;
                    ElementPos pos_new_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (SequenceId i: sequences_1_) {
                        ElementPos seq_size = (ElementPos)solution.sequences[i].sequence.size();
                        if (seq_size <= block_size)
                            break;
                        for (ElementPos pos: positions_1_) {
                            if (pos > seq_size - block_size)
                                continue;
                            compute_cost_shift(solution, i, pos, block_size, true);
                            for (ElementPos pos_new: positions_2_) {
                                if (pos == pos_new || pos_new > seq_size - block_size)
                                    continue;
                                GlobalCost c = global_costs_1d_[pos_new];
                                if (c >= c_best)
                                    continue;
                                if (pos_best != -1 && !dominates(c, c_best))
                                    continue;
                                i_best = i;
                                pos_best = pos;
                                pos_new_best = pos_new;
                                c_best = c;
                            }
                        }
                    }
                    // If an improving move has been found.
                    if (i_best != -1) {
                        //std::cout << "shift " << neighborhood
                        //    << " pos_best " << pos_best
                        //    << " pos_new_best " << pos_new_best
                        //    << std::endl;

                        // Check that the move is improving.
                        if (c_best >= global_cost(solution)) {
                            throw std::logic_error(
                                    std::to_string(block_size)
                                    + "-shift-reverse. Best cost is worse than current cost:\n"
                                    + "* Best cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Current cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        improved = true;

                        // Apply best move.
                        for (SequenceId i = 0; i < m; ++i)
                            if (i != i_best)
                                solution_tmp_.sequences[i] = solution.sequences[i];
                        const Sequence& sequence = solution.sequences[i_best];
                        Sequence& sequence_tmp = solution_tmp_.sequences[i_best];
                        sequence_tmp = empty_sequence(i_best);
                        if (pos_best > pos_new_best) {
                            for (ElementPos p = 0; p < pos_new_best; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_best + block_size - 1; p >= pos_best; --p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_new_best; p < pos_best; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_best + block_size; p < (ElementPos)sequence.sequence.size(); ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                        } else {
                            for (ElementPos p = 0; p < pos_best; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_best + block_size; p < pos_new_best + block_size; ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_best + block_size - 1; p >= pos_best; --p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                            for (ElementPos p = pos_new_best + block_size; p < (ElementPos)sequence.sequence.size(); ++p)
                                local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                        }
                        compute_global_cost(solution_tmp_);
                        solution = solution_tmp_;
                        // Check new current solution cost.
                        if (global_cost(solution) != c_best) {
                            throw std::logic_error(
                                    std::to_string(block_size)
                                    + "-shift-reverse. Costs do not match:\n"
                                    + "* Expected new cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Actual new cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        // Update statistics.
                        shift_reverse_number_of_sucesses_[block_size]++;
                    }
                    // Update statistics.
                    shift_reverse_number_of_explorations_[block_size]++;
                    break;

                } case Neighborhoods::Add: {
                    std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                    std::shuffle(positions_2_.begin(), positions_2_.end(), generator);
                    SequenceId i_best = -1;
                    ElementId j_best = -1;
                    ElementPos pos_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    elements_.clear();
                    for (SequenceId i = 0; i < m; ++i)
                        for (ElementId j: solution.sequences[i].sequence)
                            elements_.add(j);
                    elements_.shuffle_out(generator);
                    for (SequenceId i: sequences_1_) {
                        ElementPos seq_size = (ElementPos)solution.sequences[i].sequence.size();
                        for (auto it = elements_.out_begin();
                                it != elements_.out_end(); ++it) {
                            ElementId j = *it;
                            compute_cost_add(solution, i, j);
                            for (ElementPos pos: positions_2_) {
                                if (pos > seq_size)
                                    continue;
                                GlobalCost c = global_costs_1d_[pos];
                                if (c >= c_best)
                                    continue;
                                if (pos_best != -1 && !dominates(c, c_best))
                                    continue;
                                i_best = i;
                                j_best = j;
                                pos_best = pos;
                                c_best = c;
                            }
                        }
                    }
                    if (i_best != -1) {
                        //std::cout << "add "
                        //    << "j_best " << j_best
                        //    << " pos_best " << pos_best
                        //    << std::endl;

                        // Check that the move is improving.
                        if (c_best >= global_cost(solution)) {
                            throw std::logic_error(
                                    "Add. Best cost is worse than current cost:\n"
                                    "* Best cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Current cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        improved = true;

                        // Apply best move.
                        for (SequenceId i = 0; i < m; ++i)
                            if (i != i_best)
                                solution_tmp_.sequences[i] = solution.sequences[i];
                        const Sequence& sequence = solution.sequences[i_best];
                        Sequence& sequence_tmp = solution_tmp_.sequences[i_best];
                        sequence_tmp = empty_sequence(i_best);
                        for (ElementPos p = 0; p < pos_best; ++p)
                            local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                        local_scheme_0_.append(sequence_tmp, j_best);
                        for (ElementPos p = pos_best; p < (ElementPos)sequence.sequence.size(); ++p)
                            local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                        compute_global_cost(solution_tmp_);
                        solution = solution_tmp_;
                        // Check new current solution cost.
                        if (global_cost(solution) != c_best) {
                            throw std::logic_error(
                                    "Add. Costs do not match:\n"
                                    "* Expected new cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Actual new cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        // Update statistics.
                        add_number_of_sucesses_++;
                    }
                    // Update statistics.
                    add_number_of_explorations_++;
                    break;

                } case Neighborhoods::Remove: {
                    std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                    std::shuffle(positions_1_.begin(), positions_1_.end(), generator);
                    SequenceId i_best = -1;
                    ElementPos pos_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (SequenceId i: sequences_1_) {
                        compute_cost_remove(solution, i);
                        const Sequence& sequence = solution.sequences[i];
                        ElementPos seq_size = (ElementPos)sequence.sequence.size();
                        for (ElementPos pos: positions_1_) {
                            if (pos >= seq_size)
                                continue;
                            if (perturbation.type == Perturbations::ForceAdd
                                    && perturbation.force_add_j == sequence.sequence[pos])
                                continue;
                            GlobalCost c = global_costs_1d_[pos];
                            if (c >= c_best)
                                continue;
                            if (pos_best != -1 && !dominates(c, c_best))
                                continue;
                            i_best = i;
                            pos_best = pos;
                            c_best = c;
                        }
                    }
                    if (i_best != -1) {
                        //std::cout << "remove "
                        //    << " pos_best " << pos_best
                        //    << std::endl;

                        // Check that the move is improving.
                        if (c_best >= global_cost(solution)) {
                            throw std::logic_error(
                                    "Remove. Best cost is worse than current cost:\n"
                                    "* Best cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Current cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        improved = true;

                        // Apply best move.
                        for (SequenceId i = 0; i < m; ++i)
                            if (i != i_best)
                                solution_tmp_.sequences[i] = solution.sequences[i];
                        const Sequence& sequence = solution.sequences[i_best];
                        Sequence& sequence_tmp = solution_tmp_.sequences[i_best];
                        sequence_tmp = empty_sequence(i_best);
                        for (ElementPos p = 0; p < (ElementPos)sequence.sequence.size(); ++p) {
                            if (p == pos_best)
                                continue;
                            local_scheme_0_.append(sequence_tmp, sequence.sequence[p]);
                        }
                        compute_global_cost(solution_tmp_);
                        solution = solution_tmp_;
                        // Check new current solution cost.
                        if (global_cost(solution) != c_best) {
                            throw std::logic_error(
                                    "Remove. Costs do not match:\n"
                                    "* Expected new cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Actual new cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        // Update statistics.
                        remove_number_of_sucesses_++;
                    }
                    // Update statistics.
                    remove_number_of_explorations_++;
                    break;

                } case Neighborhoods::InterTwoOpt: {
                    std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                    std::shuffle(sequences_2_.begin(), sequences_2_.end(), generator);
                    std::shuffle(pairs_2_.begin(), pairs_2_.end(), generator);
                    SequenceId i1_best = -1;
                    SequenceId i2_best = -1;
                    ElementPos pos_1_best = -1;
                    ElementPos pos_2_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (SequenceId i1: sequences_1_) {
                        for (SequenceId i2: sequences_2_) {
                            if (i2 <= i1)
                                continue;
                            compute_cost_inter_two_opt(solution, i1, i2);
                            for (auto pair: pairs_2_) {
                                GlobalCost c = global_costs_2d_2_[pair.first][pair.second];
                                if (c >= c_best)
                                    continue;
                                if (i1_best != -1 && !dominates(c, c_best))
                                    continue;
                                i1_best = i1;
                                i2_best = i2;
                                pos_1_best = pair.first;
                                pos_2_best = pair.second;
                                c_best = c;
                            }
                        }
                    }
                    if (i1_best != -1) {
                        //std::cout << "inter two-opt"
                        //    << " i1_best " << i1_best
                        //    << " i2_best " << i2_best
                        //    << " pos_1_best " << pos_1_best
                        //    << " pos_2_best " << pos_2_best
                        //    << std::endl;

                        // Check that the move is improving.
                        if (c_best >= global_cost(solution)) {
                            throw std::logic_error(
                                    + "Inter 2-opt."
                                    " Best cost is worse than current cost:\n"
                                    "* Best cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Current cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        improved = true;

                        // Apply best move.
                        for (SequenceId i = 0; i < m; ++i)
                            if (i != i1_best && i != i2_best)
                                solution_tmp_.sequences[i] = solution.sequences[i];
                        const Sequence& sequence_1 = solution.sequences[i1_best];
                        const Sequence& sequence_2 = solution.sequences[i2_best];
                        // Sequence 1.
                        Sequence& sequence_tmp_1 = solution_tmp_.sequences[i1_best];
                        sequence_tmp_1 = empty_sequence(i1_best);
                        for (ElementPos p = 0; p < pos_1_best; ++p)
                            local_scheme_0_.append(sequence_tmp_1, sequence_1.sequence[p]);
                        for (ElementPos p = pos_2_best; p < (ElementPos)sequence_2.sequence.size(); ++p)
                            local_scheme_0_.append(sequence_tmp_1, sequence_2.sequence[p]);
                        // Sequence 2.
                        Sequence& sequence_tmp_2 = solution_tmp_.sequences[i2_best];
                        sequence_tmp_2 = empty_sequence(i2_best);
                        for (ElementPos p = 0; p < pos_2_best; ++p)
                            local_scheme_0_.append(sequence_tmp_2, sequence_2.sequence[p]);
                        for (ElementPos p = pos_1_best; p < (ElementPos)sequence_1.sequence.size(); ++p)
                            local_scheme_0_.append(sequence_tmp_2, sequence_1.sequence[p]);
                        compute_global_cost(solution_tmp_);
                        solution = solution_tmp_;
                        // Check new current solution cost.
                        if (global_cost(solution) != c_best) {
                            throw std::logic_error(
                                    + "Inter 2-opt."
                                    " Costs do not match:\n"
                                    "* Expected new cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Actual new cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        // Update statistics.
                        inter_two_opt_number_of_sucesses_++;
                    }
                    // Update statistics.
                    inter_two_opt_number_of_explorations_++;
                    break;

                } case Neighborhoods::InterShift: {
                    ElementPos block_size = std::get<1>(neighborhood);
                    std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                    std::shuffle(sequences_2_.begin(), sequences_2_.end(), generator);
                    std::shuffle(pairs_2_.begin(), pairs_2_.end(), generator);
                    SequenceId i1_best = -1;
                    SequenceId i2_best = -1;
                    ElementPos pos_1_best = -1;
                    ElementPos pos_2_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (SequenceId i1: sequences_1_) {
                        ElementPos seq_1_size = (ElementPos)solution.sequences[i1].sequence.size();
                        if (seq_1_size < block_size)
                            break;
                        for (SequenceId i2: sequences_2_) {
                            if (i2 == i1)
                                continue;
                            compute_cost_inter_shift(solution, i1, i2, block_size);
                            for (auto pair: pairs_2_) {
                                GlobalCost c = global_costs_2d_2_[pair.first][pair.second];
                                //std::cout
                                //    << "i1 " << i1
                                //    << " i2 " << i2
                                //    << " " << pair.first
                                //    << " " << pair.second
                                //    << " " << to_string(c)
                                //    << std::endl;
                                if (c >= c_best)
                                    continue;
                                if (i1_best != -1 && !dominates(c, c_best))
                                    continue;
                                i1_best = i1;
                                i2_best = i2;
                                pos_1_best = pair.first;
                                pos_2_best = pair.second;
                                c_best = c;
                            }
                        }
                    }
                    if (i1_best != -1) {
                        //std::cout << "inter-"
                        //    << block_size << "-shift"
                        //    << " i1_best " << i1_best
                        //    << " i2_best " << i2_best
                        //    << " pos_1_best " << pos_1_best
                        //    << " pos_2_best " << pos_2_best
                        //    << std::endl;

                        // Check that the move is improving.
                        if (c_best >= global_cost(solution)) {
                            throw std::logic_error(
                                    "Inter-" + std::to_string(block_size)
                                    + "-shift."
                                    + " Best cost is worse than current cost:\n"
                                    + "* Best cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Current cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        improved = true;

                        // Apply best move.
                        for (SequenceId i = 0; i < m; ++i)
                            if (i != i1_best && i != i2_best)
                                solution_tmp_.sequences[i] = solution.sequences[i];
                        const Sequence& sequence_1 = solution.sequences[i1_best];
                        const Sequence& sequence_2 = solution.sequences[i2_best];
                        // Sequence 1.
                        Sequence& sequence_tmp_1 = solution_tmp_.sequences[i1_best];
                        sequence_tmp_1 = empty_sequence(i1_best);
                        for (ElementPos p = 0; p < pos_1_best; ++p)
                            local_scheme_0_.append(sequence_tmp_1, sequence_1.sequence[p]);
                        for (ElementPos p = pos_1_best + block_size; p < (ElementPos)sequence_1.sequence.size(); ++p)
                            local_scheme_0_.append(sequence_tmp_1, sequence_1.sequence[p]);
                        // Sequence 2.
                        Sequence& sequence_tmp_2 = solution_tmp_.sequences[i2_best];
                        sequence_tmp_2 = empty_sequence(i2_best);
                        for (ElementPos p = 0; p < pos_2_best; ++p)
                            local_scheme_0_.append(sequence_tmp_2, sequence_2.sequence[p]);
                        for (ElementPos p = pos_1_best; p < pos_1_best + block_size; ++p)
                            local_scheme_0_.append(sequence_tmp_2, sequence_1.sequence[p]);
                        for (ElementPos p = pos_2_best; p < (ElementPos)sequence_2.sequence.size(); ++p)
                            local_scheme_0_.append(sequence_tmp_2, sequence_2.sequence[p]);
                        compute_global_cost(solution_tmp_);
                        // Check new current solution size.
                        if (sequence_tmp_1.sequence.size()
                               + sequence_tmp_2.sequence.size()
                               != sequence_1.sequence.size()
                               + sequence_2.sequence.size()) {
                            throw std::logic_error(
                                    "Inter-" + std::to_string(block_size)
                                    + "-shift."
                                    + " Sizes do not match:\n"
                                    + "* Old sequences: "
                                    + std::to_string(sequence_1.sequence.size())
                                    + " "
                                    + std::to_string(sequence_2.sequence.size())
                                    + "\n"
                                    + "* New sequences: "
                                    + std::to_string(sequence_tmp_1.sequence.size())
                                    + " "
                                    + std::to_string(sequence_tmp_2.sequence.size())
                                    + "\n");
                        }
                        solution = solution_tmp_;
                        // Check new current solution cost.
                        if (global_cost(solution) != c_best) {
                            throw std::logic_error(
                                    "Inter-" + std::to_string(block_size)
                                    + "-shift."
                                    + " Costs do not match:\n"
                                    + "* Expected new cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Actual new cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        // Update statistics.
                        inter_shift_number_of_sucesses_[block_size]++;
                    }
                    // Update statistics.
                    inter_shift_number_of_explorations_[block_size]++;
                    break;

                } case Neighborhoods::InterSwap: {
                    ElementPos block_size_1 = std::get<1>(neighborhood);
                    ElementPos block_size_2 = std::get<2>(neighborhood);
                    std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                    std::shuffle(sequences_2_.begin(), sequences_2_.end(), generator);
                    std::shuffle(pairs_2_.begin(), pairs_2_.end(), generator);
                    SequenceId i1_best = -1;
                    SequenceId i2_best = -1;
                    ElementPos pos_1_best = -1;
                    ElementPos pos_2_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (SequenceId i1: sequences_1_) {
                        ElementPos seq_1_size = (ElementPos)solution.sequences[i1].sequence.size();
                        if (seq_1_size < block_size_1)
                            break;
                        for (SequenceId i2: sequences_2_) {
                            ElementPos seq_2_size = (ElementPos)solution.sequences[i2].sequence.size();
                            if (seq_2_size < block_size_2)
                                break;
                            if (i2 == i1)
                                continue;
                            if (block_size_1 == block_size_2 && i2 <= i1)
                                continue;
                            compute_cost_inter_swap(solution, i1, i2, block_size_1, block_size_2);
                            for (auto pair: pairs_2_) {
                                GlobalCost c = global_costs_2d_2_[pair.first][pair.second];
                                if (c >= c_best)
                                    continue;
                                if (i1_best != -1 && !dominates(c, c_best))
                                    continue;
                                i1_best = i1;
                                i2_best = i2;
                                pos_1_best = pair.first;
                                pos_2_best = pair.second;
                                c_best = c;
                            }
                        }
                    }
                    if (i1_best != -1) {
                        //std::cout << "inter-("
                        //    << block_size_1 << "," << block_size_2 << ")-swap"
                        //    << " i1_best " << i1_best
                        //    << " i2_best " << i2_best
                        //    << " pos_1_best " << pos_1_best
                        //    << " pos_2_best " << pos_2_best
                        //    << std::endl;

                        // Check that the move is improving.
                        if (c_best >= global_cost(solution)) {
                            throw std::logic_error(
                                    "Inter-" + std::to_string(block_size_1)
                                    + "," + std::to_string(block_size_2)
                                    + "-swap."
                                    + " Best cost is worse than current cost:\n"
                                    + "* Best cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Current cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        improved = true;

                        // Apply best move.
                        for (SequenceId i = 0; i < m; ++i)
                            if (i != i1_best && i != i2_best)
                                solution_tmp_.sequences[i] = solution.sequences[i];
                        const Sequence& sequence_1 = solution.sequences[i1_best];
                        const Sequence& sequence_2 = solution.sequences[i2_best];
                        // Sequence 1.
                        Sequence& sequence_tmp_1 = solution_tmp_.sequences[i1_best];
                        sequence_tmp_1 = empty_sequence(i1_best);
                        for (ElementPos p = 0; p < pos_1_best; ++p)
                            local_scheme_0_.append(sequence_tmp_1, sequence_1.sequence[p]);
                        for (ElementPos p = pos_2_best; p < pos_2_best + block_size_2; ++p)
                            local_scheme_0_.append(sequence_tmp_1, sequence_2.sequence[p]);
                        for (ElementPos p = pos_1_best + block_size_1; p < (ElementPos)sequence_1.sequence.size(); ++p)
                            local_scheme_0_.append(sequence_tmp_1, sequence_1.sequence[p]);
                        // Sequence 2.
                        Sequence& sequence_tmp_2 = solution_tmp_.sequences[i2_best];
                        sequence_tmp_2 = empty_sequence(i2_best);
                        for (ElementPos p = 0; p < pos_2_best; ++p)
                            local_scheme_0_.append(sequence_tmp_2, sequence_2.sequence[p]);
                        for (ElementPos p = pos_1_best; p < pos_1_best + block_size_1; ++p)
                            local_scheme_0_.append(sequence_tmp_2, sequence_1.sequence[p]);
                        for (ElementPos p = pos_2_best + block_size_2; p < (ElementPos)sequence_2.sequence.size(); ++p)
                            local_scheme_0_.append(sequence_tmp_2, sequence_2.sequence[p]);
                        compute_global_cost(solution_tmp_);
                        // Check new current solution size.
                        if (sequence_tmp_1.sequence.size()
                               + sequence_tmp_2.sequence.size()
                               != sequence_1.sequence.size()
                               + sequence_2.sequence.size()) {
                            throw std::logic_error(
                                    "Inter-(" + std::to_string(block_size_1)
                                    + "," + std::to_string(block_size_2)
                                    + ")-swap."
                                    + " Sizes do not match:\n"
                                    + "* Old sequences: "
                                    + std::to_string(sequence_1.sequence.size())
                                    + " "
                                    + std::to_string(sequence_2.sequence.size())
                                    + "\n"
                                    + "* New sequences: "
                                    + std::to_string(sequence_tmp_1.sequence.size())
                                    + " "
                                    + std::to_string(sequence_tmp_2.sequence.size())
                                    + "\n");
                        }
                        solution = solution_tmp_;
                        // Check new current solution cost.
                        if (global_cost(solution) != c_best) {
                            throw std::logic_error(
                                    "Inter-(" + std::to_string(block_size_1)
                                    + "," + std::to_string(block_size_2)
                                    + ")-swap."
                                    + " Costs do not match:\n"
                                    + "* Expected new cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Actual new cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        // Update statistics.
                        inter_swap_number_of_sucesses_[block_size_1][block_size_2]++;

                    }
                    // Update statistics.
                    inter_swap_number_of_explorations_[block_size_1][block_size_2]++;
                    break;

                } case Neighborhoods::InterShiftReverse: {
                    ElementPos block_size = std::get<1>(neighborhood);
                    std::shuffle(sequences_1_.begin(), sequences_1_.end(), generator);
                    std::shuffle(sequences_2_.begin(), sequences_2_.end(), generator);
                    std::shuffle(pairs_2_.begin(), pairs_2_.end(), generator);
                    SequenceId i1_best = -1;
                    SequenceId i2_best = -1;
                    ElementPos pos_1_best = -1;
                    ElementPos pos_2_best = -1;
                    GlobalCost c_best = global_cost(solution);
                    for (SequenceId i1: sequences_1_) {
                        ElementPos seq_1_size = (ElementPos)solution.sequences[i1].sequence.size();
                        if (seq_1_size < block_size)
                            break;
                        for (SequenceId i2: sequences_2_) {
                            if (i2 == i1)
                                continue;
                            compute_cost_inter_shift(
                                    solution, i1, i2, block_size, true);
                            for (auto pair: pairs_2_) {
                                GlobalCost c = global_costs_2d_2_[pair.first][pair.second];
                                if (c >= c_best)
                                    continue;
                                if (i1_best != -1 && !dominates(c, c_best))
                                    continue;
                                i1_best = i1;
                                i2_best = i2;
                                pos_1_best = pair.first;
                                pos_2_best = pair.second;
                                c_best = c;
                            }
                        }
                    }
                    if (i1_best != -1) {
                        //std::cout << "inter-"
                        //    << block_size << "-shift-reverse"
                        //    << " i1_best " << i1_best
                        //    << " i2_best " << i2_best
                        //    << " pos_1_best " << pos_1_best
                        //    << " pos_2_best " << pos_2_best
                        //    << std::endl;

                        // Check that the move is improving.
                        if (c_best >= global_cost(solution)) {
                            throw std::logic_error(
                                    "Inter-" + std::to_string(block_size)
                                    + "-shift-reverse."
                                    + " Best cost is worse than current cost:\n"
                                    + "* Best cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Current cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        improved = true;

                        // Apply best move.
                        for (SequenceId i = 0; i < m; ++i)
                            if (i != i1_best && i != i2_best)
                                solution_tmp_.sequences[i] = solution.sequences[i];
                        const Sequence& sequence_1 = solution.sequences[i1_best];
                        const Sequence& sequence_2 = solution.sequences[i2_best];
                        // Sequence 1.
                        Sequence& sequence_tmp_1 = solution_tmp_.sequences[i1_best];
                        sequence_tmp_1 = empty_sequence(i1_best);
                        for (ElementPos p = 0; p < pos_1_best; ++p)
                            local_scheme_0_.append(sequence_tmp_1, sequence_1.sequence[p]);
                        for (ElementPos p = pos_1_best + block_size; p < (ElementPos)sequence_1.sequence.size(); ++p)
                            local_scheme_0_.append(sequence_tmp_1, sequence_1.sequence[p]);
                        // Sequence 2.
                        Sequence& sequence_tmp_2 = solution_tmp_.sequences[i2_best];
                        sequence_tmp_2 = empty_sequence(i2_best);
                        for (ElementPos p = 0; p < pos_2_best; ++p)
                            local_scheme_0_.append(sequence_tmp_2, sequence_2.sequence[p]);
                        for (ElementPos p = pos_1_best + block_size - 1; p >= pos_1_best; --p)
                            local_scheme_0_.append(sequence_tmp_2, sequence_1.sequence[p]);
                        for (ElementPos p = pos_2_best; p < (ElementPos)sequence_2.sequence.size(); ++p)
                            local_scheme_0_.append(sequence_tmp_2, sequence_2.sequence[p]);
                        compute_global_cost(solution_tmp_);
                        solution = solution_tmp_;
                        // Check new current solution cost.
                        if (global_cost(solution) != c_best) {
                            throw std::logic_error(
                                    "Inter-" + std::to_string(block_size)
                                    + "-shift-reverse."
                                    + " Costs do not match:\n"
                                    + "* Expected new cost: "
                                    + to_string(c_best) + "\n"
                                    + "* Actual new cost: "
                                    + to_string(global_cost(solution)) + "\n");
                        }
                        // Update statistics.
                        inter_shift_reverse_number_of_sucesses_[block_size]++;
                    }
                    // Update statistics.
                    inter_shift_reverse_number_of_explorations_[block_size]++;
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
            const Solution& solution) const
    {
        for (SequenceId i = 0; i < number_of_sequences(); ++i) {
            os << "Sequence " << i << ":";
            for (ElementId j: solution.sequences[i].sequence)
                os << " " << j;
            os << std::endl;
            os << "    Cost: " << to_string(global_cost(solution.sequences[i])) << std::endl;
        }
        os << "Total cost: " << to_string(solution.global_cost) << std::endl;
        return os;
    }

    inline void write(
            const Solution& solution,
            std::string certificate_path) const
    {
        if (certificate_path.empty())
            return;
        std::ofstream cert(certificate_path);
        if (!cert.good()) {
            throw std::runtime_error(
                    "Unable to open file \"" + certificate_path + "\".");
        }

        SequenceId m = number_of_sequences();
        if (m == 1) {
            for (ElementId j: solution.sequences[0].sequence)
                cert << j << " ";
        } else {
            for (SequenceId i = 0; i < m; ++i) {
                cert << solution.sequences[i].sequence.size() << std::endl;
                for (ElementId j: solution.sequences[i].sequence)
                    cert << " " << j;
                cert << std::endl;
            }
        }
    }

    void print_parameters(
            optimizationtools::Info& info) const
    {
        SequencePos m = number_of_sequences();
        ElementPos n = local_scheme_0_.number_of_elements();
        FFOT_VER(info, ""
                << "Number of sequences:                           " << m << std::endl
                << "Number of elements:                            " << n << std::endl
                << "Neighborhoods" << std::endl
                << "    Shift" << std::endl
                << "        Block maximum length:                  " << parameters_.shift_block_maximum_length << std::endl
                << "    Swap" << std::endl
                << "        Block maximum length:                  " << parameters_.swap_block_maximum_length << std::endl
                << "    Reverse:                                   " << parameters_.reverse << std::endl
                << "    Shift-reverse" << std::endl
                << "        Block maximum length:                  " << parameters_.shift_reverse_block_maximum_length << std::endl
                << "    Add/Remove:                                " << parameters_.add_remove << std::endl
                );
        if (m > 1) {
            FFOT_VER(info, ""
                << "    Inter-shift" << std::endl
                << "        Block maximum length:                  " << parameters_.inter_shift_block_maximum_length << std::endl
                << "    Inter-swap" << std::endl
                << "        Block maximum length:                  " << parameters_.inter_swap_block_maximum_length << std::endl
                << "    Inter-two-opt:                             " << parameters_.inter_two_opt << std::endl
                << "    Inter-shift-reverse" << std::endl
                << "        Block maximum length:                  " << parameters_.inter_shift_reverse_block_maximum_length << std::endl
                );
        }
        FFOT_VER(info, ""
                << "Perturbations" << std::endl
                << "    Double-bridge" << std::endl
                << "        Number of perturbations:               " << parameters_.double_bridge_number_of_perturbations << std::endl
                << "    Ruin-and-recreate" << std::endl
                << "        Number of perturbations:               " << parameters_.ruin_and_recreate_number_of_perturbations << std::endl
                << "        Number of elements removed:            " << parameters_.ruin_and_recreate_number_of_elements_removed << std::endl
                << "    Force-add:                                 " << parameters_.force_add << std::endl
                );
    }

    void print_statistics(
            optimizationtools::Info& info) const
    {
        for (ElementPos block_size = 1; block_size <= parameters_.shift_block_maximum_length; ++block_size) {
            FFOT_VER(info,
                    std::left << std::setw(28) << ("Shift " + std::to_string(block_size) + ":")
                    << shift_number_of_explorations_[block_size]
                    << " / " << shift_number_of_sucesses_[block_size]
                    << " / " << (double)shift_number_of_sucesses_[block_size] / shift_number_of_explorations_[block_size] * 100 << "%"
                    << std::endl);
            FFOT_PUT(info,
                    "Algorithm", ("Shift" + std::to_string(block_size) + "NumberOfExplorations"),
                    shift_number_of_explorations_[block_size]);
            FFOT_PUT(info,
                    "Algorithm", ("Shift" + std::to_string(block_size) + "NumberOfSuccesses"),
                    shift_number_of_explorations_[block_size]);
        }
        for (ElementPos block_size_1 = 1; block_size_1 <= parameters_.swap_block_maximum_length; ++block_size_1) {
            for (ElementPos block_size_2 = 1; block_size_2 <= block_size_1; ++block_size_2) {
                FFOT_VER(info,
                        std::left << std::setw(28) << ("Swap " + std::to_string(block_size_1) + "," + std::to_string(block_size_2) + ":")
                        << swap_number_of_explorations_[block_size_1][block_size_2]
                        << " / " << swap_number_of_sucesses_[block_size_1][block_size_2]
                        << " / " << (double)swap_number_of_sucesses_[block_size_1][block_size_2] / swap_number_of_explorations_[block_size_1][block_size_2] * 100 << "%"
                        << std::endl);
                FFOT_PUT(info,
                        "Algorithm", ("Swap" + std::to_string(block_size_1) + "," + std::to_string(block_size_2) + "NumberOfExplorations"),
                        swap_number_of_explorations_[block_size_1][block_size_2]);
                FFOT_PUT(info,
                        "Algorithm", ("Shift" + std::to_string(block_size_1) + "," + std::to_string(block_size_2) + "NumberOfSuccesses"),
                        swap_number_of_explorations_[block_size_1][block_size_2]);
            }
        }
        if (parameters_.reverse) {
            FFOT_VER(info,
                    std::left << std::setw(28) << ("Reverse:")
                    << reverse_number_of_explorations_
                    << " / " << reverse_number_of_sucesses_
                    << " / " << (double)reverse_number_of_sucesses_ / reverse_number_of_explorations_ * 100 << "%"
                    << std::endl);
            FFOT_PUT(info,
                    "Algorithm", ("ReverseNumberOfExplorations"),
                    reverse_number_of_explorations_);
            FFOT_PUT(info,
                    "Algorithm", ("ReverseNumberOfSuccesses"),
                    reverse_number_of_explorations_);
        }
        for (ElementPos block_size = 2; block_size <= parameters_.shift_reverse_block_maximum_length; ++block_size) {
            FFOT_VER(info,
                    std::left << std::setw(28) << ("Shift-reverse " + std::to_string(block_size) + ":")
                    << shift_reverse_number_of_explorations_[block_size]
                    << " / " << shift_reverse_number_of_sucesses_[block_size]
                    << " / " << (double)shift_reverse_number_of_sucesses_[block_size] / shift_reverse_number_of_explorations_[block_size] * 100 << "%"
                    << std::endl);
            FFOT_PUT(info,
                    "Algorithm", ("ShiftReverse" + std::to_string(block_size) + "NumberOfExplorations"),
                    shift_reverse_number_of_explorations_[block_size]);
            FFOT_PUT(info,
                    "Algorithm", ("ShiftReverse" + std::to_string(block_size) + "NumberOfSuccesses"),
                    shift_reverse_number_of_explorations_[block_size]);
        }

        if (parameters_.add_remove) {
            FFOT_VER(info,
                    std::left << std::setw(28) << ("Add:")
                    << add_number_of_explorations_
                    << " / " << add_number_of_sucesses_
                    << " / " << (double)add_number_of_sucesses_ / add_number_of_explorations_ * 100 << "%"
                    << std::endl);
            FFOT_PUT(info,
                    "Algorithm", ("ReverseNumberOfExplorations"),
                    add_number_of_explorations_);
            FFOT_PUT(info,
                    "Algorithm", ("ReverseNumberOfSuccesses"),
                    add_number_of_explorations_);
            FFOT_VER(info,
                    std::left << std::setw(28) << ("Add:")
                    << remove_number_of_explorations_
                    << " / " << remove_number_of_sucesses_
                    << " / " << (double)remove_number_of_sucesses_ / remove_number_of_explorations_ * 100 << "%"
                    << std::endl);
            FFOT_PUT(info,
                    "Algorithm", ("ReverseNumberOfExplorations"),
                    remove_number_of_explorations_);
            FFOT_PUT(info,
                    "Algorithm", ("ReverseNumberOfSuccesses"),
                    remove_number_of_explorations_);
        }

        for (ElementPos block_size = 1; block_size <= parameters_.inter_shift_block_maximum_length; ++block_size) {
            FFOT_VER(info,
                    std::left << std::setw(28) << ("Inter-shift " + std::to_string(block_size) + ":")
                    << inter_shift_number_of_explorations_[block_size]
                    << " / " << inter_shift_number_of_sucesses_[block_size]
                    << " / " << (double)inter_shift_number_of_sucesses_[block_size] / inter_shift_number_of_explorations_[block_size] * 100 << "%"
                    << std::endl);
            FFOT_PUT(info,
                    "Algorithm", ("InterShift" + std::to_string(block_size) + "NumberOfExplorations"),
                    inter_shift_number_of_explorations_[block_size]);
            FFOT_PUT(info,
                    "Algorithm", ("InterShift" + std::to_string(block_size) + "NumberOfSuccesses"),
                    inter_shift_number_of_explorations_[block_size]);
        }
        for (ElementPos block_size_1 = 1; block_size_1 <= parameters_.inter_swap_block_maximum_length; ++block_size_1) {
            for (ElementPos block_size_2 = 1; block_size_2 <= block_size_1; ++block_size_2) {
                FFOT_VER(info,
                        std::left << std::setw(28) << ("Inter-swap " + std::to_string(block_size_1) + "," + std::to_string(block_size_2) + ":")
                        << inter_swap_number_of_explorations_[block_size_1][block_size_2]
                        << " / " << inter_swap_number_of_sucesses_[block_size_1][block_size_2]
                        << " / " << (double)inter_swap_number_of_sucesses_[block_size_1][block_size_2] / inter_swap_number_of_explorations_[block_size_1][block_size_2] * 100 << "%"
                        << std::endl);
                FFOT_PUT(info,
                        "Algorithm", ("InterSwap" + std::to_string(block_size_1) + "," + std::to_string(block_size_2) + "NumberOfExplorations"),
                        inter_swap_number_of_explorations_[block_size_1][block_size_2]);
                FFOT_PUT(info,
                        "Algorithm", ("InterShift" + std::to_string(block_size_1) + "," + std::to_string(block_size_2) + "NumberOfSuccesses"),
                        inter_swap_number_of_explorations_[block_size_1][block_size_2]);
            }
        }
        if (parameters_.inter_two_opt) {
            FFOT_VER(info,
                    std::left << std::setw(28) << ("Inter-two-opt:")
                    << inter_two_opt_number_of_explorations_
                    << " / " << inter_two_opt_number_of_sucesses_
                    << " / " << (double)inter_two_opt_number_of_sucesses_ / inter_two_opt_number_of_explorations_ * 100 << "%"
                    << std::endl);
            FFOT_PUT(info,
                    "Algorithm", ("InterTwoOptNumberOfExplorations"),
                    inter_two_opt_number_of_explorations_);
            FFOT_PUT(info,
                    "Algorithm", ("InterTwoOptNumberOfSuccesses"),
                    inter_two_opt_number_of_explorations_);
        }
        for (ElementPos block_size = 2; block_size <= parameters_.inter_shift_reverse_block_maximum_length; ++block_size) {
            FFOT_VER(info,
                    std::left << std::setw(28) << ("Inter-shift-reverse " + std::to_string(block_size) + ":")
                    << inter_shift_reverse_number_of_explorations_[block_size]
                    << " / " << inter_shift_reverse_number_of_sucesses_[block_size]
                    << " / " << (double)inter_shift_reverse_number_of_sucesses_[block_size] / inter_shift_reverse_number_of_explorations_[block_size] * 100 << "%"
                    << std::endl);
            FFOT_PUT(info,
                    "Algorithm", ("InterShiftReverse" + std::to_string(block_size) + "NumberOfExplorations"),
                    inter_shift_reverse_number_of_explorations_[block_size]);
            FFOT_PUT(info,
                    "Algorithm", ("InterShiftReverse" + std::to_string(block_size) + "NumberOfSuccesses"),
                    inter_shift_reverse_number_of_explorations_[block_size]);
        }
    }

private:

    /*
     * global_cost(gc1, gc2).
     */

    template<typename, typename T>
    struct HasGlobalCostMergeMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasGlobalCostMergeMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().global_cost_merge(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    GlobalCost global_cost_merge(
            const GlobalCost& gc1,
            const GlobalCost& gc2,
            std::false_type) const
    {
        return gc1 + gc2;
    }

    GlobalCost global_cost_merge(
            const GlobalCost& gc1,
            const GlobalCost& gc2,
            std::true_type) const
    {
        return local_scheme_0_.global_cost_merge(gc1, gc2);
    }

    /**
     * Return the merged global cost of two global costs.
     *
     * If LocalScheme0 does not implement a method
     * 'GlobalCost global_cost(const GlobalCost&, const GlobalCost&)',
     * then it returns the sum of the two global costs, i.e., the global cost
     * of two sequences is the sum of the global cost of each sequence. If it
     * is not the case, for example when the objective is the makespan,
     * LocalScheme0 needs to implement a custom method.
     */
    GlobalCost global_cost_merge(
            const GlobalCost& gc1,
            const GlobalCost& gc2) const
    {
        return global_cost_merge(
                gc1,
                gc2,
                std::integral_constant<
                    bool,
                    HasGlobalCostMergeMethod<LocalScheme0, GlobalCost(
                        const GlobalCost&,
                        const GlobalCost&)>::value>());
    }

    /**
     * Compute and return the global cost of a solution.
     */
    inline void compute_global_cost(
            Solution& solution) const
    {
        solution.global_cost = global_cost(solution.sequences[0]);
        for (SequenceId i = 1; i < number_of_sequences(); ++i) {
            solution.global_cost = global_cost_merge(
                    solution.global_cost,
                    global_cost(solution.sequences[i]));
        }
    }

    /**
     * Compute and return the global cost of a solution without sequence
     * 'i_out'.
     */
    inline GlobalCost compute_global_cost(
            const Solution& solution,
            SequenceId i_out) const
    {
        if (number_of_sequences() == 1)
            return worst<GlobalCost>();
        if (i_out != 0) {
            GlobalCost gc = global_cost(solution.sequences[0]);
            for (SequenceId i = 1; i < number_of_sequences(); ++i)
                if (i != i_out)
                    gc = global_cost_merge(gc, global_cost(solution.sequences[i]));
            return gc;
        } else {
            GlobalCost gc = global_cost(solution.sequences[1]);
            for (SequenceId i = 2; i < number_of_sequences(); ++i)
                gc = global_cost_merge(gc, global_cost(solution.sequences[i]));
            return gc;
        }
    }

    /**
     * Compute and return the global cost of a solution without sequence
     * 'i_out'.
     */
    inline GlobalCost compute_global_cost(
            const Solution& solution,
            SequenceId i_out_1,
            SequenceId i_out_2) const
    {
        if (number_of_sequences() == 1) {
            throw std::logic_error(
                    "'compute_global_cost(solution, i_out_1, i_out_2)'"
                    " requires at least 2 sequences.");
        }
        if (number_of_sequences() == 2)
            return worst<GlobalCost>();
        if (i_out_1 != 0 && i_out_2 != 0) {
            GlobalCost gc = global_cost(solution.sequences[0]);
            for (SequenceId i = 1; i < number_of_sequences(); ++i)
                if (i != i_out_1 && i != i_out_2)
                    gc = global_cost_merge(gc, global_cost(solution.sequences[i]));
            return gc;
        } else if (i_out_1 != 1 && i_out_2 != 1) {
            GlobalCost gc = global_cost(solution.sequences[1]);
            for (SequenceId i = 2; i < number_of_sequences(); ++i)
                if (i != i_out_1 && i != i_out_2)
                    gc = global_cost_merge(gc, global_cost(solution.sequences[i]));
            return gc;
        } else {
            GlobalCost gc = global_cost(solution.sequences[2]);
            for (SequenceId i = 3; i < number_of_sequences(); ++i)
                gc = global_cost_merge(gc, global_cost(solution.sequences[i]));
            return gc;
        }
    }

    /*
     * Evaluate moves.
     */

    inline void compute_cost_shift(
            const Solution& solution,
            SequenceId i,
            ElementPos block_pos,
            ElementPos block_size,
            bool reverse = false)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        const Sequence& sequence = solution.sequences[i];
        // Global cost without sequence i.
        GlobalCost gc0 = compute_global_cost(solution, i);
        // Initialize sequence_cur_.
        sequence_cur_1_ = empty_sequence(i);
        // Reset global_costs_1d_.
        std::fill(
                global_costs_1d_.begin(),
                global_costs_1d_.end(),
                worst<GlobalCost>());

        // Loop through all new positions.
        ElementPos pos_min = std::max(
                (ElementPos)0,
                block_pos - parameters_.shift_maximum_distance);
        ElementPos pos_max = std::min(
                (ElementPos)sequence.sequence.size() - block_size,
                block_pos + parameters_.shift_maximum_distance);
        for (ElementPos pos_new = pos_min; pos_new <= pos_max; ++pos_new) {

            bool stop = false;

            // Initialize sequence_tmp_.
            sequence_tmp_1_ = sequence_cur_1_;

            // Add block to sequence_tmp_1_.
            if (!reverse) {
                for (ElementPos p = block_pos; p < block_pos + block_size && !stop; ++p) {
                    // Add element to sequence_tmp_1_.
                    ElementId j = sequence.sequence[p];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    if (m == 1) {
                        if (bound(sequence_tmp_1_) >= gc)
                            stop = true;
                    } else {
                        if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                            stop = true;
                    }
                }
            } else {
                for (ElementPos p = block_pos + block_size - 1; p >= block_pos && !stop; --p) {
                    // Add element to sequence_tmp_1_.
                    ElementId j = sequence.sequence[p];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    if (m == 1) {
                        if (bound(sequence_tmp_1_) >= gc)
                            stop = true;
                    } else {
                        if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                            stop = true;
                    }
                }
            }

            // Add the remaining elements to sequence_tmp.
            ElementPos p0 = (pos_new < block_pos)? pos_new: pos_new + block_size;
            for (ElementPos p = p0; p < (ElementPos)sequence.sequence.size() && !stop; ++p) {
                // Skip elements from the previously added bloc.
                if (block_pos <= p && p < block_pos + block_size)
                    continue;
                // Add element to sequence_tmp_1_.
                ElementId j = sequence.sequence[p];
                local_scheme_0_.append(sequence_tmp_1_, j);
                // Check early termination.
                if (m == 1) {
                    if (bound(sequence_tmp_1_) >= gc)
                        stop = true;
                } else {
                    if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                        stop = true;
                }
            }

            if (!stop) {
                if (m == 1) {
                    global_costs_1d_[pos_new] = global_cost(sequence_tmp_1_);
                } else {
                    global_costs_1d_[pos_new] = global_cost_merge(
                            gc0, local_scheme_0_.global_cost(sequence_tmp_1_));
                }
            }

            // Stop condition.
            if (pos_new == (ElementPos)sequence.sequence.size() - block_size)
                break;

            // Add j1 to sequence_cur_.
            assert(p0 < (ElementPos)sequence.sequence.size());
            ElementId j1 = sequence.sequence[p0];
            local_scheme_0_.append(sequence_cur_1_, j1);
            // Check early termination.
            if (m == 1) {
                if (bound(sequence_cur_1_) >= gc)
                    break;
            } else {
                if (global_cost_merge(gc0, bound(sequence_cur_1_)) >= gc)
                    break;
            }
        }
    }

    inline void compute_cost_swap(
            const Solution& solution,
            SequenceId i,
            Counter block_size)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        const Sequence& sequence = solution.sequences[i];
        // Global cost without sequence i.
        GlobalCost gc0 = compute_global_cost(solution, i);
        // Initialize sequence_cur_.
        sequence_cur_1_ = empty_sequence(i);
        // Reset global_costs_swap_.
        for (ElementPos pos_1 = 0; pos_1 < local_scheme_0_.number_of_elements(); ++pos_1) {
            std::fill(
                    global_costs_2d_1_[pos_1].begin(),
                    global_costs_2d_1_[pos_1].end(),
                    worst<GlobalCost>());
        }

        // Loop through all pairs.
        Counter pos_max = (ElementPos)sequence.sequence.size() - block_size;
        for (ElementPos pos_1 = 0; pos_1 <= pos_max; ++pos_1) {
            ElementPos pos_2_max = std::min(
                    pos_max,
                    pos_1 + parameters_.swap_maximum_distance);
            for (ElementPos pos_2 = pos_1 + block_size; pos_2 < pos_2_max; ++pos_2) {

                bool stop = false;

                // Initialize sequence_tmp_1_.
                sequence_tmp_1_ = sequence_cur_1_;
                // Add remaining elements.
                for (ElementPos pos = pos_1; pos < (ElementPos)sequence.sequence.size() && !stop; ++pos) {
                    ElementId j = sequence.sequence[pos];
                    // If j1 or j2, swap.
                    if (pos_1 <= pos && pos < pos_1 + block_size) {
                        ElementPos diff = pos - pos_1;
                        j = sequence.sequence[pos_2 + diff];
                    }
                    if (pos_2 <= pos && pos < pos_2 + block_size) {
                        ElementPos diff = pos - pos_2;
                        j = sequence.sequence[pos_1 + diff];
                    }
                    // Add element to sequence_tmp_1_.
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    if (m == 1) {
                        if (bound(sequence_tmp_1_) >= gc)
                            stop = true;
                    } else {
                        if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                            stop = true;
                    }
                }

                if (!stop) {
                    if (m == 1) {
                        global_costs_2d_1_[pos_1][pos_2 - pos_1 - 1]
                            = local_scheme_0_.global_cost(sequence_tmp_1_);
                    } else {
                        global_costs_2d_1_[pos_1][pos_2 - pos_1 - 1]
                            = global_cost_merge(gc0, local_scheme_0_.global_cost(sequence_tmp_1_));
                    }
                }
            }

            // Add j1 to sequence_cur_.
            ElementId j1 = sequence.sequence[pos_1];
            local_scheme_0_.append(sequence_cur_1_, j1);
            // Check early termination.
            if (m == 1) {
                if (bound(sequence_cur_1_) >= gc)
                    break;
            } else {
                if (global_cost_merge(gc0, bound(sequence_cur_1_)) >= gc)
                    break;
            }
        }
    }

    inline void compute_cost_swap(
            const Solution& solution,
            SequenceId i,
            Counter block_size_1,
            Counter block_size_2)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        const Sequence& sequence = solution.sequences[i];
        ElementPos n = (ElementPos)sequence.sequence.size();
        // Global cost without sequence i.
        GlobalCost gc0 = compute_global_cost(solution, i);
        // Initialize sequence_cur_.
        sequence_cur_1_ = empty_sequence(i);
        // Reset global_costs_swap_.
        for (ElementPos pos_1 = 0; pos_1 < local_scheme_0_.number_of_elements(); ++pos_1) {
            std::fill(
                    global_costs_2d_2_[pos_1].begin(),
                    global_costs_2d_2_[pos_1].end(),
                    worst<GlobalCost>());
        }

        // Loop through all pairs.
        for (ElementPos pos = 0; pos < n; ++pos) {

            // block 1 is at [pos, pos + block_size_1[.
            ElementPos pos_1 = pos;
            ElementPos pos_2_max = std::min(
                    n - block_size_2,
                    pos + block_size_1 + parameters_.swap_maximum_distance);
            for (ElementPos pos_2 = pos + block_size_1; pos_2 <= pos_2_max; ++pos_2) {

                bool stop = false;

                // Initialize sequence_tmp_1_.
                sequence_tmp_1_ = sequence_cur_1_;
                // Add block 2.
                for (ElementPos p = pos_2; p < pos_2 + block_size_2 && !stop; ++p) {
                    // Add element to sequence_tmp_1_.
                    ElementId j = sequence.sequence[p];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    if (m == 1) {
                        if (bound(sequence_tmp_1_) >= gc)
                            stop = true;
                    } else {
                        if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                            stop = true;
                    }
                }
                // Add middle elements.
                for (ElementPos p = pos_1 + block_size_1; p < pos_2 && !stop; ++p) {
                    // Add element to sequence_tmp_1_.
                    ElementId j = sequence.sequence[p];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    if (m == 1) {
                        if (bound(sequence_tmp_1_) >= gc)
                            stop = true;
                    } else {
                        if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                            stop = true;
                    }
                }
                // Add block 1.
                for (ElementPos p = pos_1; p < pos_1 + block_size_1 && !stop; ++p) {
                    // Add element to sequence_tmp_1_.
                    ElementId j = sequence.sequence[p];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    if (m == 1) {
                        if (bound(sequence_tmp_1_) >= gc)
                            stop = true;
                    } else {
                        if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                            stop = true;
                    }
                }
                // Add end elements.
                for (ElementPos p = pos_2 + block_size_2; p < n && !stop; ++p) {
                    // Add element to sequence_tmp_1_.
                    ElementId j = sequence.sequence[p];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    if (m == 1) {
                        if (bound(sequence_tmp_1_) >= gc)
                            stop = true;
                    } else {
                        if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                            stop = true;
                    }
                }

                if (!stop) {
                    if (m == 1) {
                        global_costs_2d_2_[pos_1][pos_2]
                            = local_scheme_0_.global_cost(sequence_tmp_1_);
                    } else {
                        global_costs_2d_2_[pos_1][pos_2]
                            = global_cost_merge(gc0, local_scheme_0_.global_cost(sequence_tmp_1_));
                    }
                }
            }

            // block 2 is at [pos, pos + block_size_2[.
            ElementPos pos_2 = pos;
            ElementPos pos_1_max = std::min(
                    n - block_size_1,
                    pos + block_size_2 + parameters_.swap_maximum_distance);
            for (ElementPos pos_1 = pos + block_size_2; pos_1 <= pos_1_max; ++pos_1) {

                bool stop = false;

                // Initialize sequence_tmp_1_.
                sequence_tmp_1_ = sequence_cur_1_;
                // Add block 1.
                for (ElementPos p = pos_1; p < pos_1 + block_size_1 && !stop; ++p) {
                    // Add element to sequence_tmp_1_.
                    ElementId j = sequence.sequence[p];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    if (m == 1) {
                        if (bound(sequence_tmp_1_) >= gc)
                            stop = true;
                    } else {
                        if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                            stop = true;
                    }
                }
                // Add middle elements.
                for (ElementPos p = pos_2 + block_size_2; p < pos_1 && !stop; ++p) {
                    // Add element to sequence_tmp_1_.
                    ElementId j = sequence.sequence[p];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    if (m == 1) {
                        if (bound(sequence_tmp_1_) >= gc)
                            stop = true;
                    } else {
                        if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                            stop = true;
                    }
                }
                // Add block 2.
                for (ElementPos p = pos_2; p < pos_2 + block_size_2 && !stop; ++p) {
                    // Add element to sequence_tmp_1_.
                    ElementId j = sequence.sequence[p];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    if (m == 1) {
                        if (bound(sequence_tmp_1_) >= gc)
                            stop = true;
                    } else {
                        if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                            stop = true;
                    }
                }
                // Add end elements.
                for (ElementPos p = pos_1 + block_size_1; p < n && !stop; ++p) {
                    // Add element to sequence_tmp_1_.
                    ElementId j = sequence.sequence[p];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    if (m == 1) {
                        if (bound(sequence_tmp_1_) >= gc)
                            stop = true;
                    } else {
                        if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                            stop = true;
                    }
                }
                //if (pos_1 == 46 && pos_2 == 44) {
                //    std::cout << "sequence: ";
                //    for (ElementId j: sequence.sequence)
                //        std::cout << " " << j;
                //    std::cout << std::endl;
                //    std::cout << "sequence: ";
                //    for (ElementId j: sequence_tmp_1_.sequence)
                //        std::cout << " " << j;
                //    std::cout << std::endl;
                //}

                if (!stop) {
                    if (m == 1) {
                        global_costs_2d_2_[pos_1][pos_2]
                            = local_scheme_0_.global_cost(sequence_tmp_1_);
                    } else {
                        global_costs_2d_2_[pos_1][pos_2]
                            = global_cost_merge(gc0, local_scheme_0_.global_cost(sequence_tmp_1_));
                    }
                }
            }

            // Add j to sequence_cur_.
            ElementId j = sequence.sequence[pos];
            local_scheme_0_.append(sequence_cur_1_, j);
            // Check early termination.
            if (m == 1) {
                if (bound(sequence_cur_1_) >= gc)
                    break;
            } else {
                if (global_cost_merge(gc0, bound(sequence_cur_1_)) >= gc)
                    break;
            }
        }
    }

    inline void compute_cost_reverse(
            const Solution& solution,
            SequenceId i)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        const Sequence& sequence = solution.sequences[i];
        // Global cost without sequence i.
        GlobalCost gc0 = compute_global_cost(solution, i);
        // Initialize sequence_cur_.
        sequence_cur_1_ = empty_sequence(i);
        // Reset global_costs_swap_.
        for (ElementPos pos_1 = 0; pos_1 < local_scheme_0_.number_of_elements(); ++pos_1) {
            std::fill(
                    global_costs_2d_1_[pos_1].begin(),
                    global_costs_2d_1_[pos_1].end(),
                    worst<GlobalCost>());
        }

        // Loop through all pairs.
        for (ElementPos pos_1 = 0; pos_1 < (ElementPos)sequence.sequence.size(); ++pos_1) {
            ElementPos pos_max = std::min(
                    (ElementPos)sequence.sequence.size(),
                    pos_1 + parameters_.reverse_maximum_length);
            for (ElementPos pos_2 = pos_1 + 2; pos_2 < pos_max; ++pos_2) {

                bool stop = false;

                // Initialize sequence_tmp_1_.
                sequence_tmp_1_ = sequence_cur_1_;
                // Add reverse sequence.
                for (ElementPos pos = pos_2; pos >= pos_1 && !stop; --pos) {
                    // Add element to sequence_tmp_1_.
                    ElementId j = sequence.sequence[pos];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    if (m == 1) {
                        if (bound(sequence_tmp_1_) >= gc)
                            stop = true;
                    } else {
                        if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                            stop = true;
                    }
                }
                // Add remaining elements.
                for (ElementPos pos = pos_2 + 1; pos < (ElementPos)sequence.sequence.size() && !stop; ++pos) {
                    // Add element to sequence_tmp_1_.
                    ElementId j = sequence.sequence[pos];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    if (m == 1) {
                        if (bound(sequence_tmp_1_) >= gc)
                            stop = true;
                    } else {
                        if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                            stop = true;
                    }
                }

                if (!stop) {
                    if (m == 1) {
                        global_costs_2d_1_[pos_1][pos_2 - pos_1 - 1]
                            = local_scheme_0_.global_cost(sequence_tmp_1_);
                    } else {
                        global_costs_2d_1_[pos_1][pos_2 - pos_1 - 1]
                            = global_cost_merge(gc0, local_scheme_0_.global_cost(sequence_tmp_1_));
                    }
                }
            }

            // Add j1 to sequence_cur_.
            ElementId j1 = sequence.sequence[pos_1];
            local_scheme_0_.append(sequence_cur_1_, j1);
            // Check early termination.
            if (m == 1) {
                if (bound(sequence_cur_1_) >= gc)
                    break;
            } else {
                if (global_cost_merge(gc0, bound(sequence_cur_1_)) >= gc)
                    break;
            }
        }
    }

    inline void compute_cost_add(
            const Solution& solution,
            SequenceId i,
            ElementId j0)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        const Sequence& sequence = solution.sequences[i];
        // Global cost without sequence i.
        GlobalCost gc0 = compute_global_cost(solution, i);
        // Initialize sequence_cur_.
        sequence_cur_1_ = empty_sequence(i);
        // Reset global_costs_1d_.
        std::fill(
                global_costs_1d_.begin(),
                global_costs_1d_.end(),
                worst<GlobalCost>());

        // Loop through all new positions.
        for (ElementPos pos = 0; pos <= (ElementPos)sequence.sequence.size(); ++pos) {

            bool stop = false;

            // Initialize sequence_tmp_1_.
            sequence_tmp_1_ = sequence_cur_1_;

            // Add element j0 to sequence_tmp_1_.
            local_scheme_0_.append(sequence_tmp_1_, j0);

            // Add the remaining elements to sequence_tmp.
            for (ElementPos p = pos; p < (ElementPos)sequence.sequence.size() && !stop; ++p) {
                // Add next element to sequence_tmp_1_.
                ElementId j = sequence.sequence[p];
                local_scheme_0_.append(sequence_tmp_1_, j);
                // Check early termination.
                if (m == 1) {
                    if (bound(sequence_tmp_1_) >= gc)
                        stop = true;
                } else {
                    if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                        stop = true;
                }
            }

            if (!stop) {
                if (m == 1) {
                    global_costs_1d_[pos] = global_cost(sequence_tmp_1_);
                } else {
                    global_costs_1d_[pos]
                        = global_cost_merge(gc0, global_cost(sequence_tmp_1_));
                }
            }

            // Stop condition.
            if (pos == (ElementPos)sequence.sequence.size())
                break;

            // Add next element to sequence_cur_.
            ElementId j = sequence.sequence[pos];
            local_scheme_0_.append(sequence_cur_1_, j);
        }
    }

    inline void compute_cost_remove(
            const Solution& solution,
            SequenceId i)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        const Sequence& sequence = solution.sequences[i];
        // Global cost without sequence i.
        GlobalCost gc0 = compute_global_cost(solution, i);
        // Initialize sequence_cur_.
        sequence_cur_1_ = empty_sequence(i);
        // Reset global_costs_1d_.
        std::fill(
                global_costs_1d_.begin(),
                global_costs_1d_.end(),
                worst<GlobalCost>());

        // Loop through all new positions.
        for (ElementPos pos = 0; pos < (ElementPos)sequence.sequence.size(); ++pos) {

            bool stop = false;

            // Initialize sequence_tmp_1_.
            sequence_tmp_1_ = sequence_cur_1_;

            // Add the remaining elements to sequence_tmp.
            for (ElementPos p = pos + 1; p < (ElementPos)sequence.sequence.size() && !stop; ++p) {
                // Add next element to sequence_tmp_1_.
                ElementId j = sequence.sequence[p];
                local_scheme_0_.append(sequence_tmp_1_, j);
                // Check early termination.
                if (m == 1) {
                    if (bound(sequence_tmp_1_) >= gc)
                        stop = true;
                } else {
                    if (global_cost_merge(gc0, bound(sequence_tmp_1_)) >= gc)
                        stop = true;
                }
            }

            if (!stop) {
                if (m == 1) {
                    global_costs_1d_[pos] = global_cost(sequence_tmp_1_);
                } else {
                    global_costs_1d_[pos]
                        = global_cost_merge(gc0, local_scheme_0_.global_cost(sequence_tmp_1_));
                }
            }

            // Add next element to sequence_cur_.
            ElementId j = sequence.sequence[pos];
            local_scheme_0_.append(sequence_cur_1_, j);
        }
    }

    inline void compute_cost_inter_shift(
            const Solution& solution,
            SequenceId i1,
            SequenceId i2,
            ElementPos block_size,
            bool reverse = false)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        const Sequence& sequence_1 = solution.sequences[i1];
        const Sequence& sequence_2 = solution.sequences[i2];
        ElementPos n1 = sequence_1.sequence.size();
        ElementPos n2 = sequence_2.sequence.size();
        // Global cost without sequence i.
        GlobalCost gc0 = compute_global_cost(solution, i1, i2);
        //std::cout << "i1 " << i1 << " i2 " << i2 << " gc0 " << to_string(gc0) << std::endl;
        // Initialize sequence_cur_.
        sequence_cur_2_ = empty_sequence(i2);
        // Reset global_costs_2d_2_.
        for (ElementPos pos_1 = 0; pos_1 < local_scheme_0_.number_of_elements(); ++pos_1) {
            std::fill(
                    global_costs_2d_2_[pos_1].begin(),
                    global_costs_2d_2_[pos_1].end(),
                    worst<GlobalCost>());
        }

        sequence_cur_1_ = empty_sequence(i1);

        for (ElementPos pos_1 = 0; pos_1 <= n1 - block_size; ++pos_1) {
            // pos_1 is the start position of the block of element to shift.

            // Compute first sequence, i.e. the one from which the block iof
            // elements is removed.
            sequence_tmp_1_ = sequence_cur_1_;
            for (ElementPos p = pos_1 + block_size; p < n1; ++p) {
                ElementId j = sequence_1.sequence[p];
                local_scheme_0_.append(sequence_tmp_1_, j);
            }

            sequence_cur_2_ = empty_sequence(i2);

            // Loop through all new positions.
            ElementPos pos_min = std::max(
                    (ElementPos)0,
                    pos_1 - parameters_.shift_maximum_distance);
            ElementPos pos_max = std::min(
                    n2 + 1,
                    pos_1 + parameters_.shift_maximum_distance);
            for (ElementPos pos_2 = pos_min; pos_2 <= pos_max; ++pos_2) {

                bool stop = false;

                // Initialize sequence_tmp_2_.
                sequence_tmp_2_ = sequence_cur_2_;

                // Add block to sequence_tmp_2_.
                if (!reverse) {
                    for (ElementPos p = pos_1; p < pos_1 + block_size && !stop; ++p) {
                        // Add next element to sequence_tmp_2_.
                        ElementId j = sequence_1.sequence[p];
                        local_scheme_0_.append(sequence_tmp_2_, j);
                        // Check early termination.
                        GlobalCost bnd = global_cost_merge(
                                bound(sequence_tmp_1_),
                                bound(sequence_tmp_2_));
                        if (m > 2)
                            bnd = global_cost_merge(gc0, bnd);
                        if (bnd >= gc)
                            stop = true;
                    }
                } else {
                    for (ElementPos p = pos_1 + block_size - 1; p >= pos_1 && !stop; --p) {
                        // Add next element to sequence_tmp_2_.
                        ElementId j = sequence_1.sequence[p];
                        local_scheme_0_.append(sequence_tmp_2_, j);
                        // Check early termination.
                        GlobalCost bnd = global_cost_merge(
                                bound(sequence_tmp_1_),
                                bound(sequence_tmp_2_));
                        if (m > 2)
                            bnd = global_cost_merge(gc0, bnd);
                        if (bnd >= gc)
                            stop = true;
                    }
                }

                // Add the remaining elements to sequence_tmp_2_.
                for (ElementPos p = pos_2; p < n2 && !stop; ++p) {
                    // Add next element to sequence_tmp_2_.
                    ElementId j = sequence_2.sequence[p];
                    local_scheme_0_.append(sequence_tmp_2_, j);
                    // Check early termination.
                    GlobalCost bnd = global_cost_merge(
                            bound(sequence_tmp_1_),
                            bound(sequence_tmp_2_));
                    if (m > 2)
                        bnd = global_cost_merge(gc0, bnd);
                    if (bnd >= gc)
                        stop = true;
                }

                if (!stop) {
                    GlobalCost gc_tmp = global_cost_merge(
                            global_cost(sequence_tmp_1_),
                            global_cost(sequence_tmp_2_));
                    if (m > 2)
                        gc_tmp = global_cost_merge(gc0, gc_tmp);
                    if (gc_tmp < gc) {
                        //for (ElementId j: sequence_tmp_1_.sequence)
                        //    std::cout << " " << j;
                        //std::cout << std::endl;
                        //for (ElementId j: sequence_tmp_2_.sequence)
                        //    std::cout << " " << j;
                        //std::cout << std::endl;
                        //std::cout << "gc0 " << to_string(gc0) << std::endl;
                        //std::cout << "gc1 " << to_string(global_cost(sequence_tmp_1_)) << std::endl;
                        //std::cout << "gc2 " << to_string(global_cost(sequence_tmp_2_)) << std::endl;
                        //std::cout << "gc_tmp " << to_string(gc_tmp) << std::endl;
                        global_costs_2d_2_[pos_1][pos_2] = gc_tmp;
                    }
                }

                if (pos_2 == n2)
                    break;

                // Add next element to sequence_cur_2_.
                ElementId j = sequence_2.sequence[pos_2];
                local_scheme_0_.append(sequence_cur_2_, j);
                // Check early termination.
                GlobalCost bnd = global_cost_merge(
                        bound(sequence_cur_1_),
                        bound(sequence_cur_2_));
                if (m > 2)
                    bnd = global_cost_merge(gc0, bnd);
                if (bnd >= gc)
                    break;
            }

            // Add next element to sequence_cur_1_.
            ElementId j = sequence_1.sequence[pos_1];
            local_scheme_0_.append(sequence_cur_1_, j);
            // Check early termination.
            GlobalCost bnd = bound(sequence_cur_1_);
            if (m > 2)
                bnd = global_cost_merge(gc0, bnd);
            if (bnd >= gc)
                break;
        }
    }

    inline void compute_cost_inter_two_opt(
            const Solution& solution,
            SequenceId i1,
            SequenceId i2)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        const Sequence& sequence_1 = solution.sequences[i1];
        const Sequence& sequence_2 = solution.sequences[i2];
        ElementPos n1 = sequence_1.sequence.size();
        ElementPos n2 = sequence_2.sequence.size();
        // Global cost without sequence i.
        GlobalCost gc0 = compute_global_cost(solution, i1, i2);
        // Reset global_costs_2d_2_.
        for (ElementPos pos_1 = 0; pos_1 < local_scheme_0_.number_of_elements(); ++pos_1) {
            std::fill(
                    global_costs_2d_2_[pos_1].begin(),
                    global_costs_2d_2_[pos_1].end(),
                    worst<GlobalCost>());
        }

        sequence_cur_1_ = empty_sequence(i1);

        for (ElementPos pos_1 = 0; pos_1 <= n1; ++pos_1) {

            sequence_cur_2_ = empty_sequence(i1);

            for (ElementPos pos_2 = 0; pos_2 <= n2; ++pos_2) {

                bool stop = false;

                sequence_tmp_1_ = sequence_cur_1_;
                sequence_tmp_2_ = sequence_cur_2_;

                for (ElementPos p = pos_2; p < n2 && !stop; ++p) {
                    // Add next element to sequence_tmp_1_.
                    ElementId j = sequence_2.sequence[p];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    GlobalCost bnd = global_cost_merge(
                            bound(sequence_tmp_1_),
                            bound(sequence_tmp_2_));
                    if (m > 2)
                        bnd = global_cost_merge(gc0, bnd);
                    if (bnd >= gc)
                        stop = true;
                }
                for (ElementPos p = pos_1; p < n1 && !stop; ++p) {
                    // Add next element to sequence_tmp_2_.
                    ElementId j = sequence_1.sequence[p];
                    local_scheme_0_.append(sequence_tmp_2_, j);
                    // Check early termination.
                    GlobalCost bnd = global_cost_merge(
                            bound(sequence_tmp_1_),
                            bound(sequence_tmp_2_));
                    if (m > 2)
                        bnd = global_cost_merge(gc0, bnd);
                    if (bnd >= gc)
                        stop = true;
                }

                if (!stop) {
                    GlobalCost gc_tmp = global_cost_merge(
                            global_cost(sequence_tmp_1_),
                            global_cost(sequence_tmp_2_));
                    if (m > 2)
                        gc_tmp = global_cost_merge(gc0, gc_tmp);
                    if (gc_tmp < gc)
                        global_costs_2d_2_[pos_1][pos_2] = gc_tmp;
                }

                if (pos_2 == n2)
                    break;

                // Update sequence_cur_2.
                ElementId j = sequence_2.sequence[pos_2];
                local_scheme_0_.append(sequence_cur_2_, j);
                // Check early termination.
                GlobalCost bnd = global_cost_merge(
                        bound(sequence_cur_1_),
                        bound(sequence_cur_2_));
                if (m > 2)
                    bnd = global_cost_merge(gc0, bnd);
                if (bnd >= gc)
                    break;
            }

            if (pos_1 == n1)
                break;

            // Update sequence_cur_1_.
            ElementId j = sequence_1.sequence[pos_1];
            local_scheme_0_.append(sequence_cur_1_, j);
            // Check early termination.
            GlobalCost bnd = bound(sequence_cur_1_);
            if (m > 2)
                bnd = global_cost_merge(gc0, bnd);
            if (bnd >= gc)
                break;
        }
    }

    inline void compute_cost_inter_swap(
            const Solution& solution,
            SequenceId i1,
            SequenceId i2,
            ElementPos block_size_1,
            ElementPos block_size_2)
    {
        SequenceId m = number_of_sequences();
        GlobalCost gc = global_cost(solution);
        const Sequence& sequence_1 = solution.sequences[i1];
        const Sequence& sequence_2 = solution.sequences[i2];
        ElementPos n1 = sequence_1.sequence.size();
        ElementPos n2 = sequence_2.sequence.size();
        // Global cost without sequence i.
        GlobalCost gc0 = compute_global_cost(solution, i1, i2);
        // Reset global_costs_2d_2_.
        for (ElementPos pos_1 = 0; pos_1 < local_scheme_0_.number_of_elements(); ++pos_1) {
            std::fill(
                    global_costs_2d_2_[pos_1].begin(),
                    global_costs_2d_2_[pos_1].end(),
                    worst<GlobalCost>());
        }

        sequence_cur_1_ = empty_sequence(i1);

        for (ElementPos pos_1 = 0; pos_1 <= n1 - block_size_1; ++pos_1) {

            sequence_cur_2_ = empty_sequence(i1);

            for (ElementPos pos_2 = 0; pos_2 <= n2 - block_size_2; ++pos_2) {

                bool stop = false;

                sequence_tmp_1_ = sequence_cur_1_;
                sequence_tmp_2_ = sequence_cur_2_;

                for (ElementPos p = pos_2; p < pos_2 + block_size_2 && !stop; ++p) {
                    // Add next element to sequence_tmp_1_.
                    ElementId j = sequence_2.sequence[p];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    GlobalCost bnd = global_cost_merge(
                            bound(sequence_tmp_1_),
                            bound(sequence_tmp_2_));
                    if (m > 2)
                        bnd = global_cost_merge(gc0, bnd);
                    if (bnd >= gc)
                        stop = true;
                }
                for (ElementPos p = pos_1 + block_size_1; p < n1 && !stop; ++p) {
                    // Add next element to sequence_tmp_1_.
                    ElementId j = sequence_1.sequence[p];
                    local_scheme_0_.append(sequence_tmp_1_, j);
                    // Check early termination.
                    GlobalCost bnd = global_cost_merge(
                            bound(sequence_tmp_1_),
                            bound(sequence_tmp_2_));
                    if (m > 2)
                        bnd = global_cost_merge(gc0, bnd);
                    if (bnd >= gc)
                        stop = true;
                }

                for (ElementPos p = pos_1; p < pos_1 + block_size_1 && !stop; ++p) {
                    // Add next element to sequence_tmp_2_.
                    ElementId j = sequence_1.sequence[p];
                    local_scheme_0_.append(sequence_tmp_2_, j);
                    // Check early termination.
                    GlobalCost bnd = global_cost_merge(
                            bound(sequence_tmp_1_),
                            bound(sequence_tmp_2_));
                    if (m > 2)
                        bnd = global_cost_merge(gc0, bnd);
                    if (bnd >= gc)
                        stop = true;
                }
                for (ElementPos p = pos_2 + block_size_2; p < n2 && !stop; ++p) {
                    // Add next element to sequence_tmp_2_.
                    ElementId j = sequence_2.sequence[p];
                    local_scheme_0_.append(sequence_tmp_2_, j);
                    // Check early termination.
                    GlobalCost bnd = global_cost_merge(
                            bound(sequence_tmp_1_),
                            bound(sequence_tmp_2_));
                    if (m > 2)
                        bnd = global_cost_merge(gc0, bnd);
                    if (bnd >= gc)
                        stop = true;
                }

                if (!stop) {
                    GlobalCost gc_tmp = global_cost_merge(
                            global_cost(sequence_tmp_1_),
                            global_cost(sequence_tmp_2_));
                    if (m > 2)
                        gc_tmp = global_cost_merge(gc0, gc_tmp);
                    if (gc_tmp < gc)
                        global_costs_2d_2_[pos_1][pos_2] = gc_tmp;
                }

                if (pos_2 == n2)
                    break;

                // Update sequence_cur_2.
                ElementId j = sequence_2.sequence[pos_2];
                local_scheme_0_.append(sequence_cur_2_, j);
                // Check early termination.
                GlobalCost bnd = global_cost_merge(
                        bound(sequence_cur_1_),
                        bound(sequence_cur_2_));
                if (m > 2)
                    bnd = global_cost_merge(gc0, bnd);
                if (bnd >= gc)
                    break;
            }

            if (pos_1 == n1)
                break;

            // Update sequence_cur_1_.
            ElementId j = sequence_1.sequence[pos_1];
            local_scheme_0_.append(sequence_cur_1_, j);
            // Check early termination.
            GlobalCost bnd = bound(sequence_cur_1_);
            if (m > 2)
                bnd = global_cost_merge(gc0, bnd);
            if (bnd >= gc)
                break;
        }
    }

    /*
     * Private attributes.
     */

    /** Input local scheme. */
    LocalScheme0& local_scheme_0_;
    /** Parameters. */
    Parameters parameters_;

    /*
     * Structures to iterate in random order.
     */

    /** Vector containing numbers from '0' to 'number_of_sequences - 1'. */
    std::vector<SequenceId> sequences_1_;
    /** Vector containing numbers from '0' to 'number_of_sequences - 1'. */
    std::vector<SequenceId> sequences_2_;
    /** IndexedSet of size 'number_of_elements'. */
    optimizationtools::IndexedSet elements_;
    /** Vector containing numbers from '0' to 'number_of_elements - 1'. */
    std::vector<ElementPos> positions_1_;
    /** Vector containing numbers from '0' to 'number_of_elements - 1'. */
    std::vector<ElementPos> positions_2_;
    /**
     * Vector of all non-ordered pairs from '0' to 'number_of_elements'.
     *
     * Contains 'number_of_element x (number_of_elements - 1) / 2' elements.
     */
    std::vector<std::pair<ElementPos, ElementPos>> pairs_1_;
    /**
     * Vector of all ordered pairs from '0' to 'number_of_elements'.
     *
     * Contains 'number_of_element x number_of_elements' elements.
     */
    std::vector<std::pair<ElementPos, ElementPos>> pairs_2_;

    /*
     * Structures storing neighborhood costs.
     */

    /**
     * Structure to store neighborhood global costs in a 1D vector of size
     * 'number_of_elements'.
     */
    std::vector<GlobalCost> global_costs_1d_;
    /**
     * Structure to store neighborhood global costs in a 2D triangular matrix
     * of size 'number_of_elements x (number_of_elements - 1) / 2'.
     */
    std::vector<std::vector<GlobalCost>> global_costs_2d_1_;
    /**
     * Structure to store neighborhood global costs in a 2D matrix of size
     * 'number_of_elements x number_of_elements'.
     */
    std::vector<std::vector<GlobalCost>> global_costs_2d_2_;

    /*
     * Temporary structures.
     */

    Sequence sequence_cur_1_;
    Sequence sequence_cur_2_;
    Solution solution_cur_;
    Sequence sequence_tmp_1_;
    Sequence sequence_tmp_2_;
    Solution solution_tmp_;

    /*
     * Statistics.
     */

    std::vector<Counter> shift_number_of_explorations_;
    std::vector<Counter> shift_number_of_sucesses_;
    std::vector<std::vector<Counter>> swap_number_of_explorations_;
    std::vector<std::vector<Counter>> swap_number_of_sucesses_;
    Counter reverse_number_of_explorations_ = 0;
    Counter reverse_number_of_sucesses_ = 0;
    std::vector<Counter> shift_reverse_number_of_explorations_;
    std::vector<Counter> shift_reverse_number_of_sucesses_;
    Counter add_number_of_explorations_ = 0;
    Counter add_number_of_sucesses_ = 0;
    Counter remove_number_of_explorations_ = 0;
    Counter remove_number_of_sucesses_ = 0;
    std::vector<Counter> inter_shift_number_of_explorations_;
    std::vector<Counter> inter_shift_number_of_sucesses_;
    std::vector<std::vector<Counter>> inter_swap_number_of_explorations_;
    std::vector<std::vector<Counter>> inter_swap_number_of_sucesses_;
    Counter inter_two_opt_number_of_explorations_ = 0;
    Counter inter_two_opt_number_of_sucesses_ = 0;
    std::vector<Counter> inter_shift_reverse_number_of_explorations_;
    std::vector<Counter> inter_shift_reverse_number_of_sucesses_;

};

}

}

