/**
 * Sequencing problems.
 */

#pragma once

#include "localsearchsolver/common.hpp"

#include "optimizationtools/containers/indexed_set.hpp"

namespace localsearchsolver
{

namespace sequencing
{

using SequenceId = int64_t;
using SequencePos = int64_t;
using ElementId = int64_t;
using ElementPos = int64_t;
using Mode = int64_t;

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
    Swap,
    Reverse,
    ShiftReverse,

    Add,
    Remove,

    InterTwoOpt,
    InterShift,
    InterSwap,
    InterShiftReverse,
    InterSwapStar,

    ShiftChangeMode,
    ModeSwap,
    SwapWithModes,
};

std::string neighborhood2string(
        Neighborhoods neighborhood,
        ElementPos k1,
        ElementPos k2)
{
    switch (neighborhood) {
    case Neighborhoods::Shift:
        return std::to_string(k1) + "-shift";
    case Neighborhoods::Swap:
        return "(" + std::to_string(k1) + "," + std::to_string(k2) + ")" + "-swap";
    case Neighborhoods::Reverse:
        return "Reverse";
    case Neighborhoods::ShiftReverse:
        return std::to_string(k1) + "-shift-reverse";
    case Neighborhoods::Add:
        return "Add";
    case Neighborhoods::Remove:
        return "Remove";
    case Neighborhoods::InterTwoOpt:
        return "Inter-two-opt";
    case Neighborhoods::InterShift:
        return std::to_string(k1) + "-inter-shift";
    case Neighborhoods::InterSwap:
        return "(" + std::to_string(k1) + "," + std::to_string(k2) + ")" + "-inter-swap";
    case Neighborhoods::InterShiftReverse:
        return std::to_string(k1) + "-inter-shift-reverse";
    case Neighborhoods::InterSwapStar:
        return "Inter-swap-star";
    case Neighborhoods::ShiftChangeMode:
        return "Shift-change-mode";
    case Neighborhoods::ModeSwap:
        return "Mode-swap";
    case Neighborhoods::SwapWithModes:
        return "Swap-with-modes";
    }
    return "";
}

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

    /*
     * Neighborhoods - Inter.
     */

    bool inter_two_opt = false;
    ElementPos inter_shift_block_maximum_length = 0;
    ElementPos inter_swap_block_maximum_length = 0;
    ElementPos inter_shift_reverse_block_maximum_length = 0;
    bool inter_swap_star = false;

    /*
     * Neighborhoods - Sub-sequence.
     */

    bool add_remove = false;

    /*
     * Neighborhoods - Modes.
     */

    bool shift_change_mode = false;
    bool mode_swap = false;
    bool swap_with_modes = false;

    /*
     * Perturbations.
     */

    Counter double_bridge_number_of_perturbations = 0;
    Counter ruin_and_recreate_number_of_perturbations = 0;
    ElementPos ruin_and_recreate_number_of_elements_removed = 4;
    bool force_add = false;

    /*
     * Crossovers.
     */

    double crossover_ox_weight = 0;
    double crossover_sjox_weight = 0;
    double crossover_sbox_weight = 0;
    double crossover_srex1_weight = 0;
    double crossover_srex2_weight = 0;
};

template <typename SequencingScheme>
class LocalScheme
{

public:

    using GlobalCost = typename SequencingScheme::GlobalCost;

    using SequenceData = typename SequencingScheme::SequenceData;

    struct SequenceElement
    {
        ElementId j;
        Mode mode;

        bool operator==(const SequenceElement& rhs) const
        {
            return j == rhs.j && mode == rhs.mode;
        }
    };

    struct Sequence
    {
        SequenceId i = -1;
        std::vector<SequenceElement> elements;
        SequenceData data;
    };

    struct SubSequence
    {
        SubSequence(
                const Sequence& sequence,
                ElementPos pos_1,
                ElementPos pos_2,
                bool reverse = false):
            sequence(sequence),
            pos_1(pos_1),
            pos_2(pos_2),
            reverse(reverse) { }

        const Sequence& sequence;
        ElementPos pos_1;
        ElementPos pos_2;
        bool reverse;
    };

    struct Move0
    {
        Neighborhoods type;
        ElementPos k1 = -1;
        ElementPos k2 = -1;

        SequenceId i1 = -1;
        SequenceId i2 = -1;
        ElementId j = -1;
        ElementPos pos_1 = -1;
        ElementPos pos_2 = -1;
        ElementPos pos_3 = -1;
        ElementPos pos_4 = -1;
        Mode mode = -1;

        GlobalCost global_cost = worst<GlobalCost>();
    };

    struct Neighborhood
    {
        std::vector<Move0> improving_moves = {};
        // modified_sequences[i] == true iff sequence i has changed since last
        // neighborhood exploration.
        std::vector<bool> modified_sequences;
        Counter number_of_explorations = 0;
        Counter number_of_successes = 0;
        double time = 0.0;
    };


    void append(
            Sequence& sequence,
            SequenceElement se) const
    {
        append(sequence.data, se);
        sequence.elements.push_back(se);
    }


    struct Solution
    {
        std::vector<Sequence> sequences;
        GlobalCost global_cost;

        std::vector<bool> modified_sequences;
    };

    using CompactSolution = std::vector<std::vector<SequenceElement>>;

    struct CompactSolutionHasher
    {
        std::hash<ElementId> hasher_j;
        std::hash<Mode> hasher_mode;

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
                for (const SequenceElement& se: sequence) {
                    optimizationtools::hash_combine(hash_tmp, hasher_j(se.j));
                    optimizationtools::hash_combine(hash_tmp, hasher_mode(se.mode));
                }
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
            for (const SequenceElement& se: compact_solution[i])
                append(solution.sequences[i], se);
        compute_global_cost(solution);
        return solution;
    }

    CompactSolution solution2compact(const Solution& solution)
    {
        CompactSolution compact_solution(number_of_sequences());
        for (SequenceId i = 0; i < number_of_sequences(); ++i)
            compact_solution[i] = solution.sequences[i].elements;
        return compact_solution;
    }

    /*
     * Constructors and destructor.
     */

    LocalScheme(
            SequencingScheme& sequencing_scheme,
            Parameters parameters = Parameters()):
        sequencing_scheme_(sequencing_scheme),
        parameters_(parameters)
    {
        SequencePos m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();

        // Initialize temporary structures.
        sequence_datas_cur_ = std::vector<std::vector<SequenceData>>(
                m, std::vector<SequenceData>(n + 1));
        sequence_datas_cur_2_ = std::vector<std::vector<std::vector<SequenceData>>>(m);
        modes_cur_ = std::vector<Mode>(n);
        solution_cur_ = empty_solution();
        solution_tmp_ = empty_solution();

        // Initialize neighborhoods_.

        neighborhoods_[int(Neighborhoods::Shift)]
            = std::vector<std::vector<Neighborhood>>(
                    parameters_.shift_block_maximum_length + 1,
                    {{Neighborhood()}});
        neighborhoods_[int(Neighborhoods::Swap)]
            = std::vector<std::vector<Neighborhood>>(
                    parameters_.swap_block_maximum_length + 1);
        for (ElementPos k2 = 1; k2 <= parameters_.swap_block_maximum_length; ++k2)
            neighborhoods_[int(Neighborhoods::Swap)][k2]
                = std::vector<Neighborhood>(
                        parameters_.swap_block_maximum_length + 1,
                        {Neighborhood()});
        neighborhoods_[int(Neighborhoods::Reverse)] = {{{Neighborhood()}}};
        neighborhoods_[int(Neighborhoods::ShiftReverse)]
            = std::vector<std::vector<Neighborhood>>(
                    parameters_.shift_reverse_block_maximum_length + 1,
                    {{Neighborhood()}});

        neighborhoods_[int(Neighborhoods::Add)] = {{{Neighborhood()}}};
        neighborhoods_[int(Neighborhoods::Remove)] = {{{Neighborhood()}}};

        neighborhoods_[int(Neighborhoods::InterTwoOpt)] = {{{Neighborhood()}}};
        neighborhoods_[int(Neighborhoods::InterShift)]
            = std::vector<std::vector<Neighborhood>>(
                    parameters_.inter_shift_block_maximum_length + 1,
                    {{Neighborhood()}});
        neighborhoods_[int(Neighborhoods::InterSwap)]
            = std::vector<std::vector<Neighborhood>>(
                    parameters_.swap_block_maximum_length + 1);
        for (ElementPos k2 = 1; k2 <= parameters_.swap_block_maximum_length; ++k2)
            neighborhoods_[int(Neighborhoods::InterSwap)][k2]
                = std::vector<Neighborhood>(
                        parameters_.swap_block_maximum_length + 1,
                        {Neighborhood()});
        neighborhoods_[int(Neighborhoods::InterShiftReverse)]
            = std::vector<std::vector<Neighborhood>>(
                    parameters_.inter_shift_reverse_block_maximum_length + 1,
                    {{Neighborhood()}});
        neighborhoods_[int(Neighborhoods::InterSwapStar)] = {{{Neighborhood()}}};

        neighborhoods_[int(Neighborhoods::ShiftChangeMode)] = {{{Neighborhood()}}};
        neighborhoods_[int(Neighborhoods::ModeSwap)] = {{{Neighborhood()}}};
        neighborhoods_[int(Neighborhoods::SwapWithModes)] = {{{Neighborhood()}}};

        for (int a = 0; a < (int)neighborhoods_.size(); ++a) {
            for (int k1 = 0; k1 < (int)neighborhoods_[a].size(); ++k1) {
                for (int k2 = 0; k2 < (int)neighborhoods_[a][k1].size(); ++k2) {
                    Neighborhood& neighborhood = neighborhoods_[a][k1][k2];
                    neighborhood.modified_sequences = std::vector<bool>(m, true);
                }
            }
        }
    }

    LocalScheme(const LocalScheme& sequencing_scheme):
        LocalScheme(sequencing_scheme.sequencing_scheme_, sequencing_scheme.parameters_) { }

    virtual ~LocalScheme() { }

    inline Sequence empty_sequence(SequenceId i) const
    {
        Sequence sequence;
        sequence.i = i;
        sequence.data = empty_sequence_data(i);
        return sequence;
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
        ElementPos n = sequencing_scheme_.number_of_elements();

        switch (initial_solution_id) {
        case 1: {
            // Fix a random order, find the best position.
            std::vector<ElementId> elements(n);
            std::iota(elements.begin(), elements.end(), 0);
            std::shuffle(elements.begin(), elements.end(), generator);
            Solution solution = empty_solution();
            compute_temporary_structures(solution);
            for (ElementId j: elements) {
                auto improving_moves = explore_add(solution, j);
                if (!improving_moves.empty()) {
                    std::shuffle(
                            improving_moves.begin(),
                            improving_moves.end(),
                            generator);
                    Move0 move_best;
                    move_best.global_cost = worst<GlobalCost>();
                    for (const Move0 move: improving_moves)
                        if (dominates(move.global_cost, move_best.global_cost))
                            move_best = move;
                    apply_move(solution, move_best);
                    compute_temporary_structures(solution, move_best.i1, move_best.i2);
                }
            }
            return solution;
        } case 2: {
            // Insert at the end, find the best element.
            std::vector<uint8_t> contains(sequencing_scheme_.number_of_elements(), 0);
            Solution solution = empty_solution();
            compute_temporary_structures(solution);
            for (;;) {
                auto improving_moves = explore_add_end(solution);
                if (!improving_moves.empty()) {
                    std::shuffle(
                            improving_moves.begin(),
                            improving_moves.end(),
                            generator);
                    Move0 move_best;
                    move_best.global_cost = worst<GlobalCost>();
                    for (const Move0 move: improving_moves)
                        if (dominates(move.global_cost, move_best.global_cost))
                            move_best = move;
                    apply_move(solution, move_best);
                    compute_temporary_structures(solution, move_best.i1, move_best.i2);
                } else {
                    break;
                }
            }
            return solution;
        } case 3: {
            // Find the best element x position.
            // Warning, this one can be expensive.
            std::vector<uint8_t> contains(sequencing_scheme_.number_of_elements(), 0);
            Solution solution = empty_solution();
            compute_temporary_structures(solution);
            for (;;) {
                auto improving_moves = explore_add_2(solution);
                if (!improving_moves.empty()) {
                    std::shuffle(
                            improving_moves.begin(),
                            improving_moves.end(),
                            generator);
                    Move0 move_best;
                    move_best.global_cost = worst<GlobalCost>();
                    for (const Move0 move: improving_moves)
                        if (dominates(move.global_cost, move_best.global_cost))
                            move_best = move;
                    apply_move(solution, move_best);
                    compute_temporary_structures(solution, move_best.i1, move_best.i2);
                } else {
                    break;
                }
            }
            return solution;
        } default: {
            // Random permutation.
            std::vector<ElementId> elements(n);
            std::iota(elements.begin(), elements.end(), 0);
            std::shuffle(elements.begin(), elements.end(), generator);
            std::uniform_int_distribution<SequenceId> di(0, m - 1);
            Solution solution = empty_solution();
            for (ElementId j: elements) {
                std::uniform_int_distribution<SequenceId> d_mode(0, number_of_modes(j) - 1);
                Mode mode = d_mode(generator);
                SequenceId i = di(generator);
                append(solution.sequences[i], {j, mode});
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
            return sequencing_scheme_.initial_solution(initial_solution_id, generator);
        return initial_solution(initial_solution_id, generator, false);
    }

    inline Solution initial_solution(
            Counter initial_solution_id,
            std::mt19937_64& generator)
    {
        auto begin = std::chrono::steady_clock::now();
        Solution solution = initial_solution(
                initial_solution_id,
                generator,
                std::integral_constant<
                    bool,
                    HasInitialSolutionMethod<SequencingScheme,
                    Solution(
                        Counter,
                        std::mt19937_64&)>::value>());
        auto end = std::chrono::steady_clock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
        initial_solution_time += time_span.count();
        number_of_initial_solution_calls++;
        return solution;
    }

    /*
     * Crossovers.
     */

    inline Solution crossover(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        auto begin = std::chrono::steady_clock::now();
        std::discrete_distribution<Counter> d_crossover({
                parameters_.crossover_ox_weight,
                parameters_.crossover_sjox_weight,
                parameters_.crossover_sbox_weight,
                parameters_.crossover_srex1_weight,
                parameters_.crossover_srex2_weight,
                });
        Counter x = d_crossover(generator);
        Solution solution;
        switch (x) {
        case 0: {
            solution = crossover_ox(solution_parent_1, solution_parent_2, generator);
            break;
        } case 1: {
            solution = crossover_sjox(solution_parent_1, solution_parent_2, generator);
            break;
        } case 2: {
            solution = crossover_sbox(solution_parent_1, solution_parent_2, generator);
            break;
        } case 3: {
            solution = crossover_srex1(solution_parent_1, solution_parent_2, generator);
            break;
        } case 4: {
            solution = crossover_srex2(solution_parent_1, solution_parent_2, generator);
            break;
        } default: {
            solution = crossover_ox(solution_parent_1, solution_parent_2, generator);
            break;
        }
        }
        auto end = std::chrono::steady_clock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
        crossover_time += time_span.count();
        number_of_crossover_calls++;
        return solution;
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

        const auto& elements_parent_1 = solution_parent_1.sequences[0].elements;
        const auto& elements_parent_2 = solution_parent_2.sequences[0].elements;
        Sequence& sequence = solution.sequences[0];

        ElementPos seq_1_size = elements_parent_1.size();

        std::vector<ElementPos> edges = optimizationtools::bob_floyd<ElementPos>(
                2, seq_1_size + 1, generator);
        std::sort(edges.begin(), edges.end());

        ElementPos pos_1 = edges[0];
        ElementPos pos_2 = edges[1];

        std::vector<uint8_t> in_substring(seq_1_size, false);
        for (ElementPos pos = pos_1; pos < pos_2; ++pos) {
            ElementId j = elements_parent_1[pos].j;
            in_substring[j] = true;
        }

        for (ElementPos pos = 0; pos < seq_1_size; ++pos) {
            if ((ElementPos)sequence.elements.size() == pos_1)
                for (ElementPos p = pos_1; p < pos_2; ++p)
                    append(sequence, elements_parent_1[p]);
            ElementId j = elements_parent_2[pos].j;
            if (in_substring[j])
                continue;
            append(sequence, elements_parent_2[pos]);
        }
        if ((ElementPos)sequence.elements.size() == pos_1)
            for (ElementPos p = pos_1; p < pos_2; ++p)
                append(sequence, elements_parent_1[p]);

        compute_global_cost(solution);
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

        const auto& elements_parent_1 = solution_parent_1.sequences[0].elements;
        const auto& elements_parent_2 = solution_parent_2.sequences[0].elements;
        Sequence& sequence = solution.sequences[0];

        ElementPos seq_1_size = elements_parent_1.size();
        std::vector<ElementPos> positions(seq_1_size, -1);

        // Add elements from parent_1 up to a given cut point.
        std::uniform_int_distribution<ElementPos> d_point(1, seq_1_size);
        ElementPos pos_0 = d_point(generator);
        for (ElementPos pos = 0; pos < pos_0; ++pos) {
            ElementId j = elements_parent_1[pos].j;
            positions[j] = pos;
            append(sequence, elements_parent_1[pos]);
        }

        // Add elements from parent_2 keeping the relative order.
        for (ElementPos pos = 0; pos < seq_1_size; ++pos) {
            // Add elements which have the same positions in both parents.
            for (;;) {
                ElementPos p = sequence.elements.size();
                if (p == seq_1_size)
                    break;
                ElementId j1 = elements_parent_1[p].j;
                ElementId j2 = elements_parent_2[p].j;
                if (j1 == j2) {
                    positions[j1] = p;
                    append(sequence, elements_parent_1[p]);
                    continue;
                }
                break;
            }
            ElementId j = elements_parent_2[pos].j;
            if (positions[j] != -1)
                continue;
            positions[j] = sequence.elements.size();
            append(sequence, elements_parent_2[pos]);
        }

        compute_global_cost(solution);
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

        const auto& elements_parent_1 = solution_parent_1.sequences[0].elements;
        const auto& elements_parent_2 = solution_parent_2.sequences[0].elements;
        Sequence& sequence = solution.sequences[0];

        ElementPos seq_1_size = elements_parent_1.size();
        std::vector<ElementPos> positions(seq_1_size, -1);

        // Add elements from parent_1 up to a given cut point.
        std::uniform_int_distribution<ElementPos> d_point(1, seq_1_size);
        ElementPos pos_0 = d_point(generator);
        for (ElementPos pos = 0; pos < pos_0; ++pos) {
            ElementId j = elements_parent_1[pos].j;
            positions[j] = pos;
            append(sequence, elements_parent_1[pos]);
        }

        // Add elements from parent_2 keeping the relative order.
        for (ElementPos pos = 0; pos < seq_1_size; ++pos) {
            // Add elements which have the same positions in both parents.
            for (;;) {
                ElementPos p = sequence.elements.size();
                if (p == seq_1_size)
                    break;
                ElementId j1 = elements_parent_1[p].j;
                ElementId j2 = elements_parent_2[p].j;
                if (p <= seq_1_size - 1) {
                    ElementId j1_next = elements_parent_1[p + 1].j;
                    ElementId j2_next = elements_parent_2[p + 1].j;
                    if (j1 == j2 && j1_next == j2_next) {
                        positions[j1] = p;
                        append(sequence, elements_parent_1[p + 1]);
                        continue;
                    }
                }
                if (p >= 1) {
                    ElementId j1_prev = elements_parent_1[p - 1].j;
                    ElementId j2_prev = elements_parent_2[p - 1].j;
                    if (j1 == j2 && j1_prev == j2_prev) {
                        positions[j1] = p;
                        append(sequence, elements_parent_1[p - 1]);
                        continue;
                    }
                }
                break;
            }
            ElementId j = elements_parent_2[pos].j;
            if (positions[j] != -1)
                continue;
            positions[j] = sequence.elements.size();
            append(sequence, elements_parent_2[pos]);
        }

        compute_global_cost(solution);
        return solution;
    }

    inline Solution crossover_srex1(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        //std::cout << "crossover_srex1 start" << std::endl;
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        std::uniform_real_distribution<double> d01(0, 1);

        Solution solution = empty_solution();
        solution.modified_sequences = std::vector<bool>(m, false);

        // Draw sequences to inherite from the first parent.
        std::vector<bool> sequences(m, false);
        std::vector<bool> added(n, false);
        ElementId n_cur = 0;

        // Add sequences from the first parent.
        for (SequenceId i = 0; i < m; ++i) {
            if (d01(generator) > 0.5)
                continue;
            solution.sequences[i] = solution_parent_1.sequences[i];
            sequences[i] = true;
            for (const auto& se: solution.sequences[i].elements) {
                added[se.j] = true;
                n_cur++;
            }
            //std::cout << "i " << to_string(global_cost(solution.sequences[i])) << std::endl;
        }
        //std::cout << "n_cur_1 " << n_cur << " / " << n << std::endl;

        // Add sequences from the second parent.
        for (SequenceId i = 0; i < m; ++i) {
            if (sequences[i])
                continue;
            // Find the sequence from the second parent with the highest number
            // of unscheduled elements.
            SequenceId i2_best = -1;
            ElementPos seq_size_best = -1;
            for (SequenceId i2 = 0; i2 < m; ++i2) {
                ElementPos seq_size = 0;
                for (const auto& se: solution_parent_2.sequences[i2].elements)
                    if (!added[se.j])
                        seq_size++;
                if (seq_size <= (SequencePos)solution_parent_2.sequences[i2].elements.size() / 2)
                    continue;
                if (i2_best == 0 || seq_size_best < seq_size) {
                    i2_best = i2;
                    seq_size_best = seq_size;
                }
            }
            if (i2_best == -1)
                break;
            //std::cout << "i2_best " << i2_best << " seq_size_best " << seq_size_best << std::endl;

            solution.sequences[i] = empty_sequence(i);
            for (const auto& se: solution_parent_2.sequences[i2_best].elements) {
                if (added[se.j]) {
                    solution.modified_sequences[i] = true;
                    continue;
                }
                append(solution.sequences[i], se);
                added[se.j] = true;
                n_cur++;
            }
            //std::cout << "i " << i
            //    << " gc " << to_string(global_cost(solution.sequences[i]))
            //    << std::endl;
        }
        //std::cout << "n_cur_2 " << n_cur << " / " << n << std::endl;

        // Add remaining elements at random positions.
        compute_temporary_structures(solution);
        for (ElementId j = 0; j < n; ++j) {
            if (added[j])
                continue;
            auto improving_moves = explore_add(solution, j);
            if (!improving_moves.empty()) {
                std::shuffle(
                        improving_moves.begin(),
                        improving_moves.end(),
                        generator);
                Move0 move_best;
                move_best.global_cost = worst<GlobalCost>();
                for (const Move0 move: improving_moves)
                    if (dominates(move.global_cost, move_best.global_cost))
                        move_best = move;
                apply_move(solution, move_best);
                solution.modified_sequences[move_best.i1] = true;
                compute_temporary_structures(solution, move_best.i1, move_best.i2);
            }
        }

        compute_global_cost(solution);
        //std::cout << "crossover_srex1 end" << std::endl;
        return solution;
    }

    inline Solution crossover_srex2(
            const Solution& solution_parent_1,
            const Solution& solution_parent_2,
            std::mt19937_64& generator)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        std::uniform_real_distribution<double> d01(0, 1);

        Solution solution = empty_solution();
        solution.modified_sequences = std::vector<bool>(m, false);

        // Draw sequences to inherite from the first parent.
        std::vector<bool> sequences(m, false);
        std::vector<bool> added(n, false);

        // Add sequences from the first parent.
        for (SequenceId i = 0; i < m; ++i) {
            if (d01(generator) > 0.5)
                continue;
            solution.sequences[i] = solution_parent_1.sequences[i];
            sequences[i] = true;
            for (const auto& se: solution.sequences[i].elements)
                added[se.j] = true;
        }

        // Add sequences from the second parent.
        for (SequenceId i = 0; i < m; ++i) {
            if (sequences[i])
                continue;
            solution.sequences[i] = empty_sequence(i);
            for (const auto& se: solution_parent_2.sequences[i].elements) {
                if (added[se.j]) {
                    solution.modified_sequences[i] = true;
                    continue;
                }
                append(solution.sequences[i], se);
                added[se.j] = true;
            }
        }

        // Add remaining elements at random positions.
        for (ElementId j = 0; j < n; ++j) {
            if (added[j])
                continue;
            auto improving_moves = explore_add(solution, j);
            if (!improving_moves.empty()) {
                std::shuffle(
                        improving_moves.begin(),
                        improving_moves.end(),
                        generator);
                Move0 move_best;
                move_best.global_cost = worst<GlobalCost>();
                for (const Move0 move: improving_moves)
                    if (dominates(move.global_cost, move_best.global_cost))
                        move_best = move;
                apply_move(solution, move_best);
                solution.modified_sequences[move_best.i1] = true;
                compute_temporary_structures(solution, move_best.i1, move_best.i2);
            }
        }

        compute_global_cost(solution);
        return solution;
    }

    /*
     * Solution properties.
     */

    inline GlobalCost global_cost(const Solution& solution) const
    {
        return solution.global_cost;
    }

    inline GlobalCost global_cost(const SequenceData& sequence_data) const
    {
        return sequencing_scheme_.global_cost(sequence_data);
    }

    inline GlobalCost global_cost(const Sequence& sequence) const
    {
        return sequencing_scheme_.global_cost(sequence.data);
    }

    inline GlobalCost bound(const SequenceData& sequence_data) const
    {
        return sequencing_scheme_.bound(sequence_data);
    }

    inline GlobalCost bound(const Sequence& sequence) const
    {
        return sequencing_scheme_.bound(sequence.data);
    }

    GlobalCost global_cost_goal(double value) const
    {
        return localsearchsolver::global_cost_goal(sequencing_scheme_, value);
    }

    /**
     * Return the distance between two solutions.
     *
     * The distance between two solutions (represented as element permutations)
     * is defined as the number of different edges between them (broken-pairs
     * distance).
     *
     * References:
     * - "Two memetic algorithms for heterogeneous fleet vehicle routing
     *   problems" (Prins, 2009)
     *   https://doi.org/10.1016/j.engappai.2008.10.006
     * - "A simple and effective hybrid genetic search for the job sequencing
     *   and tool switching problem" (Mecler et al., 2021)
     *   https://doi.org/10.1016/j.cor.2020.105153
     */
    inline ElementPos distance(
            const Solution& solution_1,
            const Solution& solution_2) const
    {
        SequencePos m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();

        // For each element, get its next element in the first solution.
        std::vector<ElementId> next_1(n, -1);
        for (SequenceId i = 0; i < m; ++i) {
            const auto& elements = solution_1.sequences[i].elements;
            ElementPos seq_size = elements.size();
            for (ElementPos pos = 0; pos < seq_size - 1; ++pos) {
                ElementId j = elements[pos].j;
                ElementId j_next = elements[pos + 1].j;
                next_1[j] = j_next;
            }
        }

        // Compare with second solution.
        ElementPos d = 0;
        for (SequenceId i = 0; i < m; ++i) {
            const auto& elements = solution_2.sequences[i].elements;
            ElementPos seq_size = elements.size();
            for (ElementPos pos = 0; pos < seq_size - 1; ++pos) {
                ElementId j = elements[pos].j;
                ElementId j_next = elements[pos + 1].j;
                if (j_next != next_1[j])
                    d++;
            }
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
        /** ForceAdd: element to add. */
        ElementId force_add_j = -1;
        /** Global cost of the move. */
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
            std::vector<uint8_t> contains(sequencing_scheme_.number_of_elements(), 0);
            for (SequenceId i = 0; i < number_of_sequences(); ++i)
                for (SequenceElement se: solution.sequences[i].elements)
                    contains[se.j] = 1;
            for (ElementId j = 0; j < sequencing_scheme_.number_of_elements(); ++j) {
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
            Move& move,
            std::mt19937_64& generator)
    {
        SequenceId m = number_of_sequences();

        switch (move.type) {
        case Perturbations::None: {
            break;

        } case Perturbations::DoubleBridge: {
            std::uniform_int_distribution<SequenceId> distribution(0, m - 1);
            SequenceId i = distribution(generator);
            auto positions = optimizationtools::bob_floyd<ElementPos>(
                    4,
                    solution.sequences[i].elements.size() + 1,
                    generator);
            std::sort(positions.begin(), positions.end());

            const auto& elements = solution.sequences[i].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[i];
            sequence_tmp = empty_sequence(i);
            for (ElementPos pos = 0; pos < positions[0]; ++pos)
                append(sequence_tmp, elements[pos]);
            for (ElementPos pos = positions[2]; pos < positions[3]; ++pos)
                append(sequence_tmp, elements[pos]);
            for (ElementPos pos = positions[1]; pos < positions[2]; ++pos)
                append(sequence_tmp, elements[pos]);
            for (ElementPos pos = positions[0]; pos < positions[1]; ++pos)
                append(sequence_tmp, elements[pos]);
            for (ElementPos pos = positions[3]; pos < seq_size; ++pos)
                append(sequence_tmp, elements[pos]);
            compute_global_cost(solution_tmp_);
            solution = solution_tmp_;
            break;

        } case Perturbations::RuinAndRecreate: {
            ElementId n = sequencing_scheme_.number_of_elements();
            std::vector<SequenceId> sequences(n, -1);
            std::vector<Mode> modes(n, -1);
            std::vector<ElementPos> positions(n, -1);
            optimizationtools::IndexedSet elts(n);
            for (SequenceId i = 0; i < m; ++i) {
                const auto& elements = solution.sequences[i].elements;
                ElementPos seq_size = (ElementPos)elements.size();
                for (ElementPos pos = 0; pos < seq_size; ++pos) {
                    ElementId j = elements[pos].j;
                    elts.add(j);
                    sequences[j] = i;
                    modes[j] = elements[pos].mode;
                    positions[j] = pos;
                }
            }

            Solution solution_cur_ = empty_solution();
            solution_cur_.modified_sequences = std::vector<bool>(m, false);

            ElementPos number_of_elements_removed = std::min(
                    parameters_.ruin_and_recreate_number_of_elements_removed,
                    elts.size());
            elts.shuffle_in(generator);
            while (elts.size() > number_of_elements_removed) {
                ElementId j = *(elts.begin() + (elts.size() - 1));
                elts.remove(j);
                solution_cur_.modified_sequences[sequences[j]] = true;
            }
            for (SequenceId i = 0; i < m; ++i) {
                const auto& elements = solution.sequences[i].elements;
                SequencePos seq_size = elements.size();
                Sequence& sequence_cur = solution_cur_.sequences[i];
                for (ElementPos pos = 0; pos < seq_size; ++pos) {
                    ElementId j = elements[pos].j;
                    if (!elts.contains(j))
                        append(sequence_cur, elements[pos]);
                }
            }
            // Add back removed elements.
            compute_temporary_structures(solution_cur_);
            for (ElementId j: elts) {
                ElementId j_prec_old = -1;
                if (positions[j] > 0) {
                    const auto& elements = solution.sequences[sequences[j]].elements;
                    j_prec_old = elements[positions[j] - 1].j;
                }
                auto improving_moves = explore_add(
                        solution_cur_,
                        j,
                        sequences[j],
                        j_prec_old,
                        modes[j]);
                if (!improving_moves.empty()) {
                    std::shuffle(
                            improving_moves.begin(),
                            improving_moves.end(),
                            generator);
                    Move0 move_best;
                    move_best.global_cost = worst<GlobalCost>();
                    for (const Move0 move: improving_moves)
                        if (dominates(move.global_cost, move_best.global_cost))
                            move_best = move;
                    apply_move(solution_cur_, move_best);
                    solution_cur_.modified_sequences[move_best.i1] = true;
                    compute_temporary_structures(solution_cur_, move_best.i1, move_best.i2);
                }
            }
            solution = solution_cur_;
            break;

        } case Perturbations::ForceAdd: {
            ElementId j = move.force_add_j;
            // Draw sequence.
            std::uniform_int_distribution<SequenceId> d_i(0, m - 1);
            SequenceId i = d_i(generator);
            // Draw position in sequence.
            const auto& elements = solution.sequences[i].elements;
            SequencePos seq_size = elements.size();
            std::uniform_int_distribution<SequenceId> d_pos(0, seq_size);
            ElementPos pos = d_pos(generator);
            // Draw mode.
            std::uniform_int_distribution<SequenceId> d_mode(0, number_of_modes(j) - 1);
            Mode mode = d_mode(generator);
            // Apply move.
            Move0 move;
            move.type = Neighborhoods::Add;
            move.i1 = i;
            move.j = j;
            move.mode = mode;
            move.pos_1 = pos;
            apply_move(solution, move);
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

        auto local_search_begin = std::chrono::steady_clock::now();
        SequencePos m = number_of_sequences();
        //ElementPos n = sequencing_scheme_.number_of_elements();

        // Get neighborhoods.
        std::vector<std::tuple<Neighborhoods, ElementPos, ElementPos>> neighborhoods;
        for (ElementPos block_size = 1;
                block_size <= parameters_.shift_block_maximum_length;
                ++block_size) {
            neighborhoods.push_back({Neighborhoods::Shift, block_size, 0});
        }
        for (ElementPos block_size_1 = 1;
                block_size_1 <= parameters_.swap_block_maximum_length;
                ++block_size_1) {
            for (ElementPos block_size_2 = 1;
                    block_size_2 <= block_size_1;
                    ++block_size_2) {
                neighborhoods.push_back({
                        Neighborhoods::Swap, block_size_1, block_size_2});
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
            for (ElementPos block_size = 2;
                    block_size <= parameters_.inter_shift_reverse_block_maximum_length;
                    ++block_size) {
                neighborhoods.push_back({
                        Neighborhoods::InterShiftReverse, block_size, 0});
            }
            if (parameters_.inter_swap_star)
                neighborhoods.push_back({Neighborhoods::InterSwapStar, 0, 0});
        }

        if (parameters_.shift_change_mode)
            neighborhoods.push_back({Neighborhoods::ShiftChangeMode, 0, 0});
        if (parameters_.mode_swap)
            neighborhoods.push_back({Neighborhoods::ModeSwap, 0, 0});
        if (parameters_.swap_with_modes)
            neighborhoods.push_back({Neighborhoods::SwapWithModes, 0, 0});

        for (int a = 0; a < (int)neighborhoods_.size(); ++a) {
            for (int k1 = 0; k1 < (int)neighborhoods_[a].size(); ++k1) {
                for (int k2 = 0; k2 < (int)neighborhoods_[a][k1].size(); ++k2) {
                    std::string s = neighborhood2string((Neighborhoods)a, k1, k2);
                    Neighborhood& neighborhood = neighborhoods_[a][k1][k2];
                    if (solution.modified_sequences.empty()) {
                        std::fill(
                                neighborhood.modified_sequences.begin(),
                                neighborhood.modified_sequences.end(),
                                true);
                    } else {
                        for (SequenceId i = 0; i < m; ++i)
                            neighborhood.modified_sequences[i] = solution.modified_sequences[i];
                    }
                }
            }
        }

        compute_temporary_structures(solution);

        Counter it = 0;
        for (;; ++it) {
            //std::cout << "it " << it
            //    << " c " << to_string(global_cost(solution))
            //    << std::endl;
            //print(std::cout, solution);

            std::vector<GlobalCost> gcis(m);
            for (SequenceId i = 0; i < m; ++i)
                gcis[i] = global_cost(solution.sequences[i]);

            std::shuffle(neighborhoods.begin(), neighborhoods.end(), generator);

            bool improved = false;
            // Loop through neighborhoods.
            for (auto neighborhood_id: neighborhoods) {
                auto begin = std::chrono::steady_clock::now();
                Neighborhoods type = std::get<0>(neighborhood_id);
                ElementPos k1 = std::get<1>(neighborhood_id);
                ElementPos k2 = std::get<2>(neighborhood_id);
                //std::cout << neighborhood2string(type, k1, k2) << std::endl;
                Neighborhood& neighborhood = neighborhoods_[(int)type][k1][k2];
                // Remove moves which have changed from improving_moves.
                for (auto it = neighborhood.improving_moves.begin();
                        it != neighborhood.improving_moves.end();) {
                    if (neighborhood.modified_sequences[it->i1]
                            || (it->i2 != -1 && neighborhood.modified_sequences[it->i2])) {
                        *it = neighborhood.improving_moves.back();
                        neighborhood.improving_moves.pop_back();
                    } else {
                        it++;
                    }
                }

                switch (std::get<0>(neighborhood_id)) {
                case Neighborhoods::Shift:
                    explore_shift(solution, k1);
                    break;
                case Neighborhoods::Swap:
                    explore_swap(solution, k1, k2);
                    break;
                case Neighborhoods::Reverse:
                    explore_reverse(solution);
                    break;
                case Neighborhoods::ShiftReverse:
                    explore_shift(solution, k1, true);
                    break;
                case Neighborhoods::Add:
                    explore_add(solution);
                    break;
                case Neighborhoods::Remove:
                    explore_remove(solution, perturbation);
                    break;
                case Neighborhoods::ShiftChangeMode:
                    explore_shift_change_mode(solution);
                    break;
                case Neighborhoods::ModeSwap:
                    explore_mode_swap(solution);
                    break;
                case Neighborhoods::SwapWithModes:
                    explore_swap_with_modes(solution);
                    break;
                case Neighborhoods::InterTwoOpt:
                    explore_inter_two_opt(solution);
                    break;
                case Neighborhoods::InterShift:
                    explore_inter_shift(solution, k1);
                    break;
                case Neighborhoods::InterSwap:
                    explore_inter_swap(solution, k1, k2);
                    break;
                case Neighborhoods::InterShiftReverse:
                    explore_inter_shift(solution, k1, true);
                    break;
                case Neighborhoods::InterSwapStar:
                    explore_inter_swap_star(solution);
                    break;
                }

                // Update neighborhood.
                std::fill(
                        neighborhood.modified_sequences.begin(),
                        neighborhood.modified_sequences.end(),
                        false);
                neighborhood.number_of_explorations++;

                if (!neighborhood.improving_moves.empty()) {
                    //std::cout << neighborhood2string(type, k1, k2) << std::endl;
                    improved = true;
                    neighborhood.number_of_successes++;
                    std::shuffle(
                            neighborhood.improving_moves.begin(),
                            neighborhood.improving_moves.end(),
                            generator);
                    Move0 move_best;
                    move_best.global_cost = worst<GlobalCost>();
                    for (const Move0 move: neighborhood.improving_moves)
                        if (dominates(move.global_cost, move_best.global_cost))
                            move_best = move;
                    apply_move(solution, move_best);
                    compute_temporary_structures(solution, move_best.i1, move_best.i2);
                    // Check new current solution cost.
                    GlobalCost gci = (move_best.i2 == -1)?
                        gcis[move_best.i1]:
                        global_cost_merge(gcis[move_best.i1], gcis[move_best.i2]);
                    GlobalCost gci_new = (move_best.i2 == -1)?
                        global_cost(solution.sequences[move_best.i1]):
                        global_cost_merge(
                                global_cost(solution.sequences[move_best.i1]),
                                global_cost(solution.sequences[move_best.i2]));
                    if (gci_new - gci != move_best.global_cost) {
                        //std::cout
                        //    << "k1 " << move_best.k1
                        //    << " k2 " << move_best.k2
                        //    << " i1 " << move_best.i1
                        //    << " i2 " << move_best.i2
                        //    << " pos_1 " << move_best.pos_1
                        //    << " pos_2 " << move_best.pos_2
                        //    << std::endl;
                        throw std::logic_error(
                                neighborhood2string(type, k1, k2)
                                + ", costs do not match:\n"
                                + "* Old cost: "
                                + to_string(gci) + "\n"
                                + "* Move cost: "
                                + to_string(move_best.global_cost) + "\n"
                                + "* New cost: "
                                + to_string(gci_new) + "\n");
                    }

                    // Update modified sequences.
                    for (int a = 0; a < (int)neighborhoods_.size(); ++a) {
                        for (int k1 = 0; k1 < (int)neighborhoods_[a].size(); ++k1) {
                            for (int k2 = 0; k2 < (int)neighborhoods_[a][k1].size(); ++k2) {
                                Neighborhood& neighborhood_2 = neighborhoods_[a][k1][k2];
                                neighborhood_2.modified_sequences[move_best.i1] = true;
                                if (move_best.i2 != -1)
                                    neighborhood_2.modified_sequences[move_best.i2] = true;
                            }
                        }
                    }
                }
                auto end = std::chrono::steady_clock::now();
                auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - begin);
                neighborhood.time += time_span.count();

                if (improved)
                    break;
            }

            if (!improved)
                break;
        }

        //print(std::cout, solution);
        //std::cout << to_string(global_cost(solution)) << std::endl;
        auto local_search_end = std::chrono::steady_clock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(local_search_end - local_search_begin);
        local_search_time += time_span.count();
        number_of_local_search_calls++;
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
            for (SequenceElement se: solution.sequences[i].elements) {
                if (maximum_number_of_modes_ == 1) {
                    os << " " << se.j;
                } else {
                    os << " " << se.j << " " << se.mode;
                }
            }
            os << std::endl;
            os << "    Cost: " << to_string(global_cost(solution.sequences[i])) << std::endl;
        }
        os << "Total cost: " << to_string(solution.global_cost) << std::endl;
        return os;
    }

    inline void write(
            const Solution& solution,
            std::string certificate_path,
            ElementPos offset = 0) const
    {
        if (certificate_path.empty())
            return;
        std::ofstream file(certificate_path);
        if (!file.good()) {
            throw std::runtime_error(
                    "Unable to open file \"" + certificate_path + "\".");
        }

        SequenceId m = number_of_sequences();
        if (m == 1) {
            for (SequenceElement se: solution.sequences[0].elements) {
                if (maximum_number_of_modes_ == 1) {
                    file << se.j + offset << " ";
                } else {
                    file << se.j + offset << " " << se.mode << " ";
                }
            }
        } else {
            for (SequenceId i = 0; i < m; ++i) {
                file << solution.sequences[i].elements.size() << std::endl;
                for (SequenceElement se: solution.sequences[i].elements) {
                    if (maximum_number_of_modes_ == 1) {
                        file << se.j + offset << " ";
                    } else {
                        file << se.j + offset << " " << se.mode << " ";
                    }
                }
                file << std::endl;
            }
        }
    }

    void print_parameters(
            optimizationtools::Info& info) const
    {
        SequencePos m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        info.os()
            << "Number of sequences:                           " << m << std::endl
            << "Number of elements:                            " << n << std::endl
            << "Maximum number of modes:                       " << maximum_number_of_modes_ << std::endl
            << "Neighborhoods" << std::endl
            << "    Shift" << std::endl
            << "        Block maximum length:                  " << parameters_.shift_block_maximum_length << std::endl
            << "    Swap" << std::endl
            << "        Block maximum length:                  " << parameters_.swap_block_maximum_length << std::endl
            << "    Reverse:                                   " << parameters_.reverse << std::endl
            << "    Shift-reverse" << std::endl
            << "        Block maximum length:                  " << parameters_.shift_reverse_block_maximum_length << std::endl
            << "    Add/Remove:                                " << parameters_.add_remove << std::endl;
        if (m > 1) {
            info.os()
                << "    Inter-shift" << std::endl
                << "        Block maximum length:                  " << parameters_.inter_shift_block_maximum_length << std::endl
                << "    Inter-swap" << std::endl
                << "        Block maximum length:                  " << parameters_.inter_swap_block_maximum_length << std::endl
                << "    Inter-two-opt:                             " << parameters_.inter_two_opt << std::endl
                << "    Inter-shift-reverse" << std::endl
                << "        Block maximum length:                  " << parameters_.inter_shift_reverse_block_maximum_length << std::endl
                << "    Inter-swap-star:                           " << parameters_.inter_swap_star << std::endl;
        }
        if (maximum_number_of_modes_ > 1) {
            info.os()
                << "    Shfit-change-mode:                         " << parameters_.shift_change_mode << std::endl
                << "    Mode-swap:                                 " << parameters_.mode_swap << std::endl
                << "    Swap-with-modes:                           " << parameters_.swap_with_modes << std::endl;
        }
        info.os()
            << "Perturbations" << std::endl
            << "    Double-bridge" << std::endl
            << "        Number of perturbations:               " << parameters_.double_bridge_number_of_perturbations << std::endl
            << "    Ruin-and-recreate" << std::endl
            << "        Number of perturbations:               " << parameters_.ruin_and_recreate_number_of_perturbations << std::endl
            << "        Number of elements removed:            " << parameters_.ruin_and_recreate_number_of_elements_removed << std::endl
            << "    Force-add:                                 " << parameters_.force_add << std::endl;
        info.os()
            << "Crossovers" << std::endl
            << "    OX:                                        " << parameters_.crossover_ox_weight << std::endl
            << "    SBOX:                                      " << parameters_.crossover_sbox_weight << std::endl
            << "    SJOX:                                      " << parameters_.crossover_sjox_weight << std::endl
            << "    SREX1:                                     " << parameters_.crossover_srex1_weight << std::endl
            << "    SREX2:                                     " << parameters_.crossover_srex2_weight << std::endl;
    }

    void print_statistics(
            optimizationtools::Info& info) const
    {
        info.os() << "General" << std::endl;
        info.os()
            << "    "
            << std::left << std::setw(28) << "Initial solution:"
            << number_of_initial_solution_calls
            << " / " << initial_solution_time << "s"
            << " / " << initial_solution_time / number_of_initial_solution_calls << "s"
            << std::endl;
        info.os()
            << "    "
            << std::left << std::setw(28) << "Crossover time:"
            << number_of_crossover_calls
            << " / " << crossover_time << "s"
            << " / " << crossover_time / number_of_crossover_calls << "s"
            << std::endl;
        info.os()
            << "    "
            << std::left << std::setw(28) << "Local search time:"
            << number_of_local_search_calls
            << " / " << local_search_time << "s"
            << " / " << local_search_time / number_of_local_search_calls << "s"
            << std::endl;
        info.os() << "Neighborhoods" << std::endl;
        for (int a = 0; a < (int)neighborhoods_.size(); ++a) {
            for (int k1 = 0; k1 < (int)neighborhoods_[a].size(); ++k1) {
                for (int k2 = 0; k2 < (int)neighborhoods_[a][k1].size(); ++k2) {
                    std::string s = neighborhood2string((Neighborhoods)a, k1, k2);
                    const Neighborhood& neighborhood = neighborhoods_[a][k1][k2];
                    if (neighborhood.number_of_explorations == 0)
                        continue;
                    info.os()
                        << "    "
                        << std::left << std::setw(28) << s
                        << neighborhood.number_of_explorations
                        << " / " << neighborhood.number_of_successes
                        << " / " << (double)neighborhood.number_of_successes / neighborhood.number_of_explorations * 100 << "%"
                        << " / " << neighborhood.time << "s"
                        << std::endl;
                    info.add_to_json(
                            "Algorithm",
                            s + "NumberOfExplorations",
                            neighborhood.number_of_explorations);
                    info.add_to_json(
                            "Algorithm",
                            s + "NumberOfSuccesses",
                            neighborhood.number_of_successes);

                }
            }
        }
    }

private:

    /*
     * SequencePos number_of_sequences()
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
        return sequencing_scheme_.number_of_sequences();;
    }

    SequencePos number_of_sequences() const
    {
        return number_of_sequences(
                std::integral_constant<
                    bool,
                    HasNumberOfSequencesMethod<SequencingScheme, SequencePos()>::value>());
    }

    /*
     * Mode number_of_modes(ElementId)
     */

    template<typename, typename T>
    struct HasNumberOfModesMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasNumberOfModesMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().number_of_modes(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    SequencePos number_of_modes(
            ElementId,
            std::false_type) const
    {
        return 1;
    }

    SequencePos number_of_modes(
            ElementId j,
            std::true_type) const
    {
        return sequencing_scheme_.number_of_modes(j);;
    }

    Mode number_of_modes(ElementId j) const
    {
        return number_of_modes(
                j,
                std::integral_constant<
                    bool,
                    HasNumberOfModesMethod<SequencingScheme, Mode(ElementId)>::value>());
    }

    /*
     * void append(SequenceData&, ElementId)
     * void append(SequenceData&, ElementId, ModeId)
     */

    template<typename, typename T>
    struct HasAppendMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasAppendMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().append(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    void append(
            SequenceData& sequence_data,
            SequenceElement se,
            std::false_type) const
    {
        sequencing_scheme_.append(sequence_data, se.j);;
    }

    void append(
            SequenceData& sequence_data,
            SequenceElement se,
            std::true_type) const
    {
        sequencing_scheme_.append(sequence_data, se.j, se.mode);;
    }

    void append(
            SequenceData& sequence_data,
            SequenceElement se) const
    {
        return append(
                sequence_data,
                se,
                std::integral_constant<
                    bool,
                    HasAppendMethod<SequencingScheme, void(SequenceData&, ElementId, Mode)>::value>());
    }

    /*
     * SequenceData empty_sequence_data(SequenceId)
     */

    template<typename, typename T>
    struct HasEmptySequenceDataMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasEmptySequenceDataMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().empty_sequence_data(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    SequenceData empty_sequence_data(
            SequenceId,
            std::false_type) const
    {
        return SequenceData();
    }

    SequenceData empty_sequence_data(
            SequenceId i,
            std::true_type) const
    {
        return sequencing_scheme_.empty_sequence_data(i);;
    }

    SequenceData empty_sequence_data(SequenceId i) const
    {
        return empty_sequence_data(
                i,
                std::integral_constant<
                    bool,
                    HasEmptySequenceDataMethod<SequencingScheme, SequenceData(SequenceId)>::value>());
    }

    /*
     * global_cost_merge(
     *         const GlobalCost&,
     *         const GlobalCost&)
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
        return sequencing_scheme_.global_cost_merge(gc1, gc2);
    }

    /**
     * Return the merged global cost of two global costs.
     *
     * If SequencingScheme does not implement a method
     * 'GlobalCost global_cost(const GlobalCost&, const GlobalCost&)',
     * then it returns the sum of the two global costs, i.e., the global cost
     * of two sequences is the sum of the global cost of each sequence. If it
     * is not the case, for example when the objective is the makespan,
     * SequencingScheme needs to implement a custom method.
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
                    HasGlobalCostMergeMethod<SequencingScheme, GlobalCost(
                        const GlobalCost&,
                        const GlobalCost&)>::value>());
    }

    /*
     * bool concatenate(
     *         SequenceData&,
     *         const SequenceData&)
     */

    template<typename, typename T>
    struct HasConcatenateMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasConcatenateMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().concatenate(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    bool concatenate(
            SequenceData& sequence_data,
            const SubSequence& sub_sequence,
            const GlobalCost& cutoff,
            const GlobalCost* gci1,
            std::false_type) const
    {
        if (gci1 == nullptr) {
            if (!sub_sequence.reverse) {
                auto it = sub_sequence.sequence.elements.begin() + sub_sequence.pos_1;
                auto it_end = sub_sequence.sequence.elements.begin() + sub_sequence.pos_2;
                for (;; ++it) {
                    // Add next element to sequence_data.
                    append(sequence_data, *it);
                    // Check early termination.
                    if (bound(sequence_data) >= cutoff)
                        return false;
                    // End condition.
                    if (it == it_end)
                        return true;
                }
            } else {
                auto it = sub_sequence.sequence.elements.begin() + sub_sequence.pos_2;
                auto it_end = sub_sequence.sequence.elements.begin() + sub_sequence.pos_1;
                for (;; --it) {
                    // Add next element to sequence_data.
                    append(sequence_data, *it);
                    // Check early termination.
                    if (bound(sequence_data) >= cutoff)
                        return false;
                    // End condition.
                    if (it == it_end)
                        return true;
                }
            }
        } else {
            if (!sub_sequence.reverse) {
                auto it = sub_sequence.sequence.elements.begin() + sub_sequence.pos_1;
                auto it_end = sub_sequence.sequence.elements.begin() + sub_sequence.pos_2;
                for (;; ++it) {
                    // Add next element to sequence_data.
                    append(sequence_data, *it);
                    // Check early termination.
                    if (global_cost_merge(*gci1, bound(sequence_data)) >= cutoff)
                        return false;
                    // End condition.
                    if (it == it_end)
                        return true;
                }
            } else {
                auto it = sub_sequence.sequence.elements.begin() + sub_sequence.pos_2;
                auto it_end = sub_sequence.sequence.elements.begin() + sub_sequence.pos_1;
                for (;; --it) {
                    // Add next element to sequence_data.
                    append(sequence_data, *it);
                    // Check early termination.
                    if (global_cost_merge(*gci1, bound(sequence_data)) >= cutoff)
                        return false;
                    // End condition.
                    if (it == it_end)
                        return true;
                }
            }
        }
    }

    bool concatenate(
            SequenceData& sequence_data,
            const SubSequence& sub_sequence,
            const GlobalCost&,
            const GlobalCost*,
            std::true_type) const
    {
        const SequenceData& sequence_data_2 = (!sub_sequence.reverse)?
            sequence_datas_cur_2_[sub_sequence.sequence.i][sub_sequence.pos_1][sub_sequence.pos_2]:
            sequence_datas_cur_2_[sub_sequence.sequence.i][sub_sequence.pos_2][sub_sequence.pos_1];
        return sequencing_scheme_.concatenate(
                sequence_data,
                sequence_data_2);
    }

    bool concatenate(
            SequenceData& sequence_data,
            bool sequence_empty,
            const SubSequence& sub_sequence,
            const GlobalCost& cutoff,
            const GlobalCost* gci1 = nullptr) const
    {
        if (sub_sequence.pos_1 > sub_sequence.pos_2)
            return true;
        if (sequence_empty) {
            const SequenceData& sequence_data_2 = (!sub_sequence.reverse)?
                sequence_datas_cur_2_[sub_sequence.sequence.i][sub_sequence.pos_1][sub_sequence.pos_2]:
                sequence_datas_cur_2_[sub_sequence.sequence.i][sub_sequence.pos_2][sub_sequence.pos_1];
            sequence_data = sequence_data_2;
            return true;
        }
        return concatenate(
                sequence_data,
                sub_sequence,
                cutoff,
                gci1,
                std::integral_constant<
                    bool,
                    HasConcatenateMethod<SequencingScheme, bool(SequenceData&, const SequenceData&)>::value>());
    }

    /*
     * GlobalCost global_cost_concatenate(
     *         SequenceData&,
     *         const SequenceData&)
     */

    template<typename, typename T>
    struct HasGlobalCostConcatenateMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasGlobalCostConcatenateMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().global_cost_concatenate(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    GlobalCost global_cost_concatenate(
            SequenceData& sequence_data,
            const SubSequence& sub_sequence,
            const GlobalCost& cutoff,
            const GlobalCost* gci1,
            std::false_type) const
    {
        bool ok = concatenate(sequence_data, false, sub_sequence, cutoff, gci1);
        if (!ok)
            return worst<GlobalCost>();
        return global_cost(sequence_data);
    }

    GlobalCost global_cost_concatenate(
            SequenceData& sequence_data,
            const SubSequence& sub_sequence,
            const GlobalCost&,
            const GlobalCost*,
            std::true_type) const
    {
        const SequenceData& sequence_data_2 = (!sub_sequence.reverse)?
            sequence_datas_cur_2_[sub_sequence.sequence.i][sub_sequence.pos_1][sub_sequence.pos_2]:
            sequence_datas_cur_2_[sub_sequence.sequence.i][sub_sequence.pos_2][sub_sequence.pos_1];
        return sequencing_scheme_.global_cost_concatenate(
                sequence_data,
                sequence_data_2);
    }

    GlobalCost global_cost_concatenate(
            SequenceData& sequence_data,
            bool sequence_empty,
            const SubSequence& sub_sequence,
            const GlobalCost& cutoff,
            const GlobalCost* gci1 = nullptr)
    {
        if (sub_sequence.pos_1 > sub_sequence.pos_2)
            return global_cost(sequence_data);
        if (sequence_empty) {
            const SequenceData& sequence_data_2 = (!sub_sequence.reverse)?
                sequence_datas_cur_2_[sub_sequence.sequence.i][sub_sequence.pos_1][sub_sequence.pos_2]:
                sequence_datas_cur_2_[sub_sequence.sequence.i][sub_sequence.pos_2][sub_sequence.pos_1];
            return global_cost(sequence_data_2);
        }
        return global_cost_concatenate(
                sequence_data,
                sub_sequence,
                cutoff,
                gci1,
                std::integral_constant<
                    bool,
                    HasGlobalCostConcatenateMethod<SequencingScheme, GlobalCost(SequenceData&, const SequenceData&)>::value>());
    }

    /**
     * Compute temporary structures.
     *
     * If 'i1 == -1', the temporary structures for all sequences are computed.
     * If 'i1 != -1', only the temporary structures of 'i1' and 'i2' are
     * updated.
     */
    inline void compute_temporary_structures(
            const Solution& solution,
            SequenceId i1 = -1,
            SequenceId i2 = -2)
    {
        SequencePos m = number_of_sequences();
        // Update modes_cur_.
        std::fill(modes_cur_.begin(), modes_cur_.end(), -1);
        for (SequenceId i = 0; i < m; ++i)
            for (const auto& se: solution.sequences[i].elements)
                modes_cur_[se.j] = se.mode;
        // Update sequence_datas_cur_.
        for (SequenceId i = 0; i < m; ++i) {
            if (i1 != -1 && i != i1 && i != i2)
                continue;
            const auto& elements = solution.sequences[i].elements;
            ElementPos seq_size = (ElementPos)elements.size();
            sequence_datas_cur_[i][0] = empty_sequence_data(i);
            for (ElementPos pos = 0; pos < seq_size; ++pos) {
                sequence_datas_cur_[i][pos + 1] = sequence_datas_cur_[i][pos];
                append(sequence_datas_cur_[i][pos + 1], elements[pos]);
            }
        }
        // Update sequence_datas_cur_2_.
        for (SequenceId i = 0; i < m; ++i) {
            if (i1 != -1 && i != i1 && i != i2)
                continue;
            const auto& elements = solution.sequences[i].elements;
            ElementPos seq_size = (ElementPos)elements.size();
            while ((ElementPos)sequence_datas_cur_2_[i].size()  < seq_size)
                sequence_datas_cur_2_[i].push_back({});
            for (ElementPos pos_1 = 0; pos_1 < seq_size; ++pos_1) {
                while ((ElementPos)sequence_datas_cur_2_[i][pos_1].size() < seq_size)
                    sequence_datas_cur_2_[i][pos_1].push_back({});
                sequence_datas_cur_2_[i][pos_1][pos_1] = empty_sequence_data(i);
                append(sequence_datas_cur_2_[i][pos_1][pos_1], elements[pos_1]);
            }
            for (ElementPos size = 2; size <= seq_size; ++size) {
                for (ElementPos pos_1 = 0; pos_1 + size - 1 < seq_size; ++pos_1) {
                    ElementPos pos_2 = pos_1 + size - 1;
                    sequence_datas_cur_2_[i][pos_1][pos_2] = sequence_datas_cur_2_[i][pos_1][pos_2 - 1];
                    append(sequence_datas_cur_2_[i][pos_1][pos_2], elements[pos_2]);
                    sequence_datas_cur_2_[i][pos_2][pos_1] = sequence_datas_cur_2_[i][pos_2][pos_1 + 1];
                    append(sequence_datas_cur_2_[i][pos_2][pos_1], elements[pos_1]);
                }
            }
        }
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

    /*
     * Explore neighborhoods.
     */

    inline void explore_shift(
            const Solution& solution,
            ElementPos block_size,
            bool reverse = false)
    {
        SequenceId m = number_of_sequences();
        Neighborhood& neighborhood = (!reverse)?
            neighborhoods_[int(Neighborhoods::Shift)][block_size][0]:
            neighborhoods_[int(Neighborhoods::ShiftReverse)][block_size][0];

        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();
            GlobalCost gci = global_cost(solution.sequences[i]);

            for (ElementPos block_pos = 0; block_pos <= seq_size - block_size; ++block_pos) {
                // Initialize sequence_data_cur.
                SequenceData sequence_data_cur = empty_sequence_data(i);

                // Loop through all new positions.
                ElementPos pos_min = std::max(
                        (ElementPos)0,
                        block_pos - parameters_.shift_maximum_distance);
                ElementPos pos_max = std::min(
                        seq_size - block_size,
                        block_pos + parameters_.shift_maximum_distance);
                for (ElementPos pos_new = pos_min; pos_new <= pos_max; ++pos_new) {
                    ElementPos p0 = (pos_new < block_pos)? pos_new: pos_new + block_size;
                    SequenceData sequence_data = sequence_data_cur;
                    bool ok = concatenate(sequence_data, (pos_new == 0), SubSequence(sequence, block_pos, block_pos + block_size - 1, reverse), gci);
                    if (ok)
                        ok = concatenate(sequence_data, false, SubSequence(sequence, p0, block_pos - 1), gci);
                    if (ok) {
                        GlobalCost gci_tmp = global_cost_concatenate(sequence_data, false, SubSequence(sequence, std::max(p0, block_pos + block_size), seq_size - 1), gci);
                        if (!(gci_tmp >= gci)) {
                            Move0 move;
                            move.type = (!reverse)?
                                Neighborhoods::Shift:
                                Neighborhoods::ShiftReverse;
                            move.k1 = block_size;
                            move.i1 = i;
                            move.pos_1 = block_pos;
                            move.pos_2 = pos_new;
                            move.global_cost = gci_tmp - gci;
                            neighborhood.improving_moves.push_back(move);
                        }
                    }

                    // Stop condition.
                    if (pos_new == seq_size - block_size)
                        break;

                    // Add next element to sequence_data_cur.
                    append(sequence_data_cur, sequence.elements[p0]);
                    // Check early termination.
                    if (bound(sequence_data_cur) >= gci)
                        break;
                }
            }
        }
    }

    inline void explore_swap(
            const Solution& solution,
            Counter block_size_1,
            Counter block_size_2)
    {
        SequenceId m = number_of_sequences();
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::Swap)][block_size_1][block_size_2];
        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();
            GlobalCost gci = global_cost(solution.sequences[i]);

            // Loop through all pairs.
            for (ElementPos pos = 0; pos < seq_size; ++pos) {

                // block 1 is at [pos, pos + block_size_1[.
                ElementPos pos_1 = pos;
                ElementPos pos_2_max = std::min(
                        seq_size - block_size_2,
                        pos + block_size_1 + parameters_.swap_maximum_distance);
                for (ElementPos pos_2 = pos + block_size_1; pos_2 <= pos_2_max; ++pos_2) {
                    bool ok = true;
                    SequenceData sequence_data = sequence_datas_cur_[i][pos];
                    ok = concatenate(sequence_data, (pos == 0), SubSequence(sequence, pos_2, pos_2 + block_size_2 - 1), gci);
                    if (!ok)
                        continue;
                    ok = concatenate(sequence_data, false, SubSequence(sequence, pos_1 + block_size_1, pos_2 - 1), gci);
                    if (!ok)
                        continue;
                    ok = concatenate(sequence_data, false, SubSequence(sequence, pos_1, pos_1 + block_size_1 - 1), gci);
                    if (!ok)
                        continue;
                    GlobalCost gci_tmp = global_cost_concatenate(sequence_data, false, SubSequence(sequence, pos_2 + block_size_2, seq_size - 1), gci);
                    if (!(gci_tmp >= gci)) {
                        Move0 move;
                        move.type = Neighborhoods::Swap;
                        move.k1 = block_size_1;
                        move.k2 = block_size_2;
                        move.i1 = i;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = gci_tmp - gci;
                        neighborhood.improving_moves.push_back(move);
                    }
                }

                if (block_size_1 == block_size_2)
                    continue;

                // block 2 is at [pos, pos + block_size_2[.
                ElementPos pos_2 = pos;
                ElementPos pos_1_max = std::min(
                        seq_size - block_size_1,
                        pos + block_size_2 + parameters_.swap_maximum_distance);
                for (ElementPos pos_1 = pos + block_size_2; pos_1 <= pos_1_max; ++pos_1) {
                    bool ok = true;
                    SequenceData sequence_data = sequence_datas_cur_[i][pos];
                    ok = concatenate(sequence_data, (pos == 0), SubSequence(sequence, pos_1, pos_1 + block_size_1 - 1), gci);
                    if (!ok)
                        continue;
                    ok = concatenate(sequence_data, false, SubSequence(sequence, pos_2 + block_size_2, pos_1 - 1), gci);
                    if (!ok)
                        continue;
                    ok = concatenate(sequence_data, false, SubSequence(sequence, pos_2, pos_2 + block_size_2 - 1), gci);
                    if (!ok)
                        continue;
                    GlobalCost gci_tmp = global_cost_concatenate(sequence_data, false, SubSequence(sequence, pos_1 + block_size_1, seq_size - 1), gci);
                    if (!(gci_tmp >= gci)) {
                        Move0 move;
                        move.type = Neighborhoods::Swap;
                        move.k1 = block_size_1;
                        move.k2 = block_size_2;
                        move.i1 = i;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = gci_tmp - gci;
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_reverse(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::Reverse)][0][0];
        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();
            GlobalCost gci = global_cost(solution.sequences[i]);

            // Loop through all pairs.
            for (ElementPos pos_1 = 0; pos_1 < seq_size; ++pos_1) {
                ElementPos pos_max = std::min(
                        seq_size,
                        pos_1 + parameters_.reverse_maximum_length);
                for (ElementPos pos_2 = pos_1 + 2; pos_2 < pos_max; ++pos_2) {
                    SequenceData sequence_data = sequence_datas_cur_[i][pos_1];
                    bool ok = concatenate(sequence_data, (pos_1 == 0), SubSequence(sequence, pos_1, pos_2, true), gci);
                    if (!ok)
                        continue;
                    GlobalCost gci_tmp = global_cost_concatenate(sequence_data, false, SubSequence(sequence, pos_2 + 1, seq_size - 1), gci);
                    if (!(gci_tmp >= gci)) {
                        Move0 move;
                        move.type = Neighborhoods::Reverse;
                        move.i1 = i;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = gci_tmp - gci;
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_add(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::Add)][0][0];

        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();
            GlobalCost gci = global_cost(solution.sequences[i]);

            // Loop through all new positions.
            for (ElementPos pos = 0; pos <= seq_size; ++pos) {

                for (ElementId j = 0; j < n; ++j) {
                    if (modes_cur_[j] != -1)
                        continue;

                    for (Mode mode = 0; mode < number_of_modes(j); ++mode) {
                        SequenceData sequence_data = sequence_datas_cur_[i][pos];
                        append(sequence_data, {j, mode});
                        GlobalCost gci_tmp = global_cost_concatenate(sequence_data, false, SubSequence(sequence, pos, seq_size - 1), gci);
                        if (!(gci_tmp >= gci)) {
                            Move0 move;
                            move.type = Neighborhoods::Add;
                            move.i1 = i;
                            move.j = j;
                            move.mode = mode;
                            move.pos_1 = pos;
                            move.global_cost = gci_tmp - gci;
                            neighborhood.improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
    }

    inline void explore_remove(
            const Solution& solution,
            const Move& perturbation)
    {
        SequenceId m = number_of_sequences();
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::Remove)][0][0];
        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();
            GlobalCost gci = global_cost(solution.sequences[i]);

            // Loop through all new positions.
            for (ElementPos pos = 0; pos < seq_size; ++pos) {

                ElementId j = sequence.elements[pos].j;
                if (perturbation.type == Perturbations::ForceAdd
                        && j == perturbation.force_add_j)
                    continue;

                SequenceData sequence_data = sequence_datas_cur_[i][pos];
                GlobalCost gci_tmp = global_cost_concatenate(sequence_data, (pos == 0), SubSequence(sequence, pos + 1, seq_size - 1), gci);
                if (!(gci_tmp >= gci)) {
                    Move0 move;
                    move.type = Neighborhoods::Remove;
                    move.i1 = i;
                    move.pos_1 = pos;
                    move.global_cost = gci_tmp - gci;
                    neighborhood.improving_moves.push_back(move);
                }
            }
        }
    }

    inline void explore_inter_shift(
            const Solution& solution,
            ElementPos block_size,
            bool reverse = false)
    {
        SequenceId m = number_of_sequences();
        Neighborhood& neighborhood = (!reverse)?
            neighborhoods_[int(Neighborhoods::InterShift)][block_size][0]:
            neighborhoods_[int(Neighborhoods::InterShiftReverse)][block_size][0];

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (SequenceId i2 = 0; i2 < m; ++i2) {
                if (!neighborhood.modified_sequences[i1]
                        && !neighborhood.modified_sequences[i2])
                    continue;
                if (i1 == i2)
                    continue;

                const auto& sequence_2 = solution.sequences[i2];
                SequencePos seq_2_size = sequence_2.elements.size();

                GlobalCost gci = global_cost_merge(
                        global_cost(solution.sequences[i1]),
                        global_cost(solution.sequences[i2]));

                for (ElementPos pos_1 = 0; pos_1 <= seq_1_size - block_size; ++pos_1) {
                    // pos_1 is the start position of the block of element to shift.

                    // Compute first sequence, i.e. the one from which the block iof
                    // elements is removed.
                    SequenceData sequence_data_1 = sequence_datas_cur_[i1][pos_1];
                    GlobalCost gci1_tmp = global_cost_concatenate(sequence_data_1, (pos_1 == 0), SubSequence(sequence_1, pos_1 + block_size, seq_1_size - 1), gci);
                    bool ok = !(gci1_tmp == worst<GlobalCost>());
                    if (!ok)
                        continue;

                    // Loop through all new positions.
                    ElementPos pos_min = std::max(
                            (ElementPos)0,
                            pos_1 - parameters_.shift_maximum_distance);
                    ElementPos pos_max = std::min(
                            seq_2_size,
                            pos_1 + parameters_.shift_maximum_distance);
                    for (ElementPos pos_2 = pos_min; pos_2 <= pos_max; ++pos_2) {
                        bool ok = true;
                        SequenceData sequence_data_2 = sequence_datas_cur_[i2][pos_2];
                        ok = concatenate(sequence_data_2, (pos_2 == 0), SubSequence(sequence_1, pos_1, pos_1 + block_size - 1, reverse), gci, &gci1_tmp);
                        if (!ok)
                            continue;
                        GlobalCost gci2_tmp = global_cost_concatenate(sequence_data_2, false, SubSequence(sequence_2, pos_2, seq_2_size - 1), gci, &gci1_tmp);
                        ok = !(gci2_tmp == worst<GlobalCost>());
                        if (!ok)
                            continue;
                        GlobalCost gci_tmp = global_cost_merge(gci1_tmp, gci2_tmp);
                        ok = !(gci2_tmp == worst<GlobalCost>());
                        if (!(gci_tmp >= gci)) {
                            Move0 move;
                            move.type = (!reverse)?
                                Neighborhoods::InterShift:
                                Neighborhoods::InterShiftReverse;
                            move.k1 = block_size;
                            move.i1 = i1;
                            move.i2 = i2;
                            move.pos_1 = pos_1;
                            move.pos_2 = pos_2;
                            move.global_cost = gci_tmp - gci;
                            neighborhood.improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
    }

    inline void explore_inter_two_opt(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::InterTwoOpt)][0][0];

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (SequenceId i2 = 0; i2 < m; ++i2) {
                if (!neighborhood.modified_sequences[i1]
                        && !neighborhood.modified_sequences[i2])
                    continue;
                if (i1 == i2)
                    continue;

                const auto& sequence_2 = solution.sequences[i2];
                SequencePos seq_2_size = sequence_2.elements.size();

                GlobalCost gci = global_cost_merge(
                        global_cost(solution.sequences[i1]),
                        global_cost(solution.sequences[i2]));

                for (ElementPos pos_1 = 0; pos_1 <= seq_1_size; ++pos_1) {

                    for (ElementPos pos_2 = 0; pos_2 <= seq_2_size; ++pos_2) {
                        bool ok = true;
                        SequenceData sequence_data_1 = sequence_datas_cur_[i1][pos_1];
                        SequenceData sequence_data_2 = sequence_datas_cur_[i2][pos_2];
                        GlobalCost bnd_2 = bound(sequence_data_2);
                        GlobalCost gci1_tmp = global_cost_concatenate(sequence_data_1, (pos_1 == 0), SubSequence(sequence_2, pos_2, seq_2_size - 1), gci, &bnd_2);
                        ok = !(gci1_tmp == worst<GlobalCost>());
                        if (!ok)
                            continue;
                        GlobalCost gci2_tmp = global_cost_concatenate(sequence_data_2, (pos_2 == 0), SubSequence(sequence_1, pos_1, seq_1_size - 1), gci, &gci1_tmp);
                        ok = !(gci2_tmp == worst<GlobalCost>());
                        if (!ok)
                            continue;
                        GlobalCost gci_tmp = global_cost_merge(gci1_tmp, gci2_tmp);
                        if (!(gci_tmp >= gci)) {
                            Move0 move;
                            move.type = Neighborhoods::InterTwoOpt;
                            move.i1 = i1;
                            move.i2 = i2;
                            move.pos_1 = pos_1;
                            move.pos_2 = pos_2;
                            move.global_cost = gci_tmp - gci;
                            neighborhood.improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
    }

    inline void explore_inter_swap(
            const Solution& solution,
            ElementPos block_size_1,
            ElementPos block_size_2)
    {
        SequenceId m = number_of_sequences();
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::InterSwap)][block_size_1][block_size_2];

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (SequenceId i2 = 0; i2 < m; ++i2) {
                if (!neighborhood.modified_sequences[i1]
                        && !neighborhood.modified_sequences[i2])
                    continue;
                if (i1 == i2)
                    continue;
                if (block_size_1 == block_size_2 && i2 > i1)
                    continue;

                const auto& sequence_2 = solution.sequences[i2];
                SequencePos seq_2_size = sequence_2.elements.size();

                GlobalCost gci = global_cost_merge(
                        global_cost(solution.sequences[i1]),
                        global_cost(solution.sequences[i2]));

                for (ElementPos pos_1 = 0; pos_1 <= seq_1_size - block_size_1; ++pos_1) {

                    for (ElementPos pos_2 = 0; pos_2 <= seq_2_size - block_size_2; ++pos_2) {
                        bool ok = true;
                        SequenceData sequence_data_1 = sequence_datas_cur_[i1][pos_1];
                        SequenceData sequence_data_2 = sequence_datas_cur_[i2][pos_2];
                        GlobalCost bnd_2 = bound(sequence_data_2);

                        ok = concatenate(sequence_data_1, (pos_1 == 0), SubSequence(sequence_2, pos_2, pos_2 + block_size_2 - 1), gci, &bnd_2);
                        if (!ok)
                            continue;

                        GlobalCost gci1_tmp = global_cost_concatenate(sequence_data_1, false, SubSequence(sequence_1, pos_1 + block_size_1, seq_1_size - 1), gci, &bnd_2);
                        ok = !(gci1_tmp == worst<GlobalCost>());
                        if (!ok)
                            continue;

                        ok = concatenate(sequence_data_2, (pos_2 == 0), SubSequence(sequence_1, pos_1, pos_1 + block_size_1 - 1), gci, &gci1_tmp);
                        if (!ok)
                            continue;

                        GlobalCost gci2_tmp = global_cost_concatenate(sequence_data_2, false, SubSequence(sequence_2, pos_2 + block_size_2, seq_2_size - 1), gci, &gci1_tmp);
                        ok = !(gci2_tmp == worst<GlobalCost>());
                        if (!ok)
                            continue;

                        GlobalCost gci_tmp = global_cost_merge(gci1_tmp, gci2_tmp);
                        if (!(gci_tmp >= gci)) {
                            Move0 move;
                            move.type = Neighborhoods::InterSwap;
                            move.k1 = block_size_1;
                            move.k2 = block_size_2;
                            move.i1 = i1;
                            move.i2 = i2;
                            move.pos_1 = pos_1;
                            move.pos_2 = pos_2;
                            move.global_cost = gci_tmp - gci;
                            neighborhood.improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
    }

    /**
     * 'added_element_pos' is the position in the original sequence.
     */
    inline GlobalCost compute_global_cost(
            const Sequence& sequence,
            ElementPos removed_element_pos,
            ElementPos added_element_pos,
            SequenceElement se,
            const GlobalCost& cutoff,
            const GlobalCost* gci1 = nullptr)
    {
        SequencePos seq_size = sequence.elements.size();
        //std::cout << "compute_global_cost"
        //    << " seq_size " << seq_size
        //    << " r " << removed_element_pos
        //    << " a " << added_element_pos
        //    << " j " << se.j
        //    << std::endl;
        if (removed_element_pos < added_element_pos) {
            SequenceData sequence_data = sequence_datas_cur_[sequence.i][removed_element_pos];
            bool ok = concatenate(sequence_data, (removed_element_pos == 0), SubSequence(sequence, removed_element_pos + 1, added_element_pos - 1), cutoff, gci1);
            if (!ok)
                return worst<GlobalCost>();
            append(sequence_data, se);
            return global_cost_concatenate(sequence_data, false, SubSequence(sequence, added_element_pos, seq_size - 1), cutoff, gci1);
        } else {
            SequenceData sequence_data = sequence_datas_cur_[sequence.i][added_element_pos];
            append(sequence_data, se);
            bool ok = concatenate(sequence_data, false, SubSequence(sequence, added_element_pos, removed_element_pos - 1), cutoff, gci1);
            if (!ok)
                return worst<GlobalCost>();
            return global_cost_concatenate(sequence_data, false, SubSequence(sequence, removed_element_pos + 1, seq_size - 1), cutoff, gci1);
        }
    }

    inline void explore_inter_swap_star(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::InterSwapStar)][0][0];

        for (SequenceId i1 = 0; i1 < m; ++i1) {
            const auto& sequence_1 = solution.sequences[i1];
            SequencePos seq_1_size = sequence_1.elements.size();

            for (SequenceId i2 = i1 + 1; i2 < m; ++i2) {
                if (!neighborhood.modified_sequences[i1]
                        && !neighborhood.modified_sequences[i2])
                    continue;

                const auto& sequence_2 = solution.sequences[i2];
                SequencePos seq_2_size = sequence_2.elements.size();

                // Calcul the 3 best position of the elements from the first
                // sequence into the second sequence.
                std::vector<std::vector<ElementPos>> best_positions_1(seq_1_size);
                for (ElementPos pos_1 = 0; pos_1 <= seq_1_size - 1; ++pos_1) {
                    ElementPos pos_best_1 = -1;
                    ElementPos pos_best_2 = -1;
                    ElementPos pos_best_3 = -1;
                    GlobalCost gc_best_1 = worst<GlobalCost>();
                    GlobalCost gc_best_2 = worst<GlobalCost>();
                    GlobalCost gc_best_3 = worst<GlobalCost>();
                    for (ElementPos pos_2 = 0; pos_2 <= seq_2_size; ++pos_2) {
                        SequenceData sequence_data = sequence_datas_cur_[i2][pos_2];
                        append(sequence_data, sequence_1.elements[pos_1]);
                        GlobalCost gci2_tmp = global_cost_concatenate(sequence_data, false, SubSequence(sequence_2, pos_2, seq_2_size - 1), worst<GlobalCost>());
                        if (pos_best_1 == -1 || !(gci2_tmp >= gc_best_1)) {
                            pos_best_3 = pos_best_2;
                            pos_best_2 = pos_best_1;
                            pos_best_1 = pos_2;
                            gc_best_3 = gc_best_2;
                            gc_best_2 = gc_best_1;
                            gc_best_1 = gci2_tmp;
                        } else if (pos_best_2 == -1 || !(gci2_tmp >= gc_best_2)) {
                            pos_best_3 = pos_best_2;
                            pos_best_2 = pos_2;
                            gc_best_3 = gc_best_2;
                            gc_best_2 = gci2_tmp;
                        } else if (pos_best_3 == -1 || !(gci2_tmp >= gc_best_3)) {
                            pos_best_3 = pos_2;
                            gc_best_3 = gci2_tmp;
                        }
                    }
                    //std::cout << "pos_1 " << pos_1
                    //    << " pos_best_1 " << pos_best_1
                    //    << " gc_best_1 " << to_string(gc_best_1)
                    //    << " pos_best_2 " << pos_best_2
                    //    << " gc_best_2 " << to_string(gc_best_2)
                    //    << " pos_best_3 " << pos_best_3
                    //    << " gc_best_3 " << to_string(gc_best_3)
                    //    << std::endl;
                    best_positions_1[pos_1] = {pos_best_1, pos_best_2, pos_best_3};
                }

                // Calcul the 3 best position of the elements from the second
                // sequence into the first sequence.
                std::vector<std::vector<ElementPos>> best_positions_2(seq_2_size);
                for (ElementPos pos_2 = 0; pos_2 <= seq_2_size - 1; ++pos_2) {
                    ElementPos pos_best_1 = -1;
                    ElementPos pos_best_2 = -1;
                    ElementPos pos_best_3 = -1;
                    GlobalCost gc_best_1 = worst<GlobalCost>();
                    GlobalCost gc_best_2 = worst<GlobalCost>();
                    GlobalCost gc_best_3 = worst<GlobalCost>();
                    for (ElementPos pos_1 = 0; pos_1 <= seq_1_size; ++pos_1) {
                        SequenceData sequence_data = sequence_datas_cur_[i1][pos_1];
                        append(sequence_data, sequence_2.elements[pos_2]);
                        GlobalCost gci1_tmp = global_cost_concatenate(sequence_data, false, SubSequence(sequence_1, pos_1, seq_1_size - 1), worst<GlobalCost>());
                        if (pos_best_1 == -1 || !(gci1_tmp >= gc_best_1)) {
                            pos_best_3 = pos_best_2;
                            pos_best_2 = pos_best_1;
                            pos_best_1 = pos_1;
                            gc_best_3 = gc_best_2;
                            gc_best_2 = gc_best_1;
                            gc_best_1 = gci1_tmp;
                        } else if (pos_best_2 == -1 || !(gci1_tmp >= gc_best_2)) {
                            pos_best_3 = pos_best_2;
                            pos_best_2 = pos_1;
                            gc_best_3 = gc_best_2;
                            gc_best_2 = gci1_tmp;
                        } else if (pos_best_3 == -1 || !(gci1_tmp >= gc_best_3)) {
                            pos_best_3 = pos_1;
                            gc_best_3 = gci1_tmp;
                        }
                    }
                    //std::cout << "pos_2 " << pos_2
                    //    << " pos_best_1 " << pos_best_1
                    //    << " gc_best_1 " << to_string(gc_best_1)
                    //    << " pos_best_2 " << pos_best_2
                    //    << " gc_best_2 " << to_string(gc_best_2)
                    //    << " pos_best_3 " << pos_best_3
                    //    << " gc_best_3 " << to_string(gc_best_3)
                    //    << std::endl;
                    best_positions_2[pos_2] = {pos_best_1, pos_best_2, pos_best_3};
                }

                GlobalCost gci = global_cost_merge(
                        global_cost(solution.sequences[i1]),
                        global_cost(solution.sequences[i2]));

                for (ElementPos pos_1 = 0; pos_1 <= seq_1_size - 1; ++pos_1) {

                    for (ElementPos pos_2 = 0; pos_2 <= seq_2_size - 1; ++pos_2) {

                        //std::cout
                        //    << "i1 " << i1
                        //    << " i2 " << i2
                        //    << " seq_1_size " << seq_1_size
                        //    << " seq_2_size " << seq_2_size
                        //    << " pos_1 " << pos_1
                        //    << " pos_2 " << pos_2
                        //    << std::endl;

                        // Find pos_1_new.
                        ElementPos pos_1_new = -1;
                        GlobalCost gci2_tmp = worst<GlobalCost>();
                        for (int a = 0; a < 4; ++a) {
                            ElementPos p = (a < 3)? best_positions_1[pos_1][a]: pos_2;
                            if (p != -1) {
                                GlobalCost gci2 = compute_global_cost(
                                        sequence_2,
                                        pos_2,
                                        p,
                                        sequence_1.elements[pos_1],
                                        gci);
                                //std::cout << "p " << p
                                //    << " gc " << to_string(gci2)
                                //    << std::endl;
                                if (!(gci2 >= gci2_tmp)) {
                                    gci2_tmp = gci2;
                                    pos_1_new = p;
                                }
                            }
                        }
                        if (pos_1_new == -1)
                            continue;

                        // Find pos_2_new.
                        ElementPos pos_2_new = -1;
                        GlobalCost gci1_tmp = worst<GlobalCost>();
                        for (int a = 0; a < 4; ++a) {
                            ElementPos p = (a < 3)? best_positions_2[pos_2][a]: pos_1;
                            if (p != -1) {
                                GlobalCost gci1 = compute_global_cost(
                                        sequence_1,
                                        pos_1,
                                        p,
                                        sequence_2.elements[pos_2],
                                        gci,
                                        &gci2_tmp);
                                //std::cout << "p " << p
                                //    << " gc " << to_string(gci1)
                                //    << std::endl;
                                if (!(gci1 >= gci1_tmp)) {
                                    gci1_tmp = gci1;
                                    pos_2_new = p;
                                }
                            }
                        }
                        if (pos_2_new == -1)
                            continue;

                        //std::cout << "pos_1 " << pos_2
                        //    << " pos_2 " << pos_2
                        //    << " pos_1_new " << pos_1_new
                        //    << " pos_2_new " << pos_2_new
                        //    << std::endl;

                        GlobalCost gci_tmp = global_cost_merge(gci1_tmp, gci2_tmp);
                        if (!(gci_tmp >= gci)) {
                            Move0 move;
                            move.type = Neighborhoods::InterSwapStar;
                            move.i1 = i1;
                            move.i2 = i2;
                            move.pos_1 = pos_1;
                            move.pos_2 = pos_1_new;
                            move.pos_3 = pos_2;
                            move.pos_4 = pos_2_new;
                            move.global_cost = gci_tmp - gci;
                            neighborhood.improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
    }

    inline void explore_shift_change_mode(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::ShiftChangeMode)][0][0];
        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();
            GlobalCost gci = global_cost(solution.sequences[i]);

            for (ElementPos block_pos = 0; block_pos <= seq_size - 1; ++block_pos) {
                ElementId j = sequence.elements[block_pos].j;
                // Initialize sequence_data_cur.
                SequenceData sequence_data_cur = empty_sequence_data(i);

                // Loop through all new positions.
                ElementPos pos_min = std::max(
                        (ElementPos)0,
                        block_pos - parameters_.shift_maximum_distance);
                ElementPos pos_max = std::min(
                        seq_size - 1,
                        block_pos + parameters_.shift_maximum_distance);
                for (ElementPos pos_new = pos_min; pos_new <= pos_max; ++pos_new) {
                    ElementPos p0 = (pos_new < block_pos)? pos_new: pos_new + 1;

                    for (Mode mode = 0; mode < number_of_modes(j); ++mode) {
                        SequenceData sequence_data = sequence_data_cur;
                        append(sequence_data, {j, mode});
                        bool ok = concatenate(sequence_data, false, SubSequence(sequence, p0, block_pos - 1), gci);
                        if (!ok)
                            continue;
                        GlobalCost gci_tmp = global_cost_concatenate(sequence_data, false, SubSequence(sequence, std::max(p0, block_pos + 1), seq_size - 1), gci);
                        if (!(gci_tmp >= gci)) {
                            Move0 move;
                            move.type = Neighborhoods::ShiftChangeMode;
                            move.i1 = i;
                            move.pos_1 = block_pos;
                            move.pos_2 = pos_new;
                            move.mode = mode;
                            move.global_cost = gci_tmp - gci;
                            neighborhood.improving_moves.push_back(move);
                        }
                    }

                    // Stop condition.
                    if (pos_new == seq_size - 1)
                        break;

                    // Add next element to sequence_data_cur_.
                    append(sequence_data_cur, sequence.elements[p0]);
                    // Check early termination.
                    if (bound(sequence_data_cur) >= gci)
                        break;
                }
            }
        }
    }

    inline void explore_mode_swap(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::ModeSwap)][0][0];
        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();
            GlobalCost gci = global_cost(solution.sequences[i]);

            // Loop through all pairs.
            Counter pos_max = seq_size - 1;
            for (ElementPos pos_1 = 0; pos_1 <= pos_max; ++pos_1) {
                ElementId j1 = sequence.elements[pos_1].j;
                ElementId mode_1 = sequence.elements[pos_1].mode;

                ElementPos pos_2_max = std::min(
                        pos_max,
                        pos_1 + parameters_.swap_maximum_distance);
                for (ElementPos pos_2 = pos_1 + 1; pos_2 < pos_2_max; ++pos_2) {
                    ElementId j2 = sequence.elements[pos_2].j;
                    ElementId mode_2 = sequence.elements[pos_2].mode;

                    if (mode_1 == mode_2)
                        continue;

                    bool ok = true;
                    SequenceData sequence_data = sequence_datas_cur_[i][pos_1];
                    append(sequence_data, {j1, mode_2});
                    ok = concatenate(sequence_data, false, SubSequence(sequence, pos_1 + 1, pos_2 - 1), gci);
                    if (!ok)
                        continue;
                    append(sequence_data, {j2, mode_1});
                    GlobalCost gci_tmp = global_cost_concatenate(sequence_data, false, SubSequence(sequence, pos_2 + 1, seq_size - 1), gci);
                    if (!(gci_tmp >= gci)) {
                        Move0 move;
                        move.type = Neighborhoods::ModeSwap;
                        move.i1 = i;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = gci_tmp - gci;
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline void explore_swap_with_modes(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        Neighborhood& neighborhood = neighborhoods_[int(Neighborhoods::SwapWithModes)][0][0];
        for (SequenceId i = 0; i < m; ++i) {
            if (!neighborhood.modified_sequences[i])
                continue;
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();
            GlobalCost gci = global_cost(solution.sequences[i]);

            // Loop through all pairs.
            Counter pos_max = seq_size - 1;
            for (ElementPos pos_1 = 0; pos_1 <= pos_max; ++pos_1) {
                ElementId j1 = sequence.elements[pos_1].j;
                ElementId mode_1 = sequence.elements[pos_1].mode;

                ElementPos pos_2_max = std::min(
                        pos_max,
                        pos_1 + parameters_.swap_maximum_distance);
                for (ElementPos pos_2 = pos_1 + 1; pos_2 < pos_2_max; ++pos_2) {
                    ElementId j2 = sequence.elements[pos_2].j;
                    ElementId mode_2 = sequence.elements[pos_2].mode;

                    if (mode_1 == mode_2)
                        continue;

                    bool ok = true;
                    SequenceData sequence_data = sequence_datas_cur_[i][pos_1];
                    append(sequence_data, {j2, mode_1});
                    ok = concatenate(sequence_data, false, SubSequence(sequence, pos_1 + 1, pos_2 - 1), gci);
                    if (!ok)
                        continue;
                    append(sequence_data, {j1, mode_2});
                    GlobalCost gci_tmp = global_cost_concatenate(sequence_data, false, SubSequence(sequence, pos_2 + 1, seq_size - 1), gci);
                    if (!(gci_tmp >= gci)) {
                        Move0 move;
                        move.type = Neighborhoods::SwapWithModes;
                        move.i1 = i;
                        move.pos_1 = pos_1;
                        move.pos_2 = pos_2;
                        move.global_cost = gci_tmp - gci;
                        neighborhood.improving_moves.push_back(move);
                    }
                }
            }
        }
    }

    inline std::vector<Move0> explore_add_2(
            const Solution& solution)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        std::vector<Move0> improving_moves;
        for (SequenceId i = 0; i < m; ++i) {
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();
            GlobalCost gci = global_cost(solution.sequences[i]);

            // Loop through all new positions.
            for (ElementPos pos = 0; pos <= seq_size; ++pos) {

                for (ElementId j = 0; j < n; ++j) {
                    if (modes_cur_[j] != -1)
                        continue;

                    for (Mode mode = 0; mode < number_of_modes(j); ++mode) {
                        SequenceData sequence_data = sequence_datas_cur_[i][pos];
                        append(sequence_data, {j, mode});
                        GlobalCost gci_tmp = global_cost_concatenate(sequence_data, false, SubSequence(sequence, pos, seq_size - 1), gci);
                        if (!(gci_tmp >= gci)) {
                            Move0 move;
                            move.type = Neighborhoods::Add;
                            move.i1 = i;
                            move.j = j;
                            move.mode = mode;
                            move.pos_1 = pos;
                            move.global_cost = gci_tmp - gci;
                            improving_moves.push_back(move);
                        }
                    }
                }
            }
        }
        return improving_moves;
    }

    inline std::vector<Move0> explore_add(
            const Solution& solution,
            ElementId j,
            SequenceId i_old = -1,
            ElementId j_prec_old = -2,
            Mode mode_old = -1)
    {
        std::vector<Move0> improving_moves;
        SequenceId m = number_of_sequences();
        for (SequenceId i = 0; i < m; ++i) {
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();
            GlobalCost gci = global_cost(solution.sequences[i]);

            // Loop through all new positions.
            for (ElementPos pos = 0; pos <= seq_size; ++pos) {

                for (Mode mode = 0; mode < number_of_modes(j); ++mode) {

                    if (j_prec_old == -1
                            && i == i_old
                            && (seq_size > 0 || m > 1)
                            && pos == 0
                            && mode == mode_old)
                        continue;
                    if (pos > 0
                            && j_prec_old == sequence.elements[pos - 1].j
                            && mode == mode_old)
                        continue;

                    SequenceData sequence_data = sequence_datas_cur_[i][pos];
                    append(sequence_data, {j, mode});
                    GlobalCost gci_tmp = global_cost_concatenate(sequence_data, false, SubSequence(sequence, pos, seq_size - 1), gci);
                    if (!(gci_tmp >= gci)) {
                        Move0 move;
                        move.type = Neighborhoods::Add;
                        move.i1 = i;
                        move.j = j;
                        move.mode = mode;
                        move.pos_1 = pos;
                        move.global_cost = gci_tmp - gci;
                        improving_moves.push_back(move);
                    }
                }
            }
        }
        return improving_moves;
    }

    inline std::vector<Move0> explore_add_end(
            const Solution& solution)
    {
        std::vector<Move0> improving_moves;
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();
        for (SequenceId i = 0; i < m; ++i) {
            const auto& sequence = solution.sequences[i];
            SequencePos seq_size = sequence.elements.size();
            GlobalCost gci = global_cost(solution.sequences[i]);
            ElementPos pos = seq_size;

            // Loop through all new positions.
            for (ElementId j = 0; j < n; ++j) {
                if (modes_cur_[j] != -1)
                    continue;

                for (Mode mode = 0; mode < number_of_modes(j); ++mode) {
                    SequenceData sequence_data = sequence_datas_cur_[i][pos];
                    append(sequence_data, {j, mode});
                    GlobalCost gci_tmp = global_cost_concatenate(sequence_data, false, SubSequence(sequence, pos, seq_size - 1), gci);
                    if (!(gci_tmp >= gci)) {
                        Move0 move;
                        move.type = Neighborhoods::Add;
                        move.i1 = i;
                        move.j = j;
                        move.mode = mode;
                        move.pos_1 = pos;
                        move.global_cost = gci_tmp - gci;
                        improving_moves.push_back(move);
                    }
                }
            }
        }
        return improving_moves;
    }

    inline void apply_move(
            Solution& solution,
            const Move0& move)
    {
        SequenceId m = number_of_sequences();
        ElementPos n = sequencing_scheme_.number_of_elements();

        ElementPos n_old = 0;
        for (SequenceId i = 0; i < m; ++i)
            n_old += solution.sequences[i].elements.size();
        if (n_old > n) {
            throw std::logic_error(
                    "Too many element before applying move "
                    + neighborhood2string(move.type, move.k1, move.k2)
                    + ": " + std::to_string(n_old)
                    + " / " + std::to_string(n) + ".");
        }

        switch (move.type) {
        case Neighborhoods::Shift: {
            ElementPos block_size = move.k1;
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            if (move.pos_1 > move.pos_2) {
                for (ElementPos p = 0; p < move.pos_2; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1; p < move.pos_1 + block_size; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            } else {
                for (ElementPos p = 0; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size; p < move.pos_2 + block_size; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1; p < move.pos_1 + block_size; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2 + block_size; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            }
            break;
        } case Neighborhoods::Swap: {
            ElementPos block_size_1 = move.k1;
            ElementPos block_size_2 = move.k2;
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            if (move.pos_1 < move.pos_2) {
                for (ElementPos p = 0; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2; p < move.pos_2 + block_size_2; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size_1; p < move.pos_2; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1; p < move.pos_1 + block_size_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2 + block_size_2; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            } else {
                for (ElementPos p = 0; p < move.pos_2; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1; p < move.pos_1 + block_size_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2 + block_size_2; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2; p < move.pos_2 + block_size_2; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size_1; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            }
            break;
        } case Neighborhoods::Reverse: {
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp, elements[p]);
            for (ElementPos p = move.pos_2; p >= move.pos_1; --p)
                append(sequence_tmp, elements[p]);
            for (ElementPos p = move.pos_2 + 1; p < seq_size; ++p)
                append(sequence_tmp, elements[p]);
            compute_global_cost(solution_tmp_);
            solution = solution_tmp_;
            break;
        } case Neighborhoods::ShiftReverse: {
            ElementPos block_size = move.k1;
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            if (move.pos_1 > move.pos_2) {
                for (ElementPos p = 0; p < move.pos_2; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size - 1; p >= move.pos_1; --p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            } else {
                for (ElementPos p = 0; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size; p < move.pos_2 + block_size; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + block_size - 1; p >= move.pos_1; --p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_2 + block_size; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            }
            break;
        } case Neighborhoods::Add: {
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp, elements[p]);
            append(sequence_tmp, {move.j, move.mode});
            for (ElementPos p = move.pos_1; p < seq_size; ++p)
                append(sequence_tmp, elements[p]);
            break;
        } case Neighborhoods::Remove: {
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            for (ElementPos p = 0; p < seq_size; ++p) {
                if (p == move.pos_1)
                    continue;
                append(sequence_tmp, elements[p]);
            }
            break;
        } case Neighborhoods::InterTwoOpt: {
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1 && i != move.i2)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements_1 = solution.sequences[move.i1].elements;
            const auto& elements_2 = solution.sequences[move.i2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.i1];
            sequence_tmp_1 = empty_sequence(move.i1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_2; p < seq_2_size; ++p)
                append(sequence_tmp_1, elements_2[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.i2];
            sequence_tmp_2 = empty_sequence(move.i2);
            for (ElementPos p = 0; p < move.pos_2; ++p)
                append(sequence_tmp_2, elements_2[p]);
            for (ElementPos p = move.pos_1; p < seq_1_size; ++p)
                append(sequence_tmp_2, elements_1[p]);
            break;
        } case Neighborhoods::InterShift: {
            ElementPos block_size = move.k1;
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1 && i != move.i2)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements_1 = solution.sequences[move.i1].elements;
            const auto& elements_2 = solution.sequences[move.i2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.i1];
            sequence_tmp_1 = empty_sequence(move.i1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_1 + block_size; p < seq_1_size; ++p)
                append(sequence_tmp_1, elements_1[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.i2];
            sequence_tmp_2 = empty_sequence(move.i2);
            for (ElementPos p = 0; p < move.pos_2; ++p)
                append(sequence_tmp_2, elements_2[p]);
            for (ElementPos p = move.pos_1; p < move.pos_1 + block_size; ++p)
                append(sequence_tmp_2, elements_1[p]);
            for (ElementPos p = move.pos_2; p < seq_2_size; ++p)
                append(sequence_tmp_2, elements_2[p]);
            break;
        } case Neighborhoods::InterSwap: {
            ElementPos block_size_1 = move.k1;
            ElementPos block_size_2 = move.k2;
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1 && i != move.i2)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements_1 = solution.sequences[move.i1].elements;
            const auto& elements_2 = solution.sequences[move.i2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.i1];
            sequence_tmp_1 = empty_sequence(move.i1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_2; p < move.pos_2 + block_size_2; ++p)
                append(sequence_tmp_1, elements_2[p]);
            for (ElementPos p = move.pos_1 + block_size_1; p < seq_1_size; ++p)
                append(sequence_tmp_1, elements_1[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.i2];
            sequence_tmp_2 = empty_sequence(move.i2);
            for (ElementPos p = 0; p < move.pos_2; ++p)
                append(sequence_tmp_2, elements_2[p]);
            for (ElementPos p = move.pos_1; p < move.pos_1 + block_size_1; ++p)
                append(sequence_tmp_2, elements_1[p]);
            for (ElementPos p = move.pos_2 + block_size_2; p < seq_2_size; ++p)
                append(sequence_tmp_2, elements_2[p]);
            break;
        } case Neighborhoods::InterShiftReverse: {
            ElementPos block_size = move.k1;
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1 && i != move.i2)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements_1 = solution.sequences[move.i1].elements;
            const auto& elements_2 = solution.sequences[move.i2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            // Sequence 1.
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.i1];
            sequence_tmp_1 = empty_sequence(move.i1);
            for (ElementPos p = 0; p < move.pos_1; ++p)
                append(sequence_tmp_1, elements_1[p]);
            for (ElementPos p = move.pos_1 + block_size; p < seq_1_size; ++p)
                append(sequence_tmp_1, elements_1[p]);
            // Sequence 2.
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.i2];
            sequence_tmp_2 = empty_sequence(move.i2);
            for (ElementPos p = 0; p < move.pos_2; ++p)
                append(sequence_tmp_2, elements_2[p]);
            for (ElementPos p = move.pos_1 + block_size - 1; p >= move.pos_1; --p)
                append(sequence_tmp_2, elements_1[p]);
            for (ElementPos p = move.pos_2; p < seq_2_size; ++p)
                append(sequence_tmp_2, elements_2[p]);
            break;
        } case Neighborhoods::ShiftChangeMode: {
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            if (move.pos_1 > move.pos_2) {
                for (ElementPos p = 0; p < move.pos_2; ++p)
                    append(sequence_tmp, elements[p]);
                append(sequence_tmp, {elements[move.pos_1].j, move.mode});
                for (ElementPos p = move.pos_2; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + 1; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            } else {
                for (ElementPos p = 0; p < move.pos_1; ++p)
                    append(sequence_tmp, elements[p]);
                for (ElementPos p = move.pos_1 + 1; p < move.pos_2 + 1; ++p)
                    append(sequence_tmp, elements[p]);
                append(sequence_tmp, {elements[move.pos_1].j, move.mode});
                for (ElementPos p = move.pos_2 + 1; p < seq_size; ++p)
                    append(sequence_tmp, elements[p]);
            }
            break;
        } case Neighborhoods::ModeSwap: {
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            for (ElementPos p = 0; p < seq_size; ++p) {
                SequenceElement se = elements[p];
                if (p == move.pos_1) {
                    se.mode = elements[move.pos_2].mode;
                } else if (p == move.pos_2) {
                    se.mode = elements[move.pos_1].mode;
                }
                append(sequence_tmp, se);
            }
            break;
        } case Neighborhoods::SwapWithModes: {
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements = solution.sequences[move.i1].elements;
            SequencePos seq_size = elements.size();
            Sequence& sequence_tmp = solution_tmp_.sequences[move.i1];
            sequence_tmp = empty_sequence(move.i1);
            for (ElementPos p = 0; p < seq_size; ++p) {
                SequenceElement se = elements[p];
                if (p == move.pos_1) {
                    se.j = elements[move.pos_2].j;
                } else if (p == move.pos_2) {
                    se.j = elements[move.pos_1].j;
                }
                append(sequence_tmp, se);
            }
            break;
        } case Neighborhoods::InterSwapStar: {
            //std::cout << "i1 " << move.i1
            //    << " i2 " << move.i2
            //    << " pos_1_old " << move.pos_1
            //    << " pos_2_old " << move.pos_3
            //    << " pos_1_new " << move.pos_2
            //    << " pos_2_new " << move.pos_4
            //    << " gc " << to_string(move.global_cost)
            //    << std::endl;
            for (SequenceId i = 0; i < m; ++i)
                if (i != move.i1 && i != move.i2)
                    solution_tmp_.sequences[i] = solution.sequences[i];
            const auto& elements_1 = solution.sequences[move.i1].elements;
            const auto& elements_2 = solution.sequences[move.i2].elements;
            SequencePos seq_1_size = elements_1.size();
            SequencePos seq_2_size = elements_2.size();
            Sequence& sequence_tmp_1 = solution_tmp_.sequences[move.i1];
            Sequence& sequence_tmp_2 = solution_tmp_.sequences[move.i2];
            sequence_tmp_1 = empty_sequence(move.i1);
            ElementPos pos_1_old = move.pos_1;
            ElementPos pos_1_new = move.pos_2;
            ElementPos pos_2_old = move.pos_3;
            ElementPos pos_2_new = move.pos_4;
            if (pos_1_old < pos_2_new) {
                for (ElementPos p = 0; p < pos_1_old; ++p)
                    append(sequence_tmp_1, elements_1[p]);
                for (ElementPos p = pos_1_old + 1; p < pos_2_new; ++p)
                    append(sequence_tmp_1, elements_1[p]);
                append(sequence_tmp_1, elements_2[pos_2_old]);
                for (ElementPos p = pos_2_new; p < seq_1_size; ++p)
                    append(sequence_tmp_1, elements_1[p]);
            } else {
                for (ElementPos p = 0; p < pos_2_new; ++p)
                    append(sequence_tmp_1, elements_1[p]);
                append(sequence_tmp_1, elements_2[pos_2_old]);
                for (ElementPos p = pos_2_new; p < pos_1_old; ++p)
                    append(sequence_tmp_1, elements_1[p]);
                for (ElementPos p = pos_1_old + 1; p < seq_1_size; ++p)
                    append(sequence_tmp_1, elements_1[p]);
            }
            sequence_tmp_2 = empty_sequence(move.i2);
            if (pos_2_old < pos_1_new) {
                for (ElementPos p = 0; p < pos_2_old; ++p)
                    append(sequence_tmp_2, elements_2[p]);
                for (ElementPos p = pos_2_old + 1; p < pos_1_new; ++p)
                    append(sequence_tmp_2, elements_2[p]);
                append(sequence_tmp_2, elements_1[pos_1_old]);
                for (ElementPos p = pos_1_new; p < seq_2_size; ++p)
                    append(sequence_tmp_2, elements_2[p]);
            } else {
                for (ElementPos p = 0; p < pos_1_new; ++p)
                    append(sequence_tmp_2, elements_2[p]);
                append(sequence_tmp_2, elements_1[pos_1_old]);
                for (ElementPos p = pos_1_new; p < pos_2_old; ++p)
                    append(sequence_tmp_2, elements_2[p]);
                for (ElementPos p = pos_2_old + 1; p < seq_2_size; ++p)
                    append(sequence_tmp_2, elements_2[p]);
            }
            break;
        }
        }

        ElementPos n_new = 0;
        for (SequenceId i = 0; i < m; ++i)
            n_new += solution_tmp_.sequences[i].elements.size();
        if (n_new > n) {
            throw std::logic_error(
                    "Too many elements after applying move "
                    + neighborhood2string(move.type, move.k1, move.k2)
                    + ": " + std::to_string(n_new)
                    + " / " + std::to_string(n) + ".");
        }
        if (move.type != Neighborhoods::Add
                && move.type != Neighborhoods::Remove
                && n_new != n_old) {
            throw std::logic_error(
                    "Wrong number of elements after applying move "
                    + neighborhood2string(move.type, move.k1, move.k2)
                    + "; old: " + std::to_string(n_old)
                    + "; new " + std::to_string(n_new) + ".");
        }
        if (move.type == Neighborhoods::Add && n_new != n_old + 1) {
            throw std::logic_error(
                    "Wrong number of elements after applying move "
                    + neighborhood2string(move.type, move.k1, move.k2)
                    + "; old: " + std::to_string(n_old)
                    + "; new " + std::to_string(n_new) + ".");
        }
        if (move.type == Neighborhoods::Remove && n_new != n_old - 1) {
            throw std::logic_error(
                    "Wrong number of elements after applying move "
                    + neighborhood2string(move.type, move.k1, move.k2)
                    + "; old: " + std::to_string(n_old)
                    + "; new " + std::to_string(n_new) + ".");
        }

        solution_tmp_.modified_sequences = solution.modified_sequences;
        compute_global_cost(solution_tmp_);
        solution = solution_tmp_;
    }

    /*
     * Private attributes.
     */

    /** Sequencing scheme. */
    SequencingScheme& sequencing_scheme_;
    /** Parameters. */
    Parameters parameters_;
    /**
     * Maximum number of modes.
     *
     * Comupted in the constructor.
     */
    Mode maximum_number_of_modes_ = 1;
    /** Structure storing neighborhood related information. */
    std::vector<std::vector<std::vector<Neighborhood>>> neighborhoods_
        = std::vector<std::vector<std::vector<Neighborhood>>>(15);
    Counter number_of_crossover_calls = 0;
    double crossover_time = 0.0;
    Counter number_of_local_search_calls = 0;
    double local_search_time = 0.0;
    Counter number_of_initial_solution_calls = 0;
    double initial_solution_time = 0.0;

    /*
     * Temporary structures.
     */

    std::vector<std::vector<SequenceData>> sequence_datas_cur_;
    std::vector<std::vector<std::vector<SequenceData>>> sequence_datas_cur_2_;
    std::vector<Mode> modes_cur_;
    Solution solution_cur_;
    Solution solution_tmp_;

};

}

}

