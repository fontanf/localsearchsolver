#pragma once

#include "optimizationtools/utils/output.hpp"

#include <cstdint>
#include <set>
#include <iomanip>
#include <sstream>

namespace localsearchsolver
{

using Seed = int64_t;
using Counter = int64_t;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Tuple operations ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename T, T...>
struct integer_sequence { };

template<std::size_t N, std::size_t... I>
struct gen_indices : gen_indices<(N - 1), (N - 1), I...> { };
template<std::size_t... I>
struct gen_indices<0, I...> : integer_sequence<std::size_t, I...> { };

// to_string(t) //

template<typename H>
std::stringstream& to_string_impl(std::stringstream& ss, H&& h)
{
    ss << std::forward<H>(h);
    return ss;
}

template<typename H, typename... T>
std::stringstream& to_string_impl(std::stringstream& ss, H&& h, T&&... t)
{
    to_string_impl(ss, std::forward<T>(t)...);
    ss << ", " << std::forward<H>(h);
    return ss;
}

template<typename... T, std::size_t... I>
std::stringstream to_string(
        const std::tuple<T...>& tup,
        integer_sequence<std::size_t, I...>)
{
    std::stringstream ss;
    int ctx[] = { (to_string_impl(ss, std::get<I>(tup)...), 0), 0 };
    (void)ctx;
    return ss;
}

template<typename... T>
std::string to_string(const std::tuple<T...>& tup)
{
    return to_string(tup, gen_indices<sizeof...(T)>{}).str();
}

// dominates(t1, t2) //

template<typename... T, std::size_t... I>
bool dominates(
        const std::tuple<T...>& t1,
        const std::tuple<T...>& t2,
        integer_sequence<std::size_t, I...>)
{
   using unused = bool[];

   bool b { true };

   (void)unused { b, (b = b && std::get<I>(t1) <= std::get<I>(t2))... };

   return b;
}

template <typename ... T>
bool dominates(
        const std::tuple<T...>& t1,
        const std::tuple<T...>& t2)
{
    return dominates(t1, t2, gen_indices<sizeof...(T)>{});
}

// sum_t //

template <typename ... Ts, std::size_t ... Is>
std::tuple<Ts...> sum_t(
        std::tuple<Ts...> const & t1,
        std::tuple<Ts...> const & t2,
        integer_sequence<std::size_t, Is...> const &)
{
     return { (std::get<Is>(t1) + std::get<Is>(t2))... };
}

template <typename ... Ts>
std::tuple<Ts...> operator+(
        std::tuple<Ts...> const & t1,
        std::tuple<Ts...> const & t2)
{
    return sum_t(t1, t2, gen_indices<sizeof...(Ts)>{});
}

// diff_t //

template <typename ... Ts, std::size_t ... Is>
std::tuple<Ts...> diff_t(
        std::tuple<Ts...> const & t1,
        std::tuple<Ts...> const & t2,
        integer_sequence<std::size_t, Is...> const &)
{
     return { (std::get<Is>(t1) - std::get<Is>(t2))... };
}

template <typename ... Ts>
std::tuple<Ts...> operator-(
        std::tuple<Ts...> const & t1,
        std::tuple<Ts...> const & t2)
{
    return diff_t(t1, t2, gen_indices<sizeof...(Ts)>{});
}

// worst and best //

template <typename> struct Helper;

template <typename... Ts>
struct Helper<std::tuple<Ts...>>
{
    static std::tuple<Ts...> max()
    {
        return { (std::numeric_limits<Ts>::max() / 3)... };
    }

    static std::tuple<Ts...> min()
    {
        return { (std::numeric_limits<Ts>::lowest() / 3)... };
    }
};

template <typename T>
T worst() { return Helper<T>::max(); }

template <typename T>
T best() { return Helper<T>::min(); }


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// strictly_better /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename, typename T>
struct HasStrictlyBetterMethod
{
    static_assert(
        std::integral_constant<T, false>::value,
        "Second template parameter needs to be of function type.");
};

template<typename C, typename Ret, typename... Args>
struct HasStrictlyBetterMethod<C, Ret(Args...)>
{

private:

    template<typename T>
    static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().strictly_better(std::declval<Args>()...)), Ret>::type;

    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<C>(0)) type;

public:

    static constexpr bool value = type::value;

};

template<typename LocalScheme>
bool strictly_better(
        LocalScheme&,
        const typename LocalScheme::GlobalCost& global_cost_1,
        const typename LocalScheme::GlobalCost& global_cost_2,
        std::false_type)
{
    return global_cost_1 < global_cost_2;
}

template<typename LocalScheme>
bool strictly_better(
        LocalScheme& local_scheme,
        const typename LocalScheme::GlobalCost& global_cost_1,
        const typename LocalScheme::GlobalCost& global_cost_2,
        std::true_type)
{
    return local_scheme.strictly_better(
            global_cost_1,
            global_cost_2);
}

template<typename LocalScheme>
bool strictly_better(
        LocalScheme& local_scheme,
        const typename LocalScheme::GlobalCost& global_cost_1,
        const typename LocalScheme::GlobalCost& global_cost_2)
{
    return strictly_better(
            local_scheme,
            global_cost_1,
            global_cost_2,
            std::integral_constant<
                bool,
                HasStrictlyBetterMethod<LocalScheme,
                bool(
                    const typename LocalScheme::GlobalCost&,
                    const typename LocalScheme::GlobalCost&)>::value>());
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// dominates ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename, typename T>
struct HasDominatesMethod
{
    static_assert(
        std::integral_constant<T, false>::value,
        "Second template parameter needs to be of function type.");
};

template<typename C, typename Ret, typename... Args>
struct HasDominatesMethod<C, Ret(Args...)>
{

private:

    template<typename T>
    static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().dominates(std::declval<Args>()...)), Ret>::type;

    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<C>(0)) type;

public:

    static constexpr bool value = type::value;

};

template<typename LocalScheme>
bool dominates(
        LocalScheme&,
        const typename LocalScheme::GlobalCost& global_cost_1,
        const typename LocalScheme::GlobalCost& global_cost_2,
        std::false_type)
{
    return dominates(global_cost_1, global_cost_2);
}

template<typename LocalScheme>
bool dominates(
        LocalScheme& local_scheme,
        const typename LocalScheme::GlobalCost& global_cost_1,
        const typename LocalScheme::GlobalCost& global_cost_2,
        std::true_type)
{
    return local_scheme.dominates(
            global_cost_1,
            global_cost_2);
}

template<typename LocalScheme>
bool dominates(
        LocalScheme& local_scheme,
        const typename LocalScheme::GlobalCost& global_cost_1,
        const typename LocalScheme::GlobalCost& global_cost_2)
{
    return dominates(
            local_scheme,
            global_cost_1,
            global_cost_2,
            std::integral_constant<
                bool,
                HasDominatesMethod<LocalScheme,
                bool(
                    const typename LocalScheme::GlobalCost&,
                    const typename LocalScheme::GlobalCost&)>::value>());
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// to_string ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename, typename T>
struct HasToStringMethod
{
    static_assert(
        std::integral_constant<T, false>::value,
        "Second template parameter needs to be of function type.");
};

template<typename C, typename Ret, typename... Args>
struct HasToStringMethod<C, Ret(Args...)>
{

private:

    template<typename T>
    static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().to_string(std::declval<Args>()...)), Ret>::type;

    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<C>(0)) type;

public:

    static constexpr bool value = type::value;

};

template<typename LocalScheme>
std::string to_string(
        LocalScheme&,
        const typename LocalScheme::GlobalCost& global_cost,
        std::false_type)
{
    return to_string(global_cost);
}

template<typename LocalScheme>
std::string to_string(
        LocalScheme& local_scheme,
        const typename LocalScheme::GlobalCost& global_cost,
        std::true_type)
{
    return local_scheme.to_string(global_cost);
}

template<typename LocalScheme>
std::string to_string(
        LocalScheme& local_scheme,
        const typename LocalScheme::GlobalCost& global_cost)
{
    return to_string(
            local_scheme,
            global_cost,
            std::integral_constant<
                bool,
                HasToStringMethod<LocalScheme,
                std::string(const typename LocalScheme::GlobalCost&)>::value>());
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Solution Pool /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename LocalScheme>
struct SolutionPoolComparator
{
    using Solution = typename LocalScheme::Solution;

    SolutionPoolComparator(const LocalScheme& local_scheme):
        local_scheme(local_scheme) { }

    const LocalScheme& local_scheme;

    bool operator()(
            const Solution& solution_1,
            const Solution& solution_2) const {
        return strictly_better(
                local_scheme,
                local_scheme.global_cost(solution_1),
                local_scheme.global_cost(solution_2));
    }
};

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Solution pool /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename LocalScheme>
class SolutionPool
{
    using Solution = typename LocalScheme::Solution;
    using GlobalCost = typename LocalScheme::GlobalCost;

public:

    SolutionPool(const LocalScheme& local_scheme, Counter size_max):
        local_scheme_(local_scheme),
        size_max_(size_max),
        solution_pool_comparator_(local_scheme),
        solutions_(solution_pool_comparator_),
        worst_(local_scheme.empty_solution()),
        best_(local_scheme.empty_solution()) { }

    const std::multiset<Solution, SolutionPoolComparator<LocalScheme>>& solutions() const { return solutions_; };

    /** Get the best solution of the pool. */
    const Solution& best() const { return best_; }

    /** Get the worst solution fo the pool. */
    const Solution& worst() const { return worst_; }

    /** Get the number of solutions in the pool. */
    Counter size() const { return solutions_.size(); }

    /** Get the local scheme. */
    const LocalScheme& local_scheme() const { return local_scheme_; }

    /** Add a solution to the pool. */
    int add(
            const Solution& solution)
    {
        // If the solution is worse than the worst solution of the pool, stop.
        if ((Counter)solutions_.size() >= size_max_) {
            if (!strictly_better(
                        local_scheme_,
                        local_scheme_.global_cost(solution),
                        local_scheme_.global_cost(*std::prev(solutions_.end())))) {
                return 0;
            }
        }
        //for (const auto& solution_tmp: solutions_)
        //    if (solution == solution_tmp)
        //        return 0;
        // If new best solution, display.
        bool new_best = (solutions_.size() == 0)
            || strictly_better(
                    local_scheme_,
                    local_scheme_.global_cost(solution),
                    local_scheme_.global_cost(*solutions_.begin()));
        // Add new solution to solution pool.
        solutions_.insert(solution);
        // If the pool size is now above its maximum allowed size, remove worst
        // solutions from it.
        if ((Counter)solutions_.size() > size_max_)
            solutions_.erase(std::prev(solutions_.end()));
        best_ = *solutions_.begin();
        worst_ = *std::prev(solutions_.end());
        return (new_best)? 2: 1;
    }

private:

    /** Local scheme. */
    const LocalScheme& local_scheme_;

    /** Maximum number of solutions in the solution pool. */
    Counter size_max_;

    /** Solution comparator. */
    SolutionPoolComparator<LocalScheme> solution_pool_comparator_;

    /** Solutions. */
    std::multiset<Solution, SolutionPoolComparator<LocalScheme>> solutions_;

    /** Worst solution of the pool. */
    Solution worst_;

    /** Best solution of the pool. */
    Solution best_;

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * Output structure.
 */
template <typename LocalScheme>
struct Output: optimizationtools::Output
{
    /** Constructor. */
    Output(
            const LocalScheme& local_scheme,
            Counter maximum_size_of_the_solution_pool):
        solution_pool(local_scheme, maximum_size_of_the_solution_pool) { }


    /** Solution. */
    SolutionPool<LocalScheme> solution_pool;

    /** Elapsed time. */
    double time = 0.0;


    virtual nlohmann::json to_json() const
    {
        return nlohmann::json {
            //{"Solution", to_json(local_scheme_, solution)},
            {"Value", to_string(solution_pool.local_scheme(), solution_pool.local_scheme().global_cost(solution_pool.best()))},
            {"Time", time}
        };
    }

    virtual int format_width() const { return 11; }

    virtual void format(std::ostream& os) const
    {
        int width = format_width();
        os
            << std::setw(width) << std::left << "Value: " << to_string(solution_pool.local_scheme(), solution_pool.local_scheme().global_cost(solution_pool.best())) << std::endl
            << std::setw(width) << std::left << "Time (s): " << time << std::endl
            ;
    }
};

template <typename LocalScheme>
using NewSolutionCallback = std::function<void(const Output<LocalScheme>&)>;

template <typename LocalScheme>
struct Parameters: optimizationtools::Parameters
{
    using Solution = typename LocalScheme::Solution;
    using GlobalCost = typename LocalScheme::GlobalCost;

    /** Callback function called when a new best solution is found. */
    NewSolutionCallback<LocalScheme> new_solution_callback = [](const Output<LocalScheme>&) { };

    /** Ids of generated initial solutions. */
    std::vector<Counter> initial_solution_ids = {0};

    /** User-provided initial solutions. */
    std::vector<Solution> initial_solutions;

    /** Maximum size of the solution pool. */
    Counter maximum_size_of_the_solution_pool = 1;

    /** Seed. */
    Seed seed = 0;

    /**
     * Goal.
     *
     * The alglorithm stops as soon as a solution with a better global cost is
     * found.
     */
    bool has_goal = false;

    /** Goal. */
    GlobalCost goal;


    using optimizationtools::Parameters::to_json;
    virtual nlohmann::json to_json(
            const LocalScheme& local_scheme) const
    {
        nlohmann::json json = optimizationtools::Parameters::to_json();
        json.merge_patch({
            {"MaximumSizeOfTheSolutionPool", maximum_size_of_the_solution_pool},
            {"Seed", seed},
            {"Goal", ((has_goal)? to_string(local_scheme, goal): "")}});
        return json;
    }

    virtual int format_width() const override { return 23; }

    using optimizationtools::Parameters::format;
    virtual void format(
            std::ostream& os,
            const LocalScheme& local_scheme) const
    {
        optimizationtools::Parameters::format(os);
        int width = format_width();
        os
            << std::setw(width) << std::left << "Solution pool size: " << maximum_size_of_the_solution_pool << std::endl
            << std::setw(width) << std::left << "Seed: " << seed << std::endl
            << std::setw(width) << std::left << "Goal: " << ((has_goal)? to_string(local_scheme, goal): "") << std::endl
            ;
    }
};

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// global_cost_goal ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename, typename T>
struct HasGlobalCostGoalMethod
{
    static_assert(
        std::integral_constant<T, false>::value,
        "Second template parameter needs to be of function type.");
};

template<typename C, typename Ret, typename... Args>
struct HasGlobalCostGoalMethod<C, Ret(Args...)>
{

private:

    template<typename T>
    static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().global_cost_goal(std::declval<Args>()...)), Ret>::type;

    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<C>(0)) type;

public:

    static constexpr bool value = type::value;

};

template<typename LocalScheme>
typename LocalScheme::GlobalCost global_cost_goal(
        const LocalScheme&,
        double,
        std::false_type)
{
    return best<typename LocalScheme::GlobalCost>();
}

template<typename LocalScheme>
typename LocalScheme::GlobalCost global_cost_goal(
        const LocalScheme& local_scheme,
        double value,
        std::true_type)
{
    return local_scheme.global_cost_goal(value);
}

template<typename LocalScheme>
typename LocalScheme::GlobalCost global_cost_goal(
        const LocalScheme& local_scheme,
        double value)
{
    using GlobalCost = typename LocalScheme::GlobalCost;

    return global_cost_goal(
            local_scheme,
            value,
            std::integral_constant<
                bool,
                HasGlobalCostGoalMethod<LocalScheme,
                GlobalCost(double)>::value>());
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// solution_write ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename, typename T>
struct HasSolutionWriteMethod
{
    static_assert(
        std::integral_constant<T, false>::value,
        "Second template parameter needs to be of function type.");
};

template<typename C, typename Ret, typename... Args>
struct HasSolutionWriteMethod<C, Ret(Args...)>
{

private:

    template<typename T>
    static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().solution_write(std::declval<Args>()...)), Ret>::type;

    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<C>(0)) type;

public:

    static constexpr bool value = type::value;

};

template<typename LocalScheme>
void solution_write(
        const LocalScheme&,
        const typename LocalScheme::Solution&,
        const std::string&,
        std::false_type)
{
}

template<typename LocalScheme>
void solution_write(
        const LocalScheme& local_scheme,
        const typename LocalScheme::Solution& solution,
        const std::string certificate_path,
        std::true_type)
{
    return local_scheme.solution_write(solution, certificate_path);
}

template<typename LocalScheme>
void solution_write(
        const LocalScheme& local_scheme,
        const typename LocalScheme::Solution& solution,
        const std::string certificate_path)
{
    using Solution = typename LocalScheme::Solution;

    return solution_write(
            local_scheme,
            solution,
            certificate_path,
            std::integral_constant<
                bool,
                HasSolutionWriteMethod<LocalScheme,
                void(const Solution&, const std::string&)>::value>());
}

}
