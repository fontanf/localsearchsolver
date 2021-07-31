#pragma once

#include "optimizationtools/info.hpp"
#include "optimizationtools/utils.hpp"

#include <cstdint>
#include <set>

namespace localsearchsolver
{

using Seed = int64_t;
using Counter = int64_t;

/***************************** Tuple operations ******************************/

template<typename T, T...>
struct integer_sequence { };

template<std::size_t N, std::size_t... I>
struct gen_indices : gen_indices<(N - 1), (N - 1), I...> { };
template<std::size_t... I>
struct gen_indices<0, I...> : integer_sequence<std::size_t, I...> { };

/*
 * to_string(t)
 */

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

/*
 * update_move_cost(t, t_ref)
 */

template <typename ... T, std::size_t ... I>
std::tuple<T...> update_move_cost_impl(
        const std::tuple<T...>& t1,
        const std::tuple<T...>& t2,
        integer_sequence<std::size_t, I...>)
{
    return {
        (std::min(std::get<I>(t1), std::get<I>(t2)))...,
        std::get<sizeof...(I)>(t1) };
}

template <typename ... T>
std::tuple<T...> update_move_cost(
        const std::tuple<T...>& t1,
        const std::tuple<T...>& t2)
{
    return update_move_cost_impl(t1, t2, gen_indices<sizeof...(T) - 1u>{});
}

/*
 * dominates(t1, t2)
 */

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


template <typename Scheme>
struct SolutionPoolComparator
{
    typedef typename Scheme::Solution Solution;

    SolutionPoolComparator(const Scheme& local_scheme):
        local_scheme(local_scheme) {  }

    const Scheme& local_scheme;

    bool operator()(
            const Solution& solution_1,
            const Solution& solution_2) const {
        return local_scheme.global_cost(solution_1)
            < local_scheme.global_cost(solution_2);
    }
};

template <typename Scheme>
class SolutionPool
{
    typedef typename Scheme::Solution Solution;
    typedef typename Scheme::GlobalCost GlobalCost;

public:

    SolutionPool(const Scheme& local_scheme, Counter size_max):
        scheme_(local_scheme),
        size_max_(size_max),
        solution_pool_comparator_(local_scheme),
        solutions_(solution_pool_comparator_)
    {
        Solution solution = local_scheme.empty_solution();
        solutions_.insert(solution);
        worst_ = solution;
        best_ = solution;
    }

    virtual ~SolutionPool() { }

    const std::multiset<Solution, SolutionPoolComparator<Scheme>>& solutions() const { return solutions_; };
    const Solution& best() { return best_; }
    const Solution& worst() { return worst_; }
    Counter size() const { return solutions_.size(); }

    int add(
            const Solution& solution,
            const std::stringstream& ss,
            optimizationtools::Info& info)
    {
        info.output->mutex_solutions.lock();
        // If the solution is worse than the worst solution of the pool, stop.
        if ((Counter)solutions_.size() >= size_max_) {
            if (scheme_.global_cost(solution)
                    > scheme_.global_cost(*std::prev(solutions_.end()))) {
                info.output->mutex_solutions.unlock();
                return 0;
            }
        }
        //for (const auto& solution_tmp: solutions_)
        //    if (solution == solution_tmp)
        //        return 0;
        // If new best solution, display.
        bool new_best = scheme_.global_cost(solution) < scheme_.global_cost(*solutions_.begin());
        // Add new solution to solution pool.
        solutions_.insert(solution);
        if (new_best) {
            info.output->number_of_solutions++;
            double t = info.elapsed_time();
            std::string sol_str = "Solution" + std::to_string(info.output->number_of_solutions);
            PUT(info, sol_str, "Value", to_string(scheme_.global_cost(solution)));
            PUT(info, sol_str, "Time", t);
            PUT(info, sol_str, "Comment", ss.str());
            if (!info.output->only_write_at_the_end) {
                info.write_json_output();
                scheme_.write(*solutions_.begin(), info.output->certificate_path);
            }
        }
        // If the pool size is now above its maximum allowed size, remove worst
        // solutions from it.
        if ((Counter)solutions_.size() > size_max_)
            solutions_.erase(std::prev(solutions_.end()));
        best_ = *solutions_.begin();
        worst_ = *std::prev(solutions_.end());
        info.output->mutex_solutions.unlock();
        return (new_best)? 2: 1;
    }

    void display_init(optimizationtools::Info& info)
    {
        VER(info,
                "----------------------------------------------------------------------" << std::endl
                << std::left << std::setw(16) << "Time"
                << std::left << std::setw(40) << "Value"
                << std::left << std::setw(32) << "Comment"
                << std::endl
                << "----------------------------------------------------------------------" << std::endl
                );
    }

    void display(const std::stringstream& ss, optimizationtools::Info& info)
    {
        double t = info.elapsed_time();
        VER(info,
                std::left << std::setw(16) << t
                << std::left << std::setw(40) << to_string(scheme_.global_cost(best()))
                << std::left << std::setw(32) << ss.str()
                << std::endl);
    }

    void display_end(optimizationtools::Info& info)
    {
        double t = info.elapsed_time();
        std::string sol_str = "Solution";
        VER(info,
                "---" << std::endl
                << "Time: " << t << std::endl
                << "Value: " << to_string(scheme_.global_cost(*solutions_.begin())) << std::endl);
        PUT(info, sol_str, "Time", t);
        PUT(info, sol_str, "Value", to_string(scheme_.global_cost(*solutions_.begin())));
        info.write_json_output();
        scheme_.write(best_, info.output->certificate_path);
    }

private:

    const Scheme& scheme_;
    Counter size_max_;
    SolutionPoolComparator<Scheme> solution_pool_comparator_;
    std::multiset<Solution, SolutionPoolComparator<Scheme>> solutions_;
    Solution worst_;
    Solution best_;

};

}

