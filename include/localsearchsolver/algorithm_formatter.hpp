#pragma once

#include "localsearchsolver/common.hpp"

#include <mutex>

namespace localsearchsolver
{

template <typename LocalScheme>
class AlgorithmFormatter
{
    using Solution = typename LocalScheme::Solution;

public:

    /** Constructor. */
    AlgorithmFormatter(
            const LocalScheme& local_scheme,
            const Parameters<LocalScheme>& parameters,
            Output<LocalScheme>& output):
        local_scheme_(local_scheme),
        parameters_(parameters),
        output_(output),
        os_(parameters.create_os()) { }

    /** Print the header. */
    void start(
            const std::string& algorithm_name);

    /** Print the header. */
    void print_header();

    /** Print current state. */
    void print(
            const std::stringstream& s);

    /** Update the solution. */
    void update_solution(
            const Solution& solution_new,
            const std::stringstream& s);

    /** Method to call at the end of the algorithm. */
    void end();

private:

    /*
     * Private methods
     */

    /*
     * local_scheme_instance_format
     */

    template<typename, typename T>
    struct HasLocalSchemeInstanceFormatMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasLocalSchemeInstanceFormatMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().instance_format(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    void local_scheme_instance_format(
            std::false_type)
    {
    }

    void local_scheme_instance_format(
            std::true_type)
    {
        *os_
            <<  std::endl
            << "Problem" << std::endl
            << "-------" << std::endl;
        output_.solution_pool.local_scheme().instance_format(
                *os_,
                parameters_.verbosity_level);
    }

    void local_scheme_instance_format()
    {
        local_scheme_instance_format(
                std::integral_constant<
                    bool,
                    HasLocalSchemeInstanceFormatMethod<
                        LocalScheme,
                        void(std::ostream&, int)>::value>());
    }

    /*
     * local_scheme_parameters_format
     */

    template<typename, typename T>
    struct HasLocalSchemeParametersFormatMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasLocalSchemeParametersFormatMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().parameters_format(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    void local_scheme_parameters_format(
            std::false_type)
    {
    }

    void local_scheme_parameters_format(
            std::true_type)
    {
        *os_
            <<  std::endl
            << "Local scheme parameters" << std::endl
            << "-----------------------" << std::endl;
        output_.solution_pool.local_scheme().parameters_format(
                *os_,
                parameters_.verbosity_level);
    }

    void local_scheme_parameters_format()
    {
        local_scheme_parameters_format(
                std::integral_constant<
                    bool,
                    HasLocalSchemeParametersFormatMethod<
                        LocalScheme,
                        void(std::ostream&, int)>::value>());
    }

    /*
     * local_scheme_statistics_format
     */

    template<typename, typename T>
    struct HasLocalSchemeStatisticsFormatMethod
    {
        static_assert(
                std::integral_constant<T, false>::value,
                "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasLocalSchemeStatisticsFormatMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().statistics_format(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    void local_scheme_statistics_format(
            std::false_type)
    {
    }

    void local_scheme_statistics_format(
            std::true_type)
    {
        *os_
            << std::endl
            << "Local scheme statistics" << std::endl
            << "-----------------------" << std::endl;
        output_.solution_pool.local_scheme().statistics_format(
                *os_,
                parameters_.verbosity_level);
    }

    void local_scheme_statistics_format()
    {
        local_scheme_statistics_format(
                std::integral_constant<
                    bool,
                    HasLocalSchemeStatisticsFormatMethod<
                        LocalScheme,
                        void(std::ostream&, int)>::value>());
    }

    /*
     * local_scheme_solution_format
     */

    template<typename, typename T>
    struct HasLocalSchemeSolutionFormatMethod
    {
        static_assert(
            std::integral_constant<T, false>::value,
            "Second template parameter needs to be of function type.");
    };

    template<typename C, typename Ret, typename... Args>
    struct HasLocalSchemeSolutionFormatMethod<C, Ret(Args...)>
    {

    private:

        template<typename T>
        static constexpr auto check(T*) -> typename std::is_same<decltype(std::declval<T>().solution_format(std::declval<Args>()...)), Ret>::type;

        template<typename>
        static constexpr std::false_type check(...);

        typedef decltype(check<C>(0)) type;

    public:

        static constexpr bool value = type::value;

    };

    void local_scheme_solution_format(
            const Solution&,
            std::false_type)
    {
    }

    void local_scheme_solution_format(
            const Solution& solution,
            std::true_type)
    {
        *os_
            << std::endl
            << "Local scheme solution" << std::endl
            << "---------------------" << std::endl;
        return local_scheme_.solution_format(
                solution,
                *os_,
                parameters_.verbosity_level);
    }

    void local_scheme_solution_format(
            const Solution& solution)
    {
        using Solution = typename LocalScheme::Solution;

        return local_scheme_solution_format(
                solution,
                std::integral_constant<
                    bool,
                    HasLocalSchemeSolutionFormatMethod<
                        LocalScheme,
                        void(const Solution&, std::ostream&, int)>::value>());
    }

    /*
     * Private attributes
     */

    /** Local scheme. */
    const LocalScheme& local_scheme_;

    /** Parameters. */
    const Parameters<LocalScheme>& parameters_;

    /** Output. */
    Output<LocalScheme>& output_;

    /** Messages stream. */
    std::unique_ptr<optimizationtools::ComposeStream> os_;

    /** mutex. */
    std::mutex mutex_;

};

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// Template implementations //////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename LocalScheme>
void AlgorithmFormatter<LocalScheme>::start(
        const std::string& algorithm_name)
{
    output_.json["Parameters"] = parameters_.to_json(output_.solution_pool.local_scheme());

    if (parameters_.verbosity_level == 0)
        return;
    *os_
        << "=======================================" << std::endl
        << "           LocalSearchSolver           " << std::endl
        << "=======================================" << std::endl
        ;
    local_scheme_instance_format();
    local_scheme_parameters_format();
    *os_
        << std::endl
        << "Algorithm" << std::endl
        << "---------" << std::endl
        << algorithm_name << std::endl
        << std::endl
        << "Parameters" << std::endl
        << "----------" << std::endl;
    parameters_.format(*os_, output_.solution_pool.local_scheme());
}

template <typename LocalScheme>
void AlgorithmFormatter<LocalScheme>::print_header()
{
    if (parameters_.verbosity_level == 0)
        return;
    *os_
        << std::right
        << std::endl
        << std::setw(10) << "Time"
        << std::setw(40) << "Value"
        << std::setw(40) << "Comment"
        << std::endl
        << std::setw(10) << "----"
        << std::setw(40) << "-----"
        << std::setw(40) << "-------"
        << std::endl;
}

template <typename LocalScheme>
void AlgorithmFormatter<LocalScheme>::print(
        const std::stringstream& s)
{
    if (parameters_.verbosity_level == 0)
        return;
    std::streamsize precision = std::cout.precision();
    *os_
        << std::setw(10) << std::fixed << std::setprecision(3) << output_.time << std::defaultfloat << std::setprecision(precision)
        << std::setw(40) << to_string(local_scheme_, local_scheme_.global_cost(output_.solution_pool.best()))
        << std::setw(40) << s.str()
        << std::endl;
}

template <typename LocalScheme>
void AlgorithmFormatter<LocalScheme>::update_solution(
        const typename LocalScheme::Solution& solution_new,
        const std::stringstream& s)
{
    mutex_.lock();
    if (output_.solution_pool.add(solution_new) == 2) {
        output_.time = parameters_.timer.elapsed_time();
        print(s);
        output_.json["IntermediaryOutputs"].push_back(output_.to_json());
        parameters_.new_solution_callback(output_);
    }
    mutex_.unlock();
}

template <typename LocalScheme>
void AlgorithmFormatter<LocalScheme>::end()
{
    output_.time = parameters_.timer.elapsed_time();
    output_.json["Output"] = output_.to_json();

    if (parameters_.verbosity_level == 0)
        return;
    *os_
        << std::endl
        << "Algorithm statistics" << std::endl
        << "--------------------" << std::endl;
    output_.format(*os_);
    local_scheme_statistics_format();
    local_scheme_solution_format(output_.solution_pool.best());
}

}
