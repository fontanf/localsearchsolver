cc_library(
        name = "localsearchsolver",
        hdrs = [
                "common.hpp",
                "algorithm_formatter.hpp",
                "multi_start_local_search.hpp",
                "iterated_local_search.hpp",
                "best_first_local_search.hpp",
                "genetic_local_search.hpp",
                "sequencing.hpp",
                "read_args.hpp",
        ],
        deps = [
                "@optimizationtools//optimizationtools/utils:utils",
                "@optimizationtools//optimizationtools/containers:containers",
        ],
        visibility = ["//visibility:public"],
)

cc_library(
        name = "read_args",
        hdrs = ["read_args.hpp"],
        deps = [
                ":localsearchsolver",
                "@boost//:program_options",
        ],
        visibility = ["//visibility:public"],
)
