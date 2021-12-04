cc_library(
        name = "localsearchsolver",
        hdrs = [
                "restarting_local_search.hpp",
                "iterated_local_search.hpp",
                "best_first_local_search.hpp",
                "genetic_local_search.hpp",
                "sequencing.hpp",
                "read_args.hpp",
                "common.hpp",
        ],
        deps = [
                "@optimizationtools//optimizationtools:info",
                "@optimizationtools//optimizationtools:utils",
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
