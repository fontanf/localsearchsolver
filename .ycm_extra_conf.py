def Settings(**kwargs):
    return {
            'flags': [
                '-x', 'c++',
                '-Wall', '-Wextra', '-Werror',

                '-I', '.',
                '-I', './bazel-localsearchsolver/',
                '-I', './bazel-localsearchsolver/external/'
                'json/single_include/',
                '-I', './bazel-localsearchsolver/external/'
                'googletest/googletest/include/',
                '-I', './bazel-localsearchsolver/external/'
                'boost/',
                '-I', './bazel-localsearchsolver/external/'
                'simdjson/simdjson-0.9.2/singleheader',

                # CPLEX
                '-DCPLEX_FOUND',
                '-DIL_STD',
                '-I', '/opt/ibm/ILOG/CPLEX_Studio129/concert/include/',
                '-I', '/opt/ibm/ILOG/CPLEX_Studio129/cplex/include/',
                '-I', '/opt/ibm/ILOG/CPLEX_Studio129/cpoptimizer/include/',

                # Gurobi
                '-DGUROBI_FOUND',
                '-I', '/home/florian/Programmes/gurobi811/linux64/include/',

                # optimizationtools
                '-I', './bazel-localsearchsolver/external/'
                # '-I', './../',
                'optimizationtools/',

                # orproblems
                '-I', './bazel-localsearchsolver/external/'
                # '-I', './../',
                'orproblems/',
                ],
            }
