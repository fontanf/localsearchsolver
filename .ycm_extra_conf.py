def Settings(**kwargs):
    return {
            'flags': [
                '-x', 'c++',
                '-Wall', '-Wextra', '-Werror',
                '-DGUROBI_FOUND',
                '-DCPLEX_FOUND',
                '-DIL_STD', # Cplex
                '-I', '.',
                '-I', './bazel-localsearchsolver/',
                '-I', './bazel-localsearchsolver/external/json/single_include/',
                '-I', './bazel-localsearchsolver/external/googletest/googletest/include/',
                '-I', './bazel-localsearchsolver/external/boost/',
                '-I', './bazel-localsearchsolver/external/simdjson/simdjson-0.9.2/singleheader',
                '-I', './bazel-localsearchsolver/external/optimizationtools/',
                # '-I', './../optimizationtools/',
                '-I', './bazel-localsearchsolver/external/orproblems/',
                # '-I', './../orproblems/',
                '-I', '/opt/ibm/ILOG/CPLEX_Studio129/concert/include/',
                '-I', '/opt/ibm/ILOG/CPLEX_Studio129/cplex/include/',
                '-I', '/opt/ibm/ILOG/CPLEX_Studio129/cpoptimizer/include/',
                '-I', '/home/florian/Programmes/gurobi811/linux64/include/',
                ],
            }
