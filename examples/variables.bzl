CPLEX_COPTS = select({
            "//examples:cplex_build": [
                    "-DCPLEX_FOUND",
                    "-m64",
                    "-DIL_STD"],
            "//conditions:default": []})
GUROBI_COPTS = select({
            "//examples:gurobi_build": ["-DGUROBI_FOUND"],
            "//conditions:default": []})

CPLEX_DEP = select({
            "//examples:cplex_build": ["@cplex//:cplex"],
            "//conditions:default": []})
GUROBI_DEP = select({
            "//examples:gurobi_build": ["@gurobi//:gurobi"],
            "//conditions:default": []})

