#include "examples/roadef2020.hpp"

#include "localsearchsolver/a_star.hpp"

#include <numeric>
#include <unordered_set>

#if CPLEX_FOUND
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
#endif

#if GUROBI_FOUND
#include "gurobi_c++.h"
#endif

using namespace localsearchsolver::roadef2020;

void LocalScheme::add(
        Solution& solution,
        InterventionId j,
        Time t_start) const
{
    assert(t_start >= 0);
    assert(t_start < instance_.horizon());
    assert(t_start <= instance_.start_max(j));
    assert(solution.intervention_starts[j] == -1);

    Time t_end = t_start + instance_.duration(j, t_start);
    ScenarioId lcm = instance_.least_common_multiple();

    // Update intervention_times_.
    solution.intervention_number++;
    solution.intervention_starts[j] = t_start;

    for (Time t_cur = t_start; t_cur < t_end; ++t_cur) {
        // Update solution.time_steps.
        for (ResourcePos r_pos = 0; r_pos < instance_.resource_number(j); ++r_pos) {
            ResourceId r = instance_.resource(j, r_pos);
            Workload w = instance_.workload(j, r_pos, t_cur, t_start);
            Workload w_min = instance_.workload_min(r, t_cur);
            Workload w_max = instance_.workload_max(r, t_cur);
            // Update overwork
            if (solution.time_steps[t_cur].workloads[r] >= w_max) {
                solution.overwork += w;
            } else if (solution.time_steps[t_cur].workloads[r] + w <= w_max) {
            } else {
                solution.overwork += (solution.time_steps[t_cur].workloads[r] + w - w_max);
            }
            // Update underwork
            if (solution.time_steps[t_cur].workloads[r] >= w_min) {
            } else if (solution.time_steps[t_cur].workloads[r] + w_min <= w_min) {
                solution.underwork -= w;
            } else {
                solution.underwork -= (w_min - solution.time_steps[t_cur].workloads[r]);
            }
            // Update workload
            solution.time_steps[t_cur].workloads[r] += w;
            assert(solution.time_steps[t_cur].workloads[r] >= 0);
        }

        // Update conflicts_
        SeasonId season = instance_.season(t_cur);
        if (season != -1) {
            for (ExclusionId e: instance_.exclusions(j, season)) {
                InterventionId j2 = instance_.exclusion(e).j(j);
                if (solution.intervention_starts[j2] == -1)
                    continue;
                Time t2_start = solution.intervention_starts[j2];
                if (t2_start <= t_cur && t_cur < t2_start + instance_.duration(j2, t2_start)) {
                    SolutionConflict conflict;
                    conflict.j1 = j;
                    conflict.j2 = j2;
                    conflict.t_cur = t_cur;
                    solution.conflicts.push_back(conflict);
                }
            }
        }

        // Remove previous cost
        ScenarioId s_pos = std::ceil(instance_.quantile() * instance_.scenario_number(t_cur)) - 1;
        solution.mean_cost -= lcm / instance_.scenario_number(t_cur) * solution.time_steps[t_cur].risk_sum;
        Risk ee = lcm * solution.time_steps[t_cur].risks[solution.time_steps[t_cur].sorted_scenarios[s_pos]]
            - lcm / instance_.scenario_number(t_cur) * solution.time_steps[t_cur].risk_sum;
        if (ee > 0)
            solution.expected_excess -= ee;
        assert(solution.mean_cost >= 0);
        assert(solution.expected_excess >= 0);

        // Update risks
        solution.time_steps[t_cur].risk_sum = 0;
        for (ScenarioId s = 0; s < instance_.scenario_number(t_cur); ++s) {
            solution.time_steps[t_cur].risks[s] += instance_.risk(j, t_cur, t_start, s);
            solution.time_steps[t_cur].risk_sum += solution.time_steps[t_cur].risks[s];
        }
        // Update sorted scenarios
        std::nth_element(
                solution.time_steps[t_cur].sorted_scenarios.begin(),
                solution.time_steps[t_cur].sorted_scenarios.begin() + s_pos,
                solution.time_steps[t_cur].sorted_scenarios.end(),
                [&solution, t_cur](ScenarioId s1, ScenarioId s2)
                {
                    return solution.time_steps[t_cur].risks[s1]
                        < solution.time_steps[t_cur].risks[s2];
                });
        // Update cost
        solution.mean_cost += lcm / instance_.scenario_number(t_cur) * solution.time_steps[t_cur].risk_sum;
        ee = lcm * solution.time_steps[t_cur].risks[solution.time_steps[t_cur].sorted_scenarios[s_pos]]
            - lcm / instance_.scenario_number(t_cur) * solution.time_steps[t_cur].risk_sum;
        if (ee > 0)
            solution.expected_excess += ee;
        assert(solution.mean_cost >= 0);
        assert(solution.expected_excess >= 0);
    }
}

void LocalScheme::remove(
        Solution& solution,
        InterventionId j) const
{
    assert(solution.intervention_starts[j] != -1);
    Time t_start = solution.intervention_starts[j];
    Time t_end = t_start + instance_.duration(j, t_start);
    ScenarioId lcm = instance_.least_common_multiple();

    // Update intervention_times_.
    solution.intervention_number--;
    solution.intervention_starts[j] = -1;

    // Update solution.conflicts
    for (auto it = solution.conflicts.begin(); it != solution.conflicts.end();) {
        if (it->j1 != j && it->j2 != j) {
            ++it;
        } else {
            *it = solution.conflicts.back();
            solution.conflicts.pop_back();
        }
    }

    // Update solution.time_steps.
    for (Time t_cur = t_start; t_cur < t_end; ++t_cur) {
        // Update solution.time_steps.
        for (ResourcePos r_pos = 0; r_pos < instance_.resource_number(j); ++r_pos) {
            ResourceId r = instance_.resource(j, r_pos);
            Workload w = instance_.workload(j, r_pos, t_cur, t_start);
            Workload w_min = instance_.workload_min(r, t_cur);
            Workload w_max = instance_.workload_max(r, t_cur);
            // Update overwork
            if (solution.time_steps[t_cur].workloads[r] - w >= w_max) {
                solution.overwork -= w;
            } else if (solution.time_steps[t_cur].workloads[r] <= w_max) {
            } else {
                solution.overwork -= (solution.time_steps[t_cur].workloads[r] - w_max);
            }
            // Update underwork
            if (solution.time_steps[t_cur].workloads[r] - w >= w_min) {
            } else if (solution.time_steps[t_cur].workloads[r] <= w_min) {
                solution.underwork += w;
            } else {
                solution.underwork += (solution.time_steps[t_cur].workloads[r] - w_min);
            }
            // Update workload
            solution.time_steps[t_cur].workloads[r] -= w;
        }

        // Remove previous cost
        ScenarioId s_pos = std::ceil(instance_.quantile() * instance_.scenario_number(t_cur)) - 1;
        solution.mean_cost -= lcm / instance_.scenario_number(t_cur) * solution.time_steps[t_cur].risk_sum;
        Risk ee = lcm * solution.time_steps[t_cur].risks[solution.time_steps[t_cur].sorted_scenarios[s_pos]]
            - lcm / instance_.scenario_number(t_cur) * solution.time_steps[t_cur].risk_sum;
        if (ee > 0)
            solution.expected_excess -= ee;
        assert(solution.mean_cost >= 0);
        assert(solution.expected_excess >= 0);

        // Update risks
        solution.time_steps[t_cur].risk_sum = 0;
        for (ScenarioId s = 0; s < instance_.scenario_number(t_cur); ++s) {
            solution.time_steps[t_cur].risks[s] -= instance_.risk(j, t_cur, t_start, s);
            solution.time_steps[t_cur].risk_sum += solution.time_steps[t_cur].risks[s];
        }
        // Update sorted scenarios
        std::nth_element(
                solution.time_steps[t_cur].sorted_scenarios.begin(),
                solution.time_steps[t_cur].sorted_scenarios.begin() + s_pos,
                solution.time_steps[t_cur].sorted_scenarios.end(),
                [&solution, t_cur](ScenarioId s1, ScenarioId s2)
                {
                    return solution.time_steps[t_cur].risks[s1]
                        < solution.time_steps[t_cur].risks[s2];
                });
        // Update cost
        solution.mean_cost += lcm / instance_.scenario_number(t_cur) * solution.time_steps[t_cur].risk_sum;
        ee = lcm * solution.time_steps[t_cur].risks[solution.time_steps[t_cur].sorted_scenarios[s_pos]]
            - lcm / instance_.scenario_number(t_cur) * solution.time_steps[t_cur].risk_sum;
        if (ee > 0)
            solution.expected_excess += ee;
        assert(solution.mean_cost >= 0);
        assert(solution.expected_excess >= 0);
    }
}

LocalScheme::GlobalCost LocalScheme::cost_add(
        const Solution& solution,
        InterventionId j,
        Time t_start,
        const GlobalCost& cutoff)
{
    assert(t_start >= 0);
    assert(t_start < instance_.horizon());
    assert(t_start <= instance_.start_max(j));
    assert(solution.intervention_starts[j] == -1);

    Time t_end = t_start + instance_.duration(j, t_start);
    GlobalCost c = global_cost(solution);
    intervention_number(c)--; // because -intervention_number

    // Check disjonctive constraints.
    for (Time t_cur = t_start; t_cur < t_end; ++t_cur) {
        SeasonId season = instance_.season(t_cur);
        if (season != -1) {
            for (ExclusionId e: instance_.exclusions(j, season)) {
                InterventionId j2 = instance_.exclusion(e).j(j);
                if (solution.intervention_starts[j2] == -1)
                    continue;
                Time t2_start = solution.intervention_starts[j2];
                if (t2_start <= t_cur && t_cur < t2_start + instance_.duration(j2, t2_start)) {
                    conflict_number(c)++;
                    if (c >= cutoff)
                        return cutoff;
                }
            }
        }
    }

    // Check resource constraint.
    for (Time t_cur = t_start; t_cur < t_end; ++t_cur) {
        for (ResourcePos r_pos = 0; r_pos < instance_.resource_number(j); ++r_pos) {
            Workload w = instance_.workload(j, r_pos, t_cur, t_start);
            if (w == 0)
                continue;
            ResourceId r = instance_.resource(j, r_pos);
            Workload w_cur = solution.time_steps[t_cur].workloads[r];
            Workload w_max = instance_.workload_max(r, t_cur);
            if (w_cur >= w_max) {
                overwork(c) += w;
                if (c >= cutoff)
                    return cutoff;
            } else if (w_cur + w <= w_max) {
            } else {
                overwork(c) += (w_cur + w - w_max);
                if (c >= cutoff)
                    return cutoff;
            }
        }
    }

    Cost mean_cost = solution.mean_cost;
    Cost expected_excess = solution.expected_excess;
    ScenarioId lcm = instance_.least_common_multiple();
    for (Time t_cur = t_start; t_cur < t_end; ++t_cur) {
        ScenarioId st = instance_.scenario_number(t_cur);
        ScenarioId s_pos = std::ceil(instance_.quantile() * st) - 1;

        // Remove previous cost
        mean_cost -= lcm / st * solution.time_steps[t_cur].risk_sum;
        Risk ee = lcm * solution.time_steps[t_cur].risks[solution.time_steps[t_cur].sorted_scenarios[s_pos]]
            - lcm / st * solution.time_steps[t_cur].risk_sum;
        if (ee > 0)
            expected_excess -= ee;

        // Update risks
        Risk risk_sum = 0;
        for (ScenarioId s = 0; s < st; ++s) {
            sorted_scenarios_[s] = solution.time_steps[t_cur].sorted_scenarios[s];
            risks_[s] = solution.time_steps[t_cur].risks[s] + instance_.risk(j, t_cur, t_start, s);
            risk_sum += risks_[s];
        }
        // Update sorted scenarios
        std::nth_element(
                sorted_scenarios_.begin(),
                sorted_scenarios_.begin() + s_pos,
                sorted_scenarios_.begin() + st,
                [this](ScenarioId s1, ScenarioId s2)
                {
                    return risks_[s1] < risks_[s2];
                });
        // Update cost
        mean_cost += lcm / st * risk_sum;
        ee = lcm * risks_[sorted_scenarios_[s_pos]]
            - lcm / st * risk_sum;
        if (ee > 0)
            expected_excess += ee;

        cost(c) = instance_.alpha_1() * mean_cost
            + instance_.alpha_2() * expected_excess;
        if (c >= cutoff)
            return cutoff;
    }

    cost(c) = instance_.alpha_1() * mean_cost
        + instance_.alpha_2() * expected_excess;
    assert(cost(c) >= 0);
    return c;
}

/***************************** Initial solutions *****************************/

LocalScheme::Solution LocalScheme::empty_solution() const
{
    Solution solution;
    solution.intervention_starts.resize(instance_.intervention_number(), -1);
    solution.time_steps.resize(instance_.horizon());
    for (Time t = 0; t < instance_.horizon(); ++t) {
        solution.time_steps[t].workloads.resize(instance_.resource_number(), 0);
        solution.time_steps[t].risks.resize(instance_.scenario_number(t), 0);
        solution.time_steps[t].sorted_scenarios.resize(instance_.scenario_number(t));
        for (ScenarioId s = 0; s < instance_.scenario_number(t); ++s)
            solution.time_steps[t].sorted_scenarios[s] = s;
        // Initialize underwork
        for (ResourceId r = 0; r < instance_.resource_number(); ++r)
            solution.underwork += instance_.workload_min(r, t);
    }
    return solution;
}

LocalScheme::Solution LocalScheme::initial_solution(
        Counter,
        std::mt19937_64& generator)
{
    Solution solution = empty_solution();

    // Solve the reduced instance.
    if (parameters_.reduced_instance_time > 0
            && instance_.reduced_instance() != nullptr) {
        Parameters parameters_local_scheme = parameters_;
        LocalScheme local_scheme(*instance_.reduced_instance(), parameters_local_scheme);
        AStarOptionalParameters<LocalScheme> parameters_a_star;
        //parameters_a_star.info.set_verbose(true);
        parameters_a_star.info.set_timelimit(parameters_.reduced_instance_time);
        parameters_a_star.thread_number_1 = 1;
        parameters_a_star.thread_number_2 = 1;
        parameters_a_star.initial_solution_ids = {1};
        std::uniform_int_distribution<Seed> d(0);
        parameters_a_star.seed = d(generator);
        auto output = a_star(local_scheme, parameters_a_star);
        const Solution& solr = output.solution_pool.best();
        for (InterventionId jr = 0; jr < (*instance_.reduced_instance()).intervention_number(); ++jr) {
            std::string name = (*instance_.reduced_instance()).intervention_name(jr);
            InterventionId j = instance_.intervention(name);
            add(solution, j, solr.intervention_starts[jr]);
        }
    }

    // Schedule the remaining interventions with a greedy strategy.
    std::shuffle(interventions_.begin(), interventions_.end(), generator);
    for (InterventionId j: interventions_) {
        if (solution.intervention_starts[j] != -1)
            continue;
        Cost c_best = cost(solution);
        Time t_start_best = -1;
        std::shuffle(times_.begin(), times_.end(), generator);
        for (Time t_start: times_) {
            if (t_start > instance_.start_max(j))
                continue;
            if (instance_.fixed(j, t_start) == 0)
                continue;
            GlobalCost gc = cost_add(solution, j, t_start, global_cost_worst());
            Cost c = cost(gc);
            if (t_start_best == -1 || c_best > c) {
                t_start_best = t_start;
                c_best = c;
            }
        }
        add(solution, j, t_start_best);
    }

    // Apply local search.
    local_search(solution, generator);

    // If the solution is infeasible, repair it with the ILP model.
    if (parameters_.repair == 1 && !feasible(solution)) {
        //std::cout << "ILP" << std::endl;
        auto res = repair(solution);
        if (res.first == -1) {
            std::cerr << "\033[31m" << "WARNING, instance is not feasible." << "\033[0m" << std::endl;
        } else if (res.first == -2) {
            std::cerr << "\033[31m" << "WARNING, no MILP solver available." << "\033[0m" << std::endl;
        } else {
            solution = res.second;
        }
    }
    return solution;
}

std::pair<int, LocalScheme::Solution> LocalScheme::repair(
        const Solution& infeasible_solution)
{
    Time h = instance_.horizon();
    InterventionId n = instance_.intervention_number();

#if GUROBI_FOUND
    {
        GRBEnv env;
        GRBModel model(env);

        // Variables

        // x[j][t] = 1 iff intervention j starts at time t
        //           0 otherwise
        std::vector<GRBVar*> x;
        for (InterventionId j = 0; j < n; ++j)
            x.push_back(model.addVars(instance_.start_max(j) + 1, GRB_BINARY));

        // Objective

        GRBLinExpr obj;
        for (InterventionId j = 0; j < instance_.intervention_number(); ++j)
            if (infeasible_solution.intervention_starts[j] != -1)
                obj += x[j][infeasible_solution.intervention_starts[j]];
        model.setObjective(obj);
        model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

        // Constraints

        // Resource constraints
        std::vector<std::vector<GRBLinExpr>> workloads(instance_.horizon(),
                std::vector<GRBLinExpr>(instance_.resource_number()));
        for (InterventionId j = 0; j < n; ++j) {
            for (ResourcePos r_pos = 0; r_pos < instance_.resource_number(j); ++r_pos) {
                ResourceId r = instance_.resource(j, r_pos);
                for (Time t_start = 0; t_start <= instance_.start_max(j); ++t_start)
                    for (Time t_cur = t_start; t_cur < t_start + instance_.duration(j, t_start); ++t_cur)
                        if (instance_.fixed(j, t_start) != 0
                                && !instance_.unconstrained(r, t_cur))
                            workloads[t_cur][r] += instance_.workload(j, r_pos, t_cur, t_start) * x[j][t_start];
            }
        }
        for (Time t_cur = 0; t_cur < h; ++t_cur) {
            for (ResourceId r = 0; r < instance_.resource_number(); ++r) {
                if (instance_.unconstrained(r, t_cur))
                    continue;
                model.addConstr(workloads[t_cur][r] <= instance_.workload_max(r, t_cur));
            }
        }

        // One alternative per item constraint
        for (InterventionId j = 0; j < n; j++) {
            GRBLinExpr expr;
            for (Time t_start = 0; t_start <= instance_.start_max(j); ++t_start)
                if (instance_.fixed(j, t_start) != 0)
                    expr += x[j][t_start];
            model.addConstr(expr == 1);
        }

        // Exclusions
        for (ExclusionId e = 0; e < instance_.exclusion_number(); ++e) {
            const Exclusion& exclusion = instance_.exclusion(e);
            InterventionId j1 = exclusion.j1;
            InterventionId j2 = exclusion.j2;
            for (Time t1_start = 0; t1_start <= instance_.start_max(j1); ++t1_start) {
                if (instance_.fixed(j1, t1_start) == 0)
                    continue;
                for (Time t2_start = 0; t2_start <= instance_.start_max(j2); ++t2_start) {
                    if (instance_.fixed(j2, t2_start) == 0)
                        continue;
                    bool in_season = false;
                    for (Time t_cur = std::max(t1_start, t2_start);
                            t_cur < std::min(
                                t1_start + instance_.duration(j1, t1_start),
                                t2_start + instance_.duration(j2, t2_start));
                            ++t_cur)
                        if (instance_.season(t_cur) == exclusion.season)
                            in_season = true;
                    if (in_season)
                        model.addConstr(x[j1][t1_start] + x[j2][t2_start] <= 1);
                }
            }
        }

        // Redirect standard output to log file
        model.set(GRB_StringParam_LogFile, "gurobi.log");
        model.set(GRB_IntParam_LogToConsole, 0);

        model.set(GRB_DoubleParam_MIPGap, 0); // Fix precision issue
        model.set(GRB_DoubleParam_NodefileStart, 0.5); // Avoid running out of memory
        model.set(GRB_IntParam_Threads, 1);

        // Optimize
        model.optimize();

        int ret;
        Solution solution = empty_solution();
        int optimstatus = model.get(GRB_IntAttr_Status);
        if (optimstatus == GRB_INFEASIBLE) {
            ret = -1;
        } else if (optimstatus == GRB_OPTIMAL) {
            ret = 0;
            for (InterventionId j = 0; j < n; ++j)
                for (Time t_start = 0; t_start <= instance_.start_max(j); ++t_start)
                    if (instance_.fixed(j, t_start) != 0)
                        if (x[j][t_start].get(GRB_DoubleAttr_X) > 0.5)
                            add(solution, j, t_start);
        } else if (model.get(GRB_IntAttr_SolCount) > 0) {
            ret = 1;
            for (InterventionId j = 0; j < n; ++j)
                for (Time t_start = 0; t_start <= instance_.start_max(j); ++t_start)
                    if (instance_.fixed(j, t_start) != 0)
                        if (x[j][t_start].get(GRB_DoubleAttr_X) > 0.5)
                            add(solution, j, t_start);
        } else {
            ret = 2;
        }

        return {ret, solution};
    }
#endif

#if CPLEX_FOUND
    {
        IloEnv env;
        IloModel model(env);

        // Variables

        // x[j][t] = 1 iff intervention j starts at time t
        //           0 otherwise
        std::vector<IloNumVarArray> x;
        for (InterventionId j = 0; j < n; ++j)
            x.push_back(IloNumVarArray(env, instance_.start_max(j) + 1, 0, 1, ILOBOOL));

        // Objective

        IloExpr obj(env);
        for (InterventionId j = 0; j < instance_.intervention_number(); ++j)
            if (infeasible_solution.intervention_starts[j] != -1)
                obj += x[j][infeasible_solution.intervention_starts[j]];
        IloObjective objective = IloMaximize(env, obj);
        model.add(objective);

        // Constraints

        // Resource constraints
        std::vector<std::vector<IloExpr>> workloads;
        for (Time t = 0; t < h; ++t) {
            workloads.push_back({});
            for (ResourceId r = 0; r < instance_.resource_number(); ++r)
                workloads[t].push_back(IloExpr(env));
        }
        for (InterventionId j = 0; j < n; ++j) {
            for (ResourcePos r_pos = 0; r_pos < instance_.resource_number(j); ++r_pos) {
                ResourceId r = instance_.resource(j, r_pos);
                for (Time t_start = 0; t_start <= instance_.start_max(j); ++t_start)
                    for (Time t_cur = t_start; t_cur < t_start + instance_.duration(j, t_start); ++t_cur)
                        if (instance_.fixed(j, t_start) != 0
                                && !instance_.unconstrained(r, t_cur))
                            workloads[t_cur][r] += instance_.workload(j, r_pos, t_cur, t_start) * x[j][t_start];
            }
        }
        for (Time t_cur = 0; t_cur < h; ++t_cur) {
            for (ResourceId r = 0; r < instance_.resource_number(); ++r) {
                if (instance_.unconstrained(r, t_cur))
                    continue;
                IloRange constraint(env, workloads[t_cur][r], instance_.workload_max(r, t_cur));
                model.add(constraint);
            }
        }

        // One alternative per item constraint
        for (InterventionId j = 0; j < n; j++) {
            IloExpr expr(env);
            for (Time t_start = 0; t_start <= instance_.start_max(j); ++t_start)
                if (instance_.fixed(j, t_start) != 0)
                    expr += x[j][t_start];
            model.add(expr == 1);
        }

        // Exclusions
        for (ExclusionId e = 0; e < instance_.exclusion_number(); ++e) {
            const Exclusion& exclusion = instance_.exclusion(e);
            InterventionId j1 = exclusion.j1;
            InterventionId j2 = exclusion.j2;
            for (Time t1_start = 0; t1_start <= instance_.start_max(j1); ++t1_start) {
                if (instance_.fixed(j1, t1_start) == 0)
                    continue;
                for (Time t2_start = 0; t2_start <= instance_.start_max(j2); ++t2_start) {
                    if (instance_.fixed(j2, t2_start) == 0)
                        continue;
                    bool in_season = false;
                    for (Time t_cur = std::max(t1_start, t2_start);
                            t_cur < std::min(
                                t1_start + instance_.duration(j1, t1_start),
                                t2_start + instance_.duration(j2, t2_start));
                            ++t_cur)
                        if (instance_.season(t_cur) == exclusion.season)
                            in_season = true;
                    if (in_season)
                        model.add(x[j1][t1_start] + x[j2][t2_start] <= 1);
                }
            }
        }

        IloCplex cplex(model);

        // Redirect standard output to log file
        std::ofstream logfile("cplex.log");
        cplex.setOut(logfile);

        cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.0); // Fix precision issue
        cplex.setParam(IloCplex::Param::MIP::Strategy::File, 2); // Avoid running out of memory
        cplex.setParam(IloCplex::Param::Threads, 1);

        // Time limit
        //if (parameters.info.timelimit != std::numeric_limits<double>::infinity())
        //    cplex.setParam(IloCplex::TiLim, parameters.info.remaining_time());

        // Optimize
        cplex.solve();

        int ret;
        Solution solution = empty_solution();
        if (cplex.getStatus() == IloAlgorithm::Infeasible) {
            ret = -1;
        } else if (cplex.getStatus() == IloAlgorithm::Optimal) {
            ret = 0;
            for (InterventionId j = 0; j < n; ++j)
                for (Time t_start = 0; t_start <= instance_.start_max(j); ++t_start)
                    if (instance_.fixed(j, t_start) != 0)
                        if (cplex.getValue(x[j][t_start]) > 0.5)
                            add(solution, j, t_start);
        } else if (cplex.isPrimalFeasible()) {
            ret = 1;
            for (InterventionId j = 0; j < n; ++j)
                for (Time t_start = 0; t_start <= instance_.start_max(j); ++t_start)
                    if (instance_.fixed(j, t_start) != 0)
                        if (cplex.getValue(x[j][t_start]) > 0.5)
                            add(solution, j, t_start);
        } else {
            ret = 2;
        }

        env.end();
        return {ret, solution};
    }
#endif

    (void)h;
    (void)n;
    (void)infeasible_solution;
    return {-2, empty_solution()};
}

/******************************** Local search *******************************/

std::vector<LocalScheme::Move> LocalScheme::perturbations(
        Solution& solution,
        std::mt19937_64&)
{
    std::vector<Move> moves;
    for (InterventionId j = 0; j < instance_.intervention_number(); ++j) {
        Time t_start_old = solution.intervention_starts[j];
        if (instance_.fixed(j, t_start_old) == 1)
            continue;
        remove(solution, j);
        for (Time t_start = 0; t_start <= instance_.start_max(j); ++t_start) {
            if (instance_.fixed(j, t_start) == 0)
                continue;
            auto c = cost_add(solution, j, t_start, global_cost_worst());
            Move move;
            move.j = j;
            move.t_start = t_start;
            move.global_cost = c;
            moves.push_back(move);
        }
        add(solution, j, t_start_old);
    }
    return moves;
}

void LocalScheme::local_search(
        Solution& solution,
        std::mt19937_64& generator,
        const Move& tabu)
{
    //if (tabu.j != -1)
    //    std::cout << "j " << tabu.j
    //        << " name " << instance_.intervention_name(tabu.j)
    //        << " t_start " << tabu.t_start
    //        << " d " << instance_.duration(tabu.j, tabu.t_start)
    //        << std::endl;
    //std::cout << to_string(global_cost(solution)) << std::endl;
    //std::cout << real_cost(solution) << std::endl;
    std::vector<Counter> neighborhoods = {0};
    Counter it;
    for (it = 0; ; ++it) {
        GlobalCost c_best = global_cost(solution);

        bool improved = false;
        // Loop through neighborhoods.
        for (Counter neighborhood: neighborhoods) {
            switch (neighborhood) {
            case 0: { // Shift neighborhood.
                std::shuffle(interventions_.begin(), interventions_.end(), generator);
                std::shuffle(times_.begin(), times_.end(), generator);
                InterventionId j_best = -1;
                Time t_start_best = -1;
                for (InterventionId j: interventions_) {
                    if (j == tabu.j)
                        continue;
                    Time t_start_old = solution.intervention_starts[j];
                    assert(t_start_old != -1);
                    if (instance_.fixed(j, t_start_old) == 1)
                        continue;
                    remove(solution, j);
                    for (Time t_start: times_) {
                        if (t_start > instance_.start_max(j))
                            continue;
                        if (t_start == t_start_old)
                            continue;
                        if (instance_.fixed(j, t_start) == 0)
                            continue;
                        GlobalCost c = cost_add(solution, j, t_start, c_best);
                        if (c >= c_best)
                            continue;
                        if (j_best != -1 && !dominates(c, c_best))
                            continue;
                        j_best = j;
                        t_start_best = t_start;
                        c_best = c;
                    }
                    add(solution, j, t_start_old);
                }
                if (j_best != -1) {
                    improved = true;
                    // Apply best move.
                    remove(solution, j_best);
                    add(solution, j_best, t_start_best);
                }
                break;
            } default: {
            }
            }
            if (improved)
                break;
        }
        if (!improved) {
            break;
        }
    }
    //std::cout << "it " << it << std::endl;
    //std::cout << to_string(global_cost(solution)) << std::endl;
    //std::cout << real_cost(solution) << std::endl;
}

