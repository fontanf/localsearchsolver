#include "examples/roadef2020.hpp"

#include <numeric>
#include <unordered_set>

using namespace localsearchsolver::roadef2020;

void LocalScheme::add(
        Solution& solution,
        InterventionId intervention_id,
        Time time_start) const
{
    assert(time_start >= 0);
    assert(time_start < instance_.horizon());
    assert(time_start <= instance_.start_max(intervention_id));
    assert(solution.intervention_starts[intervention_id] == -1);

    Time time_end = time_start + instance_.duration(intervention_id, time_start);
    ScenarioId lcm = instance_.least_common_multiple();

    // Update intervention_times_.
    solution.intervention_starts[intervention_id] = time_start;

    for (Time time_cur = time_start; time_cur < time_end; ++time_cur) {
        // Update solution.time_steps.
        for (ResourcePos resource_pos = 0;
                resource_pos < instance_.number_of_resources(intervention_id);
                ++resource_pos) {
            ResourceId resource_id = instance_.resource(intervention_id, resource_pos);
            Workload w = instance_.workload(
                    intervention_id,
                    resource_pos,
                    time_cur,
                    time_start);
            Workload w_min = instance_.workload_min(resource_id, time_cur);
            Workload w_max = instance_.workload_max(resource_id, time_cur);
            // Update overwork
            if (solution.time_steps[time_cur].workloads[resource_id] >= w_max) {
                solution.overwork += w;
            } else if (solution.time_steps[time_cur].workloads[resource_id] + w <= w_max) {
            } else {
                solution.overwork += (solution.time_steps[time_cur].workloads[resource_id] + w - w_max);
            }
            // Update underwork
            if (solution.time_steps[time_cur].workloads[resource_id] >= w_min) {
            } else if (solution.time_steps[time_cur].workloads[resource_id] + w_min <= w_min) {
                solution.underwork -= w;
            } else {
                solution.underwork -= (w_min - solution.time_steps[time_cur].workloads[resource_id]);
            }
            // Update workload
            solution.time_steps[time_cur].workloads[resource_id] += w;
            assert(solution.time_steps[time_cur].workloads[resource_id] >= 0);
        }

        // Update conflicts_
        SeasonId season = instance_.season(time_cur);
        if (season != -1) {
            for (ExclusionId e: instance_.exclusions(intervention_id, season)) {
                InterventionId intervention_id_2 = instance_.exclusion(e).j(intervention_id);
                if (solution.intervention_starts[intervention_id_2] == -1)
                    continue;
                Time t2_start = solution.intervention_starts[intervention_id_2];
                if (t2_start <= time_cur
                        && time_cur < t2_start
                        + instance_.duration(intervention_id_2, t2_start)) {
                    SolutionConflict conflict;
                    conflict.intervention_id_1 = intervention_id;
                    conflict.intervention_id_2 = intervention_id_2;
                    conflict.time_cur = time_cur;
                    solution.conflicts.push_back(conflict);
                }
            }
        }

        // Remove previous cost
        ScenarioId s_pos = std::ceil(instance_.quantile() * instance_.number_of_scenarios(time_cur)) - 1;
        solution.mean_cost -= lcm / instance_.number_of_scenarios(time_cur) * solution.time_steps[time_cur].risk_sum;
        Risk ee = lcm * solution.time_steps[time_cur].risks[solution.time_steps[time_cur].sorted_scenarios[s_pos]]
            - lcm / instance_.number_of_scenarios(time_cur) * solution.time_steps[time_cur].risk_sum;
        if (ee > 0)
            solution.expected_excess -= ee;
        assert(solution.mean_cost >= 0);
        assert(solution.expected_excess >= 0);

        // Update risks
        solution.time_steps[time_cur].risk_sum = 0;
        for (ScenarioId scenario_id = 0;
                scenario_id < instance_.number_of_scenarios(time_cur);
                ++scenario_id) {
            solution.time_steps[time_cur].risks[scenario_id]
                += instance_.risk(intervention_id, time_cur, time_start, scenario_id);
            solution.time_steps[time_cur].risk_sum
                += solution.time_steps[time_cur].risks[scenario_id];
        }
        // Update sorted scenarios
        std::nth_element(
                solution.time_steps[time_cur].sorted_scenarios.begin(),
                solution.time_steps[time_cur].sorted_scenarios.begin() + s_pos,
                solution.time_steps[time_cur].sorted_scenarios.end(),
                [&solution, time_cur](ScenarioId s1, ScenarioId s2)
                {
                    return solution.time_steps[time_cur].risks[s1]
                        < solution.time_steps[time_cur].risks[s2];
                });
        // Update cost
        solution.mean_cost += lcm / instance_.number_of_scenarios(time_cur) * solution.time_steps[time_cur].risk_sum;
        ee = lcm * solution.time_steps[time_cur].risks[solution.time_steps[time_cur].sorted_scenarios[s_pos]]
            - lcm / instance_.number_of_scenarios(time_cur) * solution.time_steps[time_cur].risk_sum;
        if (ee > 0)
            solution.expected_excess += ee;
        assert(solution.mean_cost >= 0);
        assert(solution.expected_excess >= 0);
    }
}

void LocalScheme::remove(
        Solution& solution,
        InterventionId intervention_id) const
{
    assert(solution.intervention_starts[intervention_id] != -1);
    Time time_start = solution.intervention_starts[intervention_id];
    Time time_end = time_start + instance_.duration(intervention_id, time_start);
    ScenarioId lcm = instance_.least_common_multiple();

    // Update intervention_times_.
    solution.intervention_starts[intervention_id] = -1;

    // Update solution.conflicts
    for (auto it = solution.conflicts.begin(); it != solution.conflicts.end();) {
        if (it->intervention_id_1 != intervention_id
                && it->intervention_id_2 != intervention_id) {
            ++it;
        } else {
            *it = solution.conflicts.back();
            solution.conflicts.pop_back();
        }
    }

    // Update solution.time_steps.
    for (Time time_cur = time_start; time_cur < time_end; ++time_cur) {
        // Update solution.time_steps.
        for (ResourcePos resource_pos = 0;
                resource_pos < instance_.number_of_resources(intervention_id);
                ++resource_pos) {
            ResourceId resource_id = instance_.resource(intervention_id, resource_pos);
            Workload w = instance_.workload(
                    intervention_id,
                    resource_pos,
                    time_cur,
                    time_start);
            Workload w_min = instance_.workload_min(resource_id, time_cur);
            Workload w_max = instance_.workload_max(resource_id, time_cur);
            // Update overwork
            if (solution.time_steps[time_cur].workloads[resource_id] - w >= w_max) {
                solution.overwork -= w;
            } else if (solution.time_steps[time_cur].workloads[resource_id] <= w_max) {
            } else {
                solution.overwork -= (solution.time_steps[time_cur].workloads[resource_id] - w_max);
            }
            // Update underwork
            if (solution.time_steps[time_cur].workloads[resource_id] - w >= w_min) {
            } else if (solution.time_steps[time_cur].workloads[resource_id] <= w_min) {
                solution.underwork += w;
            } else {
                solution.underwork += (solution.time_steps[time_cur].workloads[resource_id] - w_min);
            }
            // Update workload
            solution.time_steps[time_cur].workloads[resource_id] -= w;
        }

        // Remove previous cost
        ScenarioId s_pos = std::ceil(instance_.quantile() * instance_.number_of_scenarios(time_cur)) - 1;
        solution.mean_cost -= lcm / instance_.number_of_scenarios(time_cur) * solution.time_steps[time_cur].risk_sum;
        Risk ee = lcm * solution.time_steps[time_cur].risks[solution.time_steps[time_cur].sorted_scenarios[s_pos]]
            - lcm / instance_.number_of_scenarios(time_cur) * solution.time_steps[time_cur].risk_sum;
        if (ee > 0)
            solution.expected_excess -= ee;
        assert(solution.mean_cost >= 0);
        assert(solution.expected_excess >= 0);

        // Update risks
        solution.time_steps[time_cur].risk_sum = 0;
        for (ScenarioId scenario_id = 0;
                scenario_id < instance_.number_of_scenarios(time_cur);
                ++scenario_id) {
            solution.time_steps[time_cur].risks[scenario_id]
                -= instance_.risk(intervention_id, time_cur, time_start, scenario_id);
            solution.time_steps[time_cur].risk_sum
                += solution.time_steps[time_cur].risks[scenario_id];
        }
        // Update sorted scenarios
        std::nth_element(
                solution.time_steps[time_cur].sorted_scenarios.begin(),
                solution.time_steps[time_cur].sorted_scenarios.begin() + s_pos,
                solution.time_steps[time_cur].sorted_scenarios.end(),
                [&solution, time_cur](ScenarioId s1, ScenarioId s2)
                {
                    return solution.time_steps[time_cur].risks[s1]
                        < solution.time_steps[time_cur].risks[s2];
                });
        // Update cost
        solution.mean_cost += lcm / instance_.number_of_scenarios(time_cur) * solution.time_steps[time_cur].risk_sum;
        ee = lcm * solution.time_steps[time_cur].risks[solution.time_steps[time_cur].sorted_scenarios[s_pos]]
            - lcm / instance_.number_of_scenarios(time_cur) * solution.time_steps[time_cur].risk_sum;
        if (ee > 0)
            solution.expected_excess += ee;
        assert(solution.mean_cost >= 0);
        assert(solution.expected_excess >= 0);
    }
}

LocalScheme::GlobalCost LocalScheme::cost_add(
        const Solution& solution,
        InterventionId intervention_id,
        Time time_start,
        const GlobalCost& cutoff)
{
    assert(time_start >= 0);
    assert(time_start < instance_.horizon());
    assert(time_start <= instance_.start_max(intervention_id));
    assert(solution.intervention_starts[intervention_id] == -1);

    Time time_end = time_start + instance_.duration(intervention_id, time_start);
    GlobalCost c = global_cost(solution);

    // Check disjonctive constraints.
    for (Time time_cur = time_start; time_cur < time_end; ++time_cur) {
        SeasonId season = instance_.season(time_cur);
        if (season != -1) {
            for (ExclusionId e: instance_.exclusions(intervention_id, season)) {
                InterventionId intervention_id_2 = instance_.exclusion(e).j(intervention_id);
                if (solution.intervention_starts[intervention_id_2] == -1)
                    continue;
                Time t2_start = solution.intervention_starts[intervention_id_2];
                if (t2_start <= time_cur
                        && time_cur < t2_start
                        + instance_.duration(intervention_id_2, t2_start)) {
                    number_of_conflicts(c)++;
                    if (c >= cutoff)
                        return cutoff;
                }
            }
        }
    }

    // Check resource constraint.
    for (Time time_cur = time_start; time_cur < time_end; ++time_cur) {
        for (ResourcePos resource_pos = 0;
                resource_pos < instance_.number_of_resources(intervention_id);
                ++resource_pos) {
            Workload w = instance_.workload(
                    intervention_id,
                    resource_pos,
                    time_cur,
                    time_start);
            if (w == 0)
                continue;
            ResourceId resource_id = instance_.resource(intervention_id, resource_pos);
            Workload w_cur = solution.time_steps[time_cur].workloads[resource_id];
            Workload w_max = instance_.workload_max(resource_id, time_cur);
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
    for (Time time_cur = time_start; time_cur < time_end; ++time_cur) {
        ScenarioId st = instance_.number_of_scenarios(time_cur);
        ScenarioId s_pos = std::ceil(instance_.quantile() * st) - 1;

        // Remove previous cost
        mean_cost -= lcm / st * solution.time_steps[time_cur].risk_sum;
        Risk ee = lcm * solution.time_steps[time_cur].risks[solution.time_steps[time_cur].sorted_scenarios[s_pos]]
            - lcm / st * solution.time_steps[time_cur].risk_sum;
        if (ee > 0)
            expected_excess -= ee;

        // Update risks
        Risk risk_sum = 0;
        for (ScenarioId scenario_id = 0; scenario_id < st; ++scenario_id) {
            sorted_scenarios_[scenario_id] = solution.time_steps[time_cur].sorted_scenarios[scenario_id];
            risks_[scenario_id] = solution.time_steps[time_cur].risks[scenario_id]
                + instance_.risk(intervention_id, time_cur, time_start, scenario_id);
            risk_sum += risks_[scenario_id];
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

LocalScheme::Solution LocalScheme::empty_solution() const
{
    Solution solution;
    solution.intervention_starts.resize(instance_.number_of_interventions(), -1);
    solution.time_steps.resize(instance_.horizon());
    for (Time time = 0; time < instance_.horizon(); ++time) {
        solution.time_steps[time].workloads.resize(instance_.number_of_resources(), 0);
        solution.time_steps[time].risks.resize(instance_.number_of_scenarios(time), 0);
        solution.time_steps[time].sorted_scenarios.resize(instance_.number_of_scenarios(time));
        for (ScenarioId scenario_id = 0; scenario_id < instance_.number_of_scenarios(time); ++scenario_id)
            solution.time_steps[time].sorted_scenarios[scenario_id] = scenario_id;
        // Initialize underwork
        for (ResourceId resource_id = 0;
                resource_id < instance_.number_of_resources();
                ++resource_id) {
            solution.underwork += instance_.workload_min(resource_id, time);
        }
    }
    return solution;
}

LocalScheme::Solution LocalScheme::initial_solution(
        Counter,
        std::mt19937_64& generator)
{
    Solution solution = empty_solution();

    for (InterventionId intervention_id = 0;
            intervention_id < instance_.number_of_interventions();
            ++intervention_id) {
        std::vector<Time> starts;
        for (Time time_start = 0;
                time_start <= instance_.start_max(intervention_id);
                ++time_start) {
            if (instance_.fixed(intervention_id, time_start) == 0)
                continue;
            starts.push_back(time_start);
        }
        std::shuffle(starts.begin(), starts.end(), generator);
        add(solution, intervention_id, starts.front());
    }
    return solution;
}

std::vector<LocalScheme::Perturbation> LocalScheme::perturbations(
        Solution& solution,
        std::mt19937_64&)
{
    std::vector<Perturbation> perturbations;
    for (InterventionId intervention_id = 0;
            intervention_id < instance_.number_of_interventions();
            ++intervention_id) {
        Time time_start_old = solution.intervention_starts[intervention_id];
        if (instance_.fixed(intervention_id, time_start_old) == 1)
            continue;
        remove(solution, intervention_id);
        for (Time time_start = 0; time_start <= instance_.start_max(intervention_id); ++time_start) {
            if (instance_.fixed(intervention_id, time_start) == 0)
                continue;
            auto c = cost_add(solution, intervention_id, time_start, worst<GlobalCost>());
            Perturbation perturbation;
            perturbation.intervention_id = intervention_id;
            perturbation.time_start = time_start;
            perturbation.global_cost = c;
            number_of_conflicts(perturbation.global_cost) = number_of_conflicts(global_cost(solution));
            overwork(perturbation.global_cost) = overwork(global_cost(solution));
            underwork(perturbation.global_cost) = underwork(global_cost(solution));
            perturbations.push_back(perturbation);
        }
        add(solution, intervention_id, time_start_old);
    }
    return perturbations;
}

void LocalScheme::local_search(
        Solution& solution,
        std::mt19937_64& generator,
        const Perturbation& tabu)
{
    //if (tabu.j != -1)
    //    std::cout << "j " << tabu.j
    //        << " name " << instance_.intervention_name(tabu.j)
    //        << " time_start " << tabu.time_start
    //        << " d " << instance_.duration(tabu.j, tabu.time_start)
    //        << std::endl;
    //std::cout << to_string(global_cost(solution)) << std::endl;
    //std::cout << real_cost(solution) << std::endl;
    std::vector<Counter> neighborhoods = {0};
    Counter it;
    (void)it;
    for (it = 0; ; ++it) {
        GlobalCost c_best = global_cost(solution);

        bool improved = false;
        // Loop through neighborhoods.
        for (Counter neighborhood: neighborhoods) {
            switch (neighborhood) {
            case 0: { // Shift neighborhood.
                std::shuffle(interventions_.begin(), interventions_.end(), generator);
                std::shuffle(times_.begin(), times_.end(), generator);
                InterventionId intervention_id_best = -1;
                Time time_start_best = -1;
                for (InterventionId intervention_id: interventions_) {
                    if (intervention_id == tabu.intervention_id)
                        continue;
                    Time time_start_old = solution.intervention_starts[intervention_id];
                    assert(time_start_old != -1);
                    if (instance_.fixed(intervention_id, time_start_old) == 1)
                        continue;
                    remove(solution, intervention_id);
                    for (Time time_start: times_) {
                        if (time_start > instance_.start_max(intervention_id))
                            continue;
                        if (time_start == time_start_old)
                            continue;
                        if (instance_.fixed(intervention_id, time_start) == 0)
                            continue;
                        GlobalCost c = cost_add(
                                solution,
                                intervention_id,
                                time_start,
                                c_best);
                        if (c >= c_best)
                            continue;
                        if (intervention_id_best != -1 && !dominates(c, c_best))
                            continue;
                        intervention_id_best = intervention_id;
                        time_start_best = time_start;
                        c_best = c;
                    }
                    add(solution, intervention_id, time_start_old);
                }
                if (intervention_id_best != -1) {
                    improved = true;
                    // Apply best move.
                    remove(solution, intervention_id_best);
                    add(solution, intervention_id_best, time_start_best);
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

