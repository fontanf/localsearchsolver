#include "examples/roadef2020.hpp"

#include "localsearchsolver/a_star.hpp"

#include "simdjson.h"

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
using namespace simdjson;

/********************************** Instance *********************************/

Risk gcd(Risk a, Risk b)
{
    for (;;)
    {
        if (a == 0)
            return b;
        b %= a;
        if (b == 0)
            return a;
        a %= b;
    }
}

Risk lcm(Risk a, Risk b)
{
    Risk temp = gcd(a, b);
    return temp? (a / temp * b): 0;
}

Instance::Instance(
        std::string filepath,
        std::string)
{
    std::ifstream f(filepath);
    if (!f.good())
        std::cerr << "\033[31m" << "ERROR, unable to open file \"" << filepath << "\"" << "\033[0m" << std::endl;

    // Read json file.
    ondemand::parser parser;
    auto json = padded_string::load(filepath);
    ondemand::document doc = parser.iterate(json); // position a pointer at the beginning of the JSON data

    std::cout << "Read instance..." << std::endl;
    for (auto field: doc.get_object()) {
        auto key = std::string_view(field.unescaped_key());
        if (key == "Resources") {
            std::cout << "* Read resources..." << std::endl;
            ResourceId resource_id = 0;
            for (auto resource_field: field.value().get_object()) {
                auto resource_name = std::string_view(resource_field.unescaped_key());
                string_to_resource_[resource_name.to_string()] = resource_id;
                resources_.push_back({});
                resources_.back().name = resource_name.to_string();
                Resource& resource = resources_.back();
                for (auto m: resource_field.value().get_object()) {
                    auto m_str = std::string_view(m.unescaped_key());
                    if (m_str == "min") {
                        for (double val: m.value().get_array()) {
                            resource.min_double.push_back(val);
                            while (fabs(val * resource.multiplier
                                        - std::round(val * resource.multiplier)) > 0.000001)
                                resource.multiplier *= 10;
                        }
                    } else if (m_str == "max") {
                        for (double val: m.value().get_array()) {
                            resource.max_double.push_back(val);
                            while (fabs(val * resource.multiplier
                                        - std::round(val * resource.multiplier)) > 0.000001)
                                resource.multiplier *= 10;
                        }
                    }
                }
                resource_id++;
            }
        } else if (key == "Seasons") {
            std::cout << "* Read seasons..." << std::endl;
            SeasonId season_id = 0;
            for (auto season: field.value().get_object()) {
                auto season_name = std::string_view(season.unescaped_key());
                season_strings_.push_back(season_name.to_string());
                string_to_season_[season_name.to_string()] = season_id;
                for (auto val: season.value().get_array()) {
                    Time t = std::stol(std::string_view(val).to_string()) - 1;
                    if ((Time)seasons_.size() <= t)
                        seasons_.resize(t + 1, -1);
                    seasons_[t] = season_id;
                }
                season_id++;
            }
        } else if (key == "Interventions") {
            std::cout << "* Read interventions..." << std::endl;
            InterventionId intervention_id = 0;
            for (auto intervention_field: field.value().get_object()) {
                auto intervention_name = std::string_view(intervention_field.unescaped_key());
                string_to_interv_[intervention_name.to_string()] = intervention_id;
                interventions_.push_back({});
                Intervention& intervention = interventions_.back();
                intervention.name = intervention_name.to_string();
                intervention.exclusions = std::vector<std::vector<ExclusionId>>(season_number());
                for (auto m: intervention_field.value().get_object()) {
                    auto m_str = std::string_view(m.unescaped_key());
                    if (m_str == "tmax") {
                        intervention.t_start_max
                            = std::stol(std::string_view(m.value()).to_string()) - 1;
                    } else if (m_str == "Delta") {
                        for (double val: m.value().get_array())
                            intervention.deltas.push_back(val);
                    } else if (m_str == "workload") {
                        ResourcePos resource_pos = 0;
                        for (auto resource_field: m.value().get_object()) {
                            intervention.resources.push_back({});
                            InterventionResource& resource = intervention.resources.back();
                            auto resource_name = std::string_view(resource_field.unescaped_key());
                            ResourceId resource_id = string_to_resource_[resource_name.to_string()];
                            resource.r = resource_id;
                            Time t_pos = 0;
                            for (auto t_cur_field: resource_field.value().get_object()) {
                                Time t_cur = std::stol(std::string_view(t_cur_field.unescaped_key()).to_string()) - 1;
                                while ((Time)resource.workloads.size() <= t_cur) {
                                    resource.workloads.push_back({});
                                    resource.workloads.back().t_first = resource.workloads.size();
                                }
                                InterventionResourceWorkload& iworkload = resource.workloads[t_cur];
                                for (auto t_start_field: t_cur_field.value().get_object()) {
                                    Time t_start = std::stol(std::string_view(t_start_field.unescaped_key()).to_string()) - 1;
                                    if (iworkload.t_first > t_start)
                                        iworkload.t_first = t_start;
                                    Time t_pos = t_start - iworkload.t_first;
                                    if ((Time)iworkload.workloads_double.size() <= t_pos)
                                        iworkload.workloads_double.resize(t_pos + 1, 0);
                                    iworkload.workloads_double[t_pos] = t_start_field.value();
                                    while (fabs((double)t_start_field.value() * resources_[resource_id].multiplier
                                                - std::round((double)t_start_field.value() * resources_[resource_id].multiplier)) > 0.000001)
                                        resources_[resource_id].multiplier *= 10;
                                }
                                t_pos++;
                            }
                            resource_pos++;
                        }
                        std::sort(intervention.resources.begin(), intervention.resources.end(),
                                [](const InterventionResource& r1, const InterventionResource& r2) {
                                    return r1.r < r2.r;
                                });
                    } else if (m_str == "risk") {
                        for (auto t_cur_field: m.value().get_object()) {
                            Time t_cur = std::stol(std::string_view(t_cur_field.unescaped_key()).to_string()) - 1;
                            while ((Time)intervention.risks.size() <= t_cur) {
                                intervention.risks.push_back({});
                                intervention.risks.back().t_first = intervention.risks.size();
                            }
                            InterventionRisk& irisk = intervention.risks[t_cur];
                            for (auto t_start_field: t_cur_field.value().get_object()) {
                                Time t_start = std::stol(std::string_view(t_start_field.unescaped_key()).to_string()) - 1;
                                if (irisk.t_first > t_start)
                                    irisk.t_first = t_start;
                                Time t_pos = t_start - irisk.t_first;
                                if ((Time)irisk.risks_double.size() <= t_pos)
                                    irisk.risks_double.resize(t_pos + 1);
                                for (double risk: t_start_field.value().get_array()) {
                                    irisk.risks_double[t_pos].push_back(risk);
                                    while (fabs((double)risk * risk_multiplier_
                                                - std::round((double)risk * risk_multiplier_)) > 0.000001)
                                        risk_multiplier_ *= 10;
                                }
                            }
                        }
                    }
                }
                intervention_id++;
            }
        } else if (key == "Exclusions") {
            std::cout << "* Read exclusions..." << std::endl;
            ExclusionId exclusion_id = 0;
            for (auto exclusion_field: field.value().get_object()) {
                exclusions_.push_back({});
                Exclusion& exclusion = exclusions_.back();
                InterventionId id = 0;
                for (auto val: exclusion_field.value().get_array()) {
                    if (id == 0) {
                        std::string intervention_name = std::string_view(val).to_string();
                        exclusion.j1 = string_to_interv_[intervention_name];
                    } else if (id == 1) {
                        std::string intervention_name = std::string_view(val).to_string();
                        exclusion.j2 = string_to_interv_[intervention_name];
                    } else if (id == 2) {
                        std::string season_name = std::string_view(val).to_string();
                        exclusion.season = string_to_season_[season_name];
                    } else {

                    }
                    id++;
                }
                interventions_[exclusion.j1].exclusions[exclusion.season].push_back(exclusion_id);
                interventions_[exclusion.j1].exclusion_number++;
                interventions_[exclusion.j2].exclusions[exclusion.season].push_back(exclusion_id);
                interventions_[exclusion.j2].exclusion_number++;
                exclusion_id++;
            }
        } else if (key == "T") {
            std::cout << "* Read horizon..." << std::endl;
            horizon_ = (int64_t)field.value();
        } else if (key == "Scenarios_number") {
            std::cout << "* Read scenarios..." << std::endl;
            for (ScenarioId scenario_number: field.value().get_array()) {
                scenario_numbers_.push_back(scenario_number);
                least_common_multiple_ = lcm(least_common_multiple_, scenario_number);
                if (least_common_multiple_ >= 10000000)
                    least_common_multiple_ =  10000000;
                assert(least_common_multiple_ >= 0);
            }
        } else if (key == "Alpha") {
            std::cout << "* Read alpha..." << std::endl;
            alpha_ = field.value();
            while (alpha_ * alpha_multiplier_
                    - std::floor(alpha_ * alpha_multiplier_) != 0) {
                alpha_multiplier_ *= 10;
            }
        } else if (key == "Quantile") {
            std::cout << "* Read quantile..." << std::endl;
            quantile_ = field.value();
        } else if (key == "ComputationTime") {
        }
    }

    // Compute integer values.
    std::cout << "Convert instance..." << std::endl;

    std::cout << "* Convert alpha..." << std::endl;
    alpha_1_ = alpha_ * alpha_multiplier_;
    alpha_2_ = alpha_multiplier_ - alpha_1_;

    // Convert resources.
    std::cout << "* Convert resources..." << std::endl;
    for (ResourceId r = 0; r < resource_number(); ++r) {
        Resource& resource = resources_[r];
        resource.multiplier = std::min((ResourceId)1000000, resources_[r].multiplier);
        for (double val: resource.min_double)
            resource.min.push_back(std::round(val * resources_[r].multiplier));
        std::vector<double> min_tmp;
        resource.min_double.swap(min_tmp);
        for (double val: resource.max_double)
            resource.max.push_back(std::round(val * resources_[r].multiplier));
        std::vector<double> max_tmp;
        resource.min_double.swap(max_tmp);
    }

    // Convert risks.
    std::cout << "* Convert risks..." << std::endl;
    risk_multiplier_ = std::min((ResourceId)10000, risk_multiplier_);
    for (InterventionId j = 0; j < intervention_number(); ++j) {
        for (Time t_cur = 0; t_cur < (Time)interventions_[j].risks.size(); ++t_cur) {
            InterventionRisk& irisk = interventions_[j].risks[t_cur];
            for (const auto& vec: irisk.risks_double) {
                irisk.risks.push_back({});
                for (double val: vec)
                    irisk.risks.back().push_back(std::round(val * risk_multiplier_));
            }
            std::vector<std::vector<double>> tmp;
            irisk.risks_double.swap(tmp);
        }
    }
    //for (InterventionId j = 0; j < intervention_number(); ++j) {
    //    for (Time t_start = 0; t_start <= start_max(j); ++t_start) {
    //        Time t_end = t_start + duration(j, t_start);
    //        for (Time t_cur = t_start; t_cur < t_end; ++t_cur) {
    //            for (ScenarioId s = 0; s < scenario_number(t_cur); ++s) {
    //                std::cout << "j " << j
    //                    << " t_cur " << t_cur
    //                    << " t_start " << t_start
    //                    << " s " << s
    //                    << " d " << risk_double(j, t_cur, t_start, s)
    //                    << " m " << risk_multiplier_
    //                    << " i " << risk_int(j, t_cur, t_start, s)
    //                    << std::endl;
    //            }
    //        }
    //    }
    //}

    // Convert workloads.
    std::cout << "* Convert workloads..." << std::endl;
    for (InterventionId j = 0; j < intervention_number(); ++j) {
        for (ResourcePos r_pos = 0; r_pos < resource_number(j); ++r_pos) {
            ResourceId r = resource(j, r_pos);
            for (Time t_cur = 0; t_cur < (Time)interventions_[j].resources[r_pos].workloads.size(); ++t_cur) {
                InterventionResourceWorkload& iworkload
                    = interventions_[j].resources[r_pos].workloads[t_cur];
                for (double w: iworkload.workloads_double)
                    iworkload.workloads.push_back(
                            std::round(w * resources_[r].multiplier));
                std::vector<double> w_tmp;
                iworkload.workloads_double.swap(w_tmp);
            }
        }
    }
    //for (InterventionId j = 0; j < intervention_number(); ++j) {
    //    std::cout << "j " << j << " name " << intervention_name(j) << std::endl;
    //    for (ResourcePos r_pos = 0; r_pos < resource_number(j); ++r_pos) {
    //        ResourceId r = resource(j, r_pos);
    //        std::cout << "r_pos " << r_pos << " r " << r
    //            << " name " << resource_name(r)
    //            << std::endl;
    //        for (Time t_start = 0; t_start < start_max(j); ++t_start) {
    //            Time t_end = t_start + duration(j, t_start);
    //            std::cout << "t_start " << t_start << " t_end " << t_end << std::endl;
    //            for (Time t_cur = t_start; t_cur < t_end; ++t_cur) {
    //                Workload w = workload(j, r_pos, t_cur, t_start);
    //                std::cout << "t_cur " << t_cur
    //                    << " w " << w << std::endl;
    //            }
    //        }
    //    }
    //}

    fix_assignments();
}

void Instance::fix_assignments()
{
    std::cout << "Fix assignments..." << std::endl;
    Counter fixed_assignment_number = 0;
    Counter assignment_number = 0;
    for (InterventionId j = 0; j < intervention_number(); ++j) {
        interventions_[j].fixed_assignments = std::vector<int8_t>(start_max(j) + 1, -1);
        assignment_number += start_max(j) + 1;
    }
    Counter total_resource_number = 0;
    Counter unconstrained_resource_number = 0;
    for (ResourceId r = 0; r < resource_number(); ++r) {
        resources_[r].unconstrained = std::vector<bool>(horizon(), false);
        total_resource_number += horizon();
    }

    // Dominated starts.
    for (InterventionId j = 0; j < intervention_number(); ++j) {
        for (Time t_start = 0; t_start < start_max(j); ++t_start) {
            if (fixed(j, t_start) != -1)
                continue;
            Time t_end = t_start + duration(j, t_start);
            Time t_end_2 = t_start + 1 + duration(j, t_start + 1);
            if (t_end_2 > t_end)
                continue;
            bool dominated = true;
            for (Time t_cur = t_start + 1; t_cur < t_end; ++t_cur) {
                for (ResourcePos r_pos = 0; r_pos < resource_number(j); ++r_pos) {
                    Workload w = workload(j, r_pos, t_cur, t_start);
                    Workload w2 = workload(j, r_pos, t_cur, t_start + 1);
                    if (w < w2) {
                        dominated = false;
                        break;
                    }
                }
                if (!dominated)
                    break;
                for (ScenarioId s = 0; s < scenario_number(t_cur); ++s) {
                    Risk r = risk(j, t_cur, t_start, s);
                    Risk r2 = risk(j, t_cur, t_start + 1, s);
                    if (r < r2) {
                        dominated = false;
                        break;
                    }
                }
                if (!dominated)
                    break;
            }
            if (dominated) {
                interventions_[j].fixed_assignments[t_start] = 0;
                fixed_assignment_number++;
            }
        }
    }
    std::cout << "* Dominated starts:                "
        << fixed_assignment_number << " / " << assignment_number << std::endl;

    for (;;) {
        bool new_fixed = false;

        // Conflicts.
        for (InterventionId j = 0; j < intervention_number(); ++j) {
            std::vector<bool> mandatory(horizon(), true);
            for (Time t_start = 0; t_start <= start_max(j); ++t_start) {
                if (fixed(j, t_start) == 0)
                    continue;
                for (Time t = 0; t < t_start; ++t)
                    mandatory[t] = false;
                for (Time t = t_start + duration(j, t_start); t < horizon(); ++t)
                    mandatory[t] = false;
            }
            for (Time t_cur = 0; t_cur < horizon(); ++t_cur) {
                if (!mandatory[t_cur])
                    continue;
                SeasonId s = seasons_[t_cur];
                for (ExclusionId e: exclusions(j, s)) {
                    InterventionId j2 = exclusion(e).j(j);
                    for (Time t_start_2 = 0; t_start_2 <= start_max(j2); ++t_start_2) {
                        if (interventions_[j2].fixed_assignments[t_start_2] != -1)
                            continue;
                        if (t_start_2 <= t_cur
                                && t_cur < t_start_2 + duration(j2, t_start_2)) {
                            interventions_[j2].fixed_assignments[t_start_2] = 0;
                            fixed_assignment_number++;
                            new_fixed = true;
                        }
                    }
                }
            }
        }
        std::cout << "* Conflicts:                       "
            << fixed_assignment_number << " / " << assignment_number << std::endl;

        if (!new_fixed)
            break;
    }

    // Single starts.
    for (InterventionId j = 0; j < intervention_number(); ++j) {
        std::vector<Time> starts;
        for (Time t_start = 0; t_start <= start_max(j); ++t_start)
            if (fixed(j, t_start) != 0)
                starts.push_back(t_start);
        if (starts.size() == 1) {
            interventions_[j].fixed_assignments[starts.front()] = 1;
            fixed_assignment_number++;
        }
    }
    std::cout << "* Single starts:                   "
        << fixed_assignment_number << " / " << assignment_number << std::endl;

    // Unconstrained resources.
    std::vector<std::vector<Workload>> workloads_max(
            horizon(),
            std::vector<Workload>(resource_number(), 0));
    for (InterventionId j = 0; j < intervention_number(); ++j) {
        std::vector<std::vector<Workload>> workloads_max_j(
                horizon(),
                std::vector<Workload>(resource_number(), 0));
        for (Time t_start = 0; t_start <= start_max(j); ++t_start) {
            if (fixed(j, t_start) == 0)
                continue;
            Time t_end = t_start + duration(j, t_start);
            for (Time t_cur = t_start; t_cur < t_end; ++t_cur) {
                for (ResourcePos r_pos = 0; r_pos < resource_number(j); ++r_pos) {
                    ResourceId r = resource(j, r_pos);
                    Workload w = workload(j, r_pos, t_cur, t_start);
                    if (workloads_max_j[t_cur][r] < w)
                        workloads_max_j[t_cur][r] = w;
                }
            }
        }
        for (Time t_cur = 0; t_cur < horizon(); ++t_cur) {
            for (ResourcePos r_pos = 0; r_pos < resource_number(j); ++r_pos) {
                ResourceId r = resource(j, r_pos);
                workloads_max[t_cur][r] += workloads_max_j[t_cur][r];
            }
        }
    }
    for (Time t_cur = 0; t_cur < horizon(); ++t_cur) {
        for (ResourceId r = 0; r < resource_number(); ++r) {
            if (resources_[r].unconstrained[t_cur])
                continue;
            if (workloads_max[t_cur][r] <= workload_max(r, t_cur)) {
                resources_[r].unconstrained[t_cur] = true;
                unconstrained_resource_number++;
                //std::cout << "t " << t_cur << " r " << r
                //    << " w " << workloads_max[t_cur][r] << " / " << workload_max(r, t_cur)
                //    << std::endl;
            }
        }
    }
    std::cout << "* Unconstrained resources:         "
        << unconstrained_resource_number << " / " << total_resource_number << std::endl;

    for (InterventionId j = 0; j < intervention_number(); ++j) {

        interventions_[j].workload_mean.resize(resource_number(), 0);
        for (ResourcePos r_pos = 0; r_pos < resource_number(j); ++r_pos) {
            ResourceId r = resource(j, r_pos);
            Time start_number = 0;
            for (Time t_start = 0; t_start <= start_max(j); ++t_start) {
                double r_max = 0;
                Time t_end = t_start + duration(j, t_start);
                for (Time t_cur = t_start; t_cur < t_end; ++t_cur) {
                    Workload w = workload(j, r_pos, t_cur, t_start);
                    Workload w_max = workload_max(r, t_cur);
                    double r = (double)w / w_max;
                    if (r_max < r)
                        r_max = r;
                }
                interventions_[j].workload_mean[r] += r_max;
                start_number++;
            }
            interventions_[j].workload_mean[r] /= start_number;
        }

        Time d_sum = 0;
        Time start_number = 0;
        for (Time t_start = 0; t_start <= start_max(j); ++t_start) {
            if (fixed(j, t_start) == 0)
                continue;
            d_sum += duration(j, t_start);
            start_number++;
        }
        interventions_[j].duration_mean = (double)d_sum / start_number;
    }

    // Reduced instance.
    std::vector<InterventionId> interventions(intervention_number());
    std::iota(interventions.begin(), interventions.end(), 0);
    std::sort(interventions.begin(), interventions.end(),
            [this](InterventionId j1, InterventionId j2) {
                    return duration_mean(j1) < duration_mean(j2);
            });
    std::vector<bool> selected_interventions(intervention_number(), true);
    InterventionId n = intervention_number();
    Time d_max = 0;
    for (InterventionId j: interventions) {
        Time d = duration_mean(j);
        if (n < (double)intervention_number() / 2
                && d > d_max)
            break;
        if (d_max < ceil(d - 0.5) + 0.5) {
            //std::cout << "d " << d_max << " n " << n << std::endl;
            d_max = ceil(d - 0.5) + 0.5;
        }
        selected_interventions[j] = false;
        n--;
    }
    //std::cout << "n " << n << std::endl;
    for (InterventionId j: interventions) {
        if (selected_interventions[j])
            continue;
        if (exclusion_number(j) > 0) {
            selected_interventions[j] = true;
            n++;
        }
    }
    if (n < intervention_number()) {
        reduced_instance_ = std::unique_ptr<Instance>(new Instance(
                    reduced_instance(selected_interventions)));
        //std::cout << *reduced_instance << std::endl;
    }
}

Instance Instance::reduced_instance(const std::vector<bool>& interventions) const
{
    Instance instance_new;
    instance_new.alpha_ = alpha_;
    instance_new.quantile_ = quantile_;
    instance_new.horizon_ = horizon_;
    instance_new.resources_ = resources_;
    instance_new.scenario_numbers_ = scenario_numbers_;
    instance_new.seasons_ = seasons_;
    instance_new.string_to_resource_ = string_to_resource_;
    instance_new.string_to_season_ = string_to_season_;
    instance_new.season_strings_ = season_strings_;
    instance_new.least_common_multiple_ = least_common_multiple_;
    instance_new.alpha_1_ = alpha_1_;
    instance_new.alpha_2_ = alpha_2_;
    instance_new.alpha_multiplier_ = alpha_multiplier_;
    instance_new.risk_multiplier_ = risk_multiplier_ ;
    for (InterventionId j = 0; j < intervention_number(); ++j) {
        if (!interventions[j])
            continue;
        instance_new.string_to_interv_[interventions_[j].name] = instance_new.interventions_.size();
        instance_new.interventions_.push_back(interventions_[j]);
        for (auto& e: instance_new.interventions_.back().exclusions)
            e.clear();
        instance_new.interventions_.back().exclusion_number = 0;
    }
    for (const Exclusion exclusion: exclusions_) {
        if (!interventions[exclusion.j1])
            continue;
        if (!interventions[exclusion.j2])
            continue;
        ExclusionId exclusion_id = instance_new.exclusions_.size();
        Exclusion e;
        e.j1 = instance_new.string_to_interv_[intervention_name(exclusion.j1)];
        e.j2 = instance_new.string_to_interv_[intervention_name(exclusion.j2)];
        e.season = exclusion.season;
        instance_new.exclusions_.push_back(e);
        instance_new.interventions_[e.j1].exclusions[e.season].push_back(exclusion_id);
        instance_new.interventions_[e.j1].exclusion_number++;
        instance_new.interventions_[e.j2].exclusions[e.season].push_back(exclusion_id);
        instance_new.interventions_[e.j2].exclusion_number++;
    }
    return instance_new;
}

std::ostream& localsearchsolver::roadef2020::operator<<(
        std::ostream& os,
        const Instance& instance)
{
    os << "Instance:" << std::endl;
    os << "* Intervention number:             " << instance.intervention_number() << std::endl;
    os << "* Horizon:                         " << instance.horizon() << std::endl;
    os << "* Resource number:                 " << instance.resource_number() << std::endl;
    os << "* Season number:                   " << instance.season_number() << std::endl;
    os << "* Exclusion number:                " << instance.exclusion_number() << std::endl;
    os << "* Alpha:                           " << instance.alpha() << std::endl;
    os << "* Quantile:                        " << instance.quantile() << std::endl;
    return os;
}

/******************************** LocalScheme ********************************/

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

