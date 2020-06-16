#include "diversity_score_subset.h"
 
#include "../option_parser.h"
#include "../plugin.h"

#include "../algorithms/max_cliques.h"
//#include "../lp/lp_solver.h"

#include <list>
#include <math.h>       /* fabs */

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

//namespace diversity_score {

DiversityScoreSubset::DiversityScoreSubset(const Options &opts) : DiversityScore(opts),
        plans_subset_size(opts.get<int>("plans_subset_size")),
        exact_method(opts.get<bool>("exact_method")),
        dump_plans(opts.get<bool>("dump_plans"))
{

}


void DiversityScoreSubset::compute_metrics() {
    if (plans_subset_size > (int) _plans.size()) {
        cerr << "Number of plans in the subset cannot be greater than the set size" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }
    if (plans_subset_size == (int) _plans.size()) {
        cout << "Computing metric score for the entire set of plans" << endl;
        compute_metrics_exact_set();
        if (dump_plans) {
            print_all_plans();
        }
        return;
    }

    if (all_metrics) {
        vector<bool> tf_opts = {true, false};
        for (bool stability : tf_opts) {
            for (bool state : tf_opts) {
                for (bool uniqueness : tf_opts) {
                    if (!stability && !state && !uniqueness)
                        continue;
                    vector<size_t> selected_plan_indexes;
                    compute_metrics_greedy(stability, state, uniqueness, selected_plan_indexes);
                    if (dump_plans) {
                        cout << "Found plans for metric " << get_metric_name(stability, state, uniqueness) << endl;
                        print_plans(selected_plan_indexes);
                    }
                }
            }
        }
    } else {
        vector<size_t> selected_plan_indexes;
        if (exact_method) {
            compute_metrics_mip_external(compute_stability_metric, compute_states_metric, compute_uniqueness_metric, selected_plan_indexes);
        } else {
            compute_metrics_greedy(compute_stability_metric, compute_states_metric, compute_uniqueness_metric, selected_plan_indexes);
        }
        
        if (dump_plans) {
            cout << "Found plans for metric " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
            print_plans(selected_plan_indexes);
        }
    }
}



void DiversityScoreSubset::compute_metrics_greedy(bool stability, bool state, bool uniqueness,
            vector<size_t>& selected_plan_indexes) {

    cout << "Computing metrics " << get_metric_name(stability, state, uniqueness) << endl;
    if (plans_sets.size() == 0)
        return;

    if (plans_sets.size() == 1 || plans_subset_size <= 1) {
        float cluster_score = (plans_subset_size == 1 ? 1.0 : 0.0);
        cout << "Score after clustering " << cluster_score << ", cluster size 1, metrics " << get_metric_name(stability, state, uniqueness) << endl;
        if (plans_subset_size == 1) {
            selected_plan_indexes.push_back(0);
        }
        return;
    }

    // Applying a simple greedy algorithm for selecting plans one by one, until number_of_plans are chosen or plans_sets.size() is reached.
    // First, selecting a pair of plans with a maximal diversity
    float best_score = 0.0;
    size_t best_i = 0;
    size_t best_j = 1;

    for (size_t i = 0; i < plans_sets.size() - 1; ++i) {
        for (size_t j = i+1; j < plans_sets.size(); ++j) {
            float current_score = compute_score_for_pair(stability, state, uniqueness, i, j);
            if (current_score > best_score) {
                best_score = current_score;
                best_i = i;
                best_j = j;
            }
        }
    }
    // Choosing the two plans
    selected_plan_indexes.push_back(best_i);
    selected_plan_indexes.push_back(best_j);
    list<size_t> candidates;
    for (size_t i = 0; i < plans_sets.size(); ++i) {
        if (i == best_i || i == best_j)
            continue;
        candidates.push_back(i);
    }

    while (true) {
        // Finding best match
        if ((int) selected_plan_indexes.size() == plans_subset_size 
            || selected_plan_indexes.size() == plans_sets.size() 
            || candidates.size() == 0)
            break;
        pair<float, size_t> res = find_best_next_candidate(stability, state, uniqueness, selected_plan_indexes, candidates);
        selected_plan_indexes.push_back(res.second);
    }
    if ((int) selected_plan_indexes.size() == plans_subset_size) {
        float cluster_score = compute_score_for_set(stability, state, uniqueness, selected_plan_indexes);
        cout << "Score after clustering " << cluster_score << ", cluster size " << selected_plan_indexes.size() << ", metrics " << get_metric_name(stability, state, uniqueness) << endl;
    }
}

void DiversityScoreSubset::compute_metrics_mip_external(bool stability, bool state, bool uniqueness,
            vector<size_t>& selected_plan_indexes) {
    // For average!
    if (aggregator_metric != Aggregator::AVG) {
        cerr << "The aggregator for the selected formulation should be average!" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }
    (void) selected_plan_indexes;
    cout << "Computing metrics " << get_metric_name(stability, state, uniqueness) << endl;
    if (plans_sets.size() == 0)
        return;

    if (plans_sets.size() == 1)
        return; 

    cout << "Creating the mixed integer program" << endl;
    /*
    Variables: binary variable x_p per plan (whether the plan is in the set)
               continuous variable d_p per plan (for summed diversity per plan)
    Constraints: exactly k plans in the set: sum_p v_p = k
                 d_p >= 0
                 d_p <= M * x_p (if x_p == 0, then d_p = 0) 
                 d_p <= sum_p' x_p' * score(p, p')  (d_p cannot be more than the summed score to all chosen plans)
    Objective: maximize sum_p d_p  (twice the summed pairwise score of all chosen plans)
    */

    ofstream os("test_instance.lp");
    os << "Maximize" << endl;
    os << "obj: sum" << endl;
    os << "Subject To" << endl;

    // Constraints
    size_t constr_count = 1;
    // sum var, for simplicity
    os << "c" << constr_count << ": sum " ;
    for (size_t i=0; i < plans_sets.size(); ++i) {
        os << " - d" << i; 
    }
    os << " = 0" << endl;
    constr_count++;

    // summed score constraints, computing big M
    //float big_m = 1.0;
    for (size_t i = 0; i < plans_sets.size(); ++i) {
        // Adding a constraint
        os << "c" << constr_count << ": d" << i;
        for (size_t j = 0; j < plans_sets.size(); ++j) {
            if (i == j)
                continue;
            float current_score = compute_score_for_pair(stability, state, uniqueness, i, j);
//            big_m += current_score;
            os << " - " << current_score << " x" << j;
        }
        os << " <= 0" << endl;
        constr_count++;
    }    

    // Big M constraints
    for (size_t i = 0; i < plans_sets.size(); ++i) {
        // Adding a constraint
        os << "c" << constr_count << ": d" << i << " - " << plans_sets.size() + 1 << " x" << i << " <= 0" << endl;
        constr_count++;
    }    

    // Num plans constraint
    os << "c" << constr_count << ": x0";
    for (size_t i=1; i < plans_sets.size(); ++i) {
        os << " + x" << i; 
    }
    os << " = " << plans_subset_size << endl;
    constr_count++;

    os << "Bounds" << endl;
    for (size_t i=0; i < plans_sets.size(); ++i) {
        os << "0 <= d" << i << endl; 
    }
    os << "Binaries" << endl;
    for (size_t i=0; i < plans_sets.size(); ++i) {
        os << "x" << i << " "; 
    }
    os << endl;

    os << "End" << endl;
}

pair<float, size_t> DiversityScoreSubset::find_best_next_candidate(bool stability, bool state, bool uniqueness,
        const std::vector<size_t>& selected_plan_indexes, list<size_t>& candidates) {
    // Going over all candidates, check if not already chosen, compute the score, return the best one
    float best_score = 0.0;
    size_t best_candidate = plans_sets.size();
    std::list<size_t>::iterator best_candidate_it;
    std::list<size_t>::iterator it = candidates.begin();
    for (; it != candidates.end(); ++it) {
        size_t i = *it;
        // computing the score as a sum over scores with the existing elements
        float current_score = 0.0;
        for (size_t j : selected_plan_indexes) {
            current_score += compute_score_for_pair(stability, state, uniqueness, i, j);
        }
        //cout << "Current score: " << current_score << ", candidate " << i << endl;
        if (current_score > best_score) {
            best_score = current_score;
            best_candidate = i;
            best_candidate_it = it;
        }
    }
    if (best_candidate == plans_sets.size()) {
        // We just choose the first one
        best_candidate_it = candidates.begin();
        best_candidate = *(best_candidate_it);
    }
    // Remove the best candidate from the list
    candidates.erase(best_candidate_it);
    return make_pair(best_score, best_candidate);
}

static shared_ptr<DiversityScore> _parse_subset(OptionParser &parser) {
    add_diversity_score_options_to_parser(parser);
    add_diversity_score_subset_options_to_parser(parser);    
    Options opts = parser.parse();

    shared_ptr<DiversityScoreSubset> engine;
    if (!parser.dry_run()) {
        engine = make_shared<DiversityScoreSubset>(opts);
    }

    return engine;
}

static Plugin<DiversityScore> _plugin("subset", _parse_subset);
