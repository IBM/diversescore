#include "diversity_score_subset_optimal_long.h"
 
#include "../option_parser.h"
#include "../plugin.h"

#include "../algorithms/max_cliques.h"
#include "../lp/lp_solver.h"

#include <list>
#include <math.h>       /* fabs */

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

//namespace diversity_score {

DiversityScoreSubsetOptimalLong::DiversityScoreSubsetOptimalLong(const Options &opts) : DiversityScore(opts),
        computation_method(ComputationalMethod(opts.get_enum("metric_computation_method"))),
        plans_subset_size(opts.get<int>("plans_subset_size")),
        dump_plans(opts.get<bool>("dump_plans"))
{

}

void DiversityScoreSubsetOptimalLong::compute_metrics() {

    vector<size_t> selected_plan_indexes;

    if (computation_method == ComputationalMethod::MIP)
        compute_metric_mip(selected_plan_indexes);
    else if (computation_method == ComputationalMethod::MIP_EXTERNAL)
        compute_metric_mip_external(selected_plan_indexes);
    else {
        cerr << "Unknown option for optimal metric computation method" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }    
    if (dump_plans) {
        cout << "Found plans for metric " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
        print_plans(selected_plan_indexes);
    }
}

size_t DiversityScoreSubsetOptimalLong::get_binary_var_index(size_t plan_index) const {
    return plan_index + 1;
}

inline size_t pairing_function(size_t x, size_t y, size_t sz) {
    // Cantor pairing function, reversed on y (sz - 1 -y  instead of y)
    size_t val = (x  + sz - y )*(x + sz - y - 1) / 2 + sz - 1 - y;
    return val;
}

size_t DiversityScoreSubsetOptimalLong::get_binary_var_index(size_t plan_index1, size_t plan_index2) const {
    if (plan_index1 < plan_index2)
        return plans_sets.size() + pairing_function(plan_index1, plan_index2, plans_sets.size()) + 1;
    else 
        return plans_sets.size() + pairing_function(plan_index2, plan_index1, plans_sets.size()) + 1;
}

void DiversityScoreSubsetOptimalLong::compute_metric_mip(vector<size_t>& selected_plan_indexes) {
    cout << "Computing metrics " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
    if (plans_sets.size() == 0)
        return;

    if (ordered_plan_indexes.size() == 0 || _plans.size() == 0) {
        cout << " no plans" << endl;
    }
    
    cout << "Creating the integer program" << endl;
    /*
    Variables: binary variable per plan (whether the plan is in the clique) : X_v
               binary variable for each pair of plans (whether the weight should be considered) : X_{u,v}
               single continuous variable for bounding the pairwise diversity : d
    Constraints: the number of edges per chosen node is exactly k-1:   \forall u: \sum_{v} X_{u,v} = X_u * (k-1)    
                 exactly k plans in a clique:                          \sum_{v} X_v = k
                 d is bounded by the diversity of each chosen pair:    \forall u,v :  0 <= d + X_{u,v} <= d(u,v) + 1
    Objective: maximize d
    */

    lp::LPSolver lp_solver(lp::LPSolverType::CPLEX);

    vector<lp::LPVariable> variables;
    variables.push_back(lp::LPVariable(0.0, 1.0, 1.0, false));

    // X_v variables, X_{u,v} variables
    size_t num_binary_vars = plans_sets.size() * (plans_sets.size() + 1) / 2;
    for (size_t i=0; i<num_binary_vars; ++i)
        variables.push_back(lp::LPVariable(0.0, 1.0, 0.0, true));

    cout << "Total number of variables is " << variables.size() << endl;
    vector<lp::LPConstraint> constraints;
    for (size_t i = 0; i < plans_sets.size(); ++i) {
        for (size_t j = i + 1; j < plans_sets.size(); ++j) {
            float current_score = compute_score_for_pair(compute_stability_metric, compute_states_metric, compute_uniqueness_metric, i, j);
            //cout << " Plans " << i << ", " << j << ", score: " << current_score << endl;
            //cout << "Diff: " << (current_score - metric_bound) << endl;
            // Adding a constraint
            lp::LPConstraint d = lp::LPConstraint(0.0, current_score + 1.0 );
            d.insert(0, 1.0);
            size_t ind = get_binary_var_index(i, j);
            //cout << "Variable index " << ind << " for plan pair " << i << ", " << j << endl;
            d.insert(ind, 1.0);
            constraints.push_back(d);            
        }
    }
    for (size_t i = 0; i < plans_sets.size(); ++i) {
        lp::LPConstraint d = lp::LPConstraint(0.0, 0.0);
        size_t u_ind = get_binary_var_index(i);
        //cout << "Variable index " << u_ind << " for plan " << i << endl;
        d.insert(u_ind, plans_subset_size - 1 );
        for (size_t j = 0; j < plans_sets.size(); ++j) {
            if (i ==j) {
                continue;   
            }
            size_t uv_ind = get_binary_var_index(i,j);
            d.insert(uv_ind, -1.0 );
        }
        constraints.push_back(d);            
    }

    lp::LPConstraint exactly_k = lp::LPConstraint(plans_subset_size, plans_subset_size);

    for (size_t i = 0; i < plans_sets.size(); ++i) {
        size_t u_ind = get_binary_var_index(i);
        exactly_k.insert(u_ind, 1.0);
    }
    constraints.push_back(exactly_k);            

    cout << "Loading the problem" << endl;
    lp_solver.load_problem(lp::LPObjectiveSense::MAXIMIZE, variables, constraints);
    cout << "Solving the problem" << endl;
    lp_solver.solve();
    cout << "Done solving" <<endl;
    if (lp_solver.has_optimal_solution()) {
        cout << "Optimal solution is found" << endl;
        cout << "Found subset of size k with minimal pairwise diversity " << lp_solver.get_objective_value() << endl;
        vector<double> solution = lp_solver.extract_solution();
        for (size_t i=0; i < plans_sets.size(); ++i) {
            size_t u_ind = get_binary_var_index(i);
            if (solution[u_ind] > 0.0)
                selected_plan_indexes.push_back(i);
        }
    }

    float score = 0.0;
    if (selected_plan_indexes.size() > 0) {
        score = lp_solver.get_objective_value();
    }
    cout << "Score after clustering " << score << ", cluster size " << selected_plan_indexes.size() << ", metrics " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
}

void DiversityScoreSubsetOptimalLong::compute_metric_mip_external(vector<size_t>& selected_plan_indexes) {
    generate_mip_file(selected_plan_indexes);

}

void DiversityScoreSubsetOptimalLong::generate_mip_file(vector<size_t>& selected_plan_indexes) {

    (void) selected_plan_indexes;
    cout << "Computing metrics " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
    if (plans_sets.size() == 0)
        return;

    if (ordered_plan_indexes.size() == 0 || _plans.size() == 0) {
        cout << " no plans" << endl;
    }
    
    cout << "Creating the integer program" << endl;
    /*
    Variables: binary variable per plan (whether the plan is in the clique) : X_v
               binary variable for each pair of plans (whether the weight should be considered) : X_{u,v}
               single continuous variable for bounding the pairwise diversity : d
    Constraints: the number of edges per chosen node is exactly k-1:   \forall u: \sum_{v} X_{u,v} = X_u * (k-1)    
                 exactly k plans in a clique:                          \sum_{v} X_v = k
                 d is bounded by the diversity of each chosen pair:    \forall u,v :  0 <= d + X_{u,v} <= d(u,v) + 1
    Objective: maximize d
    */


    ofstream os("test_instance.lp");
    os << "Maximize" << endl;
    os << "obj: d" << endl;
    os << "Subject To" << endl;
    size_t constr_count = 1;
    // Constraints

    for (size_t i = 0; i < plans_sets.size(); ++i) {
        for (size_t j = i + 1; j < plans_sets.size(); ++j) {
            float current_score = compute_score_for_pair(compute_stability_metric, compute_states_metric, compute_uniqueness_metric, i, j);
            //cout << " Plans " << i << ", " << j << ", score: " << current_score << endl;
            //cout << "Diff: " << (current_score - metric_bound) << endl;
            // Adding a constraint

            size_t ind = get_binary_var_index(i, j);
            os << "c" << constr_count << ": d + x"<< ind << "<= " << current_score + 1.0 << endl; 
            constr_count++;
        }
    }
    for (size_t i = 0; i < plans_sets.size(); ++i) {
        lp::LPConstraint d = lp::LPConstraint(0.0, 0.0);
        size_t u_ind = get_binary_var_index(i);

        os << "c" << constr_count << ": " << plans_subset_size - 1 << " * x"<< u_ind ; 

        for (size_t j = 0; j < plans_sets.size(); ++j) {
            if (i ==j) {
                continue;   
            }

            size_t uv_ind = get_binary_var_index(i,j);
            os << " - x" << uv_ind; 
        }
        os << " = 0" << endl;
        constr_count++;
    }

    os << "c" << constr_count << ": " ; 
    for (size_t i=0; i < plans_sets.size(); ++i) {
        size_t u_ind = get_binary_var_index(i);
        os << " + x" << u_ind; 
    }
    os << " = " << plans_subset_size << endl;

    os << "Binaries" << endl;
    size_t num_binary_vars = plans_sets.size() * (plans_sets.size() + 1) / 2;

    for (size_t i=0; i < num_binary_vars; ++i) {
        os << "x" << i + 1 << " "; 
    }
    os << endl;
    os << "End" << endl;
}


static shared_ptr<DiversityScore> _parse_subset_optimal(OptionParser &parser) {
    add_diversity_score_options_to_parser(parser);
    add_diversity_score_subset_options_to_parser(parser);    
    add_diversity_score_subset_optimal_options_to_parser(parser);
    Options opts = parser.parse();

    shared_ptr<DiversityScoreSubsetOptimalLong> engine;
    if (!parser.dry_run()) {
        engine = make_shared<DiversityScoreSubsetOptimalLong>(opts);
    }

    return engine;
}

static Plugin<DiversityScore> _plugin("subset_optimal_long", _parse_subset_optimal);
