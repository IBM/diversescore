#include "diversity_score_subset_optimal.h"
 
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

DiversityScoreSubsetOptimal::DiversityScoreSubsetOptimal(const Options &opts) : DiversityScore(opts),
        computation_method(ComputationalMethod(opts.get_enum("metric_computation_method"))),
        plans_subset_size(opts.get<int>("plans_subset_size")),
        dump_plans(opts.get<bool>("dump_plans"))
{

}

void DiversityScoreSubsetOptimal::compute_metrics() {

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

size_t DiversityScoreSubsetOptimal::get_binary_var_index(size_t plan_index) const {
    return plan_index + 1;
}

void DiversityScoreSubsetOptimal::compute_metric_mip(vector<size_t>& selected_plan_indexes) {
    cout << "Computing metrics " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
    if (plans_sets.size() == 0)
        return;

    if (ordered_plan_indexes.size() == 0 || _plans.size() == 0) {
        cout << " no plans" << endl;
    }
    
    cout << "Creating the integer program" << endl;
    /*
    Variables: binary variable per plan (whether the plan is in the clique) : X_v
               single continuous variable for bounding the pairwise diversity : d
    Constraints: exactly k plans in a clique:                          \sum_{v} X_v = k
                 d is bounded by the diversity of each chosen pair:    \forall u,v :  0 <= d <= d(u,v) + 2 - X_u - X_v
    Objective: maximize d
    */

    lp::LPSolver lp_solver(lp::LPSolverType::CPLEX);

    vector<lp::LPVariable> variables;
    variables.push_back(lp::LPVariable(0.0, 1.0, 1.0, false));

    // X_v variables
    for (size_t i=0; i < plans_sets.size(); ++i)
        variables.push_back(lp::LPVariable(0.0, 1.0, 0.0, true));

    cout << "Total number of variables is " << variables.size() << endl;
    vector<lp::LPConstraint> constraints;
    for (size_t i = 0; i < plans_sets.size(); ++i) {
        for (size_t j = i + 1; j < plans_sets.size(); ++j) {
            float current_score = compute_score_for_pair(compute_stability_metric, compute_states_metric, compute_uniqueness_metric, i, j);
            //cout << " Plans " << i << ", " << j << ", score: " << current_score << endl;
            //cout << "Diff: " << (current_score - metric_bound) << endl;
            // Adding a constraint
            lp::LPConstraint d = lp::LPConstraint(0.0, current_score + 2.0 );
            d.insert(0, 1.0);
            d.insert(get_binary_var_index(i), 1.0);
            d.insert(get_binary_var_index(j), 1.0);
            constraints.push_back(d);            
        }
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

void DiversityScoreSubsetOptimal::compute_metric_mip_external(vector<size_t>& selected_plan_indexes) {
    generate_mip_file(selected_plan_indexes);

}

void DiversityScoreSubsetOptimal::generate_mip_file(vector<size_t>& selected_plan_indexes) {

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

            size_t u_ind = get_binary_var_index(i);
            size_t v_ind = get_binary_var_index(j);
            os << "c" << constr_count << ": d + x"<< u_ind << " + x"<< v_ind << " <= " << current_score + 2.0 << endl; 
            constr_count++;
        }
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

    shared_ptr<DiversityScoreSubsetOptimal> engine;
    if (!parser.dry_run()) {
        engine = make_shared<DiversityScoreSubsetOptimal>(opts);
    }

    return engine;
}

static Plugin<DiversityScore> _plugin("subset_optimal", _parse_subset_optimal);
