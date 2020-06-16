#include "diversity_score_subset_bounded.h"
 
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

DiversityScoreSubsetBounded::DiversityScoreSubsetBounded(const Options &opts) : DiversityScore(opts),
        metric_bound((float)opts.get<double>("metric_bound")),
        computation_method(ComputationalMethod(opts.get_enum("metric_computation_method"))),
        plans_subset_size(opts.get<int>("plans_subset_size")),
        dump_plans(opts.get<bool>("dump_plans"))
{

}

void DiversityScoreSubsetBounded::compute_metrics() {

    vector<size_t> selected_plan_indexes;

    if (computation_method == ComputationalMethod::DFS)
        compute_metric_dfs(selected_plan_indexes);
    else if (computation_method == ComputationalMethod::MAX_CLIQUES)
        compute_metric_maxclique(selected_plan_indexes);
    else if (computation_method == ComputationalMethod::MIP)
        compute_metric_mip(selected_plan_indexes);
    else if (computation_method == ComputationalMethod::MIP_EXTERNAL)
        compute_metric_mip_external(selected_plan_indexes);
    else {
        cerr << "Unknown option for bounded metric computation method" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }    
    if (dump_plans) {
        cout << "Found plans for metric " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
        print_plans(selected_plan_indexes);
    }
}

void DiversityScoreSubsetBounded::compute_metric_maxclique(vector<size_t>& selected_plan_indexes) {

    cout << "Computing metrics " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
    if (plans_sets.size() == 0)
        return;

    if (ordered_plan_indexes.size() == 0 || _plans.size() == 0) {
        cout << " no plans" << endl;
    }
    
    cout << "Computing max cliques" << endl;      
    vector<vector<int>> cgraph;
    cgraph.resize(plans_sets.size());

    cout << "Building connectivity graph" <<endl;
    int num_nodes = 0;
    for (size_t i = 0; i < plans_sets.size(); ++i) {

        for (size_t j = i + 1; j < plans_sets.size(); ++j) {
            float current_score = compute_score_for_pair(true, false, false, i, j);
            if (current_score >= metric_bound) {
                cgraph[i].push_back(j);
                cgraph[j].push_back(i);
            }
        }
        if (cgraph[i].size() > 0)
            num_nodes++;
        //cout << "---------------------------------------" << endl;
    }
    cout << "Constructed connectivity graph with " << num_nodes << " nodes" << endl;
    /*
     Going over the graph and removing nodes with the number of edges below the number_of_plans
     These nodes cannot be part of a valid clique of size >= number_of_plans
     Iterating while there is any change to the graph
    */
    bool change = true;
    while(change) {
        unordered_set<int> toremove;
        for (size_t i = 0; i < cgraph.size(); ++i) {
            if (cgraph[i].size() == 0) 
                continue;

            if ((int)cgraph[i].size() < plans_subset_size) {
                toremove.insert(i);
                vector<int> empty;
                cgraph[i].swap(empty);
            }
        }
        change = (toremove.size() > 0);
        cout << "Removing " << toremove.size() << " nodes" << endl;
        for (size_t i = 0; i < cgraph.size(); ++i) {
            if (cgraph[i].size() == 0) 
                continue;
            vector<int> remaining;
            for (int v : cgraph[i]) {
                if (toremove.find(v) == toremove.end())
                    remaining.push_back(v);
            }
            cgraph[i].swap(remaining);
        }
    }
    cout << "Done Building connectivity graph" <<endl;
    num_nodes = 0;
    for (size_t i = 0; i < cgraph.size(); ++i) {
        if (cgraph[i].size() > 0)
            num_nodes++; 
    }
    cout << "Number of nodes remaining " << num_nodes << endl;
    vector<vector<int>> max_cliques;
    max_cliques::compute_max_cliques(cgraph, max_cliques);
    cout << "Done computing max cliques" <<endl;
    int max_clique_size = 0;
    for (vector<int> clique : max_cliques) {
        if ((int) clique.size() > max_clique_size)
            max_clique_size = clique.size();
    }
    cout << "Maximal clique size " << max_clique_size << endl;

    for (vector<int> clique : max_cliques) {
        if ((int)clique.size() >= plans_subset_size) {
            selected_plan_indexes.insert(selected_plan_indexes.begin(), clique.begin(), clique.begin() + plans_subset_size);
            break;
        }
    }
    float score = 0.0;
    if (selected_plan_indexes.size() > 0) {
        score = 1.0;
    }
    cout << "Score after clustering " << score << ", cluster size " << selected_plan_indexes.size() << ", metrics " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
}

void DiversityScoreSubsetBounded::compute_metric_mip(vector<size_t>& selected_plan_indexes) {

    cout << "Computing metrics " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
    if (plans_sets.size() == 0)
        return;

    if (ordered_plan_indexes.size() == 0 || _plans.size() == 0) {
        cout << " no plans" << endl;
    }
    
    cout << "Creating the integer program" << endl;
    /*
    Variables: binary variable per plan (whether the plan is in the clique)
    Constraints: if (u,v) below the bound, then at most one of them can be in the same clique: u+v<=1
                 at least k plans in a clique: sum_v >= k
    Objective: minimize sum_v
    */

    lp::LPSolver lp_solver(lp::LPSolverType::CPLEX);

    vector<lp::LPVariable> variables;
    for (auto plan : plans_sets)
        variables.push_back(lp::LPVariable(0.0, 1.0, 1.0, true));

    vector<lp::LPConstraint> constraints;

//    cout << "Minimize" << endl;
    float EPSILON=0.00001;
    // Constraints

    lp::LPConstraint at_least_k = lp::LPConstraint(plans_subset_size, lp_solver.get_infinity());
    for (size_t i = 0; i < variables.size(); ++i )
        at_least_k.insert(i, 1.0);
    constraints.push_back(at_least_k);

    for (size_t i = 0; i < plans_sets.size(); ++i) {
        for (size_t j = i + 1; j < plans_sets.size(); ++j) {
            float current_score = compute_score_for_pair(compute_stability_metric, compute_states_metric, compute_uniqueness_metric, i, j);
            //cout << " Plans " << i << ", " << j << ", score: " << current_score << endl;
            //cout << "Diff: " << (current_score - metric_bound) << endl;
            if (current_score + EPSILON < metric_bound) {
                // Adding a constraint
                lp::LPConstraint d = lp::LPConstraint(0.0, 1.0);
                d.insert(i, 1.0);
                d.insert(j, 1.0);
                constraints.push_back(d);
            }
        }
    }


    cout << "Loading the problem" << endl;
    lp_solver.load_problem(lp::LPObjectiveSense::MINIMIZE, variables, constraints);
    cout << "Solving the problem" << endl;
    lp_solver.solve();
    cout << "Done solving" <<endl;
    if (lp_solver.has_optimal_solution()) {
        cout << "Optimal solution is found" << endl;
        cout << "Found clique of size " << lp_solver.get_objective_value() << endl;
        vector<double> solution = lp_solver.extract_solution();
        for (size_t i=0; i < solution.size(); ++i) {
            if (solution[i] > 0.0)
                selected_plan_indexes.push_back(i);

        }
    }

    float score = 0.0;
    if (selected_plan_indexes.size() > 0) {
        score = 1.0;
    }
    cout << "Score after clustering " << score << ", cluster size " << selected_plan_indexes.size() << ", metrics " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
}

void DiversityScoreSubsetBounded::compute_metric_mip_external(vector<size_t>& selected_plan_indexes) {
    generate_mip_file(selected_plan_indexes);

}

void DiversityScoreSubsetBounded::generate_mip_file(vector<size_t>& selected_plan_indexes) {
    (void) selected_plan_indexes;
    cout << "Computing metrics " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
    if (plans_sets.size() == 0)
        return;

    if (ordered_plan_indexes.size() == 0 || _plans.size() == 0) {
        cout << " no plans" << endl;
    }
    
    cout << "Creating the integer program" << endl;
    /*
    Variables: binary variable per plan (whether the plan is in the clique)
    Constraints: if (u,v) below the bound, then at most one of them can be in the same clique: u+v<=1
                 at least k plans in a clique: sum_v >= k
    Objective: minimize sum_v
    */

    ofstream os("test_instance.lp");
    os << "Minimize" << endl;
    os << "obj: sum" << endl;
    os << "Subject To" << endl;
    float EPSILON=0.00001;
    size_t constr_count = 1;
    // Constraints
    for (size_t i = 0; i < plans_sets.size(); ++i) {
        for (size_t j = i + 1; j < plans_sets.size(); ++j) {
            float current_score = compute_score_for_pair(compute_stability_metric, compute_states_metric, compute_uniqueness_metric, i, j);
            //cout << " Plans " << i << ", " << j << ", score: " << current_score << endl;
            //cout << "Diff: " << (current_score - metric_bound) << endl;
            if (current_score + EPSILON < metric_bound) {
                // Adding a constraint
                os << "c" << constr_count << ": x" << i << " + x"<< j << "<= 1" << endl; 
                constr_count++;
            }
        }
    }

    os << "c" << constr_count << ": sum " ;
    for (size_t i=0; i < plans_sets.size(); ++i) {
        os << " - x" << i; 
    }
    os << " = 0" << endl;
    constr_count++;
    os << "c" << constr_count << ": sum >= " << plans_subset_size << endl;
    os << "Binaries" << endl;
    for (size_t i=0; i < plans_sets.size(); ++i) {
        os << "x" << i << " "; 
    }
    os << endl;

    os << "End" << endl;
    
}

void DiversityScoreSubsetBounded::compute_metric_dfs(vector<size_t>& selected_plan_indexes) {

    cout << "Computing metrics " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
    if (plans_sets.size() == 0)
        return;

    if (ordered_plan_indexes.size() == 0 || _plans.size() == 0) {
        cout << " no plans" << endl;
    }
    
    cout << "Running DFS to find a clique of a given size" << endl;      
    /*
    Naive implementation: Nodes are cliques, starting with an empty clique.
    */
	vector<shared_ptr<CliqueDFSNode>> dfs_queue;

	shared_ptr<CliqueDFSNode> init_node = make_shared<CliqueDFSNode>(-1, nullptr);
	dfs_queue.push_back(init_node);
	while (!dfs_queue.empty()) {
		shared_ptr<CliqueDFSNode> curr = dfs_queue.back();
        vector<size_t> clique;
        curr->get_clique(clique);
        if ((int)clique.size() == plans_subset_size) {
            // Found the clique
            selected_plan_indexes.insert(selected_plan_indexes.begin(), clique.begin(), clique.end());
            break;
        }
		dfs_queue.pop_back();
        bool done = false;
        // Going over the nodes, if not in clique, trying to see if adding creates a clique  
        // Pruning nodes by considering only increasing order in added indicies to the clique
        // Since we are pruning nodes bases on the order, we need to add the successors in a reverse order
        size_t max_index = clique.size();
        if (clique.size() > 0)
            max_index = clique[clique.size()-1] + 1;


        std::stringstream clique_stream;
        std::copy(clique.begin(), clique.end(), std::ostream_iterator<size_t>(clique_stream, " "));

        for (size_t i = plans_sets.size(); i-- > max_index; ) {
            //if (clique.find(i) != clique.end())
            //    continue;
            bool is_clique = true;
            for (size_t j : clique) {
                float current_score = compute_score_for_pair(true, false, false, i, j);
                if (current_score < metric_bound) {
                    is_clique = false;
                    break;
                }
            }
            if (is_clique) {
                // Check if not done
                if ((int)clique.size() == plans_subset_size - 1) {
                    selected_plan_indexes.insert(selected_plan_indexes.begin(), clique.begin(), clique.end());
                    selected_plan_indexes.push_back(i);
                    done = true;
                    break;
                }
                // Adding a successor
                cout << "Generated successor of size " << clique.size() + 1 << "  " << clique_stream.str().c_str() << i << endl;

                shared_ptr<CliqueDFSNode> succ = make_shared<CliqueDFSNode>(i,curr);
				dfs_queue.push_back(succ);
            }
        }
        if (done)
            break;
	}

    float score = 0.0;
    if ((int)selected_plan_indexes.size() == plans_subset_size) {
        score = 1.0;
    }
    cout << "Selected indices: ";
    for (auto i : selected_plan_indexes)
        cout << " " << i;
    cout << endl;
    cout << "Score after clustering " << score << ", cluster size " << selected_plan_indexes.size() << ", metrics " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
}


static shared_ptr<DiversityScore> _parse_subset_bounded(OptionParser &parser) {
    add_diversity_score_options_to_parser(parser);
    add_diversity_score_subset_options_to_parser(parser);    
    add_diversity_score_subset_bounded_options_to_parser(parser);
    Options opts = parser.parse();

    shared_ptr<DiversityScoreSubsetBounded> engine;
    if (!parser.dry_run()) {
        engine = make_shared<DiversityScoreSubsetBounded>(opts);
    }

    return engine;
}

static Plugin<DiversityScore> _plugin("subset_bounded", _parse_subset_bounded);
