#include "diversity_score.h"

#include "../option_parser.h"
#include "../plugin.h"

#include "../state_registry.h"
#include "../tasks/root_task.h"

#include <list>
#include <math.h>       /* fabs */

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream> 
#include <string>

using namespace std;

//namespace diversity_score {

DiversityScore::DiversityScore(const Options &opts) :
        task(tasks::g_root_task),
        task_proxy(*task),
        state_registry(task_proxy),
        search_space(state_registry),
        cost_bound(opts.get<int>("cost_bound")),
        plans_as_multisets(opts.get<bool>("plans_as_multisets")),
        use_cache(opts.get<bool>("use_cache")),
        similarity(opts.get<bool>("similarity")),
        reduce_labels(opts.get<bool>("reduce_labels")),
        discounted_prefixes(opts.get<bool>("discounted_prefixes")),
        discount_factor((float)opts.get<double>("discount_factor")),
        plans_seed_set_size(opts.get<int>("plans_seed_set_size")),
        compute_states_metric(opts.get<bool>("compute_states_metric")),
        compute_stability_metric(opts.get<bool>("compute_stability_metric")),
        compute_uniqueness_metric(opts.get<bool>("compute_uniqueness_metric")),
        aggregator_metric(Aggregator(opts.get_enum("aggregator_metric"))),
        dump_pairs(opts.get<bool>("dump_pairs")),
        all_metrics(opts.get<bool>("all_metrics")) 
{
    // Reading the label reduction from json
    if (reduce_labels) {
        if (compute_states_metric) {
            cerr << "Label reduction is currently not implemented with states metric" << endl;
            utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);

        }
        read_label_reduction(opts.get<string>("label_reduction_file"));
    }
}


void DiversityScore::read_label_reduction(string file) {
    cout << "Reading the label reduction from file " << file.c_str() << endl;
    ifstream lr_file;
    lr_file.open(file);
    if (!lr_file.is_open()) {
        throw std::system_error(errno, std::system_category(), "failed to open file");

        cerr << "File is not open!" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }

    std::unordered_map<string, utils::HashSet<OperatorID>> ops_by_reduced_labels;
    std::unordered_map<string, OperatorID> ops_by_names;
    plan_manager.get_ops_by_names(ops_by_names, task_proxy);

    // Reading from CSV
    string line;
    while(std::getline(lr_file, line)) {
        if (line.size() == 0 || line[0] == ';')
            continue;

        //  Splitting the line into two, operator and reduced label
        vector<string> row;
        string op_name, reduced_label;
        stringstream line_to_stream(line);

        while(std::getline(line_to_stream, op_name, ',')) {
            row.push_back(op_name);
        }
        assert(row.size() == 2);
        op_name = row[0];
        reduced_label = row[1];

        //cout << op_name << endl;
        //  Getting the operator id from name
        auto it = ops_by_names.find(op_name);
        if (it == ops_by_names.end()) {
            // Trying adding a trailing space
            string op_name_trailing_space = op_name + " ";
            it = ops_by_names.find(op_name_trailing_space);
            if (it == ops_by_names.end()) {
                cerr << "#" << op_name << "#   Operator not found!!!" << endl;
                cerr << "Operator names:" << endl;
                for (auto name : ops_by_names) {
                    cerr << "#" << name.first << "#" << endl;
                }
                utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
            }
        }
        if (ops_by_reduced_labels.find(reduced_label) == ops_by_reduced_labels.end()) {
            // Adding an empty set
            ops_by_reduced_labels.insert(std::make_pair(reduced_label, utils::HashSet<OperatorID>()));
        }
        ops_by_reduced_labels[reduced_label].insert(it->second);
    }

    // Finding a representative for each class (no matter what the reduced label is)

    for (std::pair<std::string, utils::HashSet<OperatorID>> rc : ops_by_reduced_labels) {
        bool first = true;
        OperatorID op_to(-1);
        for (OperatorID op_from : rc.second) {
            if (first) {
                op_to = op_from;
                first = false;
            }
            label_reduction.insert(std::make_pair(op_from, op_to));
        }
    }
}

OperatorID DiversityScore::get_reduced_label(OperatorID op) const {
    if (!reduce_labels) {
        return op;
    }
    auto it = label_reduction.find(op);
    if (it == label_reduction.end()) {
        return op;
    }
    return it->second;
}


void DiversityScore::compute_metrics_exact_set() {
    vector<size_t> selected_plan_indexes;
    for (size_t i=0; i < _plans.size(); ++i)
        selected_plan_indexes.push_back(i);

    vector<bool> tf_opts = {true, false};
    if (all_metrics) {
        for (bool stability : tf_opts) {
            for (bool state : tf_opts) {
                for (bool uniqueness : tf_opts) {
                    if (!stability && !state && !uniqueness)
                        continue;
                    float cluster_score = compute_score_for_set(stability, state, uniqueness, selected_plan_indexes);
                    // cout << "Score: " << cluster_score << ", metrics " << get_metric_name(stability, state, uniqueness) << endl;
                    cout << "Score after clustering " << cluster_score << ", cluster size " << selected_plan_indexes.size() << ", metrics " << get_metric_name(stability, state, uniqueness) << endl;
                }
            }
        }
    } else {
        float cluster_score = compute_score_for_set(compute_stability_metric, compute_states_metric, compute_uniqueness_metric, selected_plan_indexes);
        // cout << "Score: " << cluster_score << ", metrics " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
        cout << "Score after clustering " << cluster_score << ", cluster size " << selected_plan_indexes.size() << ", metrics " << get_metric_name(compute_stability_metric, compute_states_metric, compute_uniqueness_metric) << endl;
    }
}

void DiversityScore::read_plans() {
    vector<Plan> plans;
    plan_manager.load_plans(plans, task_proxy); 

    assert(_plans.empty());
    if (cost_bound >= 0) {
        cout << "Specified cost bound, ignoring plans with cost under the bound" << endl;
    }
    // We start by removing non-unique plans    
    utils::HashSet<Plan> unique_plans;
    for (Plan plan : plans) {
        if (cost_bound >= 0) {
            int plan_cost = calculate_plan_cost(plan, task_proxy);
            if (plan_cost > cost_bound) {
                cout << "Skipping plan with cost above the bound (" << plan_cost << ")" << endl;
                continue;
            }
        }
        auto it = unique_plans.insert(plan);
        if (it.second) {
            _plans.push_back(plan);
        }
    }
    cout << "Number of unique plans read is " << _plans.size() << " out of " << plans.size() << endl;
    prepare_plans();
}

string DiversityScore::get_metric_name(bool stability, bool state, bool uniqueness) const {
    /*
     * value_attribute_names = {'value_cost' : ['--metric-cost'], 'value_stability' : ['--metric-stability'],
                         'value_state' : ['--metric-state'], 'value_uniqueness' : ['--metric-uniqueness'],
                         'value_state_stability' : ['--metric-state', '--metric-stability'],
                         'value_uniqueness_stability' : ['--metric-uniqueness', '--metric-stability'],
                         'value_state_uniqueness' : ['--metric-state', '--metric-uniqueness'],
                         'value_stability_uniqueness_state' : ['--metric-stability', '--metric-uniqueness', '--metric-state']}
     */
    string name = "";
    if (state && stability && uniqueness) {
        name = "stability_uniqueness_state";
    } else if (state && stability) {
        name = "state_stability";
    } else if (state && uniqueness) {
        name = "state_uniqueness";
    } else if (stability && uniqueness) {
        name = "uniqueness_stability";
    } else if (state) {
        name = "state";
    } else if (stability) {
        name = "stability";
    } else if (uniqueness) {
        name = "uniqueness";
    } else {
        cerr << "At least one of state, stability, uniqueness should be selected" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }
    if (aggregator_metric == Aggregator::AVG) {
        name += "_average";
    } else if (aggregator_metric == Aggregator::MIN) {
        name += "_minimum";
    } else {
        cerr << "At least one of average or minimum should be selected for aggregation method" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }
    if (discounted_prefixes) {
        name += " discounting prefixes";
    }
    return name;
}

bool cmp_plans(pair<size_t, int> a, pair<size_t, int> b) {
    return (a.second < b.second || (a.second == b.second && a.first < b.first));
}

void DiversityScore::prepare_plans() {
    // First, order by cost
    ordered_plan_indexes.clear();
    for (size_t it = 0; it < _plans.size(); ++it) {
        auto plan = _plans[it];
        int plan_cost = calculate_plan_cost(plan, task_proxy);
        if (plans_seed_set_size > 0 && it < (size_t)plans_seed_set_size)
            plan_cost = -1;  // 0 would work too...
        ordered_plan_indexes.push_back(make_pair(it,plan_cost));
    }
    sort(ordered_plan_indexes.begin(), ordered_plan_indexes.end(), cmp_plans);

    // If states metric needs to be computed, we need to compute the traces
    if (compute_states_metric || all_metrics) {
        for (auto plan_index : ordered_plan_indexes) {
            auto plan = _plans[plan_index.first];
            vector<StateID> trace;
            search_space.get_states_trace_from_path(plan, trace, task_proxy);
            plan_traces.push_back(trace);
        }
        assert(_plans.size() == plan_traces.size());
    }
    // Depending on whether we use sets or multisets, we create the plan_set in what follows.
    for (auto plan_index : ordered_plan_indexes) {
        auto plan = _plans[plan_index.first];
        plan_set set_a;
        plan_to_set(set_a, plan, plans_as_multisets);
        plans_sets.push_back(set_a);
    }
}


void DiversityScore::plan_to_set(plan_set &set_a, const Plan &plan, bool plans_as_multisets) const {
    for (auto op : plan) {
        int op_no = get_reduced_label(op).get_index();
        size_t num_elems = 1;
        if (plans_as_multisets) {
            auto e = set_a.find(op_no);
            if (e != set_a.end()) {
                num_elems += e->second;
            }
        }
        set_a[op_no] = num_elems;
    }
}


void DiversityScore::print_plan(size_t ind) {
    cout << "---------------------------------------------" << endl 
         << "Plan index: " << ind << endl;
    const Plan& plan = get_plan(ind);
    plan_manager.save_plan(plan, task_proxy, true);
}


void DiversityScore::print_plans(const std::vector<size_t>& selected_plan_indexes) {
    for (size_t ind : selected_plan_indexes) {
        print_plan(ind);
    }
}


void DiversityScore::print_all_plans() {
    for (size_t ind =0; ind < ordered_plan_indexes.size(); ++ind) {
        print_plan(ind);
    }
}


void DiversityScore::print_plans_json(const std::vector<size_t>& selected_plan_indexes, std::ostream& os) {
    os << "{ \"plans\" : [" << endl;
    bool first_dumped = false;
    for (size_t ind : selected_plan_indexes) {
        if (first_dumped)
            os << "," << endl;
        plan_manager.dump_plan_json(get_plan(ind), task_proxy, os);
        first_dumped = true;
    }
    os << "]}" << endl;
}

void DiversityScore::print_all_plans_json(std::ostream& os) {
    os << "{ \"plans\" : [" << endl;
    bool first_dumped = false;
    for (size_t ind =0; ind < ordered_plan_indexes.size(); ++ind) {
        if (first_dumped)
            os << "," << endl;
        plan_manager.dump_plan_json(get_plan(ind), task_proxy, os);
        first_dumped = true;
    }
    os << "]}" << endl;
}


size_t DiversityScore::get_num_actions(const plan_set& actions_set) const {
    size_t set_size = 0;
    for (auto ea : actions_set) {
        set_size += ea.second;
    }
    return set_size;
}

size_t DiversityScore::get_num_actions(size_t ind) const {
    return get_plan(ind).size();
}


const Plan& DiversityScore::get_plan(size_t ind) const {
    size_t plan_ind = ordered_plan_indexes[ind].first;
    return _plans[plan_ind];
}


float DiversityScore::compute_score_for_set(bool stability, bool state, bool uniqueness,
        const vector<size_t>& selected_plan_indexes) {
    if (aggregator_metric == Aggregator::AVG) {
        return compute_score_for_set_avg(stability, state, uniqueness, selected_plan_indexes);
    } else if (aggregator_metric == Aggregator::MIN) {
        return compute_score_for_set_min(stability, state, uniqueness, selected_plan_indexes);
    } else {
        cerr << "Undefined aggregation method" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }
    return -1.0;
}

float DiversityScore::compute_score_for_set_avg(bool stability, bool state, bool uniqueness,
        const vector<size_t>& selected_plan_indexes) {
    float res = 0.0;
    size_t count_pairs = 0;
    for (size_t i=0; i < selected_plan_indexes.size() - 1; ++i) {
        for (size_t j=i+1; j < selected_plan_indexes.size(); ++j) {
            count_pairs++;
            float r = compute_score_for_pair(stability, state, uniqueness, selected_plan_indexes[i], selected_plan_indexes[j]);
            size_t plan_ind1 = ordered_plan_indexes[selected_plan_indexes[i]].first;
            size_t plan_ind2 = ordered_plan_indexes[selected_plan_indexes[j]].first;
            if (dump_pairs)
                cout << "Score for pair " << plan_ind1 <<", " << plan_ind2 << ": " << r << endl;
            res += r;
        }
    }
    return res / (float) count_pairs;
}

float DiversityScore::compute_score_for_set_min(bool stability, bool state, bool uniqueness,
        const vector<size_t>& selected_plan_indexes) {
    // Not assuming any value on the elements, returning 0 when no pairs exist
    if (selected_plan_indexes.size() < 2) {
        return 0.0;
    }
    float min_res = numeric_limits<float>::max();
    for (size_t i=0; i < selected_plan_indexes.size() - 1; ++i) {
        for (size_t j=i+1; j < selected_plan_indexes.size(); ++j) {
            float res = compute_score_for_pair(stability, state, uniqueness, selected_plan_indexes[i], selected_plan_indexes[j]);
            size_t plan_ind1 = ordered_plan_indexes[selected_plan_indexes[i]].first;
            size_t plan_ind2 = ordered_plan_indexes[selected_plan_indexes[j]].first;
            if (dump_pairs)
                cout << "Score for pair " << plan_ind1 <<", " << plan_ind2 << ": " << res << endl;
            if (res < min_res)
                min_res = res;
        }
    }
    return min_res;
}


float DiversityScore::compute_discounted_prefix_similarity(bool stability, bool state, bool uniqueness,
        const Plan& plan1, const Plan& plan2, float gamma) const {

    size_t m = min(plan1.size(), plan2.size());
    float res = 0.0;
    float gamma_to_i = 1.0;
    Plan prefix1, prefix2;
    for (size_t i=0; i<m; ++i) {
        gamma_to_i *= gamma;
        prefix1.push_back(plan1[i]);
        prefix2.push_back(plan2[i]);

        // OperatorsProxy operators = task_proxy.get_operators();
        // cout <<"Prefix1:";

        // for (OperatorID a : prefix1) {
        //     cout << "  " << operators[a].get_name();
        // }
        // cout << endl << "Prefix2:";
        // for (OperatorID a : prefix2) {
        //     cout << "  " << operators[a].get_name();
        // }
        // cout << endl;

        float s = compute_similarity_for_prefix_no_cache(stability, state, uniqueness, prefix1, prefix2);
        res += s * gamma_to_i;
        // cout << "i=" << i+1 << ", gamma^i="<<gamma_to_i << ", s=" << s<< ", res=" << res << endl;
    }
    gamma_to_i *= gamma;
    float geosum = (gamma - gamma_to_i) / (1.0 - gamma);
    float ret = res / geosum;
    // cout << "S_m="<< geosum<< ", sym="<< ret << endl;
    return 1.0 - ret;
}

float DiversityScore::compute_similarity_for_prefix_no_cache(bool stability, bool state, bool uniqueness,
        const Plan& plan1, const Plan& plan2) const {

    float divide_by = 0.0;
    float res = 0.0;

    plan_set set_a;
    plan_set set_b;

    if (stability || uniqueness) {
        plan_to_set(set_a, plan1, plans_as_multisets);
        plan_to_set(set_b, plan2, plans_as_multisets);
    }

    if (stability) {
        divide_by++;
        res += compute_jaccard_similarity_score(set_a, set_b);
    }
    if (state) {
        divide_by++;
        vector<StateID> trace_a;
        vector<StateID> trace_b;

        search_space.get_states_trace_from_path(plan1, trace_a, task_proxy);
        search_space.get_states_trace_from_path(plan2, trace_b, task_proxy);
        res += compute_state_similarity_score(trace_a, trace_b);
    }
    if (uniqueness) {
        divide_by++;
        res += compute_uniqueness_similarity_score(set_a, set_b);
    }

    float ret = (res / divide_by);
    //cout << "Score for a pair: " << ret << ", indices " << plan_index1 << ", " << plan_index2 << endl;
    if (ret > 1.0) {
        cout << "Negative score between plans:" << endl;
        OperatorsProxy operators = task_proxy.get_operators();
        for (OperatorID a : plan1) {
            cout << operators[a].get_name() <<  endl;
        }
        cout << " ----------- " << endl;
        for (OperatorID a : plan2) {
            cout << operators[a].get_name() <<  endl;
        }
        cout << "-------------------------------- " << endl;
    }
    return ret;
}


float DiversityScore::compute_score_for_pair(bool stability, bool state, bool uniqueness,
        size_t plan_index1, size_t plan_index2) {

    if (discounted_prefixes) {
        float ret = compute_discounted_prefix_similarity(stability, state, uniqueness, get_plan(plan_index1), get_plan(plan_index2), discount_factor);
        return similarity ? 1.0 - ret: ret;

    }

    float divide_by = 0.0;
    float res = 0.0;
    if (stability) {
        divide_by++;
        float stab = compute_stability_similarity_score(plan_index1, plan_index2);
        res += stab;
    }
    if (state) {
        divide_by++;
        res += compute_state_similarity_score(plan_index1, plan_index2);
    }
    if (uniqueness) {
        divide_by++;
        res += compute_uniqueness_similarity_score(plan_index1, plan_index2);
    }

    float ret = 1.0 - (res / divide_by);
    //cout << "Score for a pair: " << ret << ", indices " << plan_index1 << ", " << plan_index2 << endl;
    if (ret < 0.0 || ret > 1.0) {
        cout << "Negative score between plans:" << endl;
        OperatorsProxy operators = task_proxy.get_operators();
        for (OperatorID a : get_plan(plan_index1)) {
            cout << operators[a].get_name() <<  endl;
        }
        cout << " ----------- " << endl;
        for (OperatorID a : get_plan(plan_index2)) {
            cout << operators[a].get_name() <<  endl;
        }
        cout << "-------------------------------- " << endl;
    }
    return similarity ? 1.0 - ret: ret;
}

float DiversityScore::compute_uniqueness_similarity_score(const plan_set& set_a, const plan_set& set_b) const {
    // Check if the plans are equal as sets or one contains in another
    if (set_a.size() < set_b.size())
        return compute_uniqueness_similarity_score(set_b, set_a);

    if (set_a.size() == 0 && set_b.size() == 0)
        return 1.0;

    if (set_a.size() == 0 || set_b.size() == 0)
        return 0.0;

    // Here, set_b.size() <= set_a.size(), so checking if b is subset of a
    for (auto eb : set_b) {
        auto ea = set_a.find(eb.first);
        if ( ea == set_a.end() || ea->second < eb.second) {
            // There is an element in b that is not in a, unique
            return 0.0;
        }
    }
    // All elements of b are in a
    return 1.0;
}

float DiversityScore::compute_uniqueness_similarity_score(size_t plan_index1, size_t plan_index2) {
    // Check if the plans are equal as sets or one contains in another
    const plan_set& set_a = plans_sets[plan_index1];
    const plan_set& set_b = plans_sets[plan_index2];

    if (set_a.size() < set_b.size())
        return compute_uniqueness_similarity_score(plan_index2, plan_index1);

    if (set_a.size() == 0 && set_b.size() == 0)
        return 1.0;

    if (set_a.size() == 0 || set_b.size() == 0)
        return 0.0;

    // Try to take from cache;
    float ret = get_uniqueness_from_cache(plan_index1, plan_index2);
    if (ret >= 0.0)
        return ret;

    ret = compute_uniqueness_similarity_score(set_a, set_b);
    add_uniqueness_to_cache(plan_index1, plan_index2, ret);
    return ret;
}

float DiversityScore::compute_jaccard_similarity_score(const plan_set& set_a, const plan_set& set_b) const {
    if (set_a.size() == 0 && set_b.size() == 0)
        return 1.0;

    if (set_a.size() == 0 || set_b.size() == 0)
        return 0.0;

    if (set_a.size() < set_b.size())
        return compute_jaccard_similarity_score(set_b, set_a);

    size_t set_b_minus_a_size = 0;
    for (auto eb : set_b) {
        auto ea = set_a.find(eb.first);
        if ( ea == set_a.end()) {
            set_b_minus_a_size += eb.second;
        } else if ( ea->second < eb.second) {
            set_b_minus_a_size += eb.second - ea->second;
        }         
    }

    size_t set_a_size = (plans_as_multisets) ? get_num_actions(set_a) : set_a.size();
    size_t set_b_size = (plans_as_multisets) ? get_num_actions(set_b) : set_b.size();
    
    size_t union_size = set_a_size + set_b_minus_a_size;
    size_t intersection_size = set_b_size - set_b_minus_a_size;
    //cout << "Union: " << union_size << ", intersection: " << intersection_size << endl;
    assert (union_size >= intersection_size);
    return (float) intersection_size / (float) union_size;
}

float DiversityScore::compute_stability_similarity_score(size_t plan_index1, size_t plan_index2) {
    const plan_set& set_a = plans_sets[plan_index1];
    const plan_set& set_b = plans_sets[plan_index2];

    if (set_a.size() == 0 && set_b.size() == 0)
        return 1.0;

    if (set_a.size() == 0 || set_b.size() == 0)
        return 0.0;

    if (set_a.size() < set_b.size())
        return compute_stability_similarity_score(plan_index2, plan_index1);

    // Try to take from cache;
    float ret = get_stability_from_cache(plan_index1, plan_index2);
    if (ret >= 0.0)
        return ret;

    ret = compute_jaccard_similarity_score(set_a, set_b);

    add_stability_to_cache(plan_index1, plan_index2, ret);

    return ret;
}

float DiversityScore::compute_state_similarity_score(const vector<StateID>& plan_a, const vector<StateID>& plan_b) const {
    if (plan_a.size() == 0 && plan_b.size() == 0)
        return 1.0;

    if (plan_a.size() == 0 || plan_b.size() == 0)
        return 0.0;

    if (plan_a.size() > plan_b.size())
        return compute_state_similarity_score(plan_b, plan_a);

    // Here, plan_a.size() <= plan_b.size()
    //float sum = (float) (plan_b.size() - plan_a.size());
    float sum = 0.0;
    for (size_t i = 0; i < plan_a.size(); ++i) {
        sum += compute_state_similarity_score(plan_a[i], plan_b[i]);
    }
    return sum / (float) plan_b.size();
}

float DiversityScore::compute_state_similarity_score(size_t plan_index1, size_t plan_index2) {
    const vector<StateID>& plan_a = plan_traces[plan_index1];
    const vector<StateID>& plan_b = plan_traces[plan_index2];
    if (plan_a.size() == 0 && plan_b.size() == 0)
        return 1.0;

    if (plan_a.size() == 0 || plan_b.size() == 0)
        return 0.0;

    if (plan_a.size() > plan_b.size())
        return compute_state_similarity_score(plan_index2, plan_index1);

    // Try to take from cache;
    float ret = get_state_from_cache(plan_index1, plan_index2);
    if (ret >= 0.0)
        return ret;
    ret = compute_state_similarity_score(plan_a, plan_b);
    // Add to cache
    add_state_to_cache(plan_index1, plan_index2, ret);
    return ret;
}

float DiversityScore::get_state_from_cache(size_t plan_index1, size_t plan_index2) const {
    if (!use_cache)
        return -1.0;
    auto it = state_metric_cache.find(plan_index1);
    if ( it != state_metric_cache.end()) {
        auto it2 = (*it).second.find(plan_index2);
        if (it2 != (*it).second.end())
            return (*it2).second;
    }
    return -1.0;
}

void DiversityScore::add_state_to_cache(size_t plan_index1, size_t plan_index2, float value) {
    if (!use_cache)
        return;
    auto it = state_metric_cache.find(plan_index1);

    if ( it == state_metric_cache.end()) {
        state_metric_cache.insert(std::make_pair(plan_index1, unordered_map<size_t, float>()));
        it = state_metric_cache.find(plan_index1);
        assert(it != state_metric_cache.end());
    }

    assert((*it).second.find(plan_index2) == (*it).second.end());
    (*it).second.insert(std::make_pair(plan_index2, value));
}

float DiversityScore::get_stability_from_cache(size_t plan_index1, size_t plan_index2) const {
    if (!use_cache)
        return -1.0;
    if (plan_index1 > plan_index2)
        return get_stability_from_cache(plan_index2, plan_index1);

    auto it = stability_metric_cache.find(plan_index1);
    if ( it != stability_metric_cache.end()) {
        auto it2 = (*it).second.find(plan_index2);
        if (it2 != (*it).second.end())
            return (*it2).second;
    }
    return -1.0;
}

void DiversityScore::add_stability_to_cache(size_t plan_index1, size_t plan_index2, float value) {
    if (!use_cache)
        return;
    if (plan_index1 > plan_index2) {
        add_stability_to_cache(plan_index2, plan_index1, value);
        return;
    }
    auto it = stability_metric_cache.find(plan_index1);

    if ( it == stability_metric_cache.end()) {
        stability_metric_cache.insert(std::make_pair(plan_index1, unordered_map<size_t, float>()));
        it = stability_metric_cache.find(plan_index1);
        assert(it != stability_metric_cache.end());
    }

    assert((*it).second.find(plan_index2) == (*it).second.end());
    (*it).second.insert(std::make_pair(plan_index2, value));
}


float DiversityScore::get_uniqueness_from_cache(size_t plan_index1, size_t plan_index2) const {
    if (!use_cache)
        return -1.0;
    if (plan_index1 > plan_index2)
        return get_uniqueness_from_cache(plan_index2, plan_index1);

    auto it = uniqueness_metric_cache.find(plan_index1);
    if ( it != uniqueness_metric_cache.end()) {
        auto it2 = (*it).second.find(plan_index2);
        if (it2 != (*it).second.end())
            return (*it2).second;
    }
    return -1.0;
}

void DiversityScore::add_uniqueness_to_cache(size_t plan_index1, size_t plan_index2, float value) {
    if (!use_cache)
        return;
    if (plan_index1 > plan_index2) {
        add_uniqueness_to_cache(plan_index2, plan_index1, value);
        return;
    }
    auto it = uniqueness_metric_cache.find(plan_index1);

    if ( it == uniqueness_metric_cache.end()) {
        uniqueness_metric_cache.insert(std::make_pair(plan_index1, unordered_map<size_t, float>()));
        it = uniqueness_metric_cache.find(plan_index1);
        assert(it != uniqueness_metric_cache.end());
    }

    assert((*it).second.find(plan_index2) == (*it).second.end());
    (*it).second.insert(std::make_pair(plan_index2, value));
}

bool DiversityScore::is_none_of_those(int var, int val) const {
    return task->get_fact_name(FactPair(var, val)) == "none_of_those";
}

float DiversityScore::compute_state_similarity_score(StateID state1, StateID state2) const {
    GlobalState s1 = state_registry.lookup_state(state1);
    GlobalState s2 = state_registry.lookup_state(state2);
    int num_vars = state_registry.get_num_variables();
    int same_vars = 0;
    int all_vals = 0;
    for (int var=0; var < num_vars; ++var) {
        bool s1_none_of_those = is_none_of_those(var, s1[var]);
        bool s2_none_of_those = is_none_of_those(var, s2[var]);
        if (s1_none_of_those && s2_none_of_those)
            continue;
        if (s1_none_of_those || s2_none_of_those) {
            // s1 is not the same as s2, only one of them corresponds to an actual fact
            all_vals++;
            continue;
        }
        if (s1[var] == s2[var]) {
            same_vars++;
            all_vals++;
        } else {
            all_vals += 2;
        }
    }

    // Returning the intersection size divided by union (and not the number of variables)
    return (float) same_vars / (float) all_vals;
}

void add_diversity_score_options_to_parser(OptionParser &parser) {
    parser.add_option<bool>("compute_states_metric", "Computing the metric states", "false");
    parser.add_option<bool>("compute_stability_metric", "Computing the metric stability", "false");
    parser.add_option<bool>("compute_uniqueness_metric", "Computing the metric uniqueness", "false");

    parser.add_option<bool>("similarity", "Computing similarity instead of dissimilarity", "false");
    parser.add_option<bool>("dump_pairs", "Dumping the score of each pair", "false");

    vector<string> aggregator;
    aggregator.push_back("AVG");
    aggregator.push_back("MIN");
    parser.add_enum_option(
        "aggregator_metric",
        aggregator,
        "Method for aggregating pairwise results into a metric score",
        "AVG");

    parser.add_option<bool>("all_metrics", "Computing all metrics", "false");
    parser.add_option<bool>("plans_as_multisets", "Treat plans as multisets instead of sets", "false");
    parser.add_option<bool>("use_cache", "Use cache when computing metrics", "true");

    parser.add_option<bool>("discounted_prefixes", "Computing discounted prefixes metric", "false");
    parser.add_option<double>("discount_factor", "Discount factor", "0.99");

    parser.add_option<bool>("reduce_labels", "Perform label reduction", "false");
    parser.add_option<string>("label_reduction_file", "CSV file containing the label reduction. \
                                Computing the similarity based on reduced action sequences (stability/uniqueness)", "false");

    parser.add_option<int>("cost_bound",
        "The bound on the cost of a plan",
        "-1");

    parser.add_option<int>("plans_seed_set_size",
        "The number of plans to seed the selection with",
        "0");
}

void add_diversity_score_subset_options_to_parser(OptionParser &parser) {

    parser.add_option<int>("plans_subset_size",
        "The subset size to compute a metric for",
        "0");

    parser.add_option<bool>("exact_method", "Computing the subset using exact method, generating mip", "false");

    parser.add_option<bool>("dump_plans", "Dumping (subset of) plans", "false");
    parser.add_option<string>("json_file_to_dump",
        "A path to the json file to use for dumping",
        OptionParser::NONE);
}

void add_diversity_score_subset_bounded_options_to_parser(OptionParser &parser) {

    parser.add_option<double>("metric_bound", "Metric bound", "1.0");
    
    vector<string> computation_method;
    computation_method.push_back("DFS");
    computation_method.push_back("MIP");
    computation_method.push_back("MIP_EXTERNAL");
    computation_method.push_back("MAX_CLIQUES");
    parser.add_enum_option(
        "metric_computation_method",
        computation_method,
        "Method for finding subset of plans with metric score above the bound",
        "DFS");


}

void add_diversity_score_subset_optimal_options_to_parser(OptionParser &parser) {

    vector<string> computation_method;
    computation_method.push_back("MIP");
    computation_method.push_back("MIP_EXTERNAL");
    parser.add_enum_option(
        "metric_computation_method",
        computation_method,
        "Method for finding subset of plans with optimal metric score",
        "MIP");

}


static PluginTypePlugin<DiversityScore> _type_plugin(
    "DiversityScore",
    // TODO: Replace empty string by synopsis for the wiki page.
    "");
//}

