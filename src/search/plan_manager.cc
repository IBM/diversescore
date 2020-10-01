#include "plan_manager.h"

#include "task_proxy.h"

#include "task_utils/task_properties.h"
#include "tasks/root_task.h"

#include <fstream>
#include <iostream>
#include <sstream>


using namespace std;

int calculate_plan_cost(const Plan &plan, const TaskProxy &task_proxy) {
    OperatorsProxy operators = task_proxy.get_operators();
    int plan_cost = 0;
    for (OperatorID op_id : plan) {
        plan_cost += operators[op_id].get_cost();
    }
    return plan_cost;
}

PlanManager::PlanManager()
    : plan_filename("sas_plan"),
      num_plans_to_read(0),
      num_previously_generated_plans(0),      
      is_part_of_anytime_portfolio(false) {
}

void PlanManager::set_plan_filename(const string &plan_filename_) {
    plan_filename = plan_filename_;
}

void PlanManager::set_plan_foldername(const string &plan_foldername_) {
    plan_foldername = plan_foldername_;
}

void PlanManager::set_num_previously_generated_plans(int num_previously_generated_plans_) {
    num_previously_generated_plans = num_previously_generated_plans_;
}

void PlanManager::set_is_part_of_anytime_portfolio(bool is_part_of_anytime_portfolio_) {
    is_part_of_anytime_portfolio = is_part_of_anytime_portfolio_;
}

void PlanManager::save_plan(
    const Plan &plan, const TaskProxy &task_proxy,
    bool generates_multiple_plan_files) {
    ostringstream filename;
    filename << plan_filename;
    int plan_number = num_previously_generated_plans + 1;
    if (generates_multiple_plan_files || is_part_of_anytime_portfolio) {
        filename << "." << plan_number;
    } else {
        assert(plan_number == 1);
    }
    ofstream outfile(filename.str());
    if (outfile.rdstate() & ofstream::failbit) {
        cerr << "Failed to open plan file: " << filename.str() << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }
    OperatorsProxy operators = task_proxy.get_operators();
    for (OperatorID op_id : plan) {
        cout << operators[op_id].get_name() << " (" << operators[op_id].get_cost() << ")" << endl;
        outfile << "(" << operators[op_id].get_name() << ")" << endl;
    }
    int plan_cost = calculate_plan_cost(plan, task_proxy);
    bool is_unit_cost = task_properties::is_unit_cost(task_proxy);
    outfile << "; cost = " << plan_cost << " ("
            << (is_unit_cost ? "unit cost" : "general cost") << ")" << endl;
    outfile.close();
    cout << "Plan length: " << plan.size() << " step(s)." << endl;
    cout << "Plan cost: " << plan_cost << endl;
    ++num_previously_generated_plans;
}

void PlanManager::load_plan(Plan &plan, std::string path_to_plan_file,
        const std::unordered_map<std::string, OperatorID>& ops_by_names) const {
    ifstream planfile;
    planfile.open(path_to_plan_file);

    if (!planfile.is_open()) {
        cerr << "File is not open!" << endl;
        utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
    }
    string line;
    //cout << "-----------------------------------------------------------------------" << endl;
    while(std::getline(planfile, line)) {
        if (line.size() == 0 || line[0] == ';')
            continue;
        string op_name = line.substr(1, line.size()-2);
        //cout << op_name << endl;
        auto it = ops_by_names.find(op_name);
        if (it == ops_by_names.end()) {
            // Trying adding a trailing space
            string op_name_trailing_space = op_name + " ";
            it = ops_by_names.find(op_name_trailing_space);
            if (it == ops_by_names.end()) {
                cout << "#" << op_name << "#   Operator not found!!!" << endl;
                cout << "Operator names:" << endl;
                for (auto name : ops_by_names) {
                    cout << "#" << name.first << "#" << endl;
                }
                utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
            }
        }
        //cout << "Name " << it->first << endl;
        //cout << "Found operator " << it->second->get_name() << endl;
        plan.push_back(it->second);
    }
}

void PlanManager::load_plans(std::vector<Plan> &plans, const TaskProxy &task_proxy) const {
    // Creating a hashmap from operator names to pointers
    std::unordered_map<string, OperatorID> ops_by_names;
    OperatorsProxy operators = task_proxy.get_operators();
    for (OperatorProxy op : operators) {
        ops_by_names.insert(std::make_pair<string, OperatorID>(op.get_name(), op.get_ancestor_operator_id(tasks::g_root_task.get())));
        auto it = ops_by_names.find(op.get_name());
        if (it == ops_by_names.end()) {
            cout << "Problem adding operator!!!" << endl;
            utils::exit_with(utils::ExitCode::SEARCH_INPUT_ERROR);
        }
    }

    ostringstream filename;
    filename << plan_foldername << "/" << plan_filename;
    for (int plan_no=1; plan_no <= num_plans_to_read; ++plan_no) {
        ostringstream curr_filename;
        curr_filename << filename.str() << "." << plan_no;
        Plan plan;
        load_plan(plan, curr_filename.str(), ops_by_names);
        cout << "Loaded plan " << curr_filename.str() << endl;
        plans.push_back(plan);
    }
}
