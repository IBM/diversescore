#ifndef PLAN_MANAGER_H
#define PLAN_MANAGER_H

#include <string>
#include <vector>
#include <unordered_map>

class OperatorID;
class TaskProxy;

using Plan = std::vector<OperatorID>;

class PlanManager {
    std::string plan_filename;
    std::string plan_foldername;
    int num_plans_to_read;
    int num_previously_generated_plans;
    bool is_part_of_anytime_portfolio;


    void load_plan(Plan &plan, std::string path_to_plan_file,
        const std::unordered_map<std::string, OperatorID>& ops_by_names) const;
public:
    PlanManager();

    void set_plan_filename(const std::string &plan_filename);
    void set_plan_foldername(const std::string &plan_foldername);
    void set_num_plans_to_read(int num_plans) { num_plans_to_read = num_plans; };

    void set_is_part_of_anytime_portfolio(bool is_part_of_anytime_portfolio);
    void set_num_previously_generated_plans(int num_previously_generated_plans);

    /*
      Set generates_multiple_plan_files to true if the planner can find more than
      one plan and should number the plans as FILENAME.1, ..., FILENAME.n.
    */
    void save_plan(
        const Plan &plan, const TaskProxy &task_proxy,
        bool generates_multiple_plan_files = false);

    void load_plans(std::vector<Plan> &plans, const TaskProxy &task_proxy) const;

};

extern int calculate_plan_cost(const Plan &plan, const TaskProxy &task_proxy);

#endif
