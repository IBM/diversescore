#ifndef DIVERSITY_SCORE_DIVERSITY_SCORE_SUBSET_OPTIMAL_LONG_H
#define DIVERSITY_SCORE_DIVERSITY_SCORE_SUBSET_OPTIMAL_LONG_H

#include "diversity_score.h"

#include <unordered_set>
#include <unordered_map>
#include <list>

namespace options {
class OptionParser;
class Options;
}

//namespace diversity_score {

enum class ComputationalMethod {
	MIP,
    MIP_EXTERNAL,
	MAX_CLIQUES
};

class DiversityScoreSubsetOptimalLong : public DiversityScore {

    ComputationalMethod computation_method;
    int plans_subset_size;
    // bool exact_mip_external;
    bool dump_plans;

    void generate_mip_file(std::vector<size_t>& selected_plan_indexes);

    void compute_metric_mip_external(std::vector<size_t>& selected_plan_indexes);
    void compute_metric_mip(std::vector<size_t>& selected_plan_indexes);

    size_t get_binary_var_index(size_t plan_index) const;
    size_t get_binary_var_index(size_t plan_index1, size_t plan_index2) const;

protected:

public:
    DiversityScoreSubsetOptimalLong(const options::Options &opts);
    virtual ~DiversityScoreSubsetOptimalLong() = default;

    virtual void compute_metrics();

};

//}

#endif
