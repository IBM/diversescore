#ifndef DIVERSITY_SCORE_DIVERSITY_SCORE_SUBSET_BOUNDED_H
#define DIVERSITY_SCORE_DIVERSITY_SCORE_SUBSET_BOUNDED_H

#include "diversity_score.h"

#include <unordered_set>
#include <unordered_map>
#include <list>

namespace options {
class OptionParser;
class Options;
}

//namespace diversity_score {

class CliqueDFSNode {
    int n;
    std::shared_ptr<CliqueDFSNode> parent;
public:
    CliqueDFSNode(int _n, std::shared_ptr<CliqueDFSNode> _parent) : n(_n), parent(_parent) {};
    void get_clique(std::vector<size_t>& clique) const {
        if (n < 0)
            return; 
        if (parent)
            parent->get_clique(clique);
        clique.push_back((size_t)n);
    }
};

enum class ComputationalMethod {
	DFS,
	MIP,
    MIP_EXTERNAL,
	MAX_CLIQUES
};

class DiversityScoreSubsetBounded : public DiversityScore {

    float metric_bound;
    ComputationalMethod computation_method;
    int plans_subset_size;
    // bool exact_mip_external;
    bool dump_plans;

    void generate_mip_file(std::vector<size_t>& selected_plan_indexes);

    void compute_metric_maxclique(std::vector<size_t>& selected_plan_indexes);
    void compute_metric_mip_external(std::vector<size_t>& selected_plan_indexes);
    void compute_metric_mip(std::vector<size_t>& selected_plan_indexes);
    void compute_metric_dfs(std::vector<size_t>& selected_plan_indexes);

protected:

public:
    DiversityScoreSubsetBounded(const options::Options &opts);
    virtual ~DiversityScoreSubsetBounded() = default;

    virtual void compute_metrics();

};

//}

#endif
