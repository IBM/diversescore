#ifndef DIVERSITY_SCORE_DIVERSITY_SCORE_SUBSET_H
#define DIVERSITY_SCORE_DIVERSITY_SCORE_SUBSET_H

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

enum class PairwiseStability {
	DFS,
	MIP,
    MIP_EXTERNAL,
	MAX_CLIQUES
};

class DiversityScoreSubset : public DiversityScore {

    int plans_subset_size;
    bool exact_method;
    bool dump_plans;

    void compute_metrics_greedy(bool stability, bool state, bool uniqueness,
            std::vector<size_t>& selected_plan_indexes);

    void compute_metrics_mip_external(bool stability, bool state, bool uniqueness,
            std::vector<size_t>& selected_plan_indexes);

    std::pair<float, size_t> find_best_next_candidate(bool stability, bool state, bool uniqueness,
            const std::vector<size_t>& selected_plan_indexes, std::list<size_t>& candidates);

protected:

public:
    DiversityScoreSubset(const options::Options &opts);
    virtual ~DiversityScoreSubset() = default;

    virtual void compute_metrics();

};

//}

#endif
