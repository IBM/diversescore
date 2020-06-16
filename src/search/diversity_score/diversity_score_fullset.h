#ifndef DIVERSITY_SCORE_DIVERSITY_SCORE_FULLSET_H
#define DIVERSITY_SCORE_DIVERSITY_SCORE_FULLSET_H

#include "diversity_score.h"

class DiversityScoreFullset : public DiversityScore {

public:
    DiversityScoreFullset(const options::Options &opts);
    virtual ~DiversityScoreFullset() = default;

    virtual void compute_metrics();

};

//}

#endif
