#include "diversity_score_fullset.h"

#include "../option_parser.h"
#include "../plugin.h"

using namespace std;

DiversityScoreFullset::DiversityScoreFullset(const Options &opts) : DiversityScore(opts) { }

void DiversityScoreFullset::compute_metrics() {
    compute_metrics_exact_set();
}

static shared_ptr<DiversityScore> _parse_score(OptionParser &parser) {
    add_diversity_score_options_to_parser(parser);
    Options opts = parser.parse();

    shared_ptr<DiversityScoreFullset> engine;
    if (!parser.dry_run()) {
        engine = make_shared<DiversityScoreFullset>(opts);
    }

    return engine;
}

static Plugin<DiversityScore> _plugin("score", _parse_score);
