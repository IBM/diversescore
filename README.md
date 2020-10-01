

# Diverse-score is a code that computes diversity score of a set of plans (or its subset) for a given metric.
## Implemented metrics are defined by:
1. Pairwise measure: stability, state, uniqueness (and their linear combination, obtained by specifying several measures). The parameters are (default values are *false*)
   1. *compute_stability_metric*
   2. *compute_states_metric*
   3. *compute_uniqueness_metric* 
2. Aggregation method over the pairs in the set of plans: average or minimum. The parameter is *aggregator_metric*, with values {*avg*, *min*}, default: *avg*. 
3. Treating plans as multisets or sets of actions (LPG-d planner treats plans as multisets of actions). The parameter is *plans_as_multisets* , default: *false*.

## Subset selection
Additional functionality allows for selecting a subset of plans of required size from a larger set of plans, with the aim of finding subsets with better scores. The code allows to post-process sets of plans to obtain solutions for **satisificing** and **bounded** diverse planning problems.

## Building
For building the code please use 
```
./build.py
```

**Note that the computation of bounded diversity score requires enabling CPLEX support in Fast Downward (see http://www.fast-downward.org/) and building the code with LP support.**


## Running

### Score of a given set
An example specification is provided in the script compute_stability.sh
```
# ./compute_stability.sh <domain> <problem> <plans-folder> <number-of-plans>
./compute_stability.sh domain.pddl problem.pddl found_plans 1000
```

Additional example specifications:
* *--diversity-score "score(compute_stability_metric=true, aggregator_metric=min, plans_as_multisets=true)"* <br/> for minimal pairwise stability as computed by LPG-d.
* *--diversity-score "score(compute_stability_metric=true, compute_states_metric=true, aggregator_metric=avg,plans_as_multisets=false)"* <br/> for average over linear combination of stability and state.

### Selecting subset for satisficing diverse planning
```
# ./compute_stability_subset.sh <domain> <problem> <plans-folder> <number-of-plans-in-folder> <number-of-plans-to-select> 
./compute_stability_subset.sh domain.pddl problem.pddl found_plans 1000 100
```
### Selecting subset for bounded diverse planning
```
# ./compute_stability_subset_bounded.sh <domain> <problem> <plans-folder> <number-of-plans-in-folder> <number-of-plans-to-select> <bound>
./compute_stability_subset_bounded.sh domain.pddl problem.pddl found_plans 1000 100 0.25
```

### Selecting subset for optimal diversity diverse planning
```
# ./compute_stability_subset_optimal.sh <domain> <problem> <plans-folder> <number-of-plans-in-folder> <number-of-plans-to-select> 
./compute_stability_subset_bounded.sh domain.pddl problem.pddl found_plans 1000 100 
```

## Licensing

DiverseScore is a Automated PDDL planning tool for computing the score
for a set of plans under specified diversity metrics.
Copyright (C) 2019  Michael Katz, IBM Research, USA.

Diverse score computation is built on top of Fast Downward. 
The computation code is located in the following folder:
* ./src/search/diversity_score

In addition, the code modifies the following files of Fast Downward:

* driver/run_components.py
* src/search/DownwardFiles.cmake
* src/search/command_line.{cc,h}
* src/search/search_space.{cc,h}
* src/search/plan_manager.{cc,h}
* src/search/planner.cc

The license for the extension of Fast Downward code is specified in the LICENSE file.

## Fast Downward
Fast Downward is a domain-independent planning system.

For documentation and contact information see http://www.fast-downward.org/.

The following directories are not part of Fast Downward as covered by this
license:

* ./src/search/ext

For the rest, the following license applies:

```
Fast Downward is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Fast Downward is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.
```
