#include <bogart.h>
#include <bestOverlap_DFS.h>

static bool solutionCmp(solution_t& a, solution_t& b);
static solution_t dfs(overlapGraph_t& graph, int ix);

void BOG_DFS(overlapGraph_t graph, std::set<int>& starters, std::vector<int>& result) {
	std::vector<solution_t> local_solutions;
	for (std::set<int>::iterator it = starters.begin(); it!= starters.end(); ++it) {
		local_solutions.push_back(dfs(graph,*it));
	}
	// sort local solutions
	std::sort(local_solutions.begin(),local_solutions.end(),solutionCmp);
	// store to result
	result = local_soultions[0].second;
}

static bool solutionCmp(solution_t& a, solution_t& b) {
	if (a.first != b.first) {
		return a.first > b.first;
	}
	return a.second.size() > b.second.size();
}

static solution_t dfs(overlapGraph_t& graph, int ix) {
	// all of the possible solutions from this node will be stored in this vector
	std::vector<solution_t> local_solutions;

	// iterate through all of the right overlaps of current node
	for (mhapList_t::iterator it = graph[ix].rOvlp.begin(); it!=graph[ix].rOvlp.end(); ++it) {
		mhapRead_t current_read;
		solution_t current_solution;
		int ovlp_ix;

		// get info on current overlap
		current_read = *(*it);

		// find the index of overlapping read
		if (current_read.id.first == ix) {
			ovlp_ix = current_read.id.second;
		} else {
			ovlp_ix = current_read.id.first;
		}

		if (graph[ovlp_ix].thisNodeHasBeenVisited == 0) {
			// best solution for the overlapping read has not been found yet
			graph[ovlp_ix].thisNodeHasBeenVisited = 1;
			current_solution = dfs(graph,ovlp_ix);
		} else {
			// best solution already exists
			current_solution = graph[ovlp_ix].bestSolution;
		}

		/*
			ARBITRARY METRIC
				find the best cumulative jaccard
		*/

		// multiply the jaccardScore of this overlap with cumulative jaccard score of the child node
		current_solution.first *= current_read.jaccardScore;
		// add this node into the path vector
		current_solution.second.push_back(ix);
		local_solutions.push_back(current_solution);
	}
	/*
		Leaf node
			this node must return metric-neutral value as jaccard score (1 because we are multiplying)
			this node must return only itself as the best path
	*/
	if (local_solutions.empty()) {
		solution_t tmp;
		std::vector<int> tmpv(1,ix);
		tmp.first = 1.0;
		tmp.second = tmpv;
		local_solutions.push_back(tmp);
	}
	std::sort(local_solutions.begin(),local_solutions.end(),solutionCmp);
	return local_solution[0];
}