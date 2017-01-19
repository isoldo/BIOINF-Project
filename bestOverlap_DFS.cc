#include "bogart.h"
#include "bestOverlap_DFS.h"

static bool solutionCmp(solution_t a, solution_t b);
static solution_t dfs(overlapGraph_t& graph, int ix);

void BOG_DFS(overlapGraph_t& graph, std::set<int>& starters, solution_t& result) {
	std::vector<solution_t> local_solutions;
	for (std::set<int>::iterator it = starters.begin(); it!= starters.end(); ++it) {
		local_solutions.push_back(dfs(graph,*it));
		std::cout << local_solutions[local_solutions.size()-1].second.first.size() << std::endl;
	}
	// sort local solutions
	std::sort(local_solutions.begin(),local_solutions.end(),solutionCmp);
	// store to result
	result = local_solutions[0];
}

static bool solutionCmp(solution_t a, solution_t b) {
	return a.first > b.first;
}

static solution_t dfs(overlapGraph_t& graph, int ix) {
	// all of the possible solutions from this node will be stored in this vector
	std::vector<solution_t> local_solutions;
	// std::cout << "dfs " << ix << " ";
	// iterate through all of the right overlaps of current node
	for (mhapList_t::iterator it = graph[ix].rOvlp.begin(); it!=graph[ix].rOvlp.end(); ++it) {
		mhapRead_t current_read;
		solution_t current_solution;
		int ovlp_ix;
		int tmp;

		// get info on current overlap
		current_read = *(*it);

		// find the index of overlapping read
		if (current_read.id.first == ix) {
			ovlp_ix = current_read.id.second;
			tmp = current_read.id.first;
		} else {
			ovlp_ix = current_read.id.first;
			tmp = current_read.id.second;
		}
		if (graph[ovlp_ix].thisNodeHasBeenVisited == 0) {
			// best solution for the overlapping read has not been found yet
			graph[ovlp_ix].thisNodeHasBeenVisited = 2;
			//std::cout << ovlp_ix << " to be processed" << std::endl;
			current_solution = dfs(graph,ovlp_ix);
		} else if (graph[ovlp_ix].thisNodeHasBeenVisited == 1){
			// best solution already exists
			//std::cout << ovlp_ix << " has a best solution" << std::endl;
			current_solution = graph[ovlp_ix].bestSolution;
		} else {
			//std::cout << ovlp_ix << " cycle with " << tmp << std::endl;
			continue;
		}
		// TODO break cycles !
		/*
			ARBITRARY METRIC
				find the best cumulative jaccard
		*/
		graph[ovlp_ix].thisNodeHasBeenVisited = 1;
		// multiply the jaccardScore of this overlap with cumulative jaccard score of the child node
		current_solution.first += current_read.jaccardScore;// + current_solution.first*current_solution.second.size();
		//std::cout << std::endl << ovlp_ix << " visited! Solution " << current_solution.first << std::endl;
		// add this node into the path vector
		current_solution.second.first.push_back(ix);
		current_solution.second.second.push_back(*it);
		//current_solution.first /= current_solution.second.size();
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
		std::vector<mhapRead_t*> tmpvr;
		tmp.first = 1.0;
		tmp.second.first = tmpv;
		tmp.second.second = tmpvr;
		local_solutions.push_back(tmp);
		//std::cout << "Leaf at " << ix << std::endl;
	}
	std::sort(local_solutions.begin(),local_solutions.end(),solutionCmp);
	return local_solutions[0];
}