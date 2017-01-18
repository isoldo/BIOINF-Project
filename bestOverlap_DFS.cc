void BOG_DFS(std::map<int,BOGNode_t>& graph, std::set<int>& starters, std::vector<int>& result);

void BOG_DFS(std::map<int,BOGNode_t>& graph, std::set<int>& starters, std::vector<int>& result) {
	std::vector<solution_t> local_solutions;
	for (std::set<int>::iterator it = starters.begin(); it!= starters.end(); ++it) {
		local_solutions.push_back(dfs(graph[*it],*it));
	}
	// sort local solutions
	// store to result
}

solution_t dfs(BOGNode_t& node, int ix);

solution_t dfs(BOGNode_t& node, int ix) {
	for (std::vector<mhapRead_t*>::iterator it = node.rOvlp.begin(); it!=node.rOvlp.end(); ++it) {
		// extract info from read context
		mhapRead_t currentRead = *(*it);
		if (currentRead.id.first == ix) {
			// child of this node is found on the second index
		} else {
			// child of this node is found on the second index
		}

	}
}