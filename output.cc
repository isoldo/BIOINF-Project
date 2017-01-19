#include "output.h"

void createFastaOutput(solution_t& bog, std::map<int,std::string>& input, std::string& output) {
	std::vector<int>::iterator ixIt = bog.second.first.begin();
	std::vector<mhapRead_t*>::iterator rIt = bog.second.second.begin();
	for (; ixIt != bog.second.first.end() && rIt != bog.second.second.end(); ++ixIt,++rIt) {
		int ix = *ixIt;
		mhapRead_t read = *(*rIt);
		std::string fFasta,sFasta;

		if (ix == read.id.first) {
			// the first sequence is the first element of the pair
			// this should be the left read
			if (read.end.first != (read.l.first-1) || read.start.second != 0) {
				std::cout << "Nije ti dobra logika, sinko (" << read.rc.second << ") "<< read.start.second << " " << read.end.second << " " << read.l.second << std::endl;
			}
			// ARBITRARY
			// the overlapping part will be taken from the second read
			fFasta.assign(input[read.id.first].begin(),input[read.id.first].end());
			sFasta.assign(input[read.id.second].begin()+read.end.second,input[read.id.second].end());
		} else {
			// the first sequence is the second element of the pair
			if (read.end.second != (read.l.second-1) || read.start.first != 0) {
				//std::cout << "Nije ti dobra logika 2, sinko " << (read.end.second != (read.l.second-1)) << (read.start.first != 0) << std::endl;
			}
			fFasta.assign(input[read.id.second].begin(), input[read.id.second].end());
			sFasta.assign(input[read.id.first].begin()+read.end.first, input[read.id.first].end());
		}
		if (output.empty()){
			output += fFasta;
		}
		output += sFasta;
	}
}