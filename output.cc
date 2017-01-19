#include "output.h"

void createFastaOutput(std::vector<int>& bog, std::vector<mhapRead_t>& reads, std::map<int,std::string>& input, std::string& output) {
	std::vector<int>::iterator first = bog.begin();
	std::vector<int>::iterator second = first + 1;

	while (second != bog.end()) {
		int ix1 = std::min(*first,*second);
		int ix2 = std::max(*first,*second);
		int i = 0;
		mhapRead_t mr;
		std::string firstFasta, secondFasta;
		do {
			mr = reads[i++];
		} while (mr.id.first != ix1 && mr.id.second != ix2 || i < reads.size());
		if (ix1 == mr.id.first) {
			if (mr.end.first != mr.l.first-1) {
				std::cout << "Warning: First read(" << mr.id.first << ")  does not end at its end! (ovlp: " << mr.id.first << "," << mr.id.second << ")" << std::endl;
			}
			if (mr.start.second != 0) {
				std::cout << "Warning: Second read does not start from its start! (ovlp: " << mr.id.first << "," << mr.id.second << ")" << std::endl;
			}
			firstFasta.assign(input[ix1].begin(), input[ix1].begin()+mr.start.first);
			secondFasta.assign(input[ix2].begin(),input[ix2].end());
		} else {
			if (mr.end.second != mr.l.second-1) {
				std::cout << "Warning: First read (" << mr.id.second << ") does not end at its end! (ovlp: " << mr.id.first << "," << mr.id.second << ")" << std::endl;
			}
			if (mr.start.first != 0) {
				std::cout << "Warning: Second read does not start from its start! (ovlp: " << mr.id.first << "," << mr.id.second << ")" << std::endl;
			}
			firstFasta.assign(input[ix2].begin(), input[ix2].begin()+mr.start.second);
			secondFasta.assign(input[ix1].begin(),input[ix1].end());
		}
		output += firstFasta;
		output += secondFasta;
		++first;
		++second;
	}
}