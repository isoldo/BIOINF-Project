/*
MIT License

Copyright (c) 2017 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
/*
	Program defines
	To be moved to a separate header file
*/

#define ARGC_VAL_MIN (1+2)

/*
	Custom structs
	To be moved to spearate header file
*/

typedef struct {
	std::pair<int,int> id;
	std::pair<int,int> l;
	std::pair<int,int> start;
	std::pair<int,int> end;
	std::pair<int,int> rc;
	double jaccardScore;
	int smm;
} mhapRead_t;

typedef struct {
	int id;
	std::vector<mhapRead_t*> lOvlp;
	std::vector<mhapRead_t*> rOvlp;
} BOGNode_t;

/*
	Global variables
*/

std::vector<BOGNode_t> nodes;
std::vector<mhapRead_t> mhapReads;

int main(int argc, char** argv) {
	/*
		int main() variables
	*/
	int defaultOutput = 1;
	std::string inputPath;
	std::string outputPath = "output.fasta";
	
	/*
		parse input arguments
	*/
	
	if (argc < ARGC_VAL_MIN) {
		std::cout << "Too few parameters! Expected at least " << ARGC_VAL_MIN-1 << ", got " << argc-1 << ".\r\nTerminating." << std::endl;
		return 1;
	} else {
//		std::cout << "Input parameters OK" << std::endl;
	}
	for (int i=1; i<argc; ++i) {
		if (strcmp(argv[i],"-i") == 0) {
			inputPath.assign(argv[++i]);
			// accept only *.mhap for now
			std::string inputExtension;
			std::string::size_type idx = inputPath.rfind('.');
			if (idx != std::string::npos) {
				inputExtension = inputPath.substr(idx);
				if (strcmp(inputExtension.c_str(),".mhap")) {
					std::cout << "Only *.mhap is currently supported!" << std::endl << "Terminating." << std::endl;
					return 1;
				}
			} else {
				std::cout << "Provide input file extension!" << std::endl << "Terminating." << std::endl;
				return 1;
			}
//			std::cout << "Input path: " << inputPath << std::endl;
//			std::cout << "Input extension " << inputExtension << std::endl;
		} else if (strcmp(argv[i],"-o") == 0) {
			outputPath.assign(argv[++i]);
			defaultOutput = 0;
//			std::cout << "Output path: " << outputPath << std::endl;
		} else {
			std::cout << "invalid cli options" << std::endl << "Terminating" << std::endl;
			return 1;
		}
	}

	/*
		open input file (fasta or mhap)
			mhap is better, because conversion of ecoli_corrected.fasta to mhap (graphmap) takes about 70 minutes on 4 cores and 10 GB RAM
	*/
	// only mhap for now :)
	std::ifstream inputFile(inputPath.c_str(),std::ifstream::in);
	if (!inputFile) {
		std::cout << inputPath <<": No such file" << std::endl;
		return 1;
	}
	
	/*
		parse input file
	*/
	std::string line;
	while (std::getline(inputFile,line)) {
		std::istringstream inStream(line);
		int ix1,ix2,smm,rc1,rc2,s1,e1,l1,s2,e2,l2;
		double js;
		mhapRead_t cRead;
		inStream >> ix1 >> ix2;
		inStream >> cRead.jaccardScore;
		inStream >> cRead.smm;
		inStream >> rc1 >> s1 >> e1 >> l1;
		inStream >> rc2 >> s2 >> e2 >> l2;
		cRead.id = std::make_pair(ix1,ix2);
		cRead.l = std::make_pair(l1,l2);
		cRead.start = std::make_pair(s1,s2);
		cRead.end = std::make_pair(e1,e2);
		cRead.rc = std::make_pair(rc1,rc2);
//		std::cout << ix1 << " " << ix2 << " " << js << " " << smm << " " << rc1 << " " << s1 << " " << e1 << " " << l1 << " " << rc2 << " " << s2 << " " << e2 << " " << l2 << std::endl;
		// push_back to a vector of reads
		mhapReads.push_back(cRead);
	}

	/*
		generate read and overlap contexts
	*/

	/*
		build the complete OVERLAP GRAPH
	*/

	/*
		THE ALGORITHM
		Travel through the overlap graph, remove contained reads, remove cycles, pick best overlaps
		Create the best overlap graph
			1 fasta read per chromosome is optimum
	*/

	/*
		Take the BOG and write i in FASTA format
	*/

	/*
		Notify about success, return 0
	*/
	if (defaultOutput) {
		std::cout << "Output to default: " << outputPath << std::endl;
	} else {
		std::cout << "Output to custom: " << outputPath << std::endl;
	}
	return 0;
}
