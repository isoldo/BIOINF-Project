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
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
/*
	Program defines
	To be moved to a separate header file
*/

#define ARGC_VAL_MIN (1+4)

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
	std::vector<mhapRead_t*> lOvlp;
	std::vector<mhapRead_t*> rOvlp;
} BOGNode_t;

bool cmpNodes(mhapRead_t* a, mhapRead_t* b) {
	if (a->jaccardScore != b->jaccardScore) {
		return a->jaccardScore > b->jaccardScore;
	}
	return a->smm > b->smm;
}

/*
	Global variables
*/

std::map<int,BOGNode_t> nodes;
std::map<int,BOGNode_t> nodes_filtered;
std::vector<mhapRead_t> mhapReads;
std::set<int> readIxs;
std::vector<int> result;

/*
	local function declarations
*/
int isDoveTail(mhapRead_t* mRead);
int isLeftAligned(mhapRead_t* mRead);

static void printmhapreads(void);
static void printnodes(void);

int main(int argc, char** argv) {
	/*
		int main() variables
	*/
	int defaultOutput = 1;
	std::string mhapPath;
	std::string fastaPath;
	std::string outputPath = "output.fasta";
	double minJaccardScore = 1.0;
	int minJaccardIx;
	
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
		if (strcmp(argv[i],"-mhap") == 0) {
			mhapPath.assign(argv[++i]);
			// accept only *.mhap for now
			std::string inputExtension;
			std::string::size_type idx = mhapPath.rfind('.');
			if (idx != std::string::npos) {
				inputExtension = mhapPath.substr(idx);
				if (strcmp(inputExtension.c_str(),".mhap")) {
					std::cout << "Please provide a *.mhap file. (got " <<  inputExtension << ")" << std::endl << "Terminating." << std::endl;
					return 1;
				}
			} else {
				std::cout << "Provide input file extension!" << std::endl << "Terminating." << std::endl;
				return 1;
			}
//			std::cout << "Input path: " << inputPath << std::endl;
//			std::cout << "Input extension " << inputExtension << std::endl;
		} else if (strcmp(argv[i],"-fasta") == 0) {
			fastaPath.assign(argv[++i]);
			// accept only *.mhap for now
			std::string inputExtension;
			std::string::size_type idx = fastaPath.rfind('.');
			if (idx != std::string::npos) {
				inputExtension = fastaPath.substr(idx);
				if (strcmp(inputExtension.c_str(),".fasta")) {
					std::cout << "Please provide a *.fasta file. (got " <<  inputExtension << ")" << std::endl << "Terminating." << std::endl;
					return 1;
				}
			} else {
				std::cout << "Provide input file extension!" << std::endl << "Terminating." << std::endl;
				return 1;
			}
		} else if (strcmp(argv[i],"-o") == 0) {
			outputPath.assign(argv[++i]);
			defaultOutput = 0;
			std::ofstream outputFile (outputPath.c_str());
			if (!outputFile) {
				std::cout << "Cannot generate output file!" << std::endl << "Terminating." << std::endl;
			} else {
				outputFile.close();
			}
//			std::cout << "Output path: " << outputPath << std::endl;
		} else {
			std::cout << "invalid cli options" << std::endl << "Terminating" << std::endl;
			return 1;
		}
	}
	if (fastaPath.empty() || mhapPath.empty()) {
		std::cout << "Please provide both *.fasta and *.mhap files!" << std::endl;
//		return 1;
	}
	/*
		open input file (fasta or mhap)
			mhap is better, because conversion of ecoli_corrected.fasta to mhap (graphmap) takes about 70 minutes on 4 cores and 10 GB RAM
	*/
	// only mhap for now :)
	std::ifstream inputFile(mhapPath.c_str(),std::ifstream::in);
	if (!inputFile) {
		std::cout << mhapPath <<": No such file" << std::endl;
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
	std::cout << mhapReads.size() << std::endl;
	//printmhapreads();


	/*
		generate read and overlap contexts
		build the complete OVERLAP GRAPH
	*/
	for (int i = 0; i<mhapReads.size(); ++i) {
		mhapRead_t* cRead = &mhapReads[i];
		BOGNode_t bNodeF,bNodeS;
		int ixf = cRead->id.first;
		int ixs = cRead->id.second;
		
		// load existing map record
		std::map<int,BOGNode_t>::iterator it = nodes.find(ixf);
		if (it != nodes.end()) {
			bNodeF = nodes[ixf];
		} else {
			readIxs.insert(ixf);
		}
		it = nodes.find(ixs);
		if (it != nodes.end()) {
			bNodeS = nodes[ixs];
		} else {
			readIxs.insert(ixs);
		}
		// ignore non-dovetail overlaps
		if (isDoveTail(cRead)) {
			if (isLeftAligned(cRead)) {
				// add to left overlaps
				bNodeF.lOvlp.push_back(cRead);
				bNodeS.rOvlp.push_back(cRead);
			} else {
				// add to right overlaps
				bNodeF.rOvlp.push_back(cRead);
				bNodeS.lOvlp.push_back(cRead);
			}
		}
		nodes[ixf] = bNodeF;
		nodes[ixs] = bNodeS;
	}
	std::cout << nodes.size() << std::endl;
	//printnodes();
	for (std::map<int,BOGNode_t>::iterator it=nodes.begin(); it!=nodes.end(); ++it) {
		if (it->second.lOvlp.size() != 0 || it->second.rOvlp.size() != 0) {
			// these reads have at least one overlaps
			nodes_filtered[it->first] = it->second;
		} else {
			readIxs.erase(it->first);
		}
	}
	std::cout << "Cleaned: " << nodes.size() - nodes_filtered.size() << std::endl;
	for (std::map<int,BOGNode_t>::iterator it=nodes_filtered.begin(); it!=nodes_filtered.end(); ++it) {
		// sort by Jaccard score which is not exactly Jaccard score when mhap comes from graphmap but ok
		std::sort(it->second.lOvlp.begin(),it->second.lOvlp.end(),cmpNodes);
		std::sort(it->second.rOvlp.begin(),it->second.rOvlp.end(),cmpNodes);
		if (!it->second.lOvlp.empty()) {
			if (minJaccardScore > it->second.lOvlp[0]->jaccardScore) {
				minJaccardScore = it->second.lOvlp[0]->jaccardScore;
				minJaccardIx = it->first;
			}
		}// else {
			//minJaccardScore = 0.0;
			//minJaccardIx = it->first;
		//}
	}
	//std::cout << "MJI" << minJaccardIx << std::endl;

	//for (std::vector<mhapRead_t*>::iterator it = nodes_filtered[minJaccardIx].rOvlp.begin(); it!=nodes_filtered[minJaccardIx].rOvlp.end(); ++it) {
		//std::cout << (*it)->jaccardScore << std::endl;
	//}
	//return 0;

	// thats it, we have our densely populated overlap graph

	/*
		THE ALGORITHM
		Travel through the overlap graph, remove contained reads, remove cycles, pick best overlaps
		Create the best overlap graph
			1 fasta read per chromosome is optimum
	*/

	// pick the node with the worst left overlap Jaccard score - this will be our starting node
	std::cout << "#reads: " << readIxs.size() << std::endl;
	result.push_back(minJaccardIx);
	readIxs.erase(minJaccardIx);
	while(readIxs.size()) {
		int currentReadIx = result.back();
		std::cout << "Current Read: " << currentReadIx << std::endl;
		std::pair<int,int> indices;
		if (!nodes_filtered[currentReadIx].rOvlp.empty()) {
			indices = nodes_filtered[currentReadIx].rOvlp[0]->id;
			if (indices.first == currentReadIx && readIxs.find(indices.second)!=readIxs.end()) {
				result.push_back(indices.second);
				readIxs.erase(indices.second);
			} else if (readIxs.find(indices.first)!=readIxs.end()) {
				result.push_back(indices.first);
				readIxs.erase(indices.first);
			}
			else {
				break;
			}
		} else {
			std::cout << currentReadIx << "has no rOvlp" << std::endl;
			break;
		}
	}

	std::cout << "Unused reads: " << readIxs.size() << std::endl;

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
	std::ofstream outputFile (outputPath.c_str());
	outputFile << ">The result of the BOG goes here" << std::endl;
	outputFile.close();
	return 0;
}

// not definitive
int isDoveTail(mhapRead_t* mRead) {
	if (mRead->start.first != 0 && mRead->end.first != (mRead->l.first-1)) {
		// doesnt start at the beginning nor does it end at the end, this is a pratial overlap
		return 0;
	}
	if (mRead->start.second != 0 && mRead->end.second != (mRead->l.second-1)) {
		// doesnt start at the beginning nor does it end at the end, this is a pratial overlap
		return 0;
	}
	return 1;
}

// a left aligned overlap is the one where the first index overlaps form 0 to x, x<l-1
// the end of the second read overlaps with the beginning of the first read

// not definitive
int isLeftAligned(mhapRead_t* mRead) {
	return mRead->start.first == 0;
}

static void printmhapreads(void) {
	for (int i=0; i<mhapReads.size(); ++i) {
		std::cout << mhapReads[i].id.first << " " << mhapReads[i].id.second << " " << mhapReads[i].jaccardScore << " " << mhapReads[i].smm << " " << mhapReads[i].rc.first << " " << mhapReads[i].start.first << " " << mhapReads[i].end.first << " " << mhapReads[i].l.first << " " << mhapReads[i].rc.second << " " << mhapReads[i].start.second << " " << mhapReads[i].end.second << " " << mhapReads[i].l.second << std::endl;
	}
}

static void printnodes(void) {
	for (std::map<int,BOGNode_t>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
		std::cout << "Node: " << it->first << std::endl;
		std::cout << "lOvlp: ";
		for (int j=0; j<it->second.lOvlp.size(); ++j) {
			std::cout << it->second.lOvlp[j]->id.first << "-" << it->second.lOvlp[j]->id.second << ", ";
		}
		std::cout << std::endl;
		std::cout << "rOvlp: ";
		for (int j=0; j<it->second.rOvlp.size(); ++j) {
			std::cout << it->second.rOvlp[j]->id.first << "-" << it->second.rOvlp[j]->id.second << ", ";
		}
		std::cout << std::endl << std::endl;
	}
}
