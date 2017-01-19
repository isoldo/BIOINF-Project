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

#include "bogart.h"
#include "bestOverlap_DFS.h"
#include "output.h"
#include <stdlib.h>


#include <sstream>
#include <fstream>


/*
	Program defines
	To be moved to a separate header file
*/

#define ARGC_VAL_MIN (1+4)


bool cmpNodes(mhapRead_t* a, mhapRead_t* b) {
	if (a->jaccardScore != b->jaccardScore) {
		return a->jaccardScore > b->jaccardScore;
	}
	return a->smm < b->smm;
}

/*
	Global variables
*/

/*
	local function declarations
*/
static int callMhapGenerator(CONFIG& config);
static int parseMhapInput(std::ifstream& input, std::vector<mhapRead_t>& target, CONFIG& config);

static int buildDenseGraph(std::vector<mhapRead_t>& src, std::map<int,BOGNode_t>& dest, std::set<int>& ixSet);
static int filterDeadReads(std::map<int,BOGNode_t>& graph, std::map<int,BOGNode_t>& graph_filtered, std::set<int>& ixSet, CONFIG& config);
static int findStarters(std::map<int,BOGNode_t>& graph, std::set<int>& starters, double& minJaccardScore, int& minJaccardIx);

static int isDoveTail(mhapRead_t* mRead);
static int isLeftAligned(mhapRead_t* mRead);


int main(int argc, char** argv) {
	/*
		int main() variables
	*/
	CONFIG sysConfig;

	overlapGraph_t nodes;
	overlapGraph_t nodes_filtered;

	std::vector<mhapRead_t> mhapReads;
	std::set<int> readIxs;
	std::set<int> startNodeIx;
	std::vector<std::vector<int> > results;
	std::map<int,std::string> fastaReads;
	double minimumJaccardScore = 1.0;
	int minimumJaccardScoreIndex = -1;

	sysConfig.generateMhap = 0;
	
	/*
		parse input arguments
	*/
	
	if (argc < ARGC_VAL_MIN) {
		std::cout << "Too few parameters! Expected at least " << ARGC_VAL_MIN-1 << ", got " << argc-1 << ".\r\nTerminating." << std::endl;
		return 1;
	}
	for (int i=1; i<argc; ++i) {
		if (strcmp(argv[i],"-in_mhap") == 0) {
			sysConfig.mhapPath.assign(argv[++i]);
			std::string inputExtension;
			std::string::size_type idx = sysConfig.mhapPath.rfind('.');
			if (idx != std::string::npos) {
				inputExtension = sysConfig.mhapPath.substr(idx);
				if (strcmp(inputExtension.c_str(),".mhap")) {
					std::cout << "Please set target to be a *.mhap file. (got " <<  inputExtension << ")" << std::endl << "Terminating." << std::endl;
					return 1;
				}
			} else {
				std::cout << "Provide input file extension!" << std::endl << "Terminating." << std::endl;
				return 1;
			}
		} else if (strcmp(argv[i],"-fasta") == 0) {
			sysConfig.fastaPath.assign(argv[++i]);
			std::string inputExtension;
			std::string::size_type idx = sysConfig.fastaPath.rfind('.');
			if (idx != std::string::npos) {
				inputExtension = sysConfig.fastaPath.substr(idx);
				if (strcmp(inputExtension.c_str(),".fasta")) {
					std::cout << "Please provide a *.fasta input file. (got " <<  inputExtension << ")" << std::endl << "Terminating." << std::endl;
					return 1;
				}
			} else {
				std::cout << "Provide input file extension!" << std::endl << "Terminating." << std::endl;
				return 1;
			}
		} else if (strcmp(argv[i],"-o") == 0) {
			sysConfig.outputPath.assign(argv[++i]);
			sysConfig.defaultOutput = 0;
			
		} else if (strcmp(argv[i],"-mhap_gen") == 0) {
			sysConfig.generatorPath.assign(argv[++i]);
			sysConfig.generateMhap = 1;
		} else if (strcmp(argv[i],"-v") == 0) {
			sysConfig.verbose = 1;
		} else {
			std::cout << "invalid cli options" << std::endl << "Terminating" << std::endl;
			return 1;
		}
	}

	if (sysConfig.fastaPath.empty() || sysConfig.mhapPath.empty()) {
		std::cout << "Please provide both *.fasta and *.mhap files!" << std::endl;
		return 1;
	}
	/*
		open input file (fasta or mhap)
			mhap is better, because conversion of ecoli_corrected.fasta to mhap (graphmap) takes about 70 minutes on 4 cores and 10 GB RAM
	*/
	if (sysConfig.generateMhap) {
		// TODO check for return value
		callMhapGenerator(sysConfig);
	}

	std::ifstream inputFile(sysConfig.mhapPath.c_str(),std::ifstream::in);
	if (!inputFile) {
		std::cout << sysConfig.mhapPath <<": No such file" << std::endl;
		return 1;
	}
	
	/*
		Parse input file
		Save to a vector of reads
	*/

	parseMhapInput(inputFile,mhapReads,sysConfig);

	/*
		Generate read and overlap contexts
		Build the COMPLETE OVERLAP GRAPH
	*/

	buildDenseGraph(mhapReads,nodes,readIxs);
	
	/*
		Delete nodes without any overlaps
	*/

	filterDeadReads(nodes,nodes_filtered,readIxs,sysConfig);

	/*
		Find possible starter nodes 
	*/

	findStarters(nodes_filtered,startNodeIx,minimumJaccardScore,minimumJaccardScoreIndex);

	// if there are no nodes without left overlaps, use the one with the lowest jaccard score
	if (startNodeIx.empty()) {
		startNodeIx.insert(minimumJaccardScoreIndex);
	}

	/*
		THE ALGORITHM
		Plain simple DFS. Metrics: maximize the product of jaccard score along the path.
		Room for improvement: set a low boundary for the cumulative score, create snippets with good cumulative jaccard score
	*/
	std::vector<int> bestOverlapPath;
	BOG_DFS(nodes_filtered,startNodeIx,bestOverlapPath);

	/*
		Take the BOG and write it in FASTA format
	*/

	std::ifstream fastaInput(sysConfig.fastaPath.c_str(),std::ifstream::in);
	if (!fastaInput) {
		std::cout << sysConfig.fastaPath << ": No such file" << std::endl;
		return 1;
	}
	int counter = 0;
	std::string line;
	while(std::getline(fastaInput,line)) {
		std::istringstream inStream(line);
		std::string currLine;
		inStream >> currLine;
		if (currLine[0] == '>') {
			++ counter;
		} else {
			fastaReads[counter] += currLine;
		}
	}
	std::string finalFastaResult;

	createFastaOutput(bestOverlapPath,mhapReads,fastaReads,finalFastaResult);

	/*
		Notify about success, return 0
	*/
	std::cout << "Output to: " << sysConfig.outputPath << std::endl;
	
	std::ofstream outputFile (sysConfig.outputPath.c_str());
	outputFile << ">The result of the BOG goes here" << std::endl;
	outputFile << finalFastaResult << std::endl;
	outputFile.close();
	return 0;
}

static int callMhapGenerator(CONFIG& config) {
	std::string local_g_path;
	local_g_path.assign(config.generatorPath);
	local_g_path += " owler -r ";
	local_g_path += config.fastaPath;
	local_g_path += " -d ";
	local_g_path += config.fastaPath;
	local_g_path += " -o ";
	local_g_path += config.mhapPath;
	if (config.verbose) {
		std::cout << "Calling " << local_g_path << std::endl;
	}
	return system(local_g_path.c_str());
}

static int parseMhapInput(std::ifstream& input, std::vector<mhapRead_t>& target, CONFIG& config) {
	std::string line;
	while (std::getline(input,line)) {
		std::istringstream inStream(line);
		int ix1,ix2,rc1,rc2,s1,e1,l1,s2,e2,l2;
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
		// push_back to a vector of reads
		target.push_back(cRead);
	}
	int tSize = target.size();
	if (config.verbose) {
		std::cout << "Parsed " << tSize << " reads." << std::endl;
	}
	return tSize;
}

static int buildDenseGraph(std::vector<mhapRead_t>& src, std::map<int,BOGNode_t>& dest, std::set<int>& ixSet) {
	for (int i = 0; i<src.size(); ++i) {
		mhapRead_t* cRead = &src[i];
		BOGNode_t bNodeF,bNodeS;
		int ixf = cRead->id.first;
		int ixs = cRead->id.second;

		// load existing map record
		std::map<int,BOGNode_t>::iterator it;

		it = dest.find(ixf);
		if (it != dest.end()) {
			bNodeF = dest[ixf];
		} else {
			// kill this with constructor
			bNodeF.thisNodeHasBeenVisited = 0;
			bNodeF.bestCumulativeJaccardScore = 0.0;
			ixSet.insert(ixf);
		}

		it = dest.find(ixs);
		if (it != dest.end()) {
			bNodeS = dest[ixs];
		} else {
			// kill this with constructor
			bNodeS.thisNodeHasBeenVisited = 0;
			bNodeS.bestCumulativeJaccardScore = 0.0;
			ixSet.insert(ixs);
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
		dest[ixf] = bNodeF;
		dest[ixs] = bNodeS;
	}
	return dest.size();
}

static int filterDeadReads(std::map<int,BOGNode_t>& graph, std::map<int,BOGNode_t>& graph_filtered, std::set<int>& ixSet, CONFIG& config) {
	for (std::map<int,BOGNode_t>::iterator it=graph.begin(); it!=graph.end(); ++it) {
		if (it->second.lOvlp.size() != 0 || it->second.rOvlp.size() != 0) {
			// these reads have at least one overlaps
			graph_filtered[it->first] = it->second;
		} else {
			ixSet.erase(it->first);
			if (config.verbose) {
				std::cout << "Read " << it->first << " is dead" << std::endl;
			}
		}
	}
}

static int findStarters(std::map<int,BOGNode_t>& graph, std::set<int>& starters, double& minJaccardScore, int& minJaccardIx) {
	for (std::map<int,BOGNode_t>::iterator it=graph.begin(); it!=graph.end(); ++it) {
		// sort by Jaccard score which is not exactly Jaccard score when mhap comes from graphmap but ok
		std::sort(it->second.lOvlp.begin(),it->second.lOvlp.end(),cmpNodes);
		std::sort(it->second.rOvlp.begin(),it->second.rOvlp.end(),cmpNodes);
		if (!it->second.lOvlp.empty()) {
			if (minJaccardScore > it->second.lOvlp[0]->jaccardScore) {
				minJaccardScore = it->second.lOvlp[0]->jaccardScore;
				minJaccardIx = it->first;
			}
		}else {
			minJaccardScore = 0.0;
			minJaccardIx = it->first;
			starters.insert(it->first);
		}
	}
}


// not definitive
static int isDoveTail(mhapRead_t* mRead) {
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
static int isLeftAligned(mhapRead_t* mRead) {
	return mRead->start.first == 0;
}
