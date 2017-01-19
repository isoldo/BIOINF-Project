#ifndef _BOGART_H
#define _BOGART_H

#include <vector>
#include <utility>
#include <string>
#include <cstring>
#include <map>
#include <algorithm>
#include <iostream>

typedef struct {
	int 		verbose;
	int 		defaultOutput;
	int 		generateMhap;
	std::string mhapPath;
	std::string fastaPath;
	std::string generatorPath;
	std::string outputPath;
} CONFIG;

typedef struct {
	std::pair<int,int> 	id;
	std::pair<int,int> 	l;
	std::pair<int,int> 	start;
	std::pair<int,int> 	end;
	std::pair<int,int> 	rc;
	double 				jaccardScore;
	int 				smm;
} mhapRead_t;

typedef std::vector<mhapRead_t*> mhapList_t;
typedef std::pair<double,std::vector<int> > solution_t;

typedef struct {
	int 		thisNodeHasBeenVisited;
	double 		bestCumulativeJaccardScore;
	mhapList_t 	lOvlp;
	mhapList_t 	rOvlp;
	solution_t	bestSolution;
} BOGNode_t;

typedef std::map<int, BOGNode_t> overlapGraph_t;

#endif