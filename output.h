#ifndef _OUTPUT_H
#define _OUTPUT_H

#include "bogart.h"

void createFastaOutput(std::vector<int>& bog, std::vector<mhapRead_t>& reads, std::map<int,std::string>& input, std::string& output);

#endif