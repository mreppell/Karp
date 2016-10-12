#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

#include "fastaIndex.h"
#include "common_k.h"

void fastaIndex::ProcessFIndex(std::vector<std::string>& transfasta) {

  std::string delim = "\t";

  int vec_index = 0;
  for (int ii=0;ii<transfasta.size();++ii) {
    std::string fai = transfasta[ii] + ".fai";
    const char* fn = fai.c_str();

    //Check if fasta file has been modified since index was created    
     struct stat fai_stat;
     struct stat fa_stat;
     stat(fai.c_str(), &fai_stat);
     stat(transfasta[ii].c_str(), &fa_stat);
      if (fai_stat.st_mtime < fa_stat.st_mtime) {
       std::cerr << "\nWARNING, reference fasta " << transfasta[ii] << " modified after creation of index " << transfasta[ii] << ".fai, if the fasta file has been modified such that the index is no longer valid segmentation faults and unstable behavior may occur when the program is run.\n\n";
       }

    std::ifstream myfile (fn);
    std::string line;
     
    //Read fasta index file into findex
    if (myfile.is_open()) {
      while (getline(myfile,line)) {
	std::vector<std::string> inter_line;
	split(line,delim,inter_line);
	std::unordered_map<std::string,int>::const_iterator got = ref_names.find(inter_line[0]);
	if (got == ref_names.end() ) {

	  findex_entry fe;
	  fe.offset = std::strtoull(inter_line[2].c_str(),NULL,10);
	  fe.linelength = std::strtoull(inter_line[3].c_str(),NULL,10);
	  fe.file = ii;
	  values.push_back(fe);
	  std::pair<std::string, int> curr_entry(inter_line[0],vec_index);
	  ref_names.insert(curr_entry);
	  ++vec_index;
	}
      }
      myfile.close();
    } else {
      std::cerr << "Unable to open fasta reference, " << fai << std::endl;
    }
  }
}
