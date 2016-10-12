#ifndef KARP_FASTAINDEX
#define KARP_FASTAINDEX

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>

struct findex_entry {
  uint64_t offset;
  uint64_t linelength;
  int file;
};

struct fastaIndex {

  std::unordered_map<std::string,int> ref_names;
  std::vector<findex_entry> values;
  std::vector<std::string>& transfasta;

fastaIndex(std::vector<std::string>& transfasta) : transfasta(transfasta) {
    ProcessFIndex(transfasta);
  }

  void ProcessFIndex(std::vector<std::string>& transfasta);
 
};

#endif
  
  
  
  
