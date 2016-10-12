#ifndef TABULATE_HPP
#define TABULATE_HPP

#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <zlib.h>
#include <algorithm>
#include <fstream>
#include <cmath>

#include <stdlib.h>
#include <stdio.h>
#include "common_k.h"

struct ResTable {

  std::vector<std::string> sfiles;
  std::unordered_map<std::string, std::unordered_map<unsigned int, double> > otu_table;
  std::unordered_map<unsigned int, std::string> tax_table;
  std::string& out_base;
  double min_f;

  ResTable(ProgramOptions& opt) : out_base(opt.out), min_f(opt.minimum_frequency_cutoff) {} 
  void tabulate(std::string& all_files);
  void output();

};

#endif
